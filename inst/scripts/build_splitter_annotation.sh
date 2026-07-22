#!/usr/bin/env bash
# =============================================================================
# build_splitter_annotation.sh  —  reditR pipeline template
# =============================================================================
# Purpose:
#   Build the barcode-to-cluster annotation TSV required by
#   bam_extract_barcode_reads_commandline_chr_V2.py (Step 1 of
#   scrna_preprocessing.sh).
#
#   Joins CellRanger filtered barcodes with cluster labels from an external
#   cluster assignment file (e.g. from SCEA or Seurat), normalises barcode
#   format, and drops clusters below a minimum cell count to avoid pseudobulks
#   too sparse for reliable editing detection.
#
# Output format (required by the barcode-splitter):
#   <barcode>-1<TAB>cluster_<N>
#
# Usage:
#   bash build_splitter_annotation.sh
#
#   Run on a login node — fast, no compute needed.
#
# Pipeline position:
#   CellRanger output          ─┐
#   Cluster labels (SCEA etc.) ─┴─► build_splitter_annotation.sh
#                                         │
#                                         ▼  annotation.tsv
#                              scrna_preprocessing.sh
#                                         │
#                                         ▼
#                              extracting_read_counts.sh
#                                         │
#                                         ▼
#                              reditR::filter_editing_sites()
# =============================================================================

set -euo pipefail

# ─── CONFIGURE THESE PATHS ────────────────────────────────────────────────────

# Path to CellRanger's filtered barcodes file (gzipped)
CR_BARCODES=<REPLACE: /path/to/cellranger_out/sample/outs/filtered_feature_bc_matrix/barcodes.tsv.gz>

# Tab-separated cluster label file (no header):
#   column 1: cell ID in format <SAMPLE_PREFIX>-<barcode>
#   column 2: cluster assignment at the chosen resolution (K)
# This can be exported from Seurat, Scanpy, or downloaded from SCEA.
CLUSTER_FILE=<REPLACE: /path/to/cluster_labels.tsv>

# Sample prefix used in the cluster file cell IDs (e.g. SAMEA12345678)
SAMPLE_PREFIX=<REPLACE: YOUR_SAMPLE_PREFIX>

# Cluster resolution column to use (the K value in your cluster file)
SELECTED_K=20

# Minimum cells per cluster — clusters below this are dropped
# (too few cells → sparse pseudobulk → unreliable editing detection)
MIN_CELLS=30

# Where to write the annotation TSV
OUT_DIR=<REPLACE: /path/to/output/directory>

# ──────────────────────────────────────────────────────────────────────────────

mkdir -p "$OUT_DIR"
TMPDIR=$(mktemp -d -t splitter_ann_XXXX)
trap "rm -rf $TMPDIR" EXIT

echo "============================================================"
echo "Building splitter annotation"
echo "============================================================"
echo "CellRanger barcodes: $CR_BARCODES"
echo "Cluster file:        $CLUSTER_FILE"
echo "Sample prefix:       $SAMPLE_PREFIX"
echo "Cluster resolution:  K=$SELECTED_K"
echo "Min cells/cluster:   $MIN_CELLS"
echo "Output dir:          $OUT_DIR"
echo "============================================================"

# ----- Step 1: extract CellRanger filtered barcodes --------------------------
echo
echo "[1/5] Extracting CellRanger barcodes..."
zcat "$CR_BARCODES" > "$TMPDIR/cellranger_barcodes.txt"
N_CR=$(wc -l < "$TMPDIR/cellranger_barcodes.txt")
echo "      $N_CR CellRanger barcodes"

# ----- Step 2: find which columns belong to this sample ----------------------
echo
echo "[2/5] Locating sample columns in cluster file..."
SAMPLE_COLS=$(head -1 "$CLUSTER_FILE" | tr '\t' '\n' | \
              awk -v sp="$SAMPLE_PREFIX" '$0 ~ "^"sp"-" {print NR}' | \
              paste -sd ',')
N_SAMPLE_COLS=$(echo "$SAMPLE_COLS" | tr ',' '\n' | wc -l)
echo "      $N_SAMPLE_COLS cells found for sample $SAMPLE_PREFIX"

# ----- Step 3: extract barcode + cluster pairs at chosen resolution ----------
echo
echo "[3/5] Extracting (barcode, cluster) pairs at K=$SELECTED_K..."

# Header row gives cell IDs; the K=SELECTED_K row gives cluster IDs.
# Both restricted to this sample's columns, then transposed.
head -1 "$CLUSTER_FILE" | cut -f $SAMPLE_COLS | tr '\t' '\n' \
    > "$TMPDIR/cell_ids.txt"
awk -F'\t' -v K=$SELECTED_K '$1=="TRUE" && $2==K' "$CLUSTER_FILE" | \
    cut -f $SAMPLE_COLS | tr '\t' '\n' \
    > "$TMPDIR/clusters.txt"

paste "$TMPDIR/cell_ids.txt" "$TMPDIR/clusters.txt" > "$TMPDIR/pairs.tsv"

N_PAIRS=$(wc -l < "$TMPDIR/pairs.tsv")
echo "      $N_PAIRS pairs extracted"
echo "      First 3 pairs (sanity check):"
head -3 "$TMPDIR/pairs.tsv" | sed 's/^/        /'

# ----- Step 4: normalise barcodes to CellRanger format -----------------------
echo
echo "[4/5] Normalising barcodes to CellRanger format..."

# Cluster file cell IDs:   <SAMPLE_PREFIX>-<barcode>
# CellRanger barcodes:     <barcode>-1
# Transformation:          strip prefix, add -1 suffix
awk -F'\t' -v sp="$SAMPLE_PREFIX" '
{
    bc = $1
    sub("^"sp"-", "", bc)
    bc = bc "-1"
    print bc "\t" "cluster_" $2
}
' "$TMPDIR/pairs.tsv" > "$TMPDIR/normalised.tsv"

echo "      First 3 normalised pairs:"
head -3 "$TMPDIR/normalised.tsv" | sed 's/^/        /'

# ----- Step 5: intersect with CellRanger barcodes and filter small clusters --
echo
echo "[5/5] Intersecting with CellRanger barcodes and filtering small clusters..."

awk -F'\t' 'NR==FNR {keep[$1]=1; next} ($1 in keep)' \
    "$TMPDIR/cellranger_barcodes.txt" \
    "$TMPDIR/normalised.tsv" > "$TMPDIR/intersected.tsv"

N_INTERSECT=$(wc -l < "$TMPDIR/intersected.tsv")
echo "      $N_INTERSECT barcodes after intersection"

echo
echo "      Cluster size distribution before filtering:"
cut -f2 "$TMPDIR/intersected.tsv" | sort | uniq -c | sort -k1,1nr | \
    awk -v mc=$MIN_CELLS '{
        flag = ($1 >= mc) ? "KEEP" : "DROP"
        printf "        %6d  %s  %s\n", $1, $2, flag
    }'

KEEP_CLUSTERS=$(cut -f2 "$TMPDIR/intersected.tsv" | sort | uniq -c | \
                awk -v mc=$MIN_CELLS '$1 >= mc {print $2}' | paste -sd '|')

ANNOTATION="$OUT_DIR/annotation.tsv"
awk -F'\t' -v keep="$KEEP_CLUSTERS" '
BEGIN { split(keep, k, "|"); for (i in k) ok[k[i]]=1 }
$2 in ok
' "$TMPDIR/intersected.tsv" > "$ANNOTATION"

N_FINAL=$(wc -l < "$ANNOTATION")
N_CLUSTERS=$(cut -f2 "$ANNOTATION" | sort -u | wc -l)

echo
echo "============================================================"
echo "Done!"
echo "  Annotation file: $ANNOTATION"
echo "  Cells retained:  $N_FINAL"
echo "  Clusters kept:   $N_CLUSTERS"
echo "============================================================"
echo
echo "Cluster distribution:"
cut -f2 "$ANNOTATION" | sort | uniq -c | sort -k1,1nr | sed 's/^/  /'
echo
echo "First 5 rows of annotation file:"
head -5 "$ANNOTATION" | sed 's/^/  /'
echo
echo "Next step: pass $ANNOTATION to scrna_preprocessing.sh via qsub -v ANNOTATION=..."
