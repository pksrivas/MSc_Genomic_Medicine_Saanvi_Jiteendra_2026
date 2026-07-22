#!/usr/bin/env bash
# =============================================================================
# bulk_bam_to_fastq.sh  —  reditR pipeline template
# =============================================================================
# Purpose:
#   Convert a bulk RNA-seq BAM file to FASTQ for input to SPRINT.
#   Only needed if you start from a BAM file. If you already have FASTQs,
#   skip this and run SPRINT directly.
#
# Usage (PBS/HPC — submit once per sample):
#   qsub -v BAM=/path/to/sample.bam,FASTQ=/path/to/output.fastq \
#        bulk_bam_to_fastq.sh
#
# Dependencies:
#   - samtools  (tested with 1.19)
#   - bedtools  (for bamToFastq; tested with 2.31.1)
#
# After this completes:
#   Run SPRINT on the output FASTQ, then pipe results through
#   extracting_read_counts.sh to produce the count table for reditR.
# =============================================================================

#PBS -N bulk_bam_to_fastq
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -l walltime=4:00:00
#PBS -j oe

set -euo pipefail

# ─── LOAD MODULES ─────────────────────────────────────────────────────────────
# Adjust module names for your HPC environment.
# module load SAMtools/<VERSION>
# module load BEDTools/<VERSION>
# ──────────────────────────────────────────────────────────────────────────────

# ─── CONFIGURE ────────────────────────────────────────────────────────────────
SAM=samtools    # full path if samtools is not on $PATH
CPU=4           # should match PBS ncpus above
# ──────────────────────────────────────────────────────────────────────────────

# BAM and FASTQ are passed at submission time via qsub -v
: "${BAM:?BAM variable not set — pass via qsub -v BAM=...}"
: "${FASTQ:?FASTQ variable not set — pass via qsub -v FASTQ=...}"

SORTED_BAM="${BAM%.bam}_namesorted.bam"

echo "[$(date)] Sorting BAM by read name..."
${SAM} sort -n -@ "${CPU}" "${BAM}" -o "${SORTED_BAM}"

echo "[$(date)] Converting to FASTQ..."
bamToFastq -i "${SORTED_BAM}" -fq "${FASTQ}"

rm "${SORTED_BAM}"

echo "[$(date)] Done. FASTQ written to ${FASTQ}"
echo "Next step: run SPRINT on ${FASTQ}"
