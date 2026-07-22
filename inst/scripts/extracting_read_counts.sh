#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   Bulk RNA-seq (default):
#     bash extracting_read_counts.sh
#
#   scRNA-seq or custom samples — pass a file listing sample names (one per line):
#     bash extracting_read_counts.sh samples.txt
#
# Each sample name must match a directory named {sample}_ in OUT_DIR
# containing SPRINT_identified_regular.res

OUT_DIR="${OUT_DIR:-/rds/general/user/sj1825/home/diabetes_output}"

# Build sample list: from file argument or default bulk RNA-seq names
if [[ $# -gt 0 && -f "$1" ]]; then
    mapfile -t SAMPLES < "$1"
else
    SAMPLES=(AS001 AS002 AS003 AS004 AS005 AS006 RK001 RK002 RK003 RK004 RK005 RK006)
fi

for SAMPLE in "${SAMPLES[@]}"; do
    DIR="${SAMPLE}_"
    grep -v "^#" "${OUT_DIR}/${DIR}/SPRINT_identified_regular.res" | \
    awk -v s="${SAMPLE}" '{
        split($7,a,":");
        site=$1":"$3;
        print site"\t"s"\t"a[1]"\t"a[2]
    }' > "${OUT_DIR}/${SAMPLE}_editing.txt"
done

# -------------------------
# Merge all per-sample files into one table
# -------------------------
echo -e "site\tsample\tedited\ttotal" > "${OUT_DIR}/all_samples_editing.txt"
for SAMPLE in "${SAMPLES[@]}"; do
    cat "${OUT_DIR}/${SAMPLE}_editing.txt" >> "${OUT_DIR}/all_samples_editing.txt"
done