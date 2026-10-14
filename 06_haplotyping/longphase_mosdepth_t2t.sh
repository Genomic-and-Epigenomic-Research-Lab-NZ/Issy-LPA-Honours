#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=longphase_mosdepth
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/longphase_mosdepth_%j.out
#SBATCH --error=/path/to/project/logs/longphase_mosdepth_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

SPLIT_DIR=/path/to/project/longphase_split_bams
OUT_DIR=/path/to/project/longphase_mosdepth
mkdir -p "$OUT_DIR"

BEDS_DIR=/path/to/project/t2t_qc/beds
LPA_BED="${BEDS_DIR}/LPA_expanded.bed"
KIV2_BED="${BEDS_DIR}/KIV2_unpadded.bed"

REGIONS_BED="${OUT_DIR}/haplotype_regions_named.bed"
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "LPA"}' "$LPA_BED" > "$REGIONS_BED"
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "KIV2"}' "$KIV2_BED" >> "$REGIONS_BED"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)
HAPLOTYPES=(HP1 HP2)

SUMMARY=/path/to/project/longphase_mosdepth_summary.txt
> "$SUMMARY"

for SAMPLE in "${SAMPLES[@]}"; do
  for HAP in "${HAPLOTYPES[@]}"; do
    BAM="${SPLIT_DIR}/${SAMPLE}.LPA.longphase.${HAP}.bam"
    PREFIX="${OUT_DIR}/${SAMPLE}.${HAP}.longphase"

    if [ ! -f "$BAM" ]; then
      echo "MISSING: $BAM" >> "$SUMMARY"
      continue
    fi

    if [ -f "${PREFIX}.regions.bed.gz" ]; then
      echo "Skipping ${SAMPLE} ${HAP}, already done"
    else
      echo "Running mosdepth for ${SAMPLE} ${HAP}..."
      mosdepth --no-per-base --by "$REGIONS_BED" "$PREFIX" "$BAM"
    fi

    echo "=== ${SAMPLE} / ${HAP} ===" >> "$SUMMARY"
    zcat "${PREFIX}.regions.bed.gz" >> "$SUMMARY"
    echo "" >> "$SUMMARY"
  done
done

echo "Done. Summary at $SUMMARY"
