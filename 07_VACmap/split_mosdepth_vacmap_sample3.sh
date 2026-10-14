#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=vacmap_split_mosdepth_sample3
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/vacmap_split_mosdepth_sample3_%j.out
#SBATCH --error=/path/to/project/logs/vacmap_split_mosdepth_sample3_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

HAPLOTAGGED_DIR=/path/to/project/vacmap_haplotagged
SPLIT_DIR=/path/to/project/vacmap_split_bams
MOSDEPTH_DIR=/path/to/project/vacmap_mosdepth
mkdir -p "$SPLIT_DIR" "$MOSDEPTH_DIR"

BEDS_DIR=/path/to/project/t2t_qc/beds
LPA_BED="${BEDS_DIR}/LPA_expanded.bed"
KIV2_BED="${BEDS_DIR}/KIV2_unpadded.bed"

REGIONS_BED="${MOSDEPTH_DIR}/regions_named.bed"
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "LPA"}' "$LPA_BED" > "$REGIONS_BED"
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "KIV2"}' "$KIV2_BED" >> "$REGIONS_BED"

TOOLS=(whatshap longphase)

SUMMARY="${MOSDEPTH_DIR}/SAMPLE_3_vacmap_mosdepth_summary.txt"
> "$SUMMARY"

for TOOL in "${TOOLS[@]}"; do
  IN_BAM="${HAPLOTAGGED_DIR}/SAMPLE_3.vacmap.${TOOL}.bam"
  HP1_BAM="${SPLIT_DIR}/SAMPLE_3.vacmap.${TOOL}.HP1.bam"
  HP2_BAM="${SPLIT_DIR}/SAMPLE_3.vacmap.${TOOL}.HP2.bam"

  if [ ! -f "$IN_BAM" ]; then
    echo "MISSING: $IN_BAM" >> "$SUMMARY"
    continue
  fi

  if [ -f "$HP1_BAM" ] && [ -f "$HP2_BAM" ]; then
    echo "Skipping split for ${TOOL}, already done"
  else
    echo "Splitting ${TOOL}..."
    samtools view -b -d HP:1 "$IN_BAM" > "$HP1_BAM"
    samtools view -b -d HP:2 "$IN_BAM" > "$HP2_BAM"
    samtools index "$HP1_BAM"
    samtools index "$HP2_BAM"
  fi

  N1=$(samtools view -c "$HP1_BAM")
  N2=$(samtools view -c "$HP2_BAM")
  echo "${TOOL}: HP1=${N1} reads, HP2=${N2} reads" >> "$SUMMARY"

  for HAP in HP1 HP2; do
    BAM="${SPLIT_DIR}/SAMPLE_3.vacmap.${TOOL}.${HAP}.bam"
    PREFIX="${MOSDEPTH_DIR}/SAMPLE_3.${TOOL}.${HAP}"

    if [ -f "${PREFIX}.regions.bed.gz" ]; then
      echo "Skipping mosdepth for ${TOOL} ${HAP}, already done"
    else
      mosdepth --no-per-base --by "$REGIONS_BED" "$PREFIX" "$BAM"
    fi

    echo "=== ${TOOL} / ${HAP} ===" >> "$SUMMARY"
    zcat "${PREFIX}.regions.bed.gz" >> "$SUMMARY"
  done
  echo "" >> "$SUMMARY"
done

echo "Done. Summary at $SUMMARY"
