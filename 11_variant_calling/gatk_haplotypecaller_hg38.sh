#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=gatk_hg38_fixed
#SBATCH --partition=aoraki
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/gatk_hg38_fixed_%j.out
#SBATCH --error=/path/to/project/logs/gatk_hg38_fixed_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

REF=/path/to/refs/hg38.primary_assembly.fa
BAM_BASE=/path/to/project/PROJECT_hg38

OUT_DIR=/path/to/project/gatk_hg38
FIXED_BAM_DIR=/path/to/project/gatk_hg38/read_group_fixed_bams
mkdir -p "$OUT_DIR" "$FIXED_BAM_DIR"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

CHR="chr6"
LPA_START=160300000
LPA_END=160800000
REGION="${CHR}:${LPA_START}-${LPA_END}"

WIDE_START=160050000
WIDE_END=161050000
WIDE_REGION="${CHR}:${WIDE_START}-${WIDE_END}"

for SAMPLE in "${SAMPLES[@]}"; do
  IN_BAM="${BAM_BASE}/${SAMPLE}/filtered/${SAMPLE}_hg38_q10.bam"
  FIXED_BAM="${FIXED_BAM_DIR}/${SAMPLE}_hg38_q10_rg.bam"
  OUT_VCF="${OUT_DIR}/${SAMPLE}_hg38_LPA.vcf.gz"
  OUT_VCF_WIDE="${OUT_DIR}/${SAMPLE}_hg38_LPA_wide.vcf.gz"

  if [ ! -f "$IN_BAM" ]; then
    echo "[${SAMPLE}] MISSING BAM: ${IN_BAM}, skipping" >&2
    continue
  fi

  if [ ! -f "$FIXED_BAM" ]; then
    echo "[${SAMPLE}] Adding read group (SM=${SAMPLE})..."
    gatk AddOrReplaceReadGroups \
      -I "$IN_BAM" \
      -O "$FIXED_BAM" \
      -RGID "${SAMPLE}" \
      -RGLB "lib1" \
      -RGPL "ONT" \
      -RGPU "unit1" \
      -RGSM "${SAMPLE}"
    samtools index "$FIXED_BAM"
  fi

  set +o pipefail
  samtools quickcheck "$FIXED_BAM" 2>/dev/null
  QC_STATUS=$?
  set -o pipefail
  if [ "$QC_STATUS" -ne 0 ]; then
    echo "[${SAMPLE}] quickcheck failed on fixed BAM, skipping" >&2
    continue
  fi

  if [ ! -f "$OUT_VCF" ]; then
    echo "[${SAMPLE}] Running GATK HaplotypeCaller (standard region)..."
    gatk HaplotypeCaller \
      -R "$REF" \
      -I "$FIXED_BAM" \
      -L "$REGION" \
      -O "$OUT_VCF"
  else
    echo "[${SAMPLE}] Standard region already done, skipping"
  fi

  if [ ! -f "$OUT_VCF_WIDE" ]; then
    echo "[${SAMPLE}] Running GATK HaplotypeCaller (wide region)..."
    gatk HaplotypeCaller \
      -R "$REF" \
      -I "$FIXED_BAM" \
      -L "$WIDE_REGION" \
      -O "$OUT_VCF_WIDE"
  else
    echo "[${SAMPLE}] Wide region already done, skipping"
  fi
done

echo "Done."
