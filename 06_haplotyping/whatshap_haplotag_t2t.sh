#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=whatshap_haplotag
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/whatshap_haplotag_%j.out
#SBATCH --error=/path/to/project/logs/whatshap_haplotag_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

LPA_DIR=/path/to/project/mini_bams
REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
VCF_DIR=/path/to/project/longshot_phased
OUT_DIR=/path/to/project/haplotagged_bams
mkdir -p "$OUT_DIR"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

for SAMPLE in "${SAMPLES[@]}"; do
  BAM="${LPA_DIR}/${SAMPLE}.LPA.bam"
  VCF="${VCF_DIR}/${SAMPLE}.LPA.phased.vcf"
  OUT_BAM="${OUT_DIR}/${SAMPLE}.LPA.haplotagged.bam"

  if [ -f "$OUT_BAM" ]; then
    echo "Skipping ${SAMPLE}, already done"
    continue
  fi

  if [ ! -f "$VCF" ]; then
    echo "MISSING phased VCF for ${SAMPLE}: ${VCF}" >&2
    exit 1
  fi

  VCF_GZ="${VCF}.gz"
  if [ ! -f "$VCF_GZ" ] || [ ! -f "${VCF_GZ}.tbi" ]; then
    echo "MISSING bgzip/tabix index for ${SAMPLE}: ${VCF_GZ} (.tbi)" >&2
    exit 1
  fi

  set +o pipefail
  samtools quickcheck "$BAM" 2>/dev/null
  QC_STATUS=$?
  set -o pipefail
  if [ "$QC_STATUS" -ne 0 ]; then
    echo "quickcheck failed for ${SAMPLE}: ${BAM}" >&2
    exit 1
  fi

  echo "Haplotagging ${SAMPLE}..."
  whatshap haplotag \
    --reference "$REF" \
    --ignore-read-groups \
    -o "$OUT_BAM" \
    "$VCF_GZ" \
    "$BAM"

  samtools index "$OUT_BAM"
done

echo "Done."
