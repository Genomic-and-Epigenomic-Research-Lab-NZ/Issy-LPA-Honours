#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=longshot_phase
#SBATCH --partition=aoraki
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/longshot_phase_%j.out
#SBATCH --error=/path/to/project/logs/longshot_phase_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

LPA_DIR=/path/to/project/mini_bams
REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
OUT_DIR=/path/to/project/longshot_phased
mkdir -p "$OUT_DIR"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

for SAMPLE in "${SAMPLES[@]}"; do
  BAM="${LPA_DIR}/${SAMPLE}.LPA.bam"
  OUT_VCF="${OUT_DIR}/${SAMPLE}.LPA.phased.vcf"

  if [ -f "$OUT_VCF" ]; then
    echo "Skipping ${SAMPLE}, already done"
    continue
  fi

  set +o pipefail
  samtools quickcheck "$BAM" 2>/dev/null
  QC_STATUS=$?
  set -o pipefail
  if [ "$QC_STATUS" -ne 0 ]; then
    echo "quickcheck failed for ${SAMPLE}: ${BAM}" >&2
    exit 1
  fi

  echo "Phasing ${SAMPLE}..."
  longshot \
    --bam "$BAM" \
    --ref "$REF" \
    --out "$OUT_VCF"
done

echo "Done."
