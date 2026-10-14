#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=longshot_vacmap_sample3
#SBATCH --partition=aoraki
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/longshot_vacmap_sample3_%j.out
#SBATCH --error=/path/to/project/logs/longshot_vacmap_sample3_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

BAM=/path/to/project/vacmap_out/SAMPLE_3_piped.sorted.bam
REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
OUT_DIR=/path/to/project/vacmap_longshot_phased
mkdir -p "$OUT_DIR"

OUT_VCF="${OUT_DIR}/SAMPLE_3.vacmap.LPA.phased.vcf"

if [ -f "$OUT_VCF" ]; then
  echo "Already done, skipping: $OUT_VCF"
  exit 0
fi

set +o pipefail
samtools quickcheck "$BAM" 2>/dev/null
QC_STATUS=$?
set -o pipefail
if [ "$QC_STATUS" -ne 0 ]; then
  echo "quickcheck failed for ${BAM}" >&2
  exit 1
fi

echo "Phasing SAMPLE_3 (VACmap alignment)..."
longshot \
  --bam "$BAM" \
  --ref "$REF" \
  --out "$OUT_VCF"

echo "Done."
