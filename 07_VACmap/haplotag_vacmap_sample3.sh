#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=haplotag_vacmap_sample3
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/haplotag_vacmap_sample3_%j.out
#SBATCH --error=/path/to/project/logs/haplotag_vacmap_sample3_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

BAM=/path/to/project/vacmap_out/SAMPLE_3_piped.sorted.bam
REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
VCF_GZ=/path/to/project/vacmap_longshot_phased/SAMPLE_3.vacmap.LPA.phased.vcf.gz

OUT_DIR=/path/to/project/vacmap_haplotagged
mkdir -p "$OUT_DIR"

WHATSHAP_BAM="${OUT_DIR}/SAMPLE_3.vacmap.whatshap.bam"
LONGPHASE_PREFIX="${OUT_DIR}/SAMPLE_3.vacmap.longphase"
LONGPHASE_BAM="${LONGPHASE_PREFIX}.bam"

if [ ! -f "$VCF_GZ" ] || [ ! -f "${VCF_GZ}.tbi" ]; then
  echo "MISSING bgzip/tabix index for phased VCF: ${VCF_GZ} (.tbi)" >&2
  exit 1
fi

set +o pipefail
samtools quickcheck "$BAM" 2>/dev/null
QC_STATUS=$?
set -o pipefail
if [ "$QC_STATUS" -ne 0 ]; then
  echo "quickcheck failed for ${BAM}" >&2
  exit 1
fi

if [ -f "$WHATSHAP_BAM" ]; then
  echo "whatshap output already exists, skipping"
else
  echo "Haplotagging with whatshap..."
  whatshap haplotag \
    --reference "$REF" \
    --ignore-read-groups \
    -o "$WHATSHAP_BAM" \
    "$VCF_GZ" \
    "$BAM"
  samtools index "$WHATSHAP_BAM"
fi

if [ -f "$LONGPHASE_BAM" ]; then
  echo "longphase output already exists, skipping"
else
  echo "Haplotagging with longphase..."
  longphase haplotag \
    -s "$VCF_GZ" \
    -b "$BAM" \
    -r "$REF" \
    -o "$LONGPHASE_PREFIX" \
    -q 0 \
    --tagSupplementary \
    -t "$SLURM_CPUS_PER_TASK"
  samtools index "$LONGPHASE_BAM"
fi

echo "Done."
