#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=longphase_haplotag
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/longphase_haplotag_%j.out
#SBATCH --error=/path/to/project/logs/longphase_haplotag_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

LPA_DIR=/path/to/project/mini_bams
REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
VCF_DIR=/path/to/project/longshot_phased
OUT_DIR=/path/to/project/longphase_haplotagged
mkdir -p "$OUT_DIR"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

for SAMPLE in "${SAMPLES[@]}"; do
  BAM="${LPA_DIR}/${SAMPLE}.LPA.bam"
  VCF_GZ="${VCF_DIR}/${SAMPLE}.LPA.phased.vcf.gz"
  OUT_PREFIX="${OUT_DIR}/${SAMPLE}.LPA.longphase"
  OUT_BAM="${OUT_PREFIX}.bam"

  if [ -f "$OUT_BAM" ]; then
    echo "Skipping ${SAMPLE}, already done"
    continue
  fi

  if [ ! -f "$VCF_GZ" ]; then
    echo "MISSING bgzip/tabix-indexed VCF for ${SAMPLE}: ${VCF_GZ}" >&2
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

  echo "Haplotagging ${SAMPLE} with longphase..."
  longphase haplotag \
    -s "$VCF_GZ" \
    -b "$BAM" \
    -r "$REF" \
    -o "$OUT_PREFIX" \
    -q 0 \
    --tagSupplementary \
    -t "$SLURM_CPUS_PER_TASK"

  samtools index "$OUT_BAM"

  echo "${SAMPLE} KIV-2 tagged reads (longphase):"
  samtools view "$OUT_BAM" NC_060930.1:161862121-161998520 | grep -o "HP:i:[0-9]*" | sort | uniq -c
done

echo "Done."
