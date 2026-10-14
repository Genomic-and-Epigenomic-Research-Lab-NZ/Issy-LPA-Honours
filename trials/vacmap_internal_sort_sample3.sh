#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=vacmap_sample3
#SBATCH --partition=aoraki
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/path/to/project/logs/vacmap_sample3_%j.out
#SBATCH --error=/path/to/project/logs/vacmap_sample3_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate vacmap_env
set -u

LPA_BAM=/path/to/project/mini_bams/SAMPLE_3.LPA.bam
REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna

OUT_DIR=/path/to/project/vacmap_out
mkdir -p "$OUT_DIR"

FASTQ="${OUT_DIR}/SAMPLE_3.LPA.fastq"
VACMAP_BAM="${OUT_DIR}/SAMPLE_3_vacmap.sorted.bam"

if [ -f "$FASTQ" ]; then
  echo "FASTQ already exists, skipping extraction: $FASTQ"
else
  echo "Extracting reads to FASTQ..."
  set +o pipefail
  samtools quickcheck "$LPA_BAM" 2>/dev/null
  QC_STATUS=$?
  set -o pipefail
  if [ "$QC_STATUS" -ne 0 ]; then
    echo "quickcheck failed for ${LPA_BAM}" >&2
    exit 1
  fi

  samtools fastq -F 0x900 "$LPA_BAM" > "$FASTQ"
fi

READ_COUNT=$(( $(wc -l < "$FASTQ") / 4 ))
echo "FASTQ contains ${READ_COUNT} reads."

if [ -f "$VACMAP_BAM" ]; then
  echo "VACmap output already exists, skipping: $VACMAP_BAM"
else
  echo "Running VACmap (mode H, ONT)..."
  vacmap \
    -ref "$REF" \
    -read "$FASTQ" \
    -mode H \
    -o "$VACMAP_BAM" \
    -t "$SLURM_CPUS_PER_TASK" \
    --MD \
    --force
fi

samtools index "$VACMAP_BAM"

echo "Done. Output: $VACMAP_BAM"
