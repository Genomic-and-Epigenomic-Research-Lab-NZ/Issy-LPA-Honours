#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=modkit_pileup_all
#SBATCH --partition=aoraki
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/modkit_pileup_all_%j.out
#SBATCH --error=/path/to/project/logs/modkit_pileup_all_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

if ! command -v modkit &> /dev/null; then
  echo "ERROR: modkit not found in custom_env. Check installation before running." >&2
  exit 1
fi

REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
LPA_DIR=/path/to/project/mini_bams

OUT_DIR=/path/to/project/modkit_pileup
mkdir -p "$OUT_DIR"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

for SAMPLE in "${SAMPLES[@]}"; do
  BAM="${LPA_DIR}/${SAMPLE}.LPA.bam"
  OUT_BED="${OUT_DIR}/${SAMPLE}.modkit_pileup.bed"

  set +o pipefail
  samtools quickcheck "$BAM" 2>/dev/null
  QC_STATUS=$?
  set -o pipefail
  if [ "$QC_STATUS" -ne 0 ]; then
    echo "[${SAMPLE}] quickcheck failed, skipping" >&2
    continue
  fi

  MM_COUNT=$(samtools view "$BAM" | head -1000 | grep -c "MM:Z:" || true)
  if [ "$MM_COUNT" -eq 0 ]; then
    echo "[${SAMPLE}] WARNING: no MM:Z: tags found in first 1000 reads -- methylation info may be missing from this BAM. Running anyway, but check output." >&2
  fi

  if [ -f "$OUT_BED" ]; then
    echo "[${SAMPLE}] Already done, skipping"
    continue
  fi

  echo "[${SAMPLE}] Running modkit pileup..."
  modkit pileup \
    "$BAM" \
    "$OUT_BED" \
    --ref "$REF" \
    --threads "$SLURM_CPUS_PER_TASK"
done

echo "Done."
