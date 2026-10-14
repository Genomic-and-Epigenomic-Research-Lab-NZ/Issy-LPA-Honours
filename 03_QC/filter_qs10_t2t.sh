#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=filter_qs10
#SBATCH --account=ACCOUNT
#SBATCH --partition=aoraki
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/path/to/project/logs/filter_qs10_%j.out
#SBATCH --error=/path/to/project/logs/filter_qs10_%j.err

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load miniconda3

eval "$(conda shell.bash hook)"

set +u
conda activate custom_env
set -u

echo "Threads available: $OMP_NUM_THREADS"

USER_TAG="PROJECT"
REF="T2T"
TAG_REF="${USER_TAG}_${REF}"

SAMPLES=(SAMPLE_1 SAMPLE_3 SAMPLE_2 SAMPLE_4)

for SAMPLE in "${SAMPLES[@]}"; do

  echo "Processing $SAMPLE"

  MERGE_DIR="/path/to/project/${TAG_REF}/${SAMPLE}/merged"
  MERGED_BAM="${MERGE_DIR}/${SAMPLE}_${REF}_merged.bam"

  FILTERED_DIR="/path/to/project/${TAG_REF}/${SAMPLE}/filtered"
  Q10_BAM="${FILTERED_DIR}/${SAMPLE}_${REF}_q10.bam"

  if [ ! -f "$MERGED_BAM" ]; then
    echo "SKIPPING $SAMPLE - merged bam not found at $MERGED_BAM"
    echo "(merge this sample first, then rerun this script)"
    continue
  fi

  if [ -f "$Q10_BAM" ]; then
    echo "SKIPPING $SAMPLE - q10 bam already exists at $Q10_BAM"
    continue
  fi

  mkdir -p "$FILTERED_DIR"

  echo "Filtering $SAMPLE merged BAM to qs >= 10..."
  samtools view -@ $OMP_NUM_THREADS -e '[qs]>=10' -b -o "$Q10_BAM" "$MERGED_BAM"
  samtools index -@ $OMP_NUM_THREADS "$Q10_BAM"

   echo "Checking MM/ML tags survived filtering for $SAMPLE..."
  set +o pipefail
  MM_CHECK=$(samtools view "$Q10_BAM" | head -5 | grep -o 'MM:Z:[^[:space:]]*' | wc -l)
  ML_CHECK=$(samtools view "$Q10_BAM" | head -5 | grep -o 'ML:B:C[^[:space:]]*' | wc -l)
  set -o pipefail

  if [ "$MM_CHECK" -eq 0 ] || [ "$ML_CHECK" -eq 0 ]; then
    echo "ERROR: MM or ML tags missing after qs10 filter for $SAMPLE"
    echo "STOPPING this sample - do not proceed to QC or modkit with this bam"
    continue
  else
    echo "MM/ML tags present for $SAMPLE, q10 bam ready"
  fi

done

echo "Filtering run complete for all samples (see above for any skipped/failed)."
