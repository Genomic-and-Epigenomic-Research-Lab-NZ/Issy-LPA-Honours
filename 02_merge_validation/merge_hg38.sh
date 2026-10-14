#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=merge_bams
#SBATCH --account=ACCOUNT
#SBATCH --partition=aoraki
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/path/to/project/logs/merge_hg38_%j.out
#SBATCH --error=/path/to/project/logs/merge_hg38_%j.err

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load miniconda3

eval "$(conda shell.bash hook)"

set +u
conda activate custom_env
set -u

echo "Threads available: $OMP_NUM_THREADS"

SAMPLE="SAMPLE_3"
USER_TAG="PROJECT"
REF="hg38"
TAG_REF="${USER_TAG}_${REF}"

OUTPUT_DIR="/path/to/scratch/hg38_bams/"
MERGE_DIR="/path/to/project/${TAG_REF}/${SAMPLE}/merged"
T2T_LIST="/path/to/project/PROJECT_T2T/${SAMPLE}/merged/${SAMPLE}_list.txt"

mkdir -p $MERGE_DIR
mkdir -p /path/to/project/${TAG_REF}/${SAMPLE}

find $OUTPUT_DIR -path "*${SAMPLE}*/bam_pass/*.bam" > ${MERGE_DIR}/${SAMPLE}_list.txt

COUNT=$(wc -l < ${MERGE_DIR}/${SAMPLE}_list.txt)
EXPECTED_COUNT=$(wc -l < "$T2T_LIST")
echo "$COUNT files to merge for $SAMPLE (expected $EXPECTED_COUNT, from T2T shard count)"

if [ "$COUNT" -eq 0 ]; then
  echo "ERROR: no passed bams found"
  exit 1
fi

if [ "$COUNT" -gt "$EXPECTED_COUNT" ]; then
  echo "ERROR: found $COUNT files, more than expected $EXPECTED_COUNT files"
  exit 1
fi

echo "Running samtools quickcheck on all shards..."

QUICKCHECK_LOG="${MERGE_DIR}/${SAMPLE}_quickcheck.log"

if samtools quickcheck -v $(cat ${MERGE_DIR}/${SAMPLE}_list.txt) > "$QUICKCHECK_LOG" 2>/dev/null; then
  echo "All $COUNT shards passed quickcheck"
else
  echo "ERROR: one or more BAM shards failed integrity check:"
  cat "$QUICKCHECK_LOG"
  exit 1
fi

samtools merge -@ $OMP_NUM_THREADS -c -p -f -b ${MERGE_DIR}/${SAMPLE}_list.txt ${MERGE_DIR}/${SAMPLE}_${REF}_merged.bam
samtools index -@ $OMP_NUM_THREADS ${MERGE_DIR}/${SAMPLE}_${REF}_merged.bam

echo "$SAMPLE merge complete: ${SAMPLE}_${REF}_merged.bam"
