#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=flye_pooled_all
#SBATCH --partition=aoraki
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/path/to/project/logs/flye_pooled_all_%j.out
#SBATCH --error=/path/to/project/logs/flye_pooled_all_%j.err

set -uo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

if ! command -v flye &> /dev/null; then
  echo "ERROR: flye not found in custom_env." >&2
  exit 1
fi

FASTQ_DIR=/path/to/project/flye_assembly/fastq
ASSEMBLY_DIR=/path/to/project/flye_assembly/assemblies

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

for SAMPLE in "${SAMPLES[@]}"; do
  echo "################ ${SAMPLE} ################"

  HP1_FASTQ="${FASTQ_DIR}/${SAMPLE}.HP1.fastq"
  HP2_FASTQ="${FASTQ_DIR}/${SAMPLE}.HP2.fastq"
  UNASSIGNED_FASTQ="${FASTQ_DIR}/${SAMPLE}.unassigned.fastq"
  POOLED_FASTQ="${FASTQ_DIR}/${SAMPLE}.pooled_all.fastq"
  POOLED_DIR="${ASSEMBLY_DIR}/${SAMPLE}_pooled_all"

  MISSING=0
  for f in "$HP1_FASTQ" "$HP2_FASTQ" "$UNASSIGNED_FASTQ"; do
    if [ ! -f "$f" ]; then
      echo "[${SAMPLE}] MISSING ${f}" >&2
      MISSING=1
    fi
  done
  if [ "$MISSING" -eq 1 ]; then
    echo "[${SAMPLE}] Skipping due to missing input(s)"
    continue
  fi

  if [ ! -f "$POOLED_FASTQ" ]; then
    echo "[${SAMPLE}] Pooling HP1 + HP2 + unassigned..."
    cat "$HP1_FASTQ" "$HP2_FASTQ" "$UNASSIGNED_FASTQ" > "$POOLED_FASTQ"
  fi

  N_POOLED=$(( $(wc -l < "$POOLED_FASTQ") / 4 ))
  echo "[${SAMPLE}] Pooled FASTQ contains ${N_POOLED} reads"

  if [ -f "${POOLED_DIR}/assembly.fasta" ]; then
    echo "[${SAMPLE}] Pooled assembly already exists, skipping Flye"
  else
    echo "[${SAMPLE}] Running Flye on pooled reads..."
    flye --nano-hq "$POOLED_FASTQ" --out-dir "$POOLED_DIR" --threads "$SLURM_CPUS_PER_TASK"
  fi

  echo "[${SAMPLE}] Done."
  echo ""
done

echo "ALL SAMPLES PROCESSED."
echo "NOTE: pooled assemblies mix both haplotypes together (no phasing)"
echo "-- treat as a diagnostic only, not a replacement for the"
echo "haplotype-resolved assemblies."
