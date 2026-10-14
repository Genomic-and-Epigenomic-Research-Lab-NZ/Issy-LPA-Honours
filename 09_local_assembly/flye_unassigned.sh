#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=flye_unassigned_variants
#SBATCH --partition=aoraki
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/path/to/project/logs/flye_unassigned_variants_%j.out
#SBATCH --error=/path/to/project/logs/flye_unassigned_variants_%j.err

set -uo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

FULL_BAM_DIR=/path/to/project/PROJECT_T2T
READNAMES_DIR=/path/to/project/flye_assembly/read_names
FASTQ_DIR=/path/to/project/flye_assembly/fastq
ASSEMBLY_DIR=/path/to/project/flye_assembly/assemblies

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

for SAMPLE in "${SAMPLES[@]}"; do
  echo "################ ${SAMPLE} ################"

  FULL_BAM="${FULL_BAM_DIR}/${SAMPLE}/filtered/${SAMPLE}_T2T_q10.bam"
  UNASSIGNED_NAMES="${READNAMES_DIR}/${SAMPLE}.unassigned_read_names.txt"

  if [ ! -f "$FULL_BAM" ] || [ ! -f "$UNASSIGNED_NAMES" ]; then
    echo "[${SAMPLE}] MISSING required input files, skipping" >&2
    continue
  fi

  UNASSIGNED_FASTQ="${FASTQ_DIR}/${SAMPLE}.unassigned.fastq"
  if [ ! -f "$UNASSIGNED_FASTQ" ]; then
    echo "[${SAMPLE}] Recovering unassigned FASTQ..."
    samtools view -N "$UNASSIGNED_NAMES" -b "$FULL_BAM" | samtools fastq - > "$UNASSIGNED_FASTQ" 2>/dev/null
  fi
  N_UNASSIGNED=$(( $(wc -l < "$UNASSIGNED_FASTQ") / 4 ))
  echo "[${SAMPLE}] Unassigned FASTQ contains ${N_UNASSIGNED} reads"

  HP1_FASTQ="${FASTQ_DIR}/${SAMPLE}.HP1.fastq"
  HP2_FASTQ="${FASTQ_DIR}/${SAMPLE}.HP2.fastq"
  HP1_COMBINED="${FASTQ_DIR}/${SAMPLE}.HP1_plus_unassigned.fastq"
  HP2_COMBINED="${FASTQ_DIR}/${SAMPLE}.HP2_plus_unassigned.fastq"
  HP1_DEDUP="${FASTQ_DIR}/${SAMPLE}.HP1_plus_unassigned.dedup.fastq"
  HP2_DEDUP="${FASTQ_DIR}/${SAMPLE}.HP2_plus_unassigned.dedup.fastq"

  if [ ! -f "$HP1_DEDUP" ]; then
    echo "[${SAMPLE}] Building and deduping HP1+unassigned..."
    cat "$HP1_FASTQ" "$UNASSIGNED_FASTQ" > "$HP1_COMBINED"
    seqkit rmdup -n "$HP1_COMBINED" -o "$HP1_DEDUP" 2>>"${FASTQ_DIR}/${SAMPLE}.HP1_rmdup.log"
    N_BEFORE=$(( $(wc -l < "$HP1_COMBINED") / 4 ))
    N_AFTER=$(( $(wc -l < "$HP1_DEDUP") / 4 ))
    echo "[${SAMPLE}] HP1+unassigned: ${N_BEFORE} reads before dedup, ${N_AFTER} after"
  fi

  if [ ! -f "$HP2_DEDUP" ]; then
    echo "[${SAMPLE}] Building and deduping HP2+unassigned..."
    cat "$HP2_FASTQ" "$UNASSIGNED_FASTQ" > "$HP2_COMBINED"
    seqkit rmdup -n "$HP2_COMBINED" -o "$HP2_DEDUP" 2>>"${FASTQ_DIR}/${SAMPLE}.HP2_rmdup.log"
    N_BEFORE=$(( $(wc -l < "$HP2_COMBINED") / 4 ))
    N_AFTER=$(( $(wc -l < "$HP2_DEDUP") / 4 ))
    echo "[${SAMPLE}] HP2+unassigned: ${N_BEFORE} reads before dedup, ${N_AFTER} after"
  fi

  UNASSIGNED_DIR="${ASSEMBLY_DIR}/${SAMPLE}_unassigned_only"
  if [ ! -f "${UNASSIGNED_DIR}/assembly.fasta" ]; then
    echo "[${SAMPLE}] Running Flye for unassigned-only..."
    flye --nano-hq "$UNASSIGNED_FASTQ" --out-dir "$UNASSIGNED_DIR" --threads "$SLURM_CPUS_PER_TASK"
  fi

  HP1_COMBINED_DIR="${ASSEMBLY_DIR}/${SAMPLE}_HP1_plus_unassigned"
  if [ ! -f "${HP1_COMBINED_DIR}/assembly.fasta" ]; then
    echo "[${SAMPLE}] Running Flye for HP1+unassigned (deduped)..."
    flye --nano-hq "$HP1_DEDUP" --out-dir "$HP1_COMBINED_DIR" --threads "$SLURM_CPUS_PER_TASK"
  fi

  HP2_COMBINED_DIR="${ASSEMBLY_DIR}/${SAMPLE}_HP2_plus_unassigned"
  if [ ! -f "${HP2_COMBINED_DIR}/assembly.fasta" ]; then
    echo "[${SAMPLE}] Running Flye for HP2+unassigned (deduped)..."
    flye --nano-hq "$HP2_DEDUP" --out-dir "$HP2_COMBINED_DIR" --threads "$SLURM_CPUS_PER_TASK"
  fi

  echo "[${SAMPLE}] Done."
  echo ""
done

echo "ALL SAMPLES PROCESSED."
