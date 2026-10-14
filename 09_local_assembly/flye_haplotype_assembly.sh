#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=flye_haplotype_assembly
#SBATCH --partition=aoraki
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/path/to/project/logs/flye_haplotype_assembly_%j.out
#SBATCH --error=/path/to/project/logs/flye_haplotype_assembly_%j.err

set -uo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

LONGPHASE_DIR=/path/to/project/longphase_haplotagged

FULL_BAM_DIR=/path/to/project/PROJECT_T2T

OUT_DIR=/path/to/project/flye_assembly
READNAMES_DIR="${OUT_DIR}/read_names"
FASTQ_DIR="${OUT_DIR}/fastq"
ASSEMBLY_DIR="${OUT_DIR}/assemblies"
mkdir -p "$READNAMES_DIR" "$FASTQ_DIR" "$ASSEMBLY_DIR"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

for SAMPLE in "${SAMPLES[@]}"; do
  echo "################ ${SAMPLE} ################"

  LONGPHASE_BAM="${LONGPHASE_DIR}/${SAMPLE}.LPA.longphase.bam"
  FULL_BAM="${FULL_BAM_DIR}/${SAMPLE}/filtered/${SAMPLE}_T2T_q10.bam"

  if [ ! -f "$LONGPHASE_BAM" ]; then
    echo "[${SAMPLE}] MISSING: ${LONGPHASE_BAM}, skipping" >&2
    continue
  fi
  if [ ! -f "$FULL_BAM" ]; then
    echo "[${SAMPLE}] MISSING: ${FULL_BAM}, skipping" >&2
    continue
  fi

  set +o pipefail
  samtools quickcheck "$LONGPHASE_BAM" 2>/dev/null; QC1=$?
  samtools quickcheck "$FULL_BAM" 2>/dev/null; QC2=$?
  set -o pipefail
  if [ "$QC1" -ne 0 ] || [ "$QC2" -ne 0 ]; then
    echo "[${SAMPLE}] quickcheck failed on one or both BAMs, skipping" >&2
    continue
  fi

  HP1_NAMES="${READNAMES_DIR}/${SAMPLE}.HP1_read_names.txt"
  HP2_NAMES="${READNAMES_DIR}/${SAMPLE}.HP2_read_names.txt"
  UNASSIGNED_NAMES="${READNAMES_DIR}/${SAMPLE}.unassigned_read_names.txt"

  if [ ! -f "$HP1_NAMES" ]; then
    echo "[${SAMPLE}] Extracting HP1 read names..."
    samtools view -d HP:1 "$LONGPHASE_BAM" | cut -f1 | sort -u > "$HP1_NAMES"
  fi
  if [ ! -f "$HP2_NAMES" ]; then
    echo "[${SAMPLE}] Extracting HP2 read names..."
    samtools view -d HP:2 "$LONGPHASE_BAM" | cut -f1 | sort -u > "$HP2_NAMES"
  fi
  if [ ! -f "$UNASSIGNED_NAMES" ]; then
    echo "[${SAMPLE}] Extracting unassigned read names..."
    samtools view "$LONGPHASE_BAM" | awk '!/HP:i:/' | cut -f1 | sort -u > "$UNASSIGNED_NAMES"
  fi

  N_HP1=$(wc -l < "$HP1_NAMES")
  N_HP2=$(wc -l < "$HP2_NAMES")
  N_UNASSIGNED=$(wc -l < "$UNASSIGNED_NAMES")
  echo "[${SAMPLE}] HP1=${N_HP1} reads, HP2=${N_HP2} reads, unassigned=${N_UNASSIGNED} reads"

  for HAP in HP1 HP2; do
    NAMES_FILE="${READNAMES_DIR}/${SAMPLE}.${HAP}_read_names.txt"
    FASTQ_OUT="${FASTQ_DIR}/${SAMPLE}.${HAP}.fastq"

    if [ ! -f "$FASTQ_OUT" ]; then
      echo "[${SAMPLE}] Recovering FASTQ for ${HAP} from full Q10 BAM..."
      samtools view -N "$NAMES_FILE" -b "$FULL_BAM" | samtools fastq - > "$FASTQ_OUT" 2>/dev/null
    fi
    N_FASTQ=$(( $(wc -l < "$FASTQ_OUT") / 4 ))
    echo "[${SAMPLE}] ${HAP} FASTQ contains ${N_FASTQ} reads"
  done

  for HAP in HP1 HP2; do
    FASTQ_IN="${FASTQ_DIR}/${SAMPLE}.${HAP}.fastq"
    HAP_ASSEMBLY_DIR="${ASSEMBLY_DIR}/${SAMPLE}_${HAP}"

    if [ -f "${HAP_ASSEMBLY_DIR}/assembly.fasta" ]; then
      echo "[${SAMPLE}] ${HAP} assembly already exists, skipping"
      continue
    fi

    echo "[${SAMPLE}] Running Flye for ${HAP}..."
    flye --nano-hq "$FASTQ_IN" \
      --out-dir "$HAP_ASSEMBLY_DIR" \
      --threads "$SLURM_CPUS_PER_TASK"
  done

  echo "[${SAMPLE}] Done."
  echo ""
done

echo "ALL SAMPLES PROCESSED."
