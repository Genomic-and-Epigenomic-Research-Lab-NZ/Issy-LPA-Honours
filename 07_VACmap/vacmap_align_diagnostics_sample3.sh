#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=vacmap_diagnostics
#SBATCH --partition=aoraki
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/path/to/project/logs/vacmap_diagnostics_%j.out
#SBATCH --error=/path/to/project/logs/vacmap_diagnostics_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate vacmap_env
set -u

REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
FASTQ=/path/to/project/vacmap_out/SAMPLE_3.LPA.fastq
OUT_DIR=/path/to/project/vacmap_out
mkdir -p "$OUT_DIR"

SUMMARY="${OUT_DIR}/vacmap_diagnostics_summary.txt"
> "$SUMMARY"

echo "Input FASTQ read count: $(( $(wc -l < "$FASTQ") / 4 ))" >> "$SUMMARY"
echo "" >> "$SUMMARY"

DEBUG_BAM="${OUT_DIR}/SAMPLE_3_debug.sorted.bam"
DEBUG_LOG="${OUT_DIR}/vacmap_debug_log.txt"

echo "=== TEST 1: --debug mode (Mode H) ===" >> "$SUMMARY"
if vacmap -ref "$REF" -read "$FASTQ" -mode H -o "$DEBUG_BAM" \
    -t "$SLURM_CPUS_PER_TASK" --debug --force > "$DEBUG_LOG" 2>&1; then
  echo "Run completed." >> "$SUMMARY"
  samtools flagstat "$DEBUG_BAM" >> "$SUMMARY" 2>&1 || echo "flagstat failed (file may not exist)" >> "$SUMMARY"
  echo "--- last 50 lines of debug log ---" >> "$SUMMARY"
  tail -50 "$DEBUG_LOG" >> "$SUMMARY"
else
  echo "VACmap --debug run FAILED (exit code $?)" >> "$SUMMARY"
  echo "--- last 50 lines of debug log ---" >> "$SUMMARY"
  tail -50 "$DEBUG_LOG" >> "$SUMMARY" 2>&1 || true
fi
echo "" >> "$SUMMARY"

MODES_BAM="${OUT_DIR}/SAMPLE_3_modeS.sorted.bam"

echo "=== TEST 2: Mode S (Sensitive) ===" >> "$SUMMARY"
vacmap -ref "$REF" -read "$FASTQ" -mode S -o "$MODES_BAM" \
  -t "$SLURM_CPUS_PER_TASK" --force >> "$SUMMARY" 2>&1
echo "Run completed." >> "$SUMMARY"
samtools flagstat "$MODES_BAM" >> "$SUMMARY" 2>&1 || echo "flagstat failed" >> "$SUMMARY"
echo "" >> "$SUMMARY"

PIPED_BAM="${OUT_DIR}/SAMPLE_3_piped.sorted.bam"

echo "=== TEST 3: Mode H, piped through samtools sort (per README example) ===" >> "$SUMMARY"
vacmap -ref "$REF" -read "$FASTQ" -mode H -t "$SLURM_CPUS_PER_TASK" --force \
  | samtools sort -@4 > "$PIPED_BAM"
echo "Run completed." >> "$SUMMARY"
samtools flagstat "$PIPED_BAM" >> "$SUMMARY" 2>&1 || echo "flagstat failed" >> "$SUMMARY"
echo "" >> "$SUMMARY"

echo "=== ALL TESTS DONE ===" >> "$SUMMARY"
echo "Summary at: $SUMMARY"
