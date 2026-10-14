#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=sv_calling_minimap2_all
#SBATCH --partition=aoraki
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/sv_calling_minimap2_all_%j.out
#SBATCH --error=/path/to/project/logs/sv_calling_minimap2_all_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
LPA_DIR=/path/to/project/mini_bams

OUT_DIR=/path/to/project/sv_calling_minimap2
mkdir -p "$OUT_DIR"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

CHR="NC_060930.1"
KIV2_START=161862121
KIV2_END=161998520

SUMMARY="${OUT_DIR}/sv_calling_summary.txt"
> "$SUMMARY"

for SAMPLE in "${SAMPLES[@]}"; do
  BAM="${LPA_DIR}/${SAMPLE}.LPA.bam"

  SNIFFLES_VCF="${OUT_DIR}/${SAMPLE}_sniffles2.vcf"
  CUTESV_VCF="${OUT_DIR}/${SAMPLE}_cutesv.vcf"
  CUTESV_WORKDIR="${OUT_DIR}/${SAMPLE}_cutesv_work"
  mkdir -p "$CUTESV_WORKDIR"

  set +o pipefail
  samtools quickcheck "$BAM" 2>/dev/null
  QC_STATUS=$?
  set -o pipefail
  if [ "$QC_STATUS" -ne 0 ]; then
    echo "[${SAMPLE}] quickcheck failed, skipping" >> "$SUMMARY"
    continue
  fi

  echo "[${SAMPLE}] Running Sniffles2..."
  if [ ! -f "$SNIFFLES_VCF" ]; then
    sniffles \
      --input "$BAM" \
      --vcf "$SNIFFLES_VCF" \
      --reference "$REF" \
      --threads "$SLURM_CPUS_PER_TASK"
  fi

  echo "[${SAMPLE}] Running cuteSV..."
  if [ ! -f "$CUTESV_VCF" ]; then
    cuteSV \
      "$BAM" \
      "$REF" \
      "$CUTESV_VCF" \
      "$CUTESV_WORKDIR" \
      --threads "$SLURM_CPUS_PER_TASK" \
      --max_cluster_bias_INS 100 \
      --diff_ratio_merging_INS 0.3 \
      --max_cluster_bias_DEL 100 \
      --diff_ratio_merging_DEL 0.3
  fi

  echo "=== ${SAMPLE} ===" >> "$SUMMARY"
  echo "--- Sniffles2 calls overlapping KIV-2 ---" >> "$SUMMARY"
  grep -v "^#" "$SNIFFLES_VCF" | awk -v chr="$CHR" -v s="$KIV2_START" -v e="$KIV2_END" \
    '$1==chr && $2>=s-50000 && $2<=e+50000 {print $1, $2, $3, $5, $8}' >> "$SUMMARY" || echo "(none found)" >> "$SUMMARY"

  echo "--- cuteSV calls overlapping KIV-2 ---" >> "$SUMMARY"
  grep -v "^#" "$CUTESV_VCF" | awk -v chr="$CHR" -v s="$KIV2_START" -v e="$KIV2_END" \
    '$1==chr && $2>=s-50000 && $2<=e+50000 {print $1, $2, $3, $5, $8}' >> "$SUMMARY" || echo "(none found)" >> "$SUMMARY"
  echo "" >> "$SUMMARY"
done

echo "Done. Summary at $SUMMARY"
