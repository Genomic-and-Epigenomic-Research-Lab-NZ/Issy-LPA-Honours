#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=pooled_kiv2_check
#SBATCH --partition=aoraki
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/pooled_kiv2_check_%j.out
#SBATCH --error=/path/to/project/logs/pooled_kiv2_check_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

REF_INDEX=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.w10_k15.mmi
ASSEMBLY=/path/to/project/flye_assembly/assemblies/SAMPLE_3_pooled_all/assembly.fasta
OUT_DIR=/path/to/project/flye_assembly/assembly_vs_t2t
SAM_OUT="${OUT_DIR}/SAMPLE_3_pooled_all_vs_T2T.sam"
CSV_OUT="${OUT_DIR}/SAMPLE_3_pooled_all_kiv2_copies.csv"
SCRIPT_PATH="${OUT_DIR}/estimate_kiv2_copies.py"

echo "Aligning pooled assembly to T2T..."
minimap2 -ax map-ont "$REF_INDEX" "$ASSEMBLY" > "$SAM_OUT"

echo "=== SAM file check ==="
head -5 "$SAM_OUT"
echo "Total lines: $(wc -l < "$SAM_OUT")"

echo ""
echo "=== All alignment records (incl. secondary/supplementary) ==="
samtools view "$SAM_OUT" | awk '{print $1, $2, $3, $4, length($10)}'

echo ""
echo "=== KIV-2 overlap quantification ==="
if [ -f "$SCRIPT_PATH" ]; then
  python3 "$SCRIPT_PATH" "$SAM_OUT" "$CSV_OUT"
  cat "$CSV_OUT"
else
  echo "estimate_kiv2_copies.py not found at ${SCRIPT_PATH} -- run kiv2_quantification.sh first, or copy it here."
fi

echo ""
echo "Done."
