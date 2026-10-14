#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=vacmap_k9_t1_test
#SBATCH --partition=aoraki
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=/path/to/project/logs/vacmap_k9_t1_test_%j.out
#SBATCH --error=/path/to/project/logs/vacmap_k9_t1_test_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate vacmap_env
set -u

REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
FASTQ=/path/to/project/vacmap_out/SAMPLE_3.LPA.fastq
OUT_DIR=/path/to/project/vacmap_out/k9_t1_test
mkdir -p "$OUT_DIR"
BAM="${OUT_DIR}/SAMPLE_3_vacmap_k9_t1.sorted.bam"

echo "Running VACmap with k=9 (paper's local-index k-mer for repeat resolution)..."
vacmap -ref "$REF" -read "$FASTQ" -mode H -k 9 -t 1 --force \
  | samtools sort -@4 > "$BAM"
samtools index "$BAM"

echo "=== Read counts ==="
samtools flagstat "$BAM"

CHR="NC_060930.1"
KIV2_START=161862121
KIV2_END=161998520

echo ""
echo "=== KIV-2 region ==="
TOTAL=$(samtools view -c "$BAM" "${CHR}:${KIV2_START}-${KIV2_END}")
SA_COUNT=$(samtools view "$BAM" "${CHR}:${KIV2_START}-${KIV2_END}" | awk '
  { for(i=12;i<=NF;i++) if($i ~ /^SA:Z:/) c++ } END{print c+0}')
echo "total_alignments: $TOTAL"
echo "sa_count: $SA_COUNT"
if [ "$TOTAL" -gt 0 ]; then
  echo "sa_pct: $(awk -v a="$SA_COUNT" -v b="$TOTAL" 'BEGIN{printf "%.2f", (a/b)*100}')"
fi

MAPQ_TMP="${OUT_DIR}/sample3_k9_mapq.tmp"
samtools view "$BAM" "${CHR}:${KIV2_START}-${KIV2_END}" | awk '{print $5}' | sort -n > "$MAPQ_TMP"
awk '{a[NR]=$1; sum+=$1} END{
  n=NR
  if(n==0){print "no reads"; exit}
  mean=sum/n
  if(n%2==1) med=a[int((n+1)/2)]; else med=(a[n/2]+a[n/2+1])/2
  q1=a[int(n*0.25)+1]; q3=a[int(n*0.75)+1]
  printf "mapq_mean: %.2f\nmapq_median: %.2f\nmapq_q1: %.2f\nmapq_q3: %.2f\n", mean, med, q1, q3
}' "$MAPQ_TMP"
rm -f "$MAPQ_TMP"

echo ""
echo "Done."
