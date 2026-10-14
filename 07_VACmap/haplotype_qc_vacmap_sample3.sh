#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=vacmap_sample3_hapqc_gapfill
#SBATCH --partition=aoraki
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/vacmap_sample3_hapqc_gapfill_%j.out
#SBATCH --error=/path/to/project/logs/vacmap_sample3_hapqc_gapfill_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

SPLIT_DIR=/path/to/project/vacmap_split_bams
OUT=/path/to/project/vacmap_mosdepth/vacmap_haplotype_qc_summary.txt

CHR="NC_060930.1"
LPA_START=161500000
LPA_END=162300000
KIV2_START=161862121
KIV2_END=161998520

get_region_args () {
  case "$1" in
    LPA)           echo "${CHR}:${LPA_START}-${LPA_END}" ;;
    KIV2_unpadded) echo "${CHR}:${KIV2_START}-${KIV2_END}" ;;
  esac
}

qc_block () {
  local bam="$1" region_args="$2" label="$3" out="$4"
  echo "=== ${label} ===" >> "$out"
  local total sa_count softclip_count mapq_tmp
  total=$(samtools view -c "$bam" $region_args)
  echo "total_alignments: $total" >> "$out"
  sa_count=$(samtools view "$bam" $region_args | awk '
    { for(i=12;i<=NF;i++) if($i ~ /^SA:Z:/) c++ } END{print c+0}')
  echo "sa_count: $sa_count" >> "$out"
  softclip_count=$(samtools view "$bam" $region_args | awk '$6 ~ /S/' | wc -l)
  echo "softclip_count: $softclip_count" >> "$out"
  mapq_tmp="${out}.mapq.tmp"
  samtools view "$bam" $region_args | awk '{print $5}' | sort -n > "$mapq_tmp"
  awk '{a[NR]=$1; sum+=$1} END{
    n=NR
    if(n==0){print "mapq: no reads"; exit}
    mean=sum/n
    if(n%2==1) med=a[int((n+1)/2)]; else med=(a[n/2]+a[n/2+1])/2
    q1=a[int(n*0.25)+1]; q3=a[int(n*0.75)+1]
    printf "mapq_mean: %.2f\nmapq_median: %.2f\nmapq_q1: %.2f\nmapq_q3: %.2f\n", mean, med, q1, q3
  }' "$mapq_tmp" >> "$out"
  rm -f "$mapq_tmp"
  if [ "$total" -gt 0 ]; then
    echo "sa_pct: $(awk -v a="$sa_count" -v b="$total" 'BEGIN{printf "%.2f", (a/b)*100}')" >> "$out"
    echo "softclip_pct: $(awk -v a="$softclip_count" -v b="$total" 'BEGIN{printf "%.2f", (a/b)*100}')" >> "$out"
  fi
  echo "" >> "$out"
}

for TOOL in whatshap longphase; do
  for HAP in HP1 HP2; do
    BAM="${SPLIT_DIR}/SAMPLE_3.vacmap.${TOOL}.${HAP}.bam"
    for REGION in LPA KIV2_unpadded; do
      REGION_ARGS=$(get_region_args "$REGION")
      qc_block "$BAM" "$REGION_ARGS" "SAMPLE_3 / ${TOOL} / ${HAP} / ${REGION}" "$OUT"
    done
  done
done

echo "Done. Appended to $OUT"
