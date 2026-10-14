#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=vacmap_qc_sample3
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/vacmap_qc_sample3_%j.out
#SBATCH --error=/path/to/project/logs/vacmap_qc_sample3_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate vacmap_env
set -u

BAM=/path/to/project/vacmap_out/SAMPLE_3_piped.sorted.bam

if [ ! -f "${BAM}.bai" ]; then
  echo "Index not found, indexing BAM..."
  samtools index "$BAM"
fi

CHR="NC_060930.1"
LPA_START=161500000
LPA_END=162300000
KIV2_START=161862121
KIV2_END=161998520
KIV2_PAD_START=161812121
KIV2_PAD_END=162048520

OUT=/path/to/project/vacmap_out/SAMPLE_3_vacmap_qc.txt
> "$OUT"

get_region_args () {
  local region="$1"
  case "$region" in
    LPA)           echo "${CHR}:${LPA_START}-${LPA_END}" ;;
    KIV2_unpadded) echo "${CHR}:${KIV2_START}-${KIV2_END}" ;;
    KIV2_padded)   echo "${CHR}:${KIV2_PAD_START}-${KIV2_PAD_END}" ;;
    FLANK)         echo "${CHR}:${KIV2_PAD_START}-${KIV2_START} ${CHR}:${KIV2_END}-${KIV2_PAD_END}" ;;
  esac
}

REGIONS=(LPA KIV2_unpadded KIV2_padded FLANK)

for REGION in "${REGIONS[@]}"; do
  REGION_ARGS=$(get_region_args "$REGION")

  echo "=== SAMPLE_3 / VACmap / ${REGION} ===" >> "$OUT"

  TOTAL=$(samtools view -c "$BAM" $REGION_ARGS)
  echo "total_alignments: $TOTAL" >> "$OUT"

  SA_COUNT=$(samtools view "$BAM" $REGION_ARGS | awk '
    {
      for(i=12;i<=NF;i++)
        if($i ~ /^SA:Z:/)
          c++
    }
    END{print c+0}')
  echo "sa_count: $SA_COUNT" >> "$OUT"

  SOFTCLIP_COUNT=$(samtools view "$BAM" $REGION_ARGS | awk '$6 ~ /S/' | wc -l)
  echo "softclip_count: $SOFTCLIP_COUNT" >> "$OUT"

  MAPQ_TMP="$(dirname "$OUT")/SAMPLE_3.vacmap.${REGION}.mapq.txt"
  samtools view "$BAM" $REGION_ARGS | awk '{print $5}' | sort -n > "$MAPQ_TMP"
  awk '{a[NR]=$1; sum+=$1} END{
    n=NR
    if(n==0){print "mapq: no reads"; exit}
    mean=sum/n
    if(n%2==1) med=a[int((n+1)/2)]; else med=(a[n/2]+a[n/2+1])/2
    q1=a[int(n*0.25)+1]
    q3=a[int(n*0.75)+1]
    printf "mapq_mean: %.2f\nmapq_median: %.2f\nmapq_q1: %.2f\nmapq_q3: %.2f\n", mean, med, q1, q3
  }' "$MAPQ_TMP" >> "$OUT"
  rm -f "$MAPQ_TMP"

  if [ "$TOTAL" -gt 0 ]; then
    SA_PCT=$(awk -v a="$SA_COUNT" -v b="$TOTAL" 'BEGIN{printf "%.2f", (a/b)*100}')
    SOFTCLIP_PCT=$(awk -v a="$SOFTCLIP_COUNT" -v b="$TOTAL" 'BEGIN{printf "%.2f", (a/b)*100}')
  else
    SA_PCT="NA"; SOFTCLIP_PCT="NA"
  fi
  echo "sa_pct: $SA_PCT" >> "$OUT"
  echo "softclip_pct: $SOFTCLIP_PCT" >> "$OUT"

  echo "" >> "$OUT"
done

echo "Done. Results at $OUT"
