#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=softclip_T2T
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/softclip_T2T_%j.out
#SBATCH --error=/path/to/project/logs/softclip_T2T_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

LPA_DIR=/path/to/project/mini_bams
KIV2_DIR=/path/to/project/t2t_qc

CHR="NC_060930.1"
KIV2_UNPAD_START=161862121
KIV2_UNPAD_END=161998520
KIV2_PAD_START=161812121
KIV2_PAD_END=162048520

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)
REGIONS=(LPA KIV2_unpadded KIV2_padded KIV2_flank)

OUT=/path/to/project/alignment_assessment_summary.txt
> "$OUT"

get_bam_path () {
  local sample="$1"
  local region="$2"
  case "$region" in
    LPA)           echo "${LPA_DIR}/${sample}.LPA.bam" ;;
    KIV2_unpadded) echo "${KIV2_DIR}/${sample}/filtered/${sample}_T2T_KIV2.bam" ;;
    KIV2_padded)   echo "${KIV2_DIR}/${sample}/filtered/${sample}_T2T_KIV2padded.bam" ;;
    KIV2_flank)    echo "${KIV2_DIR}/${sample}/filtered/${sample}_T2T_KIV2padded.bam" ;;
  esac
}

get_region_args () {
  local region="$1"
  case "$region" in
    KIV2_flank)
      echo "${CHR}:${KIV2_PAD_START}-${KIV2_UNPAD_START} ${CHR}:${KIV2_UNPAD_END}-${KIV2_PAD_END}"
      ;;
    *)
      echo ""
      ;;
  esac
}

for SAMPLE in "${SAMPLES[@]}"; do
  for REGION in "${REGIONS[@]}"; do
    BAM=$(get_bam_path "$SAMPLE" "$REGION")
    REGION_ARGS=$(get_region_args "$REGION")

    if [ ! -f "$BAM" ]; then
      echo "MISSING: $BAM" >> "$OUT"
      continue
    fi

    if [ -n "$REGION_ARGS" ] && [ ! -f "${BAM}.bai" ]; then
      samtools index "$BAM"
    fi

    echo "=== ${SAMPLE} / ${REGION} ===" >> "$OUT"

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

    MAPQ_TMP="$(dirname "$OUT")/${SAMPLE}.${REGION}.mapq.txt"
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

    echo "" >> "$OUT"
  done
done

echo "Done. Results at $OUT"
