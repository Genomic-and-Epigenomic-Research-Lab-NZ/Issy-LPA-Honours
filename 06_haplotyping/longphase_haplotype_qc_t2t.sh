#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=longphase_qc
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/longphase_qc_%j.out
#SBATCH --error=/path/to/project/logs/longphase_qc_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

SPLIT_DIR=/path/to/project/longphase_split_bams

CHR="NC_060930.1"
LPA_START=161500000
LPA_END=162300000
KIV2_START=161862121
KIV2_END=161998520

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)
HAPLOTYPES=(HP1 HP2)
REGIONS=(LPA KIV2)

OUT=/path/to/project/longphase_qc_summary.txt
> "$OUT"

get_region_args () {
  local region="$1"
  case "$region" in
    LPA)  echo "${CHR}:${LPA_START}-${LPA_END}" ;;
    KIV2) echo "${CHR}:${KIV2_START}-${KIV2_END}" ;;
  esac
}

for SAMPLE in "${SAMPLES[@]}"; do
  for HAP in "${HAPLOTYPES[@]}"; do
    BAM="${SPLIT_DIR}/${SAMPLE}.LPA.longphase.${HAP}.bam"

    if [ ! -f "$BAM" ]; then
      echo "MISSING: $BAM" >> "$OUT"
      continue
    fi

    for REGION in "${REGIONS[@]}"; do
      REGION_ARGS=$(get_region_args "$REGION")

      echo "=== ${SAMPLE} / ${HAP} / ${REGION} ===" >> "$OUT"

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

      MAPQ_TMP="$(dirname "$OUT")/${SAMPLE}.${HAP}.${REGION}.longphase.mapq.txt"
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

      MISMATCH_TMP="$(dirname "$OUT")/${SAMPLE}.${HAP}.${REGION}.longphase.mismatch.txt"
      samtools view "$BAM" $REGION_ARGS | awk '
        {
          nm = -1
          for (i = 12; i <= NF; i++) {
            if ($i ~ /^NM:i:/) {
              split($i, a, ":")
              nm = a[3]
            }
          }
          if (nm == -1) next

          cigar = $6
          num = ""
          alnlen = 0
          for (j = 1; j <= length(cigar); j++) {
            c = substr(cigar, j, 1)
            if (c ~ /[0-9]/) {
              num = num c
            } else {
              if (c == "M" || c == "I" || c == "D") {
                alnlen += num
              }
              num = ""
            }
          }
          if (alnlen > 0) {
            print nm, alnlen
          }
        }' > "$MISMATCH_TMP"

      awk '{
        nm_sum += $1
        aln_sum += $2
        n++
      } END{
        if (n == 0) { print "mismatch: no reads"; exit }
        overall_rate = (aln_sum > 0) ? (nm_sum / aln_sum * 100) : 0
        printf "mismatch_reads: %d\nmismatch_mean_nm: %.2f\nmismatch_overall_rate_pct: %.3f\n", n, nm_sum/n, overall_rate
      }' "$MISMATCH_TMP" >> "$OUT"
      rm -f "$MISMATCH_TMP"

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
  done
done

echo "Done. Results at $OUT"
