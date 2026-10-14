#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=vacmap_mismatch_rate
#SBATCH --partition=aoraki
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/vacmap_mismatch_rate_%j.out
#SBATCH --error=/path/to/project/logs/vacmap_mismatch_rate_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

SPLIT_DIR=/path/to/project/vacmap_split_bams
OUT=/path/to/project/vacmap_mosdepth/vacmap_mismatch_summary.txt
> "$OUT"

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

mismatch_block () {
  local bam="$1" region_args="$2" label="$3" out="$4"
  echo "=== ${label} ===" >> "$out"

  local mismatch_tmp="${out}.mismatch.tmp"
  samtools view "$bam" $region_args | awk '
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
    }' > "$mismatch_tmp"

  awk '{
    nm_sum += $1
    aln_sum += $2
    n++
  } END{
    if (n == 0) { print "mismatch: no reads"; exit }
    overall_rate = (aln_sum > 0) ? (nm_sum / aln_sum * 100) : 0
    printf "mismatch_reads: %d\nmismatch_mean_nm: %.2f\nmismatch_overall_rate_pct: %.3f\n", n, nm_sum/n, overall_rate
  }' "$mismatch_tmp" >> "$out"
  rm -f "$mismatch_tmp"
  echo "" >> "$out"
}

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)
TOOLS=(whatshap longphase)
HAPLOTYPES=(HP1 HP2)
REGIONS=(LPA KIV2_unpadded)

for SAMPLE in "${SAMPLES[@]}"; do
  for TOOL in "${TOOLS[@]}"; do
    for HAP in "${HAPLOTYPES[@]}"; do
      BAM="${SPLIT_DIR}/${SAMPLE}.vacmap.${TOOL}.${HAP}.bam"
      if [ ! -f "$BAM" ]; then
        echo "MISSING: $BAM" >> "$OUT"
        continue
      fi
      for REGION in "${REGIONS[@]}"; do
        REGION_ARGS=$(get_region_args "$REGION")
        mismatch_block "$BAM" "$REGION_ARGS" "${SAMPLE} / ${TOOL} / ${HAP} / ${REGION}" "$OUT"
      done
    done
  done
done

echo "Done. Results at $OUT"
