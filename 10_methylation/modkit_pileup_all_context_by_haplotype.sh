#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=methylation_by_haplotype
#SBATCH --partition=aoraki
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/methylation_by_haplotype_%j.out
#SBATCH --error=/path/to/project/logs/methylation_by_haplotype_%j.err

set -uo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

if ! command -v modkit &> /dev/null; then
  echo "ERROR: modkit not found in custom_env." >&2
  exit 1
fi

REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
WHATSHAP_SPLIT_DIR=/path/to/project/haplotype_split_bams
LONGPHASE_SPLIT_DIR=/path/to/project/longphase_split_bams
OUT_DIR=/path/to/project/methylation_by_haplotype
mkdir -p "$OUT_DIR"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

LPA_START=161500000
LPA_END=162300000
KIV2_START=161862121
KIV2_END=161998520

SUMMARY="${OUT_DIR}/methylation_by_haplotype_summary.txt"
> "$SUMMARY"

get_bed_path () {
  local sample="$1" tool="$2" hap="$3"
  if [ "$tool" == "whatshap" ]; then
    echo "${WHATSHAP_SPLIT_DIR}/${sample}.LPA.${hap}.bam"
  else
    echo "${LONGPHASE_SPLIT_DIR}/${sample}.LPA.longphase.${hap}.bam"
  fi
}

for SAMPLE in "${SAMPLES[@]}"; do
  for TOOL in whatshap longphase; do
    for HAP in HP1 HP2; do
      BAM=$(get_bed_path "$SAMPLE" "$TOOL" "$HAP")
      OUT_BED="${OUT_DIR}/${SAMPLE}.${TOOL}.${HAP}.modkit_pileup.bed"

      if [ ! -f "$BAM" ]; then
        echo "[${SAMPLE} / ${TOOL} / ${HAP}] MISSING BAM: ${BAM}" >> "$SUMMARY"
        continue
      fi

      if [ ! -f "$OUT_BED" ]; then
        echo "[${SAMPLE} / ${TOOL} / ${HAP}] Running modkit pileup..."
        modkit pileup "$BAM" "$OUT_BED" --ref "$REF" --threads "$SLURM_CPUS_PER_TASK" 2>/dev/null
      fi

      echo "=== ${SAMPLE} / ${TOOL} / ${HAP} ===" >> "$SUMMARY"

      for CODE in a h m; do
        case "$CODE" in
          a) LABEL="6mA" ;;
          h) LABEL="5hmC" ;;
          m) LABEL="5mC" ;;
        esac

        LPA_MEAN=$(awk -v code="$CODE" -v s="$LPA_START" -v e="$LPA_END" \
          '$2>=s && $2<=e && $4==code {sum+=$11; n++} END{if(n>0) printf "%.4f", sum/n; else print "NA"}' "$OUT_BED")

        KIV2_MEAN=$(awk -v code="$CODE" -v s="$KIV2_START" -v e="$KIV2_END" \
          '$2>=s && $2<=e && $4==code {sum+=$11; n++} END{if(n>0) printf "%.4f", sum/n; else print "NA"}' "$OUT_BED")

        echo "${LABEL}: LPA=${LPA_MEAN}%, KIV2=${KIV2_MEAN}%" >> "$SUMMARY"
      done
      echo "" >> "$SUMMARY"
    done
  done
done

echo "Done. Summary at $SUMMARY"
