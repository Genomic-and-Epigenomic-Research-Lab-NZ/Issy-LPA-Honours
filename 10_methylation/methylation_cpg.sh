#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=methylation_cpg_rerun
#SBATCH --partition=aoraki
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/methylation_cpg_rerun_%j.out
#SBATCH --error=/path/to/project/logs/methylation_cpg_rerun_%j.err

set -uo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
WHOLE_DIR=/path/to/project/mini_bams
WHATSHAP_SPLIT_DIR=/path/to/project/haplotype_split_bams
LONGPHASE_SPLIT_DIR=/path/to/project/longphase_split_bams
OLD_ALLCONTEXT_DIR=/path/to/project/methylation_by_haplotype
OUT_DIR=/path/to/project/methylation_cpg
mkdir -p "$OUT_DIR"

CHR="NC_060930.1"
LPA_S=161500000;  LPA_E=162300000
KIV_S=161862121;  KIV_E=161998520
FL_LS=161812121;  FL_RE=162048520

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)
MINCOVS=(1 3 5)

SUMMARY="${OUT_DIR}/methylation_cpg_summary.tsv"
printf "sample\ttool\thaplotype\tregion\tmod\tmin_cov\tn_sites\tmean_site_pct\tpooled_pct\n" > "$SUMMARY"

summarise () {
  local bed="$1" sample="$2" tool="$3" hap="$4" mods="$5"
  for MOD in $mods; do
    for MC in "${MINCOVS[@]}"; do
      for REGION in LPA KIV2 FLANK; do
        awk -v chr="$CHR" -v mod="$MOD" -v mc="$MC" -v reg="$REGION" \
            -v ls="$LPA_S" -v le="$LPA_E" -v ks="$KIV_S" -v ke="$KIV_E" \
            -v fls="$FL_LS" -v fre="$FL_RE" \
            -v sample="$sample" -v tool="$tool" -v hap="$hap" '
          BEGIN{OFS="\t"}
          $1==chr && $4==mod && $10>=mc {
            p=$2; keep=0
            if (reg=="LPA"   && p>=ls && p<=le) keep=1
            if (reg=="KIV2"  && p>=ks && p<=ke) keep=1
            if (reg=="FLANK" && ((p>=fls && p<ks) || (p>ke && p<=fre))) keep=1
            if (keep) { n++; sum+=$11; nm+=$12; nv+=$10 }
          }
          END{
            mean=(n>0)?sprintf("%.4f",sum/n):"NA"
            pooled=(nv>0)?sprintf("%.4f",nm/nv*100):"NA"
            print sample,tool,hap,reg,mod,mc,n+0,mean,pooled
          }' "$bed" >> "$SUMMARY"
      done
    done
  done
}

run_pileup () {
  local bam="$1" out="$2"
  if [ ! -f "$bam" ]; then echo "MISSING BAM: $bam" >&2; return 1; fi
  if [ ! -s "$out" ]; then
    modkit pileup "$bam" "$out" --ref "$REF" --cpg --combine-strands \
      --threads "$SLURM_CPUS_PER_TASK" 2>> "${OUT_DIR}/modkit_errors.log"
    if [ $? -ne 0 ] || [ ! -s "$out" ]; then
      echo "FAILED: modkit pileup on $bam (see ${OUT_DIR}/modkit_errors.log)" >&2
      rm -f "$out"
      return 1
    fi
  fi
  return 0
}

for SAMPLE in "${SAMPLES[@]}"; do
  OUT="${OUT_DIR}/${SAMPLE}.whole.cpg.bed"
  run_pileup "${WHOLE_DIR}/${SAMPLE}.LPA.bam" "$OUT" && summarise "$OUT" "$SAMPLE" "none" "all" "m h"

  for TOOL in longphase; do
    for HAP in HP1 HP2; do
      if [ "$TOOL" == "whatshap" ]; then
        BAM="${WHATSHAP_SPLIT_DIR}/${SAMPLE}.LPA.${HAP}.bam"
      else
        BAM="${LONGPHASE_SPLIT_DIR}/${SAMPLE}.LPA.longphase.${HAP}.bam"
      fi
      OUT="${OUT_DIR}/${SAMPLE}.${TOOL}.${HAP}.cpg.bed"
      run_pileup "$BAM" "$OUT" && summarise "$OUT" "$SAMPLE" "$TOOL" "$HAP" "m h"
    done
  done
done

for SAMPLE in "${SAMPLES[@]}"; do
  OLD_WHOLE="/path/to/project/modkit_pileup/${SAMPLE}.modkit_pileup.bed"
  [ -f "$OLD_WHOLE" ] && summarise "$OLD_WHOLE" "$SAMPLE" "none" "all" "a"
  for TOOL in longphase; do
    for HAP in HP1 HP2; do
      OLD="${OLD_ALLCONTEXT_DIR}/${SAMPLE}.${TOOL}.${HAP}.modkit_pileup.bed"
      [ -f "$OLD" ] && summarise "$OLD" "$SAMPLE" "$TOOL" "$HAP" "a"
    done
  done
done

echo "Done. Summary: $SUMMARY"
