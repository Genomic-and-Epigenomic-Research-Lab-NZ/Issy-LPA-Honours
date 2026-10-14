#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=gatk_summary
#SBATCH --account=ACCOUNT
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/gatk_summary_%j.out
#SBATCH --error=/path/to/project/logs/gatk_summary_%j.err

set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

VCF_DIR=/path/to/project/gatk_hg38

HG38_ROOT=/path/to/project/PROJECT_hg38
BAM_PATTERN='${HG38_ROOT}/${SAMPLE}/filtered/${SAMPLE}_hg38_q10.bam'

SR_VCF_DIR=""

REF=/path/to/refs/hg38.primary_assembly.fa
OUT=/path/to/project/gatk_summary
mkdir -p ${OUT}/filtered

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)
declare -A FHS=([SAMPLE_1]=SAMPLE_1 [SAMPLE_2]=SAMPLE_2 [SAMPLE_3]=SAMPLE_3 [SAMPLE_4]=SAMPLE_4)

LPA_WIN=chr6:160300000-160800000
KIV2=chr6:160612753-160645586
RS10455872=160589086
RS3798220=160540105

printf "participant\tsample\traw_total\traw_SNV\traw_indel\traw_in_KIV2\tfilt_total\tfilt_SNV\tfilt_indel\tfilt_het\tfilt_homalt\tfilt_in_KIV2\n" > ${OUT}/variant_counts.tsv
printf "participant\tsample\trsID\tposition\tVCF_genotype\tread_depth_MAPQ20\n" > ${OUT}/rsid_genotypes.tsv
printf "participant\tposition\ttype\tgenotype\n" > ${OUT}/variant_positions.tsv

FILT_LIST=()

for SAMPLE in "${SAMPLES[@]}"; do
  F=${FHS[$SAMPLE]}
  echo "=== ${SAMPLE} (${F}) | $(date) ==="

  VCF=${VCF_DIR}/${SAMPLE}_hg38_LPA.vcf.gz
  if [ ! -f "$VCF" ]; then
    echo "[${SAMPLE}] no VCF found in ${VCF_DIR}, skipping" >&2
    continue
  fi
  echo "[${SAMPLE}] using ${VCF}"

  bcftools norm -m -any -f ${REF} ${VCF} -Ou \
    | bcftools view -t ${LPA_WIN} -Oz -o ${OUT}/${SAMPLE}_norm.vcf.gz
  bcftools index -f ${OUT}/${SAMPLE}_norm.vcf.gz

  bcftools view -i 'QUAL>=30 && FORMAT/DP>=10' ${OUT}/${SAMPLE}_norm.vcf.gz \
    -Oz -o ${OUT}/filtered/${SAMPLE}_filt.vcf.gz
  bcftools index -f ${OUT}/filtered/${SAMPLE}_filt.vcf.gz
  FILT_LIST+=(${OUT}/filtered/${SAMPLE}_filt.vcf.gz)

  N=${OUT}/${SAMPLE}_norm.vcf.gz
  P=${OUT}/filtered/${SAMPLE}_filt.vcf.gz
  RT=$(bcftools view -H $N | wc -l)
  RS=$(bcftools view -H -v snps $N | wc -l)
  RI=$(bcftools view -H -v indels $N | wc -l)
  RK=$(bcftools view -H -r ${KIV2} $N | wc -l)
  FT=$(bcftools view -H $P | wc -l)
  FS=$(bcftools view -H -v snps $P | wc -l)
  FI=$(bcftools view -H -v indels $P | wc -l)
  FH=$(bcftools view -H -g het $P | wc -l)
  FA=$(bcftools view -H -g hom $P | wc -l)
  FK=$(bcftools view -H -r ${KIV2} $P | wc -l)
  printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" \
    $F $SAMPLE $RT $RS $RI $RK $FT $FS $FI $FH $FA $FK >> ${OUT}/variant_counts.tsv

  bcftools query -f "${F}\t%POS\t%TYPE\t[%GT]\n" $P >> ${OUT}/variant_positions.tsv

  BAM=$(eval echo ${BAM_PATTERN})
  [ -f "$BAM" ] || BAM=${HG38_ROOT}/${SAMPLE}/merged/${SAMPLE}_hg38_merged.bam
  echo "[${SAMPLE}] depth check BAM: ${BAM}"
  for RS in rs10455872:${RS10455872} rs3798220:${RS3798220}; do
    ID=${RS%%:*}; POS=${RS##*:}
    GT=$(bcftools query -r chr6:${POS}-${POS} -f '[%GT]' $N || true)
    [ -z "$GT" ] && GT="no call (reference)"
    if [ -f "$BAM" ]; then
      DP=$(samtools depth -a -Q 20 -r chr6:${POS}-${POS} ${BAM} | awk '{print $3}')
    else
      DP="BAM not found"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" $F $SAMPLE $ID $POS "$GT" "$DP" >> ${OUT}/rsid_genotypes.tsv
  done
done

if [ ${#FILT_LIST[@]} -gt 1 ]; then
  bcftools isec -n+1 -c none ${FILT_LIST[@]} -p ${OUT}/isec_all > /dev/null
  awk -F'\t' '{n=gsub(/1/,"1",$5); c[n]++} END {for (k in c) printf "carried_by_%s_participants\t%d\n", k, c[k]}' \
    ${OUT}/isec_all/sites.txt | sort > ${OUT}/shared_counts.tsv

  bcftools merge --force-samples -m none ${FILT_LIST[@]} -Ou \
    | bcftools view -G -Ou \
    | bcftools annotate -x INFO,ID,QUAL,FILTER -Oz -o ${OUT}/LPA_sites_for_VEP.vcf.gz
  bcftools index -f ${OUT}/LPA_sites_for_VEP.vcf.gz
fi

if [ -n "$SR_VCF_DIR" ]; then
  printf "participant\tnanopore_only\tshort_read_only\tboth\n" > ${OUT}/concordance.tsv
  for SAMPLE in "${SAMPLES[@]}"; do
    F=${FHS[$SAMPLE]}
    SR=${SR_VCF_DIR}/${F}.LPA.vcf.gz
    LR=${OUT}/filtered/${SAMPLE}_filt.vcf.gz
    [ -f "$SR" ] && [ -f "$LR" ] || { echo "[${F}] short-read or nanopore VCF missing, skipping" >&2; continue; }
    bcftools norm -m -any -f ${REF} ${SR} -Ou \
      | bcftools view -t ${LPA_WIN} -i 'QUAL>=30 && FORMAT/DP>=10' -Ou \
      | bcftools view -t ^${KIV2} -Oz -o ${OUT}/${F}_SR_filt.vcf.gz
    bcftools index -f ${OUT}/${F}_SR_filt.vcf.gz
    bcftools view -t ^${KIV2} ${LR} -Oz -o ${OUT}/${F}_LR_noKIV2.vcf.gz
    bcftools index -f ${OUT}/${F}_LR_noKIV2.vcf.gz
    bcftools isec -c none ${OUT}/${F}_LR_noKIV2.vcf.gz ${OUT}/${F}_SR_filt.vcf.gz -p ${OUT}/isec_${F} > /dev/null
    LO=$(grep -vc '^#' ${OUT}/isec_${F}/0000.vcf || true)
    SO=$(grep -vc '^#' ${OUT}/isec_${F}/0001.vcf || true)
    BO=$(grep -vc '^#' ${OUT}/isec_${F}/0002.vcf || true)
    printf "%s\t%s\t%s\t%s\n" $F $LO $SO $BO >> ${OUT}/concordance.tsv
  done
fi

echo; echo "=== Variant counts (LPA window) ==="; column -t ${OUT}/variant_counts.tsv
echo; echo "=== rs10455872 / rs3798220 ==="; column -t -s $'\t' ${OUT}/rsid_genotypes.tsv
[ -f ${OUT}/shared_counts.tsv ] && { echo; echo "=== Shared variants ==="; column -t ${OUT}/shared_counts.tsv; }
[ -f ${OUT}/concordance.tsv ] && { echo; echo "=== Nanopore vs short read (outside KIV-2) ==="; column -t ${OUT}/concordance.tsv; }
echo; echo "VEP upload file: ${OUT}/LPA_sites_for_VEP.vcf.gz"
echo "Positions for figure: ${OUT}/variant_positions.tsv"
echo "Done."
