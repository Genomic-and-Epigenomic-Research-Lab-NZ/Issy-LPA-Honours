#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
eval "$(conda shell.bash hook)"
set +u; conda activate custom_env; set -u

OUT=/path/to/project/gatk_summary
HG38_ROOT=/path/to/project/PROJECT_hg38

for RS in rs10455872:160589086 rs3798220:160540105 rs41272114:160585045 rs41267813:160577167; do
  ID=${RS%%:*}; POS=${RS##*:}
  for SAMPLE in SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4; do
    GT=$(bcftools query -r chr6:${POS}-${POS} -f '[%GT]' ${OUT}/filtered/${SAMPLE}_filt.vcf.gz || true)
    [ -z "$GT" ] && GT="no call (reference)"
    DP=$(samtools depth -a -Q 20 -r chr6:${POS}-${POS} ${HG38_ROOT}/${SAMPLE}/filtered/${SAMPLE}_hg38_q10.bam | awk '{print $3}')
    printf "%s\t%s\t%s\t%s\tdepth=%s\n" $SAMPLE $ID $POS "$GT" "$DP"
  done
done
