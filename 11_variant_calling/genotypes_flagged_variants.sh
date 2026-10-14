#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=assoc_gt
#SBATCH --account=ACCOUNT
#SBATCH --partition=aoraki
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/path/to/project/logs/assoc_gt_%j.out
#SBATCH --error=/path/to/project/logs/assoc_gt_%j.err

set -euo pipefail
eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

SITES=/path/to/project/gatk_summary/lpa_assoc_sites.tsv
VCF_DIR=/path/to/project/gatk_summary/filtered
OUT=/path/to/project/gatk_summary
SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)
declare -A FHS=([SAMPLE_1]=SAMPLE_1 [SAMPLE_2]=SAMPLE_2 [SAMPLE_3]=SAMPLE_3 [SAMPLE_4]=SAMPLE_4)

cut -f1,2 ${SITES} > ${OUT}/assoc_regions.tsv

for S in "${SAMPLES[@]}"; do
  bcftools query -R ${OUT}/assoc_regions.tsv -f '%POS\t[%GT]\t[%AD]\n' ${VCF_DIR}/${S}_filt.vcf.gz \
    > ${OUT}/assoc_${S}.tsv
done

printf "pos\trsID\tgene\tconsequence\tgnomAD_AF\tClinVar\ttraits" > ${OUT}/assoc_genotypes.tsv
for S in "${SAMPLES[@]}"; do printf "\t%s" ${FHS[$S]} >> ${OUT}/assoc_genotypes.tsv; done
printf "\tn_carriers\n" >> ${OUT}/assoc_genotypes.tsv

awk -F'\t' -v OFS='\t' '
  FILENAME ~ /assoc_SAMPLE_1/ {a[$1]=$2" ("$3")"; next}
  FILENAME ~ /assoc_SAMPLE_2/ {b[$1]=$2" ("$3")"; next}
  FILENAME ~ /assoc_SAMPLE_3/ {c[$1]=$2" ("$3")"; next}
  FILENAME ~ /assoc_SAMPLE_4/ {d[$1]=$2" ("$3")"; next}
  {
    p=$2; n=0
    ga=(p in a)?a[p]:"ref"; gb=(p in b)?b[p]:"ref"; gc=(p in c)?c[p]:"ref"; gd=(p in d)?d[p]:"ref"
    n=(ga!="ref")+(gb!="ref")+(gc!="ref")+(gd!="ref")
    gene=($4=="")?"intergenic":$4
    print p,$3,gene,$5,$6,$7,$8,ga,gb,gc,gd,n
  }' ${OUT}/assoc_SAMPLE_1.tsv ${OUT}/assoc_SAMPLE_2.tsv ${OUT}/assoc_SAMPLE_3.tsv ${OUT}/assoc_SAMPLE_4.tsv ${SITES} \
  >> ${OUT}/assoc_genotypes.tsv

echo "=== Coding variants ==="
awk -F'\t' 'NR==1 || $4 ~ /missense|frameshift/' ${OUT}/assoc_genotypes.tsv | cut -f1-6,8-12 | column -t -s $'\t'
echo
echo "=== Associated variants: number carried per participant ==="
awk -F'\t' 'NR>1 && $7!="" {for(i=8;i<=11;i++) if($i!="ref") n[i]++; tot++}
  END {print "sites with association:", tot; print "SAMPLE_1", n[8]+0; print "SAMPLE_2", n[9]+0; print "SAMPLE_3", n[10]+0; print "SAMPLE_4", n[11]+0}' ${OUT}/assoc_genotypes.tsv
echo
echo "=== Associated variants with gnomAD AF < 0.05 ==="
awk -F'\t' 'NR==1 || ($7!="" && $5!="" && $5<0.05 && $12>0)' ${OUT}/assoc_genotypes.tsv | cut -f1-3,5,7-12 | column -t -s $'\t'
echo "Full table: ${OUT}/assoc_genotypes.tsv"
