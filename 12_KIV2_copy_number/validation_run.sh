#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=kilda_validation
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/path/to/project/logs/kilda/%x_%j.out
#SBATCH --error=/path/to/project/logs/kilda/%x_%j.err
set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
eval "$(conda shell.bash hook)"
set +u; conda activate custom_env; set -u

OUT=/path/to/project/kilda_validation
REF=${OUT}/ref/1KG_ONT_VIENNA_hg38.fa
KILDA=/path/to/project/KILDA
K=31
KIV2_KMERS=${KILDA}/data/kmers_GRCh38_LPA_k31/KIV2_hg38_kmers_6copies_specific.tsv
NORM_KMERS=${KILDA}/data/kmers_GRCh38_LPA_k31/LPA_hg38_kmers_1copies_specific.tsv
cd ${OUT}
mkdir -p fastq counts phasing
mapfile -t SAMPLES < available_samples.txt
[ ${#SAMPLES[@]} -gt 0 ] || { echo "STOP: available_samples.txt is empty"; exit 1; }
echo "Samples to validate (${#SAMPLES[@]}): ${SAMPLES[*]}"

cat <<'EOF' | tr -s ' ' '\t' > truth.tsv
HG03836 STU SAS 22.05 15.15 37.20 37.61 21.11 16.51
HG00657 CHS EAS 23.14 7.03 30.17 32.12 25.04 7.08
HG01302 PUR AMR 22.14 8.05 30.19 29.85 21.43 8.41
NA10851 CEU EUR 17.07 10.05 27.13 27.42 17.81 9.62
HG00254 GBR EUR 13.04 11.03 24.07 24.99 13.01 11.97
EOF

: > counts.list
for S in "${SAMPLES[@]}"; do
  samtools fastq -F 0x900 bam/${S}_LPA.bam > fastq/${S}.fastq 2>/dev/null
  jellyfish count -t ${SLURM_CPUS_PER_TASK} -m ${K} -s 100M -C -o counts/${S}_${K}mer fastq/${S}.fastq
  jellyfish dump -c -t counts/${S}_${K}mer > counts/${S}.counts
  printf "%s\t%s\n" ${S} ${OUT}/counts/${S}.counts >> counts.list
done
KOUT=kilda_results_${SLURM_JOB_ID}
python3 ${KILDA}/bin/kilda.py -c counts.list -k ${KIV2_KMERS} -l ${NORM_KMERS} -o ${KOUT}/ -v

: > phasing/regions.tsv; : > phasing/labels.tsv
for P in 160613786 160619338 160624884 160630428 160635977 160641520; do
  printf "chr6\t%s\n" ${P} >> phasing/regions.tsv; printf "%s\tSNP296\tT\n" ${P} >> phasing/labels.tsv; done
for P in 160614754 160620306 160625852 160631396 160636945 160642488; do
  printf "chr6\t%s\n" ${P} >> phasing/regions.tsv; printf "%s\tSNP1264\tC\n" ${P} >> phasing/labels.tsv; done
sort -k2,2n phasing/regions.tsv -o phasing/regions.tsv
: > phasing/fractions.tsv
for S in "${SAMPLES[@]}"; do
  bcftools mpileup -f ${REF} -R phasing/regions.tsv -q 0 -Q 10 -B -d 100000 -a FORMAT/AD --ignore-RG -Ou bam/${S}_LPA.bam 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' > phasing/${S}_sites.tsv
  awk -v S=${S} 'BEGIN{FS=OFS="\t"}
    FNR==NR {lab[$1]=$2; er[$1]=$3; next}
    ($2 in lab) { if ($3!=er[$2]) bad++
      na=split($4,al,","); split($5,ad,","); g=0
      for (i=1;i<=na;i++) if (al[i]=="G") g=ad[i+1]
      r+=ad[1]; a+=g }
    END { t=r+a; printf "%s\t%d\t%d\t%.3f\t%s\n", S, r, a, (t>0? a/t:0), (bad? "REF_MISMATCH":"ok") }' \
    phasing/labels.tsv phasing/${S}_sites.tsv >> phasing/fractions.tsv
done
echo "=== Phasing read counts (sample, ref reads, G reads, G fraction, ref check) ==="
cat phasing/fractions.tsv

printf "sample\tpop\tKILDA\tBionano_total\tDRAGEN\tKILDA/Bionano\tG_frac\tphased\tour_larger\tour_smaller\tBN_larger\tBN_smaller\tref_check\n" > validation_results.tsv
awk 'BEGIN{FS=OFS="\t"}
  FILENAME==ARGV[1] {pop[$1]=$2; bl[$1]=$4; bs[$1]=$5; bt[$1]=$6; dt[$1]=$7; next}
  FILENAME==ARGV[2] {f[$1]=$4; chk[$1]=$5; next}
  FNR>1 { s=$1; k=$2; g=f[s]; ph=(g>=0.05 && g<=0.95)? "yes":"no"
          a=k*g; b=k*(1-g); L=(a>b)?a:b; Sm=(a>b)?b:a
          printf "%s\t%s\t%.2f\t%.2f\t%.2f\t%.3f\t%.3f\t%s\t%s\t%s\t%.2f\t%.2f\t%s\n", s, pop[s], k, bt[s], dt[s], k/bt[s], g, ph,
                 (ph=="yes")? sprintf("%.1f",L):"-", (ph=="yes")? sprintf("%.1f",Sm):"-", bl[s], bs[s], chk[s] }' \
  truth.tsv phasing/fractions.tsv ${KOUT}/kilda_kiv2.CNs >> validation_results.tsv
echo "=== Validation: KILDA + phasing on 1KG ONT vs Bionano/DRAGEN ==="
column -t validation_results.tsv
awk -F'\t' 'NR>1 {n++; r+=$6; d=$3-$4; e+=(d<0?-d:d)} END {printf "\nn=%d | mean KILDA/Bionano = %.3f | mean abs difference = %.2f copies\n", n, r/n, e/n}' validation_results.tsv
