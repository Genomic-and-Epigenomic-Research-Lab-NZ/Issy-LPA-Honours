#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=kiv2_finalise
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/kiv2_finalise_%j.out
#SBATCH --error=/path/to/project/logs/kiv2_finalise_%j.err

set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

REF=/path/to/refs/hg38.primary_assembly.fa
HG38_ROOT=/path/to/project/PROJECT_hg38
ARRAY_TSV=/path/to/project/array_depth/results_array_depth.tsv
KILDA_DIR=/path/to/project/KILDA
VALID_DIR=/path/to/project/kilda_validation
SAMTOOLS_REMOTE=/path/to/home/.conda/envs/samtools_remote/bin/samtools
OUT=/path/to/project/kiv2_final
mkdir -p ${OUT}

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)
declare -A FHS=([SAMPLE_1]=SAMPLE_1 [SAMPLE_2]=SAMPLE_2 [SAMPLE_3]=SAMPLE_3 [SAMPLE_4]=SAMPLE_4)

{
  echo "Recorded: $(date)"
  echo "samtools (PATH):  $(which samtools) | $(samtools --version | sed -n 1p)"
  echo "samtools_remote:  $(${SAMTOOLS_REMOTE} --version | sed -n 1p)"
  echo "bcftools:         $(bcftools --version | sed -n 1p)"
  echo "mosdepth:         $(mosdepth --version 2>&1)"
  echo "minimap2:         $(minimap2 --version 2>&1 || echo 'not in custom_env')"
  echo "jellyfish:        $(jellyfish --version 2>&1)"
  echo "KILDA commit:     $(git -C ${KILDA_DIR} log -1 --format='%H (%cd)' --date=short 2>/dev/null || echo 'not a git repo - check how KILDA was installed')"
} > ${OUT}/versions.txt

cat > ${OUT}/sites.tsv <<EOF
chr6	160613786	296	1	T
chr6	160619338	296	2	T
chr6	160624884	296	3	T
chr6	160630428	296	4	T
chr6	160635977	296	5	T
chr6	160641520	296	6	T
chr6	160614754	1264	1	C
chr6	160620306	1264	2	C
chr6	160625852	1264	3	C
chr6	160631396	1264	4	C
chr6	160636945	1264	5	C
chr6	160642488	1264	6	C
EOF

cut -f1,2 ${OUT}/sites.tsv > ${OUT}/regions.tsv

: > ${OUT}/sites_checked.tsv
while IFS=$'\t' read -r CHR POS SNP COPY EXP; do
  OBS=$(samtools faidx ${REF} ${CHR}:${POS}-${POS} | tail -n 1 | tr 'a-z' 'A-Z')
  if [ "$OBS" == "$EXP" ]; then CHK=ok; else CHK=REF_MISMATCH; fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$CHR" "$POS" "$SNP" "$COPY" "$EXP" "$OBS" "$CHK" >> ${OUT}/sites_checked.tsv
done < ${OUT}/sites.tsv

printf "participant\tsample\tsnp\tcopy\tpos\tref_check\tref_reads\tG_reads\tG_fraction\n" > ${OUT}/table_A5_sites.tsv
printf "participant\tsample\tarray_depth_total\tG_frac_296\tG_frac_1264\tsnp_difference\tG_frac_pooled\tref_reads\tG_reads\tphaseable\tallele_G\tallele_other\n" > ${OUT}/allele_sizes.tsv

for SAMPLE in "${SAMPLES[@]}"; do
  echo "=== ${SAMPLE} (${FHS[$SAMPLE]}) | $(date) ==="
  BAM=${HG38_ROOT}/${SAMPLE}/merged/${SAMPLE}_hg38_merged.bam

  if [ ! -f "$BAM" ]; then
    echo "[${SAMPLE}] MISSING BAM: ${BAM}, skipping" >&2
    continue
  fi

  bcftools mpileup -f ${REF} -R ${OUT}/regions.tsv -q 0 -Q 10 -B -I -d 100000 -a FORMAT/AD ${BAM} \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' > ${OUT}/${SAMPLE}_pileup.tsv

  CN=$(awk -F'\t' -v S=${SAMPLE} '$1==S {print $4}' ${ARRAY_TSV})
  if [ -z "$CN" ]; then
    echo "[${SAMPLE}] no array-depth total found in ${ARRAY_TSV}" >&2
    CN=NA
  fi

  awk -F'\t' -v S=${SAMPLE} -v F=${FHS[$SAMPLE]} -v CN=${CN} \
      -v A5=${OUT}/table_A5_sites.tsv -v AS=${OUT}/allele_sizes.tsv '
    FNR==NR {
      n=split($3","$4, A, ","); split($5, D, ",")
      r[$2]=D[1]+0; g[$2]=0
      for (i=2; i<=n; i++) if (A[i]=="G") g[$2]=D[i]+0
      next
    }
    {
      pos=$2; snp=$3; cp=$4; chk=$7
      rr=(pos in r)? r[pos] : 0
      gg=(pos in g)? g[pos] : 0
      f=(rr+gg>0)? gg/(rr+gg) : 0
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.4f\n", F, S, snp, cp, pos, chk, rr, gg, f >> A5
      R[snp]+=rr; G[snp]+=gg; RT+=rr; GT+=gg
    }
    END {
      for (s in R) {
        fp=(R[s]+G[s]>0)? G[s]/(R[s]+G[s]) : 0
        printf "%s\t%s\t%s\tpooled\t-\t-\t%d\t%d\t%.4f\n", F, S, s, R[s], G[s], fp >> A5
        fs[s]=fp
      }
      fpool=(RT+GT>0)? GT/(RT+GT) : 0
      diff=fs["296"]-fs["1264"]; if (diff<0) diff=-diff
      ph=(RT>=10 && GT>=10 && fpool>0.05 && fpool<0.95)? "yes" : "no"
      if (ph=="yes" && CN!="NA") { a1=sprintf("%.2f", CN*fpool); a2=sprintf("%.2f", CN*(1-fpool)) }
      else { a1="-"; a2="-" }
      printf "%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\t%s\t%s\t%s\n", F, S, CN, fs["296"], fs["1264"], diff, fpool, RT, GT, ph, a1, a2 >> AS
    }' ${OUT}/${SAMPLE}_pileup.tsv ${OUT}/sites_checked.tsv
done

{
  echo "Validation result files in ${VALID_DIR} (tsv/txt under 50 kb):"
  find ${VALID_DIR} -maxdepth 3 -type f \( -name "*.tsv" -o -name "*.txt" -o -name "*.CNs" \) -size -50k | sort | while read -r f; do
    echo
    echo "== ${f}"
    cat "${f}"
  done
} > ${OUT}/validation_files.txt

echo
echo "=== Tool versions ==="
cat ${OUT}/versions.txt
echo
echo "=== Table 37: per-copy phasing SNP counts ==="
column -t ${OUT}/table_A5_sites.tsv
echo
echo "=== Allele sizes (exact) ==="
column -t ${OUT}/allele_sizes.tsv
echo
echo "Any REF_MISMATCH rows above mean coordinates or naming are wrong - do not use those rows."
echo "Validation files written to ${OUT}/validation_files.txt"
echo "Done."
