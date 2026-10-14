#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=depth_ratio
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/path/to/project/logs/depth_ratio/%x_%j.out
#SBATCH --error=/path/to/project/logs/depth_ratio/%x_%j.err

set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
T=$SLURM_CPUS_PER_TASK

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

for tool in samtools minimap2 mosdepth; do
  command -v ${tool} >/dev/null || { echo "MISSING TOOL: ${tool} (not in custom_env)"; exit 1; }
done

REF_HG38=/path/to/refs/hg38.primary_assembly.fa
HG38_ROOT=/path/to/project/PROJECT_hg38
KILDA_CN=/path/to/project/kilda_out_all/results/kilda_kiv2.CNs
SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

WIN_CHR=chr6
WIN_START=160300000
WIN_END=160800000
WIN_LEN=$(( WIN_END - WIN_START + 1 ))

UNIT_START=160624589
UNIT_END=160630132

OUT=/path/to/project/depth_ratio
mkdir -p ${OUT}/ref ${OUT}/fastq ${OUT}/bam ${OUT}/depth ${OUT}/stats
cd ${OUT}

samtools faidx ${REF_HG38} ${WIN_CHR}:${WIN_START}-${WIN_END} | sed "1s/.*/>window/" > ref/window.fa
samtools faidx ${REF_HG38} ${WIN_CHR}:${UNIT_START}-${UNIT_END} | sed "1s/.*/>unit/" > ref/unit.fa
samtools faidx ref/window.fa
ULEN=$(( UNIT_END - UNIT_START + 1 ))

minimap2 -c -x asm5 -N 50 --secondary=yes -t ${T} ref/window.fa ref/unit.fa > ref/unit_hits.paf
awk 'BEGIN{OFS="\t"} ($11 >= 5000 && $10/$11 >= 0.95) {print $8, $9, $10/$11}' ref/unit_hits.paf \
  | sort -k1,1n > ref/kiv2_copies.tsv
echo "KIV-2 copies found in hg38 window (0-based start, end, identity):"
cat ref/kiv2_copies.tsv
NCOPIES=$(wc -l < ref/kiv2_copies.tsv)
FIRST=$(head -1 ref/kiv2_copies.tsv | cut -f1)
LAST=$(tail -1 ref/kiv2_copies.tsv | cut -f2)
SPAN=$(( LAST - FIRST ))
echo "Hits passing filter: ${NCOPIES} | array span: ${SPAN} bp (expect ~33,281 = 6 units; Behera et al. 2024)"
if [ "${NCOPIES}" -lt 4 ] || [ "${SPAN}" -lt 32800 ] || [ "${SPAN}" -gt 33800 ]; then
  echo "STOP: array span does not match 6 KIV-2 units. Check ref/unit_hits.paf."
  exit 1
fi
REMOVED=$(( LAST - FIRST - ULEN ))
echo "Array in window: ${FIRST}-${LAST} | unit length ${ULEN} | removed ${REMOVED} bp"

echo ">collapsed" > ref/collapsed.fa
{ samtools faidx ref/window.fa window:1-${FIRST} | grep -v "^>"
  grep -v "^>" ref/unit.fa
  samtools faidx ref/window.fa window:$(( LAST + 1 ))-${WIN_LEN} | grep -v "^>"
} | tr -d '\n' | fold -w 60 >> ref/collapsed.fa
echo "" >> ref/collapsed.fa
samtools faidx ref/collapsed.fa

to_collapsed () {
  local g=$1
  local w=$(( g - WIN_START ))
  if [ ${w} -lt ${FIRST} ]; then echo ${w}; else echo $(( w - REMOVED )); fi
}
{
  printf "collapsed\t%d\t%d\tKIV2_unit\n" $(( FIRST + 500 )) $(( FIRST + ULEN - 500 ))
  for R in "160320000 160330000 CTRL_upstream1" \
           "160535000 160545000 CTRL_LPA3prime" \
           "160680000 160690000 CTRL_LPA5prime" \
           "160770000 160780000 CTRL_downstream"; do
    set -- ${R}
    printf "collapsed\t%d\t%d\t%s\n" $(to_collapsed $1) $(to_collapsed $2) $3
  done
} > ref/regions.bed
echo "Depth regions (collapsed coordinates):"
cat ref/regions.bed
awk '$4=="KIV2_unit"' ref/regions.bed > ref/unit.bed
awk '$4!="KIV2_unit"' ref/regions.bed > ref/controls.bed

printf "sample\tKIV2_unit_depth\tctrl_upstream1\tctrl_LPA3prime\tctrl_LPA5prime\tctrl_downstream\tctrl_median\tDEPTH_RATIO_CN\tKILDA_CN\terr_KIV2\terr_ctrl\n" \
  > results_depth_ratio.tsv

for SAMPLE in "${SAMPLES[@]}"; do
  echo "=== ${SAMPLE} | $(date) ==="
  BAM=${HG38_ROOT}/${SAMPLE}/merged/${SAMPLE}_hg38_merged.bam

  samtools view -@ ${T} -b -F 0x900 ${BAM} ${WIN_CHR}:${WIN_START}-${WIN_END} \
    | samtools fastq -@ ${T} - > fastq/${SAMPLE}.fq
  echo "Reads extracted: $(( $(wc -l < fastq/${SAMPLE}.fq) / 4 ))"

  minimap2 -ax map-ont --secondary=no -t ${T} ref/collapsed.fa fastq/${SAMPLE}.fq \
    | samtools sort -@ ${T} -o bam/${SAMPLE}.bam
  samtools index bam/${SAMPLE}.bam

  mosdepth -t ${T} -n --by ref/regions.bed depth/${SAMPLE} bam/${SAMPLE}.bam

  samtools stats -@ ${T} -t ref/unit.bed bam/${SAMPLE}.bam > stats/${SAMPLE}_unit.stats
  samtools stats -@ ${T} -t ref/controls.bed bam/${SAMPLE}.bam > stats/${SAMPLE}_ctrl.stats
  ERR_U=$(grep "^SN" stats/${SAMPLE}_unit.stats | grep "error rate:" | cut -f3)
  ERR_C=$(grep "^SN" stats/${SAMPLE}_ctrl.stats | grep "error rate:" | cut -f3)

  KCN=$(awk -v S=${SAMPLE} '$1==S {print $2}' ${KILDA_CN})

  zcat depth/${SAMPLE}.regions.bed.gz | awk -v S=${SAMPLE} -v K=${KCN} -v EU=${ERR_U} -v EC=${ERR_C} '
    BEGIN{OFS="\t"}
    { d[$4]=$5 }
    END {
      c[1]=d["CTRL_upstream1"]; c[2]=d["CTRL_LPA3prime"]; c[3]=d["CTRL_LPA5prime"]; c[4]=d["CTRL_downstream"]
      for (i=1;i<=4;i++) for (j=i+1;j<=4;j++) if (c[j]<c[i]) {t=c[i];c[i]=c[j];c[j]=t}
      med=(c[2]+c[3])/2
      cn=(med>0)? 2*d["KIV2_unit"]/med : 0
      printf "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\n", S, d["KIV2_unit"],
        d["CTRL_upstream1"], d["CTRL_LPA3prime"], d["CTRL_LPA5prime"], d["CTRL_downstream"], med, cn, K, EU, EC
    }' >> results_depth_ratio.tsv
done

echo "=== Depth ratio vs KILDA ==="
column -t results_depth_ratio.tsv
echo "=== Done | $(date) ==="
