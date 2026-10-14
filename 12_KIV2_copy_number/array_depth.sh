#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=array_depth
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/path/to/project/logs/depth_ratio/%x_%j.out
#SBATCH --error=/path/to/project/logs/depth_ratio/%x_%j.err

set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
T=$SLURM_CPUS_PER_TASK

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

HG38_ROOT=/path/to/project/PROJECT_hg38
SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)
KILDA_CN=/path/to/project/kilda_out_hg38/results/kilda_kiv2.CNs
DR=/path/to/project/depth_ratio
OUT=/path/to/project/array_depth
mkdir -p ${OUT}/bam ${OUT}/depth
cd ${OUT}

cat > regions.bed <<EOF
chr6	160613491	160646773	KIV2_array
chr6	160320000	160330000	CTRL_upstream1
chr6	160535000	160545000	CTRL_LPA3prime
chr6	160680000	160690000	CTRL_LPA5prime
chr6	160770000	160780000	CTRL_downstream
EOF

printf "sample\tarray_depth\tctrl_median\tARRAY_DEPTH_CN\tKILDA_hg38_CN\tcollapsed_softclip_pct\n" > results_array_depth.tsv

for SAMPLE in "${SAMPLES[@]}"; do
  echo "=== ${SAMPLE} | $(date) ==="
  BAM=${HG38_ROOT}/${SAMPLE}/merged/${SAMPLE}_hg38_merged.bam

  samtools view -@ ${T} -b ${BAM} chr6:160300000-160800000 > bam/${SAMPLE}_win.bam
  samtools index bam/${SAMPLE}_win.bam
  mosdepth -t ${T} -n --by regions.bed depth/${SAMPLE} bam/${SAMPLE}_win.bam

  SC=$(samtools view ${DR}/bam/${SAMPLE}.bam collapsed:313992-318535 | awk '
    { c=$6; while (match(c, /[0-9]+[MIDNSHP=X]/)) { n=substr(c,RSTART,RLENGTH-1); op=substr(c,RSTART+RLENGTH-1,1);
        if (op=="S") s+=n; if (op=="M"||op=="I"||op=="S"||op=="="||op=="X") q+=n; c=substr(c,RSTART+RLENGTH) } }
    END { printf "%.1f", (q>0)? 100*s/q : 0 }')

  KCN=$(awk -v S=${SAMPLE} '$1==S {print $2}' ${KILDA_CN})

  zcat depth/${SAMPLE}.regions.bed.gz | awk -v S=${SAMPLE} -v K=${KCN} -v SC=${SC} '
    { d[$4]=$5 }
    END {
      c[1]=d["CTRL_upstream1"]; c[2]=d["CTRL_LPA3prime"]; c[3]=d["CTRL_LPA5prime"]; c[4]=d["CTRL_downstream"]
      for (i=1;i<=4;i++) for (j=i+1;j<=4;j++) if (c[j]<c[i]) {t=c[i];c[i]=c[j];c[j]=t}
      med=(c[2]+c[3])/2
      cn=(med>0)? 2*6*d["KIV2_array"]/med : 0
      printf "%s\t%.2f\t%.2f\t%.2f\t%s\t%s\n", S, d["KIV2_array"], med, cn, K, SC
    }' >> results_array_depth.tsv
done

echo "=== Array depth vs KILDA ==="
column -t results_array_depth.tsv
