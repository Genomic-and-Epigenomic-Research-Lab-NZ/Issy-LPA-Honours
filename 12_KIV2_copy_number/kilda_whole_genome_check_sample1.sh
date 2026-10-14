#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=kilda_wgs
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/path/to/project/logs/kilda/%x_%j.out
#SBATCH --error=/path/to/project/logs/kilda/%x_%j.err

set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
T=$SLURM_CPUS_PER_TASK

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

SAMPLE=SAMPLE_1
BAM=/path/to/project/PROJECT_T2T/${SAMPLE}/filtered/${SAMPLE}_T2T_q10.bam

KILDA=/path/to/project/KILDA
K=31
KIV2_KMERS=${KILDA}/data/kmers_GRCh38_LPA_k31/KIV2_hg38_kmers_6copies_specific.tsv
NORM_KMERS=${KILDA}/data/kmers_GRCh38_LPA_k31/LPA_hg38_kmers_1copies_specific.tsv
RSIDS=${KILDA}/data/lpa_three_rsids.tsv

OUT=/path/to/project/kilda_out_wgs
mkdir -p ${OUT}/counts

echo "=== ${SAMPLE} whole genome | $(date) ==="
echo "Primary reads in BAM:  $(samtools view -@ ${T} -c -F 0x900 ${BAM})"
echo "Unmapped reads in BAM: $(samtools view -@ ${T} -c -f 4 ${BAM})"

awk 'NF {print ">k" NR; print $1}' ${KIV2_KMERS} ${NORM_KMERS} > ${OUT}/kilda_kmers.fa
awk 'NF {print ">r" NR "ref"; print $2; print ">r" NR "alt"; print $3}' ${RSIDS} >> ${OUT}/kilda_kmers.fa
echo "Query k-mers: $(grep -c '^>' ${OUT}/kilda_kmers.fa)"

jellyfish count -t ${T} -m ${K} -s 10M -C \
  --if ${OUT}/kilda_kmers.fa \
  -o ${OUT}/counts/${SAMPLE}_wgs_${K}mer \
  <(samtools fastq -@ 4 -F 0x900 ${BAM} 2>/dev/null)
jellyfish dump -c -t ${OUT}/counts/${SAMPLE}_wgs_${K}mer > ${OUT}/counts/${SAMPLE}_wgs.counts

printf "%s_WGS\t%s\n" ${SAMPLE} ${OUT}/counts/${SAMPLE}_wgs.counts > ${OUT}/counts.list

python3 ${KILDA}/bin/kilda.py \
  -c ${OUT}/counts.list \
  -k ${KIV2_KMERS} \
  -l ${NORM_KMERS} \
  -r ${RSIDS} \
  -o ${OUT}/results/ \
  -p -v

echo "=== Whole-genome result (window result was 36.41) ==="
cat ${OUT}/results/kilda_kiv2.CNs
echo "=== Done | $(date) ==="
