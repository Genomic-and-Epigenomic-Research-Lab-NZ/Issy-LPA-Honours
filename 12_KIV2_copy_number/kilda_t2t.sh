#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=kilda_all
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/path/to/project/logs/kilda/%x_%j.out
#SBATCH --error=/path/to/project/logs/kilda/%x_%j.err

set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)
REGION="NC_060930.1:161500000-162300000"

KILDA=/path/to/project/KILDA
K=31
KIV2_KMERS=${KILDA}/data/kmers_GRCh38_LPA_k31/KIV2_hg38_kmers_6copies_specific.tsv
NORM_KMERS=${KILDA}/data/kmers_GRCh38_LPA_k31/LPA_hg38_kmers_1copies_specific.tsv
RSIDS=${KILDA}/data/lpa_three_rsids.tsv

OUT=/path/to/project/kilda_out_all
mkdir -p ${OUT}/fastq ${OUT}/counts
: > ${OUT}/counts.list

for SAMPLE in "${SAMPLES[@]}"; do
  echo "=== ${SAMPLE} | $(date) ==="
  BAM=/path/to/project/PROJECT_T2T/${SAMPLE}/filtered/${SAMPLE}_T2T_q10.bam

  samtools view -b -F 0x900 ${BAM} ${REGION} \
    | samtools fastq - > ${OUT}/fastq/${SAMPLE}_LPA.fastq
  echo "Reads extracted: $(( $(wc -l < ${OUT}/fastq/${SAMPLE}_LPA.fastq) / 4 ))"

  jellyfish count -t ${SLURM_CPUS_PER_TASK} -m ${K} -s 100M -C \
    -o ${OUT}/counts/${SAMPLE}_${K}mer ${OUT}/fastq/${SAMPLE}_LPA.fastq
  jellyfish dump -c -t ${OUT}/counts/${SAMPLE}_${K}mer > ${OUT}/counts/${SAMPLE}.counts

  printf "%s\t%s\n" ${SAMPLE} ${OUT}/counts/${SAMPLE}.counts >> ${OUT}/counts.list
done

python3 ${KILDA}/bin/kilda.py \
  -c ${OUT}/counts.list \
  -k ${KIV2_KMERS} \
  -l ${NORM_KMERS} \
  -r ${RSIDS} \
  -o ${OUT}/results/ \
  -p -v

echo "=== Results ==="
cat ${OUT}/results/kilda_kiv2.CNs
echo "=== Done | $(date) ==="
