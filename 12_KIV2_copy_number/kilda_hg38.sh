#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=kilda_hg38
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
HG38_ROOT=/path/to/project/PROJECT_hg38
REGION="chr6:160300000-160800000"

KILDA=/path/to/project/KILDA
K=31
KIV2_KMERS=${KILDA}/data/kmers_GRCh38_LPA_k31/KIV2_hg38_kmers_6copies_specific.tsv
NORM_KMERS=${KILDA}/data/kmers_GRCh38_LPA_k31/LPA_hg38_kmers_1copies_specific.tsv
RSIDS=${KILDA}/data/lpa_three_rsids.tsv

OUT=/path/to/project/kilda_out_hg38
mkdir -p ${OUT}/fastq ${OUT}/counts
: > ${OUT}/counts.list

for SAMPLE in "${SAMPLES[@]}"; do
  echo "=== ${SAMPLE} | $(date) ==="
  BAM=${HG38_ROOT}/${SAMPLE}/merged/${SAMPLE}_hg38_merged.bam

  samtools view -b -F 0x900 ${BAM} ${REGION} \
    | samtools fastq - > ${OUT}/fastq/${SAMPLE}_LPA_hg38.fastq
  echo "Reads extracted: $(( $(wc -l < ${OUT}/fastq/${SAMPLE}_LPA_hg38.fastq) / 4 ))"

  jellyfish count -t ${SLURM_CPUS_PER_TASK} -m ${K} -s 100M -C \
    -o ${OUT}/counts/${SAMPLE}_${K}mer ${OUT}/fastq/${SAMPLE}_LPA_hg38.fastq
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

echo "=== hg38 results (T2T window: 36.41, 36.40, 32.57, 32.53) ==="
cat ${OUT}/results/kilda_kiv2.CNs
echo "=== Done | $(date) ==="
