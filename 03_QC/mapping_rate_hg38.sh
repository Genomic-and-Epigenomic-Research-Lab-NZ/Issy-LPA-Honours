#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=hg38_flagstat
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/path/to/project/logs/%x_%j.out

eval "$(conda shell.bash hook)"
conda activate custom_env
for S in SAMPLE_3 SAMPLE_1 SAMPLE_2 SAMPLE_4; do
  echo "=== ${S} ==="
  samtools flagstat -@ 8 /path/to/project/PROJECT_hg38/${S}/merged/${S}_hg38_merged.bam
done
