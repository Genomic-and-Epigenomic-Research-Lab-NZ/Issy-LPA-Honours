#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=split_longphase
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/split_longphase_%j.out
#SBATCH --error=/path/to/project/logs/split_longphase_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

HAPLOTAGGED_DIR=/path/to/project/longphase_haplotagged
OUT_DIR=/path/to/project/longphase_split_bams
mkdir -p "$OUT_DIR"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

for SAMPLE in "${SAMPLES[@]}"; do
  IN_BAM="${HAPLOTAGGED_DIR}/${SAMPLE}.LPA.longphase.bam"
  HP1_BAM="${OUT_DIR}/${SAMPLE}.LPA.longphase.HP1.bam"
  HP2_BAM="${OUT_DIR}/${SAMPLE}.LPA.longphase.HP2.bam"

  if [ -f "$HP1_BAM" ] && [ -f "$HP2_BAM" ]; then
    echo "Skipping ${SAMPLE}, already done"
    continue
  fi

  if [ ! -f "$IN_BAM" ]; then
    echo "MISSING longphase-tagged BAM for ${SAMPLE}: ${IN_BAM}" >&2
    exit 1
  fi

  set +o pipefail
  samtools quickcheck "$IN_BAM" 2>/dev/null
  QC_STATUS=$?
  set -o pipefail
  if [ "$QC_STATUS" -ne 0 ]; then
    echo "quickcheck failed for ${SAMPLE}: ${IN_BAM}" >&2
    exit 1
  fi

  echo "Splitting ${SAMPLE}..."
  samtools view -b -d HP:1 "$IN_BAM" > "$HP1_BAM"
  samtools view -b -d HP:2 "$IN_BAM" > "$HP2_BAM"

  samtools index "$HP1_BAM"
  samtools index "$HP2_BAM"

  N1=$(samtools view -c "$HP1_BAM")
  N2=$(samtools view -c "$HP2_BAM")
  echo "${SAMPLE}: HP1=${N1} reads, HP2=${N2} reads"
done

echo "Done."
