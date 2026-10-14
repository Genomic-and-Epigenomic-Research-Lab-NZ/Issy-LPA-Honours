#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=qc_x2
#SBATCH --account=ACCOUNT
#SBATCH --partition=aoraki
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-3
#SBATCH --output=/path/to/project/logs/qc_x2_%A_%a.out
#SBATCH --error=/path/to/project/logs/qc_x2_%A_%a.err

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load miniconda3

eval "$(conda shell.bash hook)"

set +u
conda activate custom_env
set -u

echo "Threads available: $OMP_NUM_THREADS"

USER_TAG="PROJECT"
REF="T2T"
TAG_REF="${USER_TAG}_${REF}"
REFERENCE="/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

LPA_CHR="NC_060930.1"
LPA_START="161500000"
LPA_END="162300000"

KIV2_CHR="NC_060930.1"
KIV2_START="161862121"
KIV2_END="161998520"

KIV2_PAD_START="161812121"
KIV2_PAD_END="162048520"

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "Array task $SLURM_ARRAY_TASK_ID processing sample: $SAMPLE"

INPUT_MERGED_DIR="/path/to/shared/merged_bams"
INPUT_Q10_DIR="/path/to/shared/merged_bams/qc_reports"
INPUT_LPA_DIR="/path/to/shared/lpa_analysis"
INPUT_LPA_BAM="${INPUT_LPA_DIR}/${SAMPLE}.LPA.bam"

OUT_BASE="/path/to/project/t2t_qc"
BED_DIR="${OUT_BASE}/beds"

MERGED_BAM="${INPUT_MERGED_DIR}/${SAMPLE}.merged.bam"
Q10_BAM="${INPUT_Q10_DIR}/${SAMPLE}.q10.bam"

SAMPLE_DIR="${OUT_BASE}/${SAMPLE}"
QC_NO_Q10_DIR="${SAMPLE_DIR}/qc_no_q10"
QC_DIR="${SAMPLE_DIR}/qc"
FILTERED_DIR="${SAMPLE_DIR}/filtered"

mkdir -p "$OUT_BASE"
mkdir -p "$BED_DIR"
mkdir -p "$QC_NO_Q10_DIR"
mkdir -p "$QC_DIR"
mkdir -p "$FILTERED_DIR"

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
  if [ ! -f ${BED_DIR}/LPA_expanded.bed ]; then
    echo "Writing LPA bed file..."
    printf "%s\t%s\t%s\n" "$LPA_CHR" "$LPA_START" "$LPA_END" > ${BED_DIR}/LPA_expanded.bed
  else
    echo "LPA bed file already exists, skipping"
  fi

  if [ ! -f ${BED_DIR}/KIV2_unpadded.bed ]; then
    echo "Writing unpadded KIV-2 bed file..."
    printf "%s\t%s\t%s\n" "$KIV2_CHR" "$KIV2_START" "$KIV2_END" > ${BED_DIR}/KIV2_unpadded.bed
  else
    echo "Unpadded KIV-2 bed file already exists, skipping"
  fi

  if [ ! -f ${BED_DIR}/KIV2_padded.bed ]; then
    echo "Writing padded KIV-2 bed file..."
    printf "%s\t%s\t%s\n" "$KIV2_CHR" "$KIV2_PAD_START" "$KIV2_PAD_END" > ${BED_DIR}/KIV2_padded.bed
  else
    echo "Padded KIV-2 bed file already exists, skipping"
  fi
else
  sleep 10
fi

if [ ! -f "$MERGED_BAM" ]; then
  echo "SKIPPING $SAMPLE stage 1 - merged bam not found at $MERGED_BAM"
elif [ -f "${QC_NO_Q10_DIR}/${SAMPLE}_${REF}_flagstat.tsv" ]; then
  echo "SKIPPING $SAMPLE stage 1 - qc_no_q10 already exists"
else
  echo "Running flagstat on raw merged bam..."
  samtools flagstat -@ $OMP_NUM_THREADS -O tsv "$MERGED_BAM" > ${QC_NO_Q10_DIR}/${SAMPLE}_${REF}_flagstat.tsv

  echo "Running samtools stats on raw merged bam..."
  samtools stats -@ $OMP_NUM_THREADS -r ${REFERENCE} "$MERGED_BAM" > ${QC_NO_Q10_DIR}/${SAMPLE}_${REF}_stats.txt
  grep ^SN ${QC_NO_Q10_DIR}/${SAMPLE}_${REF}_stats.txt | cut -f 2- > ${QC_NO_Q10_DIR}/${SAMPLE}_${REF}_stats_summary.tsv

  echo "Running genome-wide mosdepth on raw merged bam..."
  mosdepth -t $OMP_NUM_THREADS ${QC_NO_Q10_DIR}/${SAMPLE}_${REF}_genome "$MERGED_BAM"
  mv ${QC_NO_Q10_DIR}/${SAMPLE}_${REF}_genome.mosdepth.summary.txt ${QC_NO_Q10_DIR}/${SAMPLE}_${REF}_genome.mosdepth.summary.tsv

  echo "Running NanoPlot on raw merged bam..."
  NanoPlot --bam "$MERGED_BAM" \
    --outdir ${QC_NO_Q10_DIR}/${SAMPLE}_${REF}_nanoplot \
    --prefix ${SAMPLE}_${REF}_ \
    --threads $OMP_NUM_THREADS

  echo "$SAMPLE stage 1 (raw merged QC) complete."
fi

if [ ! -f "$Q10_BAM" ]; then
  echo "SKIPPING $SAMPLE stage 2 - q10 bam not found at $Q10_BAM"
else
  if [ -f "${QC_DIR}/${SAMPLE}_${REF}_flagstat.tsv" ]; then
    echo "SKIPPING $SAMPLE stage 2 - qc already exists"
  else
    echo "Running flagstat on q10 bam..."
    samtools flagstat -@ $OMP_NUM_THREADS -O tsv "$Q10_BAM" > ${QC_DIR}/${SAMPLE}_${REF}_flagstat.tsv

    echo "Running samtools stats on q10 bam..."
    samtools stats -@ $OMP_NUM_THREADS -r ${REFERENCE} "$Q10_BAM" > ${QC_DIR}/${SAMPLE}_${REF}_stats.txt
    grep ^SN ${QC_DIR}/${SAMPLE}_${REF}_stats.txt | cut -f 2- > ${QC_DIR}/${SAMPLE}_${REF}_stats_summary.tsv

    echo "Running genome-wide mosdepth on q10 bam..."
    mosdepth -t $OMP_NUM_THREADS ${QC_DIR}/${SAMPLE}_${REF}_genome "$Q10_BAM"
    mv ${QC_DIR}/${SAMPLE}_${REF}_genome.mosdepth.summary.txt ${QC_DIR}/${SAMPLE}_${REF}_genome.mosdepth.summary.tsv

    echo "Running NanoPlot on q10 bam..."
    NanoPlot --bam "$Q10_BAM" \
      --outdir ${QC_DIR}/${SAMPLE}_${REF}_nanoplot \
      --prefix ${SAMPLE}_${REF}_ \
      --threads $OMP_NUM_THREADS

    echo "$SAMPLE stage 2 (q10 QC) complete."
  fi
fi

if [ ! -f "$INPUT_LPA_BAM" ]; then
  echo "SKIPPING $SAMPLE stage 3 - LPA bam not found at $INPUT_LPA_BAM"
elif [ -f "${QC_DIR}/${SAMPLE}_${REF}_LPA_flagstat.tsv" ]; then
  echo "SKIPPING $SAMPLE stage 3 - LPA qc already exists"
else
  echo "Running flagstat on LPA bam..."
  samtools flagstat -@ $OMP_NUM_THREADS -O tsv "$INPUT_LPA_BAM" > ${QC_DIR}/${SAMPLE}_${REF}_LPA_flagstat.tsv

  echo "Running stats on LPA bam..."
  samtools stats -@ $OMP_NUM_THREADS -r ${REFERENCE} "$INPUT_LPA_BAM" > ${QC_DIR}/${SAMPLE}_${REF}_LPA_stats.txt
  grep ^SN ${QC_DIR}/${SAMPLE}_${REF}_LPA_stats.txt | cut -f 2- > ${QC_DIR}/${SAMPLE}_${REF}_LPA_stats_summary.tsv

      echo "Running mosdepth on LPA bam..."
  mosdepth -t $OMP_NUM_THREADS --by ${BED_DIR}/LPA_expanded.bed ${QC_DIR}/${SAMPLE}_${REF}_LPA "$INPUT_LPA_BAM"
  mv ${QC_DIR}/${SAMPLE}_${REF}_LPA.mosdepth.summary.txt ${QC_DIR}/${SAMPLE}_${REF}_LPA.mosdepth.summary.tsv

  echo "Running NanoPlot on LPA bam..."
  NanoPlot --bam "$INPUT_LPA_BAM" \
    --outdir ${QC_DIR}/${SAMPLE}_${REF}_nanoplot_LPA \
    --prefix ${SAMPLE}_${REF}_LPA_ \
    --threads $OMP_NUM_THREADS

  echo "$SAMPLE stage 3 (LPA QC) complete."
fi

if [ ! -f "$Q10_BAM" ]; then
  echo "SKIPPING $SAMPLE stage 4 - q10 bam not found, cannot subset KIV-2"
elif [ -f "${QC_DIR}/${SAMPLE}_${REF}_KIV2_flagstat.tsv" ]; then
  echo "SKIPPING $SAMPLE stage 4 - KIV-2 unpadded qc already exists"
else
  echo "Subsetting $SAMPLE to unpadded KIV-2 region..."
  samtools view -@ $OMP_NUM_THREADS -b -L ${BED_DIR}/KIV2_unpadded.bed "$Q10_BAM" > ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2.bam
  samtools index -@ $OMP_NUM_THREADS ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2.bam

  echo "Running flagstat on unpadded KIV-2..."
  samtools flagstat -@ $OMP_NUM_THREADS -O tsv ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2.bam > ${QC_DIR}/${SAMPLE}_${REF}_KIV2_flagstat.tsv

  echo "Running stats on unpadded KIV-2..."
  samtools stats -@ $OMP_NUM_THREADS -r ${REFERENCE} ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2.bam > ${QC_DIR}/${SAMPLE}_${REF}_KIV2_stats.txt
  grep ^SN ${QC_DIR}/${SAMPLE}_${REF}_KIV2_stats.txt | cut -f 2- > ${QC_DIR}/${SAMPLE}_${REF}_KIV2_stats_summary.tsv

  echo "Running mosdepth on unpadded KIV-2 (region-restricted)..."
  mosdepth -t $OMP_NUM_THREADS --by ${BED_DIR}/KIV2_unpadded.bed ${QC_DIR}/${SAMPLE}_${REF}_KIV2 ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2.bam
  mv ${QC_DIR}/${SAMPLE}_${REF}_KIV2.mosdepth.summary.txt ${QC_DIR}/${SAMPLE}_${REF}_KIV2.mosdepth.summary.tsv

  echo "Running NanoPlot on unpadded KIV-2..."
  NanoPlot --bam ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2.bam \
    --outdir ${QC_DIR}/${SAMPLE}_${REF}_nanoplot_KIV2 \
    --prefix ${SAMPLE}_${REF}_KIV2_ \
    --threads $OMP_NUM_THREADS

  echo "$SAMPLE stage 4 (KIV-2 unpadded QC) complete."
fi

if [ ! -f "$Q10_BAM" ]; then
  echo "SKIPPING $SAMPLE stage 5 - q10 bam not found, cannot subset KIV-2 padded"
elif [ -f "${QC_DIR}/${SAMPLE}_${REF}_KIV2padded_flagstat.tsv" ]; then
  echo "SKIPPING $SAMPLE stage 5 - KIV-2 padded qc already exists"
else
  echo "Subsetting $SAMPLE to padded KIV-2 region..."
  samtools view -@ $OMP_NUM_THREADS -b -L ${BED_DIR}/KIV2_padded.bed "$Q10_BAM" > ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2padded.bam
  samtools index -@ $OMP_NUM_THREADS ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2padded.bam

  echo "Running flagstat on padded KIV-2..."
  samtools flagstat -@ $OMP_NUM_THREADS -O tsv ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2padded.bam > ${QC_DIR}/${SAMPLE}_${REF}_KIV2padded_flagstat.tsv

  echo "Running stats on padded KIV-2..."
  samtools stats -@ $OMP_NUM_THREADS -r ${REFERENCE} ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2padded.bam > ${QC_DIR}/${SAMPLE}_${REF}_KIV2padded_stats.txt
  grep ^SN ${QC_DIR}/${SAMPLE}_${REF}_KIV2padded_stats.txt | cut -f 2- > ${QC_DIR}/${SAMPLE}_${REF}_KIV2padded_stats_summary.tsv

  echo "Running mosdepth on padded KIV-2 (region-restricted)..."
  mosdepth -t $OMP_NUM_THREADS --by ${BED_DIR}/KIV2_padded.bed ${QC_DIR}/${SAMPLE}_${REF}_KIV2padded ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2padded.bam
  mv ${QC_DIR}/${SAMPLE}_${REF}_KIV2padded.mosdepth.summary.txt ${QC_DIR}/${SAMPLE}_${REF}_KIV2padded.mosdepth.summary.tsv

  echo "Running NanoPlot on padded KIV-2..."
  NanoPlot --bam ${FILTERED_DIR}/${SAMPLE}_${REF}_KIV2padded.bam \
    --outdir ${QC_DIR}/${SAMPLE}_${REF}_nanoplot_KIV2padded \
    --prefix ${SAMPLE}_${REF}_KIV2padded_ \
    --threads $OMP_NUM_THREADS

  echo "$SAMPLE stage 5 (KIV-2 padded QC) complete."
fi

echo "$SAMPLE all stages complete."
