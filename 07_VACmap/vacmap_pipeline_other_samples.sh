#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=vacmap_full_pipeline
#SBATCH --partition=aoraki
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/path/to/project/logs/vacmap_full_pipeline_%j.out
#SBATCH --error=/path/to/project/logs/vacmap_full_pipeline_%j.err

set -uo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"

LPA_DIR=/path/to/project/mini_bams
REF=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
BEDS_DIR=/path/to/project/t2t_qc/beds
LPA_BED="${BEDS_DIR}/LPA_expanded.bed"
KIV2_BED="${BEDS_DIR}/KIV2_unpadded.bed"

VACMAP_OUT=/path/to/project/vacmap_out
VACMAP_PHASED=/path/to/project/vacmap_longshot_phased
VACMAP_HAPLOTAGGED=/path/to/project/vacmap_haplotagged
VACMAP_SPLIT=/path/to/project/vacmap_split_bams
VACMAP_MOSDEPTH=/path/to/project/vacmap_mosdepth
mkdir -p "$VACMAP_OUT" "$VACMAP_PHASED" "$VACMAP_HAPLOTAGGED" "$VACMAP_SPLIT" "$VACMAP_MOSDEPTH"

CHR="NC_060930.1"
LPA_START=161500000
LPA_END=162300000
KIV2_START=161862121
KIV2_END=161998520
KIV2_PAD_START=161812121
KIV2_PAD_END=162048520

REGIONS_BED="${VACMAP_MOSDEPTH}/regions_named.bed"
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "LPA"}' "$LPA_BED" > "$REGIONS_BED"
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "KIV2"}' "$KIV2_BED" >> "$REGIONS_BED"

get_region_args () {
  case "$1" in
    LPA)           echo "${CHR}:${LPA_START}-${LPA_END}" ;;
    KIV2_unpadded) echo "${CHR}:${KIV2_START}-${KIV2_END}" ;;
    KIV2_padded)   echo "${CHR}:${KIV2_PAD_START}-${KIV2_PAD_END}" ;;
    FLANK)         echo "${CHR}:${KIV2_PAD_START}-${KIV2_START} ${CHR}:${KIV2_END}-${KIV2_PAD_END}" ;;
  esac
}

qc_block () {
  local bam="$1" region_args="$2" label="$3" out="$4"
  echo "=== ${label} ===" >> "$out"
  local total sa_count softclip_count mapq_tmp
  total=$(samtools view -c "$bam" $region_args)
  echo "total_alignments: $total" >> "$out"
  sa_count=$(samtools view "$bam" $region_args | awk '
    { for(i=12;i<=NF;i++) if($i ~ /^SA:Z:/) c++ } END{print c+0}')
  echo "sa_count: $sa_count" >> "$out"
  softclip_count=$(samtools view "$bam" $region_args | awk '$6 ~ /S/' | wc -l)
  echo "softclip_count: $softclip_count" >> "$out"
  mapq_tmp="${out}.mapq.tmp"
  samtools view "$bam" $region_args | awk '{print $5}' | sort -n > "$mapq_tmp"
  awk '{a[NR]=$1; sum+=$1} END{
    n=NR
    if(n==0){print "mapq: no reads"; exit}
    mean=sum/n
    if(n%2==1) med=a[int((n+1)/2)]; else med=(a[n/2]+a[n/2+1])/2
    q1=a[int(n*0.25)+1]; q3=a[int(n*0.75)+1]
    printf "mapq_mean: %.2f\nmapq_median: %.2f\nmapq_q1: %.2f\nmapq_q3: %.2f\n", mean, med, q1, q3
  }' "$mapq_tmp" >> "$out"
  rm -f "$mapq_tmp"
  if [ "$total" -gt 0 ]; then
    echo "sa_pct: $(awk -v a="$sa_count" -v b="$total" 'BEGIN{printf "%.2f", (a/b)*100}')" >> "$out"
    echo "softclip_pct: $(awk -v a="$softclip_count" -v b="$total" 'BEGIN{printf "%.2f", (a/b)*100}')" >> "$out"
  fi
  echo "" >> "$out"
}

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_4)
ALIGN_QC_SUMMARY="${VACMAP_OUT}/vacmap_alignment_qc_summary.txt"
HAP_QC_SUMMARY="${VACMAP_MOSDEPTH}/vacmap_haplotype_qc_summary.txt"
MOSDEPTH_SUMMARY="${VACMAP_MOSDEPTH}/vacmap_mosdepth_summary.txt"
PIPELINE_LOG="${VACMAP_OUT}/vacmap_pipeline_run_log.txt"
> "$PIPELINE_LOG"

for SAMPLE in "${SAMPLES[@]}"; do
  echo "################ ${SAMPLE} ################" | tee -a "$PIPELINE_LOG"

  conda activate vacmap_env
  MINI_BAM="${LPA_DIR}/${SAMPLE}.LPA.bam"
  FASTQ="${VACMAP_OUT}/${SAMPLE}.LPA.fastq"

  if [ ! -f "$FASTQ" ]; then
    echo "[${SAMPLE}] Extracting FASTQ..." | tee -a "$PIPELINE_LOG"
    samtools quickcheck "$MINI_BAM" 2>/dev/null
    if [ $? -ne 0 ]; then
      echo "[${SAMPLE}] ERROR: quickcheck failed on ${MINI_BAM}, skipping sample." | tee -a "$PIPELINE_LOG"
      continue
    fi
    samtools fastq -F 0x900 "$MINI_BAM" > "$FASTQ"
  fi

  VACMAP_BAM="${VACMAP_OUT}/${SAMPLE}_vacmap_piped.sorted.bam"

  if [ ! -f "$VACMAP_BAM" ]; then
    echo "[${SAMPLE}] Running VACmap..." | tee -a "$PIPELINE_LOG"
    vacmap -ref "$REF" -read "$FASTQ" -mode H -t "$SLURM_CPUS_PER_TASK" --force \
      | samtools sort -@4 > "$VACMAP_BAM"
    if [ ! -s "$VACMAP_BAM" ]; then
      echo "[${SAMPLE}] ERROR: VACmap output empty/missing, skipping sample." | tee -a "$PIPELINE_LOG"
      continue
    fi
  fi
  samtools index "$VACMAP_BAM"

  echo "[${SAMPLE}] Running alignment QC..." | tee -a "$PIPELINE_LOG"
  for REGION in LPA KIV2_unpadded KIV2_padded FLANK; do
    REGION_ARGS=$(get_region_args "$REGION")
    qc_block "$VACMAP_BAM" "$REGION_ARGS" "${SAMPLE} / VACmap / ${REGION}" "$ALIGN_QC_SUMMARY"
  done

  conda activate custom_env
  PHASED_VCF="${VACMAP_PHASED}/${SAMPLE}.vacmap.LPA.phased.vcf"

  if [ ! -f "$PHASED_VCF" ]; then
    echo "[${SAMPLE}] Running Longshot..." | tee -a "$PIPELINE_LOG"
    longshot --bam "$VACMAP_BAM" --ref "$REF" --out "$PHASED_VCF"
    if [ ! -s "$PHASED_VCF" ]; then
      echo "[${SAMPLE}] ERROR: Longshot output empty/missing, skipping sample." | tee -a "$PIPELINE_LOG"
      continue
    fi
  fi

  VCF_GZ="${PHASED_VCF}.gz"
  if [ ! -f "$VCF_GZ" ]; then
    bgzip -c "$PHASED_VCF" > "$VCF_GZ"
  fi
  if [ ! -f "${VCF_GZ}.tbi" ]; then
    tabix -p vcf "$VCF_GZ"
  fi
  if [ ! -f "$VCF_GZ" ] || [ ! -f "${VCF_GZ}.tbi" ]; then
    echo "[${SAMPLE}] ERROR: bgzip/tabix failed to produce a valid index. Skipping haplotagging for this sample -- do this step manually and rerun from here." | tee -a "$PIPELINE_LOG"
    continue
  fi

  WHATSHAP_BAM="${VACMAP_HAPLOTAGGED}/${SAMPLE}.vacmap.whatshap.bam"
  LONGPHASE_PREFIX="${VACMAP_HAPLOTAGGED}/${SAMPLE}.vacmap.longphase"
  LONGPHASE_BAM="${LONGPHASE_PREFIX}.bam"

  if [ ! -f "$WHATSHAP_BAM" ]; then
    echo "[${SAMPLE}] Haplotagging with whatshap..." | tee -a "$PIPELINE_LOG"
    whatshap haplotag --reference "$REF" --ignore-read-groups \
      -o "$WHATSHAP_BAM" "$VCF_GZ" "$VACMAP_BAM"
    if [ -s "$WHATSHAP_BAM" ]; then samtools index "$WHATSHAP_BAM"; fi
  fi

  if [ ! -f "$LONGPHASE_BAM" ]; then
    echo "[${SAMPLE}] Haplotagging with longphase..." | tee -a "$PIPELINE_LOG"
    longphase haplotag -s "$VCF_GZ" -b "$VACMAP_BAM" -r "$REF" \
      -o "$LONGPHASE_PREFIX" -q 0 --tagSupplementary -t "$SLURM_CPUS_PER_TASK"
    if [ -s "$LONGPHASE_BAM" ]; then samtools index "$LONGPHASE_BAM"; fi
  fi

  TOOLS=(whatshap longphase)
  for TOOL in "${TOOLS[@]}"; do
    IN_BAM="${VACMAP_HAPLOTAGGED}/${SAMPLE}.vacmap.${TOOL}.bam"
    HP1_BAM="${VACMAP_SPLIT}/${SAMPLE}.vacmap.${TOOL}.HP1.bam"
    HP2_BAM="${VACMAP_SPLIT}/${SAMPLE}.vacmap.${TOOL}.HP2.bam"

    if [ ! -f "$IN_BAM" ]; then
      echo "[${SAMPLE}] WARNING: ${TOOL} haplotagged BAM missing, skipping split for this tool." | tee -a "$PIPELINE_LOG"
      continue
    fi

    if [ ! -f "$HP1_BAM" ] || [ ! -f "$HP2_BAM" ]; then
      echo "[${SAMPLE}] Splitting ${TOOL}..." | tee -a "$PIPELINE_LOG"
      samtools view -b -d HP:1 "$IN_BAM" > "$HP1_BAM"
      samtools view -b -d HP:2 "$IN_BAM" > "$HP2_BAM"
      samtools index "$HP1_BAM"
      samtools index "$HP2_BAM"
    fi

    N1=$(samtools view -c "$HP1_BAM" 2>/dev/null || echo 0)
    N2=$(samtools view -c "$HP2_BAM" 2>/dev/null || echo 0)
    echo "${SAMPLE} ${TOOL}: HP1=${N1} reads, HP2=${N2} reads" >> "$MOSDEPTH_SUMMARY"

    for HAP in HP1 HP2; do
      BAM="${VACMAP_SPLIT}/${SAMPLE}.vacmap.${TOOL}.${HAP}.bam"
      PREFIX="${VACMAP_MOSDEPTH}/${SAMPLE}.${TOOL}.${HAP}"
      if [ ! -f "${PREFIX}.regions.bed.gz" ]; then
        mosdepth --no-per-base --by "$REGIONS_BED" "$PREFIX" "$BAM"
      fi
      echo "=== ${SAMPLE} / ${TOOL} / ${HAP} ===" >> "$MOSDEPTH_SUMMARY"
      zcat "${PREFIX}.regions.bed.gz" >> "$MOSDEPTH_SUMMARY" 2>/dev/null || echo "mosdepth output missing" >> "$MOSDEPTH_SUMMARY"

      for REGION in LPA KIV2_unpadded; do
        REGION_ARGS=$(get_region_args "$REGION")
        qc_block "$BAM" "$REGION_ARGS" "${SAMPLE} / ${TOOL} / ${HAP} / ${REGION}" "$HAP_QC_SUMMARY"
      done
    done
  done

  echo "[${SAMPLE}] Done." | tee -a "$PIPELINE_LOG"
  echo "" | tee -a "$PIPELINE_LOG"
done

echo "ALL SAMPLES PROCESSED." | tee -a "$PIPELINE_LOG"
echo "Alignment QC:  $ALIGN_QC_SUMMARY"
echo "Mosdepth:      $MOSDEPTH_SUMMARY"
echo "Haplotype QC:  $HAP_QC_SUMMARY"
echo "Run log:       $PIPELINE_LOG"
