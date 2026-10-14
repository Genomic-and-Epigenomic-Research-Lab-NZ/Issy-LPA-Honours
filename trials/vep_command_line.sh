#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=vep_annotation
#SBATCH --partition=aoraki
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/vep_annotation_%j.out
#SBATCH --error=/path/to/project/logs/vep_annotation_%j.err

set -uo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

VEP_DIR=/path/to/project/vep
CACHE_DIR="${VEP_DIR}/cache"
CLINVAR_DIR="${VEP_DIR}/clinvar"
mkdir -p "$CACHE_DIR" "$CLINVAR_DIR"

if ! command -v vep &> /dev/null; then
  echo "ERROR: vep not found. Install manually first with:" >&2
  echo "  conda install -c bioconda -c conda-forge ensembl-vep" >&2
  exit 1
fi

if [ ! -d "${CACHE_DIR}/homo_sapiens" ]; then
  echo "Downloading VEP GRCh38 cache (this is the slow step, likely several GB)..."
  vep_install -a cf -s homo_sapiens -y GRCh38 -c "$CACHE_DIR" --NO_HTSLIB
fi

CLINVAR_VCF="${CLINVAR_DIR}/clinvar_GRCh38.vcf.gz"
if [ ! -f "$CLINVAR_VCF" ]; then
  echo "Downloading ClinVar GRCh38 VCF..."
  wget -O "$CLINVAR_VCF" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
  wget -O "${CLINVAR_VCF}.tbi" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
fi

IN_DIR=/path/to/project/gatk_hg38
OUT_DIR=/path/to/project/vep_annotated
mkdir -p "$OUT_DIR"

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)

for SAMPLE in "${SAMPLES[@]}"; do
  IN_VCF="${IN_DIR}/${SAMPLE}_hg38_LPA.vcf.gz"
  OUT_VCF="${OUT_DIR}/${SAMPLE}_hg38_LPA_annotated.vcf"
  OUT_SUMMARY="${OUT_DIR}/${SAMPLE}_hg38_LPA_annotated_summary.html"

  if [ ! -f "$IN_VCF" ]; then
    echo "[${SAMPLE}] MISSING: ${IN_VCF}, skipping" >&2
    continue
  fi

  if [ -f "$OUT_VCF" ]; then
    echo "[${SAMPLE}] Already annotated, skipping"
    continue
  fi

  echo "[${SAMPLE}] Running VEP..."
  vep \
    --input_file "$IN_VCF" \
    --output_file "$OUT_VCF" \
    --stats_file "$OUT_SUMMARY" \
    --cache \
    --offline \
    --dir_cache "$CACHE_DIR" \
    --assembly GRCh38 \
    --fasta /path/to/refs/hg38.primary_assembly.fa \
    --custom file="${CLINVAR_VCF}",short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
    --af_gnomade \
    --af_gnomadg \
    --everything \
    --vcf \
    --force_overwrite \
    --fork "$SLURM_CPUS_PER_TASK"

  echo "[${SAMPLE}] Done."
done

echo "ALL SAMPLES PROCESSED."
