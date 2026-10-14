#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=assembly_kiv2_check
#SBATCH --partition=aoraki
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/assembly_kiv2_check_%j.out
#SBATCH --error=/path/to/project/logs/assembly_kiv2_check_%j.err

set -euo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

REF_INDEX=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.w10_k15.mmi
ASSEMBLY_DIR=/path/to/project/flye_assembly/assemblies
OUT_DIR=/path/to/project/flye_assembly/assembly_vs_t2t
mkdir -p "$OUT_DIR"

CHR="NC_060930.1"
KIV2_START=161862121
KIV2_END=161998520

SAMPLES=(SAMPLE_1 SAMPLE_2 SAMPLE_3 SAMPLE_4)
VARIANTS=(HP1 HP2 unassigned_only HP1_plus_unassigned HP2_plus_unassigned)

SUMMARY="${OUT_DIR}/kiv2_overlap_summary.txt"
> "$SUMMARY"

for SAMPLE in "${SAMPLES[@]}"; do
  for VARIANT in "${VARIANTS[@]}"; do
    ASSEMBLY="${ASSEMBLY_DIR}/${SAMPLE}_${VARIANT}/assembly.fasta"
    SAM_OUT="${OUT_DIR}/${SAMPLE}_${VARIANT}_vs_T2T.sam"

    if [ ! -f "$ASSEMBLY" ]; then
      echo "[${SAMPLE}_${VARIANT}] Assembly not found, skipping" >> "$SUMMARY"
      continue
    fi

    minimap2 -ax map-ont "$REF_INDEX" "$ASSEMBLY" > "$SAM_OUT" 2>/dev/null

    echo "=== ${SAMPLE} / ${VARIANT} ===" >> "$SUMMARY"

    samtools view "$SAM_OUT" | grep -v "^@" | awk -v chr="$CHR" -v s="$KIV2_START" -v e="$KIV2_END" '
      {
        contig=$1; rname=$3; pos=$4; seqlen=length($10)
        end=pos+seqlen
        overlap = (rname==chr && pos<e && end>s) ? "YES - OVERLAPS KIV-2" : "no"
        print "  contig="contig, "chrom="rname, "pos="pos, "length="seqlen, "KIV2_overlap="overlap
      }' >> "$SUMMARY"

    echo "" >> "$SUMMARY"
  done
done

echo "Done. Summary at $SUMMARY"
