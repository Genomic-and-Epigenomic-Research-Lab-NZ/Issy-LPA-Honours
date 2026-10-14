#!/bin/bash -e
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running

#SBATCH --account=ACCOUNT
#SBATCH --job-name=prep_basecall_PROJECT_hg38
#SBATCH --partition=aoraki
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/path/to/project/SLURMs/PROJECT/logs/prep_%j.out
#SBATCH --error=/path/to/project/SLURMs/PROJECT/logs/prep_%j.err

set -euo pipefail

module load samtools

USER_TAG="PROJECT"

POD5_DIR="/path/to/shared/pod5s"

SCRATCH_BASE="/path/to/scratch"
OUTBASE="${SCRATCH_BASE}/array_output"

BASE="/path/to/project"

LIST="${BASE}/pod5_list_${USER_TAG}_hg38.txt"

REFERENCE="/path/to/refs/hg38.primary_assembly.fa"
SAMPLE_SHEET="/path/to/sample_sheet.csv"

MODELS="/path/to/shared/models_used"
MODEL_BASE="${MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.2.0"

DORADO="/path/to/shared/dorado-1.4.0-linux-x64/bin/dorado"

echo "Checking sample sheet:"
sed -n '1,20p' "$SAMPLE_SHEET" | cat -A

echo "Preparing directories..."

mkdir -p "$OUTBASE"

echo "Building POD5 list..."

find "$POD5_DIR" -type f -name "*.pod5" \
    | grep -v "no_sample_id" \
    | grep -v "fail_mixed" \
    | sort \
    > "$LIST"

echo "Validating POD5 file sizes..."

CLEAN_LIST="${BASE}/pod5_list_clean_${USER_TAG}_hg38.txt"
BAD_LIST="${BASE}/bad_pod5_${USER_TAG}_hg38.txt"

> "$CLEAN_LIST"
> "$BAD_LIST"

TOTAL_SIZE=0

while read -r f; do
    SIZE=$(stat -c%s "$f")
    TOTAL_SIZE=$((TOTAL_SIZE + SIZE))

    if [[ "$SIZE" -lt 200000 ]]; then
        echo "Flagging small POD5: $f"
        echo "$f" >> "$BAD_LIST"
    else
        echo "$f" >> "$CLEAN_LIST"
    fi
done < "$LIST"

mv "$CLEAN_LIST" "$LIST"
COUNT=$(wc -l < "$LIST")
BAD_COUNT=$(wc -l < "$BAD_LIST")

TOTAL_GB=$(echo "$TOTAL_SIZE / (1024^3)" | bc)

echo "Total dataset size: ${TOTAL_GB} GB"
echo "Valid POD5 files: $COUNT"
echo "Flagged (too small) POD5 files: $BAD_COUNT"

echo "Validating Dorado input..."

if grep -v "^$POD5_DIR/" "$LIST" | grep -q .; then
    echo "ERROR: Files outside the shared pod5 directory detected"
    exit 1
fi

EST_GPU_HOURS=$(echo "$TOTAL_GB / 10" | bc)

if [[ "$EST_GPU_HOURS" -lt 4 ]]; then
    EST_GPU_HOURS=4
fi

if [[ "$EST_GPU_HOURS" -gt 120 ]]; then
    EST_GPU_HOURS=120
fi

echo "Estimated basecalling time: ${EST_GPU_HOURS} hours"

echo "Submitting basecalling array job..."

sbatch <<EOF
#!/bin/bash -e
#SBATCH --account=ACCOUNT
#SBATCH --job-name=basecall_${USER_TAG}_hg38
#SBATCH --partition=aoraki_gpu_A100_80GB
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=${EST_GPU_HOURS}:00:00
#SBATCH --array=1-${COUNT}%20
#SBATCH --output=/path/to/project/SLURMs/${USER_TAG}/logs/%A_%a.out
#SBATCH --error=/path/to/project/SLURMs/${USER_TAG}/logs/%A_%a.err

module load cuda

FILE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [[ -z "\${FILE:-}" ]]; then
    echo "No file for task ID \$SLURM_ARRAY_TASK_ID"
    exit 0
fi

BASENAME=\$(basename "\$FILE" .pod5)
OUTDIR="${OUTBASE}/\${BASENAME}"

mkdir -p "\$OUTDIR"

echo "Processing: \$FILE"

"$DORADO" basecaller \
    --verbose \
    --models-directory "$MODELS" \
    --min-qscore 8 \
    --output-dir "\$OUTDIR" \
    --reference "$REFERENCE" \
    --sample-sheet "$SAMPLE_SHEET" \
    --mm2-opts "--secondary=no" \
    "$MODEL_BASE" \
    "\$FILE"

echo "Completed: \$FILE"
EOF
