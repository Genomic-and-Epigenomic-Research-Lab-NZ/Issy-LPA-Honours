#!/bin/bash -e
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running

#SBATCH --account=ACCOUNT
#SBATCH --job-name=prep_basecall_PROJECT
#SBATCH --partition=aoraki
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/path/to/project/SLURMs/PROJECT/logs/prep_%j.out
#SBATCH --error=/path/to/project/SLURMs/PROJECT/logs/prep_%j.err

set -euo pipefail

module load samtools

USER_TAG="PROJECT"

PROJECT_BASE="/path/to/project/${USER_TAG}_sequencing"
BASE="/path/to/scratch/${USER_TAG}"

POD5_DIR="${BASE}/working_pod5s"
OUTBASE="${BASE}/array_output"

LIST="${POD5_DIR}/pod5_list_${USER_TAG}.txt"

REFERENCE="/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
SAMPLE_SHEET="/path/to/sample_sheet.csv"

MODELS="/path/to/shared/models_used"
MODEL_BASE="${MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
MODEL_5mC="${MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mC_5hmC@v2"
MODEL_6mA="${MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1"

DORADO="/path/to/shared/dorado-1.4.0-linux-x64/bin/dorado"

echo "Checking sample sheet:"
sed -n '1,20p' "$SAMPLE_SHEET" | cat -A

echo "Estimating dataset size before rsync..."

find "$PROJECT_BASE" -type f -name "*.pod5" \
    > "${BASE}/all_pod5_paths.txt"

TOTAL_FILES=$(wc -l < "${BASE}/all_pod5_paths.txt")

TOTAL_SIZE=0

while read -r f; do
    SIZE=$(stat -c%s "$f")
    TOTAL_SIZE=$((TOTAL_SIZE + SIZE))
done < "${BASE}/all_pod5_paths.txt"

TOTAL_GB=$(echo "$TOTAL_SIZE / (1024^3)" | bc)

echo "Total POD5 files: $TOTAL_FILES"
echo "Total dataset size: ${TOTAL_GB} GB"

EST_RSYNC_HOURS=$(echo "$TOTAL_GB / 540" | bc)

if [[ "$EST_RSYNC_HOURS" -lt 1 ]]; then
    EST_RSYNC_HOURS=1
fi

if [[ "$EST_RSYNC_HOURS" -gt 24 ]]; then
    EST_RSYNC_HOURS=24
fi

echo "Estimated rsync time: ${EST_RSYNC_HOURS} hours"

echo "Preparing directories..."

mkdir -p "$POD5_DIR"
mkdir -p "$OUTBASE"

echo "Running rsync (flat)..."

rsync -av \
    --files-from="${BASE}/all_pod5_paths.txt" \
    --no-relative \
    / "$POD5_DIR/"

echo "Building POD5 list..."

find "$POD5_DIR" -type f -name "*.pod5" \
    | grep -v "no_sample_id" \
    | grep -v "fail_mixed" \
    | sort \
    > "$LIST"

echo "Validating and filtering invalid POD5 files..."

CLEAN_LIST="${POD5_DIR}/pod5_list_clean_${USER_TAG}.txt"
BAD_LIST="${POD5_DIR}/bad_pod5_${USER_TAG}.txt"

> "$CLEAN_LIST"
> "$BAD_LIST"

while read -r f; do
    SIZE=$(stat -c%s "$f")

    if [[ "$SIZE" -lt 200000 ]]; then
        echo "Removing small POD5: $f"
        echo "$f" >> "$BAD_LIST"
    else
        echo "$f" >> "$CLEAN_LIST"
    fi
done < "$LIST"

mv "$CLEAN_LIST" "$LIST"
COUNT=$(wc -l < "$LIST")
echo "Total POD5 files ready: $COUNT"

GOOD_COUNT=$(wc -l < "$LIST")
BAD_COUNT=$(wc -l < "$BAD_LIST")

echo "Valid POD5 files: $GOOD_COUNT"
echo "Removed small POD5 files: $BAD_COUNT"

echo "Validating Dorado input..."

if grep -v "^$POD5_DIR/" "$LIST" | grep -q .; then
    echo "ERROR: Files outside working_pod5s detected"
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
#SBATCH --job-name=basecall_${USER_TAG}
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
    --modified-bases-models "${MODEL_5mC},${MODEL_6mA}" \
    --min-qscore 8 \
    --output-dir "\$OUTDIR" \
    --reference "$REFERENCE" \
    --sample-sheet "$SAMPLE_SHEET" \
    --mm2-opts "--secondary=no" \
    "$MODEL_BASE" \
    "\$FILE"

echo "Completed: \$FILE"
EOF
