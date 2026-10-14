#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=methylation_plots
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/methylation_plots_%j.out
#SBATCH --error=/path/to/project/logs/methylation_plots_%j.err

set -uo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

if ! python3 -c "import matplotlib" 2>/dev/null; then
  echo "Installing matplotlib..."
  pip install matplotlib --quiet
fi

BED_DIR=/path/to/project/modkit_pileup
OUT_DIR=/path/to/project/methylation_plots
mkdir -p "$OUT_DIR"

SCRIPT_PATH="${OUT_DIR}/plot_methylation.py"

cat > "$SCRIPT_PATH" << 'PYEOF'
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LPA_START = 161500000
LPA_END = 162300000
KIV2_START = 161862121
KIV2_END = 161998520
BIN_SIZE = 2000

samples = ["SAMPLE_1", "SAMPLE_2", "SAMPLE_3", "SAMPLE_4"]
bed_dir = sys.argv[1]
out_dir = sys.argv[2]

def load_5mc_positions(bed_path):
    """Return (positions, pct_values) for 5mC ('m') calls in the LPA window."""
    positions, values = [], []
    with open(bed_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            chrom, start = fields[0], int(fields[1])
            code = fields[3]
            pct = float(fields[10])
            if chrom != "NC_060930.1" or code != "m":
                continue
            if LPA_START <= start <= LPA_END:
                positions.append(start)
                values.append(pct)
    return np.array(positions), np.array(values)

fig, axes = plt.subplots(2, 2, figsize=(14, 8))
axes = axes.flatten()

for i, sample in enumerate(samples):
    bed_path = f"{bed_dir}/{sample}.modkit_pileup.bed"
    try:
        pos, val = load_5mc_positions(bed_path)
    except FileNotFoundError:
        axes[i].set_title(f"{sample} (file not found)")
        continue

    ax = axes[i]
    ax.scatter(pos, val, s=2, alpha=0.4, color="#6a3d9a")
    ax.axvspan(KIV2_START, KIV2_END, color="orange", alpha=0.15, label="KIV-2")
    ax.set_title(sample, fontsize=11, fontweight="bold")
    ax.set_xlabel("Position (T2T-CHM13v2.0, NC_060930.1)", fontsize=9)
    ax.set_ylabel("5mC (%)", fontsize=9)
    ax.set_ylim(0, 100)
    if i == 0:
        ax.legend(fontsize=8, loc="upper right")

fig.suptitle("CpG (5mC) Methylation Across LPA", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{out_dir}/cpg_methylation_across_LPA.png", dpi=200, bbox_inches="tight")
plt.close()
print("Saved: cpg_methylation_across_LPA.png")

n_bins = (LPA_END - LPA_START) // BIN_SIZE
heatmap_data = np.full((len(samples), n_bins), np.nan)

for i, sample in enumerate(samples):
    bed_path = f"{bed_dir}/{sample}.modkit_pileup.bed"
    try:
        pos, val = load_5mc_positions(bed_path)
    except FileNotFoundError:
        continue
    bin_idx = ((pos - LPA_START) // BIN_SIZE).astype(int)
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)
    for b in range(n_bins):
        mask = bin_idx == b
        if mask.sum() > 0:
            heatmap_data[i, b] = val[mask].mean()

fig, ax = plt.subplots(figsize=(14, 4))
im = ax.imshow(heatmap_data, aspect="auto", cmap="YlOrRd", vmin=0, vmax=100)
ax.set_yticks(range(len(samples)))
ax.set_yticklabels(samples)
ax.set_xlabel(f"Position bin ({BIN_SIZE} bp each, LPA window)", fontsize=10)
ax.set_title("Mean 5mC (%) Across LPA, Binned", fontsize=13, fontweight="bold")

kiv2_start_bin = (KIV2_START - LPA_START) // BIN_SIZE
kiv2_end_bin = (KIV2_END - LPA_START) // BIN_SIZE
ax.axvline(kiv2_start_bin, color="blue", linestyle="--", linewidth=1)
ax.axvline(kiv2_end_bin, color="blue", linestyle="--", linewidth=1, label="KIV-2 boundary")

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Mean 5mC (%)", fontsize=9)
ax.legend(fontsize=8, loc="upper right")

plt.tight_layout()
plt.savefig(f"{out_dir}/cpg_methylation_heatmap.png", dpi=200, bbox_inches="tight")
plt.close()
print("Saved: cpg_methylation_heatmap.png")
PYEOF

python3 "$SCRIPT_PATH" "$BED_DIR" "$OUT_DIR"

echo "Done. Plots at ${OUT_DIR}"
