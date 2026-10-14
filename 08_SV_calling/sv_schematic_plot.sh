#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=sv_schematic
#SBATCH --partition=aoraki
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/project/logs/sv_schematic_%j.out
#SBATCH --error=/path/to/project/logs/sv_schematic_%j.err

set -uo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

if ! python3 -c "import matplotlib" 2>/dev/null; then
  pip install matplotlib --quiet
fi

VCF_DIR=/path/to/project/sv_calling_minimap2
OUT_DIR=/path/to/project/sv_schematic
mkdir -p "$OUT_DIR"
SCRIPT_PATH="${OUT_DIR}/sv_schematic.py"

cat > "$SCRIPT_PATH" << 'PYEOF'
import sys
import gzip
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

LPA_START = 161500000
LPA_END = 162300000
KIV2_START = 161862121
KIV2_END = 161998520

samples = ["SAMPLE_1", "SAMPLE_2", "SAMPLE_3", "SAMPLE_4"]
vcf_dir = sys.argv[1]
out_dir = sys.argv[2]

def parse_sv_vcf(path):
    calls = []
    opener = gzip.open if path.endswith(".gz") else open
    try:
        with opener(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                chrom, pos, alt, info = fields[0], int(fields[1]), fields[4], fields[7]
                if chrom != "NC_060930.1":
                    continue
                svtype = "?"
                svlen = 0
                for kv in info.split(";"):
                    if kv.startswith("SVTYPE="):
                        svtype = kv.split("=")[1]
                    if kv.startswith("SVLEN="):
                        try:
                            svlen = abs(int(kv.split("=")[1]))
                        except ValueError:
                            pass
                if LPA_START <= pos <= LPA_END:
                    calls.append((pos, svtype, svlen))
    except FileNotFoundError:
        pass
    return calls

fig, ax = plt.subplots(figsize=(14, 6))

color_map = {"DEL": "#e41a1c", "INS": "#377eb8", "DUP": "#4daf4a",
             "INV": "#984ea3", "BND": "#ff7f00", "?": "#999999"}

y_positions = {s: i for i, s in enumerate(samples)}

for sample in samples:
    vcf_path = f"{vcf_dir}/{sample}_sniffles2.vcf"
    calls = parse_sv_vcf(vcf_path)
    y = y_positions[sample]
    for pos, svtype, svlen in calls:
        color = color_map.get(svtype, "#999999")
        size = max(20, min(200, svlen / 20))
        ax.scatter(pos, y, s=size, color=color, alpha=0.7, edgecolor="black", linewidth=0.5)

ax.axvspan(KIV2_START, KIV2_END, color="orange", alpha=0.15, label="KIV-2")
ax.set_yticks(list(y_positions.values()))
ax.set_yticklabels(samples)
ax.set_xlabel("Position (T2T-CHM13v2.0, NC_060930.1)", fontsize=11)
ax.set_title("Structural Variant Calls Across LPA (Sniffles2)", fontsize=13, fontweight="bold")
ax.set_xlim(LPA_START, LPA_END)

legend_handles = [mpatches.Patch(color=c, label=t) for t, c in color_map.items() if t != "?"]
legend_handles.append(mpatches.Patch(color="orange", alpha=0.3, label="KIV-2 interval"))
ax.legend(handles=legend_handles, fontsize=9, loc="upper right")

plt.tight_layout()
plt.savefig(f"{out_dir}/sv_schematic.png", dpi=200, bbox_inches="tight")
print("Saved: sv_schematic.png")
PYEOF

python3 "$SCRIPT_PATH" "$VCF_DIR" "$OUT_DIR"

echo "Done. Schematic at ${OUT_DIR}/sv_schematic.png"
