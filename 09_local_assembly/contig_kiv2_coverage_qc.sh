#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=kiv2_quantification
#SBATCH --partition=aoraki
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/logs/kiv2_quantification_%j.out
#SBATCH --error=/path/to/project/logs/kiv2_quantification_%j.err

set -uo pipefail
mkdir -p /path/to/project/logs

eval "$(conda shell.bash hook)"
set +u
conda activate custom_env
set -u

SAM_DIR=/path/to/project/flye_assembly/assembly_vs_t2t
SCRIPT_PATH="${SAM_DIR}/estimate_kiv2_copies.py"

cat > "$SCRIPT_PATH" << 'PYEOF'
#!/usr/bin/env python3
import re
import sys
import csv

KIV2_START = 161862121
KIV2_END = 161998520
KIV2_REF_LENGTH = KIV2_END - KIV2_START
ASSUMED_UNIT_SIZE = 5600

def cigar_query_span_in_ref_window(cigar, ref_start, window_start, window_end):
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    ref_pos = ref_start
    query_bases_in_window = 0
    ref_bases_in_window = 0

    for length, op in ops:
        length = int(length)
        if op in ('M', '=', 'X'):
            new_ref_pos = ref_pos + length
            overlap_start = max(ref_pos, window_start)
            overlap_end = min(new_ref_pos, window_end)
            if overlap_end > overlap_start:
                overlap_len = overlap_end - overlap_start
                query_bases_in_window += overlap_len
                ref_bases_in_window += overlap_len
            ref_pos = new_ref_pos
        elif op == 'D':
            new_ref_pos = ref_pos + length
            overlap_start = max(ref_pos, window_start)
            overlap_end = min(new_ref_pos, window_end)
            if overlap_end > overlap_start:
                ref_bases_in_window += (overlap_end - overlap_start)
            ref_pos = new_ref_pos
        elif op == 'I':
            if window_start <= ref_pos < window_end:
                query_bases_in_window += length
    return query_bases_in_window, ref_bases_in_window

def main(sam_path, out_csv):
    rows = []
    with open(sam_path) as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 6:
                continue
            qname, flag, rname, pos, mapq, cigar = fields[0:6]
            flag = int(flag)
            if rname != "NC_060930.1":
                continue
            if cigar == "*":
                continue
            if flag & 0x100 or flag & 0x800:
                continue

            ref_start = int(pos) - 1
            q_span, r_span = cigar_query_span_in_ref_window(
                cigar, ref_start, KIV2_START, KIV2_END
            )
            if r_span == 0:
                continue

            net_diff_bp = q_span - r_span
            est_extra_units = net_diff_bp / ASSUMED_UNIT_SIZE

            rows.append({
                "contig": qname,
                "ref_span_in_KIV2_bp": r_span,
                "contig_bases_in_KIV2_bp": q_span,
                "net_diff_bp": net_diff_bp,
                "est_copy_number_diff_vs_ref": round(est_extra_units, 2),
                "fraction_of_KIV2_window_covered": round(r_span / KIV2_REF_LENGTH, 3),
            })

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "contig", "ref_span_in_KIV2_bp", "contig_bases_in_KIV2_bp",
            "net_diff_bp", "est_copy_number_diff_vs_ref",
            "fraction_of_KIV2_window_covered"
        ])
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"Wrote {len(rows)} contig(s) overlapping KIV-2 to {out_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 estimate_kiv2_copies.py <input.sam> <output.csv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
PYEOF

MASTER_CSV="${SAM_DIR}/all_samples_kiv2_copies.csv"
echo "sample_variant,contig,ref_span_in_KIV2_bp,contig_bases_in_KIV2_bp,net_diff_bp,est_copy_number_diff_vs_ref,fraction_of_KIV2_window_covered" > "$MASTER_CSV"

for SAM in "${SAM_DIR}"/*_vs_T2T.sam; do
  BASE=$(basename "$SAM" .sam)
  echo "################ ${BASE} ################"

  SORTED_BAM="${SAM_DIR}/${BASE}.sorted.bam"
  if [ ! -f "$SORTED_BAM" ]; then
    echo "[${BASE}] Sorting and indexing for IGV..."
    samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$SORTED_BAM" "$SAM"
    samtools index "$SORTED_BAM"
  fi

  OUT_CSV="${SAM_DIR}/${BASE}_kiv2_copies.csv"
  if [ ! -f "$OUT_CSV" ]; then
    echo "[${BASE}] Running quantification..."
    python3 "$SCRIPT_PATH" "$SAM" "$OUT_CSV"
  fi

  if [ -f "$OUT_CSV" ]; then
    tail -n +2 "$OUT_CSV" | awk -v tag="$BASE" -F',' 'BEGIN{OFS=","} {print tag, $0}' >> "$MASTER_CSV"
  fi

  CONTIG_QC_CSV="${SAM_DIR}/${BASE}_contig_qc.csv"
  if [ ! -f "$CONTIG_QC_CSV" ]; then
    echo "[${BASE}] Running contig alignment QC..."
    echo "contig,total_contig_length,aligned_length,softclip_length,softclip_pct,mapq,has_SA" > "$CONTIG_QC_CSV"

    grep -v "^@" "$SAM" | awk -F'\t' 'BEGIN{OFS=","}
      {
        flag=$2
        is_secondary = (int(flag/256) % 2 == 1)
        is_supplementary = (int(flag/2048) % 2 == 1)
        if (is_secondary || is_supplementary) next

        contig=$1; mapq=$5; cigar=$6

        num=""; aligned=0; softclip=0
        for (i=1; i<=length(cigar); i++) {
          c=substr(cigar,i,1)
          if (c ~ /[0-9]/) { num = num c }
          else {
            if (c=="M" || c=="I" || c=="=" || c=="X") aligned += num
            if (c=="S") softclip += num
            num=""
          }
        }
        total_len = aligned + softclip

        has_sa = "no"
        for (i=12; i<=NF; i++) {
          if ($i ~ /^SA:Z:/) has_sa = "yes"
        }

        softclip_pct = (total_len > 0) ? (softclip/total_len*100) : 0
        printf "%s,%d,%d,%d,%.2f,%s,%s\n", contig, total_len, aligned, softclip, softclip_pct, mapq, has_sa
      }' >> "$CONTIG_QC_CSV"
  fi

  echo ""
done

echo "Done. Combined results at $MASTER_CSV"
echo "Contig alignment QC files: ${SAM_DIR}/*_contig_qc.csv"
echo "Sorted/indexed BAMs (for IGV) at ${SAM_DIR}/*.sorted.bam"
