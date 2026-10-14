#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
#SBATCH --job-name=sample3_flye2
#SBATCH --account=ACCOUNT
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/path/to/project/flye_table/sample3_rerun2_%j.log

set -eo pipefail
F=/path/to/project/flye_assembly
MMI=/path/to/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.w10_k15.mmi
S=SAMPLE_3

eval "$(conda shell.bash hook)"; set +u; conda activate custom_env; set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

dedup() { awk 'NR%4==1{split($1,a," "); id=a[1]; keep=!(id in s); s[id]=1} keep' "$1"; }

for h in HP1 HP2; do
  cat $F/fastq/$S.$h.fastq $F/fastq/$S.unassigned.fastq | dedup /dev/stdin > $F/fastq/$S.${h}_plus_unassigned.dedup.fastq
  echo "$h + unassigned (dedup) reads: $(awk 'NR%4==2' $F/fastq/$S.${h}_plus_unassigned.dedup.fastq | wc -l)"
done
dedup $F/fastq/$S.unassigned.fastq > $F/fastq/$S.unassigned.dedup.fastq
echo "Unassigned (dedup) reads: $(awk 'NR%4==2' $F/fastq/$S.unassigned.dedup.fastq | wc -l)"

for rs in HP1_plus_unassigned HP2_plus_unassigned unassigned_only; do
  d=$F/assemblies/${S}_$rs
  if [ -d ${d}_old ]; then rm -rf $d; else [ -d $d ] && mv $d ${d}_old; fi
  case $rs in
    unassigned_only) fq=$F/fastq/$S.unassigned.dedup.fastq ;;
    *) fq=$F/fastq/$S.${rs}.dedup.fastq ;;
  esac
  flye --nano-hq $fq --out-dir $d --threads 8
  minimap2 -ax map-ont $MMI $d/assembly.fasta > $F/assembly_vs_t2t/${S}_${rs}_vs_T2T.sam
  samtools sort -@ 4 -o $F/assembly_vs_t2t/${S}_${rs}_vs_T2T.sorted.bam $F/assembly_vs_t2t/${S}_${rs}_vs_T2T.sam
  samtools index $F/assembly_vs_t2t/${S}_${rs}_vs_T2T.sorted.bam
done

python3 - << 'EOF'
import re, subprocess
F="/path/to/project/flye_assembly"; S,E=161862121,161998520; L=E-S+1
print("Sample\tRead set\tContigs\tMAPQ60\tMax softclip %\tBest KIV-2 %\tLow-MAPQ contigs")
for rs in ["HP1_plus_unassigned","HP2_plus_unassigned","unassigned_only"]:
    b=f"{F}/assembly_vs_t2t/SAMPLE_3_{rs}_vs_T2T.sorted.bam"
    n=sum(1 for l in open(f"{F}/assemblies/SAMPLE_3_{rs}/assembly.fasta") if l.startswith(">"))
    sam=subprocess.run(["samtools","view","-F","0x904",b],capture_output=True,text=True).stdout.splitlines()
    m60=set(); low=[]; clipmax=0; cov={}
    for l in sam:
        x=l.split("\t"); q,pos,mq,cig=x[0],int(x[3]),int(x[4]),x[5]
        ops=re.findall(r"(\d+)([MIDNSHP=X])",cig)
        qlen=sum(int(k) for k,o in ops if o in "MIS=X"); sc=sum(int(k) for k,o in ops if o=="S")
        clipmax=max(clipmax,100*sc/qlen if qlen else 0)
        (m60.add(q) if mq==60 else low.append(f"{q}(MAPQ {mq})"))
        if x[2]=="NC_060930.1":
            span=sum(int(k) for k,o in ops if o in "MDN=X")
            cov[q]=cov.get(q,0)+max(0,min(pos+span-1,E)-max(pos,S)+1)
    print(f"SAMPLE_3\t{rs}\t{n}\t{len(m60)}\t{clipmax:.2f}\t{100*max(cov.values())/L if cov else 0:.1f}\t{','.join(low)}")
EOF
