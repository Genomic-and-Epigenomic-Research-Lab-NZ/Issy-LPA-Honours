#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
cd /path/to/project/kilda_validation
eval "$(conda shell.bash hook)"
conda activate custom_env

printf "chr6\t160613491\t160646773\tKIV2_array\nchr6\t160320000\t160330000\tC1\nchr6\t160535000\t160545000\tC2\nchr6\t160680000\t160690000\tC3\nchr6\t160770000\t160780000\tC4\n" > array_regions.bed
mkdir -p depth
for S in HG03836 HG00657 HG01302 NA10851 HG00254; do
  mosdepth -n --by array_regions.bed depth/${S} bam/${S}_LPA.bam
  zcat depth/${S}.regions.bed.gz | awk -v S=$S '{d[$4]=$5} END{n=split(d["C1"]" "d["C2"]" "d["C3"]" "d["C4"],c," "); asort(c); m=(c[2]+c[3])/2; printf "%s\tarray=%.1f\tctrl=%.1f\tCN=%.1f\n", S, d["KIV2_array"], m, 12*d["KIV2_array"]/m}'
done
