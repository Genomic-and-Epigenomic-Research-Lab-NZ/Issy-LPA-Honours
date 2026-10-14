#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running

set -euo pipefail

eval "$(conda shell.bash hook)"
set +u
conda activate samtools_remote
set -u
samtools --version | grep -q "libcurl=yes" || { echo "STOP: this samtools cannot read remote files (libcurl=no)"; exit 1; }

BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA
OUT=/path/to/project/kilda_validation
REGION="chr6:160300000-160800000"
mkdir -p ${OUT}/ref ${OUT}/bam
cd ${OUT}

SAMPLES=(HG02731 HG03492 HG03857 HG03836 HG04214 HG04002 HG03785 HG04152 HG03022
         NA19316 HG00657 HG03279 HG01302 NA18990 NA10851 HG00254 HG01516)

if [ ! -s ref/1KG_ONT_VIENNA_hg38.fa.fai ]; then
  echo "Downloading reference (~1 GB compressed)..."
  curl -s ${BASE}/reference/1KG_ONT_VIENNA_hg38.fa.gz -o ref/1KG_ONT_VIENNA_hg38.fa.gz
  gunzip -c ref/1KG_ONT_VIENNA_hg38.fa.gz > ref/1KG_ONT_VIENNA_hg38.fa
  samtools faidx ref/1KG_ONT_VIENNA_hg38.fa
fi
grep -P "^chr6\t" ref/1KG_ONT_VIENNA_hg38.fa.fai || { echo "STOP: no 'chr6' in reference"; exit 1; }

: > available_samples.txt
for S in "${SAMPLES[@]}"; do
  URL=${BASE}/hg38/${S}.hg38.cram
  CODE=$(curl -s -o /dev/null -w "%{http_code}" -I ${URL})
  if [ "${CODE}" != "200" ]; then
    echo "${S}: not in Vienna set (HTTP ${CODE}) - skipped"
    continue
  fi
  if [ ! -s bam/${S}_LPA.bam ]; then
    echo "${S}: pulling ${REGION} ..."
    samtools view -b -T ref/1KG_ONT_VIENNA_hg38.fa ${URL} ${REGION} > bam/${S}_LPA.bam
    samtools index bam/${S}_LPA.bam
  fi
  N=$(samtools view -c -F 0x900 bam/${S}_LPA.bam)
  echo "${S}: ${N} primary reads"
  echo ${S} >> available_samples.txt
done

rm -f *.crai
echo "=== $(wc -l < available_samples.txt) samples ready ==="
cat available_samples.txt
