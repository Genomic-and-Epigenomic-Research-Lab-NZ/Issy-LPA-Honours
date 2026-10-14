#!/bin/bash
# Fill in the /path/to/ placeholders, ACCOUNT and SAMPLE names before running
VCF=/path/to/short_read_LPA.vcf.gz

bcftools query -l $VCF
bcftools view -h $VCF | grep '^##contig' | grep -i -E 'chr6|NC_000006|ID=6,'

[ -f ${VCF}.tbi ] || [ -f ${VCF}.csi ] || bcftools index -t $VCF

for pos in 160589086 160540105 160585045 160577167; do
  echo "== $pos"
  bcftools query -r chr6:$((pos-50))-$((pos+50)) -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GT]\n' $VCF
done

echo "== records in chr6:160500000-160700000"
bcftools view -H -r chr6:160500000-160700000 $VCF | wc -l

echo "== rsID search"
bcftools query -f '%ID\n' $VCF | grep -c '^rs'
bcftools query -f '%CHROM\t%POS\t%ID\n' $VCF | grep -wE 'rs10455872|rs3798220|rs41272114|rs41267813'
