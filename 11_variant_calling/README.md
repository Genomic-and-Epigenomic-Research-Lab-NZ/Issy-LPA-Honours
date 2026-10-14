# 11_variant_calling

Variant calling, filtering and genotyping on GRCh38 (Sections 2.5.5.2 and
2.5.5.3). Annotation was done on the Ensembl VEP web interface using the
sites-only VCF from `variant_filtering_summary.sh`.

| Script | What it does |
|---|---|
| `gatk_haplotypecaller_hg38.sh` | Add read groups with Picard, then GATK HaplotypeCaller over the LPA window |
| `variant_filtering_summary.sh` | Split multi-allelics, filter (QUAL ≥ 30, DP ≥ 10), count variants, check depth at known variants, write the sites-only VCF for VEP |
| `known_variant_genotypes.sh` | Genotype and read depth (MAPQ ≥ 20) at rs10455872, rs3798220, rs41272114 and rs41267813 (Table 24). Run interactively |
| `short_read_genotype_check.sh` | Looks up the same four variants by position (±50 bp) in the short-read VCF. Run interactively |
| `genotypes_flagged_variants.sh` | Genotypes and read counts at flagged coding and Lp(a)-associated variants (Tables 25 and 26) |
