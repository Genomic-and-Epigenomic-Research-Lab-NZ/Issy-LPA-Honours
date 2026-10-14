# 10_methylation

Methylation analysis with modkit (Section 2.5.4.6).

| Script | What it does |
|---|---|
| `methylation_cpg.sh` | modkit pileup with `--cpg --combine-strands` on all reads and on LongPhase HP1/HP2, then summarises KIV-2, the rest of LPA and flanking sequence |
| `modkit_pileup_all_context.sh` | Earlier all-context pileup on all reads. Used only for the 6mA summary |
| `modkit_pileup_all_context_by_haplotype.sh` | The same by haplotype. Used only for the 6mA summary |
| `methylation_plots.sh` | Methylation plots across the LPA locus |
