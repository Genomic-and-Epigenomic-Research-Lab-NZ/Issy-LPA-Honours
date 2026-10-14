# 06_haplotyping

Phasing, haplotagging and haplotype-resolved metrics on minimap2
T2T-CHM13v2.0 alignments (Section 2.5.4.3).

| Script | What it does |
|---|---|
| `longshot_phase_t2t.sh` | Longshot phasing of heterozygous SNPs in the LPA window |
| `whatshap_haplotag_t2t.sh` | WhatsHap haplotagging |
| `whatshap_split_t2t.sh` | Split WhatsHap BAMs into HP1 and HP2 |
| `whatshap_haplotype_qc_t2t.sh` | SA, soft-clip, MAPQ and edit rate per haplotype (Tables 15 and 17) |
| `whatshap_haplotype_mosdepth_t2t.sh` | Depth per haplotype |
| `longphase_haplotag_t2t.sh` | LongPhase haplotagging (`-q 0 --tagSupplementary`) |
| `longphase_split_t2t.sh` | Split LongPhase BAMs into HP1 and HP2 |
| `longphase_haplotype_qc_t2t.sh` | SA, soft-clip, MAPQ and edit rate per haplotype (Tables 16 and 18) |
| `longphase_mosdepth_t2t.sh` | Depth per haplotype |
| `haplotype_qc_all_aligners.sh` | Haplotype QC across minimap2 and VACmap, WhatsHap and LongPhase |
