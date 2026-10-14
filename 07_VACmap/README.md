# 07_VACmap

VACmap realignment of the LPA mini BAM reads (Section 2.5.4.4).

SAMPLE_3 was run first, one step at a time. The other three samples were then
run through the same steps in a single script.

| Script | What it does |
|---|---|
| `vacmap_align_diagnostics_sample3.sh` | Diagnostic runs on SAMPLE_3. Test 3 (VACmap piped through `samtools sort`) produced the final SAMPLE_3 alignment |
| `alignment_assessment_vacmap_sample3.sh` | SA, soft-clip and MAPQ by region for SAMPLE_3 |
| `longshot_phase_vacmap_sample3.sh` | Longshot phasing on the SAMPLE_3 VACmap alignment |
| `haplotag_vacmap_sample3.sh` | WhatsHap and LongPhase haplotagging for SAMPLE_3 |
| `split_mosdepth_vacmap_sample3.sh` | HP1/HP2 splitting and haplotype depth for SAMPLE_3 |
| `haplotype_qc_vacmap_sample3.sh` | Haplotype-resolved alignment metrics for SAMPLE_3 |
| `vacmap_pipeline_other_samples.sh` | All of the above for SAMPLE_1, SAMPLE_2 and SAMPLE_4 |
| `mismatch_rate_vacmap.sh` | Edit rate per haplotype on VACmap alignments, all samples |

Runs that were tested and not used are in `trials/`.
