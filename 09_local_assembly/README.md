# 09_local_assembly

Haplotype-stratified local assembly and assembly validation (Sections
2.5.4.5.1 and 2.5.4.5.2).

| Script | What it does |
|---|---|
| `flye_haplotype_assembly.sh` | Recover HP1 and HP2 reads by name from the whole-genome qs ≥ 10 BAM and assemble each with Flye |
| `flye_unassigned.sh` | Assemble unassigned reads alone and HP1/HP2 + unassigned (duplicates removed) |
| `rerun_unassigned_assemblies_sample3.sh` | Re-run of the SAMPLE_3 + unassigned assemblies with duplicate read IDs removed |
| `contig_realign_t2t.sh` | Realign contigs to T2T-CHM13v2.0 with minimap2 |
| `contig_kiv2_coverage_qc.sh` | KIV-2 coverage per contig from the CIGAR string, and contig MAPQ, soft-clip and SA (Table 21) |

Bandage was run with its `image` command in `bandage_env` (Table 3).
