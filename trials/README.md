# trials

Methods tested but not used in the final analysis.

| Script | Trial |
|---|---|
| `vacmap_internal_sort_sample3.sh` | VACmap with internal sorted output (`-o`), which kept only a fraction of reads (Appendix D) |
| `vacmap_relaxed_filters_sample3.sh` | VACmap with relaxed filtering (`-maxdivergence`, `--nodiscard`), which did not fix the read loss (Appendix D) |
| `vacmap_k9_sample3.sh` | VACmap with `-k 9` genome-wide (Section 2.5.4.4) |
| `vacmap_k9_1thread_sample3.sh` | The same `-k 9` test on one thread |
| `collapsed_reference_copy_number.sh` | Copy number from reads realigned to a reference with KIV-2 collapsed to one unit. Underestimated copy number (Appendix G, Table 36) |
| `flye_pooled_haplotypes.sh` | Flye on HP1 + HP2 + unassigned reads together, as a diagnostic |
| `pooled_assembly_check_sample3.sh` | KIV-2 coverage of the pooled SAMPLE_3 assembly |
| `vep_command_line.sh` | Command-line VEP. Annotation was done on the VEP web interface instead |
