# reference

BED files for the regions used to subset BAMs and calculate region coverage.
These are the same coordinates written by `03_QC/qc_t2t.sh` and
`03_QC/qc_hg38.sh` (Tables 5 and 6).

| File | Region |
|---|---|
| `T2T_LPA_expanded.bed` | LPA window, T2T-CHM13v2.0 |
| `T2T_KIV2_unpadded.bed` | KIV-2, T2T-CHM13v2.0 |
| `T2T_KIV2_padded.bed` | KIV-2 ± 50 kb, T2T-CHM13v2.0 |
| `GRCh38_LPA_expanded.bed` | LPA window, GRCh38 |
| `GRCh38_KIV2_unpadded.bed` | KIV-2, GRCh38 |
| `GRCh38_KIV2_padded.bed` | KIV-2 ± 50 kb, GRCh38 |
