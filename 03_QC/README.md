# 03_QC

Quality filtering and QC (Section 2.5.3), and subsetting to the LPA and KIV-2
regions (Sections 2.5.4.1 and 2.5.5.1).

| Script | What it does |
|---|---|
| `filter_qs10_t2t.sh` | Keep reads with qs ≥ 10 (T2T-CHM13v2.0) and check MM/ML tags are kept |
| `filter_qs10_hg38.sh` | Keep reads with qs ≥ 10 (GRCh38) |
| `qc_t2t.sh` | QC before and after filtering, then subset and QC the KIV-2 regions (T2T-CHM13v2.0) |
| `qc_hg38.sh` | QC before and after filtering, then subset and QC the LPA and KIV-2 regions (GRCh38) |
| `mapping_rate_hg38.sh` | flagstat on merged GRCh38 BAMs for mapping rate (Table 10) |
