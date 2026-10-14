# 12_KIV2_copy_number

KIV-2 copy number estimation (Section 2.5.6).

| Script | What it does |
|---|---|
| `kilda_t2t.sh` | KILDA on reads from the T2T-CHM13v2.0 LPA window |
| `kilda_hg38.sh` | KILDA on reads from the GRCh38 LPA window |
| `kilda_whole_genome_check_sample1.sh` | KILDA on all reads for one sample, to check the LPA window did not miss KIV-2 reads |
| `array_depth.sh` | Total copy number from depth across the six GRCh38 KIV-2 copies and four control regions |
| `kiv2_allele_sizes_final.sh` | Read counts at the Behera et al. phasing SNPs and allele sizes from array depth (Tables 28 and 37) |
| `validation_fetch.sh` | Pulls the LPA window from public 1000 Genomes ONT Vienna CRAMs. Runs on the login node, not as a SLURM job |
| `validation_run.sh` | KILDA and phasing SNPs on the validation samples, compared with optical mapping |
| `validation_array_depth.sh` | Array depth on the validation samples (Table 27). Run interactively |
