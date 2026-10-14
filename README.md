# Investigating LPA Gene Variants Associated with Premature Coronary Artery Disease in Fijian Populations

Scripts written for my BBiomedSc(Hons) thesis at the University of Otago,
Department of Pathology and Molecular Medicine, Wellington (2026).

**Author:** Isabella Yock

**Supervisors:** Dr Pritika Narayan, Dr Aaron Stevens, Professor Klaus Lehnert


## Overview

This project used Oxford Nanopore long-read sequencing to characterise the
LPA gene and the KIV-2 repeat region in four participants from the Fiji Heart
Study. Reads were aligned to two reference genomes:

- **T2T-CHM13v2.0** for alignment assessment, haplotype phasing, realignment
  with VACmap, structural variant calling, local assembly and methylation
- **GRCh38** for variant calling, annotation and KIV-2 copy number estimation

KIV-2 copy number was estimated with KILDA (k-mer based) and from read depth
across the KIV-2 array, with array depth used as the main estimate. Both
methods were validated against optical mapping in five public 1000 Genomes
samples with Oxford Nanopore data.


## Data

No sequencing data, participant data, or sample identifiers are included in
this repository. Data are held under the Fiji Heart Study ethics approval
(HDEC 2025 EXP 12014) and, in line with the study's data and tissue
management plan, individual-level data will not be uploaded to public
repositories.

The validation samples (HG03836, HG00657, HG01302, NA10851, HG00254) are
publicly available from the 1000 Genomes ONT Vienna dataset (Schloissnig et
al., 2025). Their KIV-2 allele sizes from optical mapping are from Behera et
al. (2025).


## Computing environment

All analyses were run on the University of Otago Aoraki HPC cluster. Scripts
were written in Bash and Python and submitted as SLURM batch jobs.

File paths, account names and sample identifiers have been replaced with
placeholders (e.g. `/path/to/project/`, `SAMPLE_ID`). Update these before
running. SLURM log directories need to be created before submitting jobs, as
SLURM opens log files before the script runs.

### Conda environments

| File | Used for | Thesis table |
|---|---|---|
| `envs/custom_env.yml` | Main analysis environment, including Jellyfish for KILDA | Table 1 |
| `envs/vacmap_env.yml` | VACmap realignment | Table 2 |
| `envs/bandage_env.yml` | Assembly graph visualisation | Table 3 |

Recreate an environment with:

    conda env create -f envs/custom_env.yml

VACmap needed its own environment because it requires minimap2 v2.29, which
conflicts with minimap2 v2.31 in the main environment. Bandage was kept
separate because its Qt graphics libraries share low-level dependencies with
tools in the main environment.

### Software outside the Conda environments

| Tool | Version | Use |
|---|---|---|
| Dorado | 1.4.0 | Basecalling and alignment (standalone binary) |
| KILDA | commit 3c221c2 (29 June 2026) | KIV-2 copy number from k-mers |
| BLASTn (NCBI web) | - | Deriving KIV-2 coordinates |
| Ensembl VEP (web) | - | Variant annotation |
| IGV (web) | 2.4.10 | Read visualisation |

Full tool versions are listed in Tables 1-4 of the thesis.

### Dorado models

| Model | Use |
|---|---|
| `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` | Canonical basecalling |
| `dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mC_5hmC@v2` | 5mC and 5hmC calls (T2T-CHM13v2.0 alignment only) |
| `dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1` | 6mA calls (T2T-CHM13v2.0 alignment only) |


## Repository structure

| Folder | Contents | Thesis section |
|---|---|---|
| `reference/` | BED files for LPA, KIV-2, flanking, promoter and array-depth control regions | Tables 5, 6 and 33 |
| `01_basecalling_alignment/` | POD5 filtering, Dorado basecalling with modified base models, and alignment to T2T-CHM13v2.0 and GRCh38 | 2.5.1 |
| `02_merge_validation/` | BAM integrity checks, shard count validation and per-sample merging | 2.5.2 |
| `03_QC/` | Whole-genome QC before and after qs≥10 filtering | 2.5.3 |
| `04_region_filtering/` | Subsetting BAMs to LPA and KIV-2 (padded and unpadded) on both references, and region QC | 2.5.4.1, 2.5.5.1 |
| `05_alignment_metrics/` | Supplementary alignment, soft-clipping, MAPQ and mismatch rate by region | 2.5.4.2 |
| `06_haplotyping/` | Longshot phasing, WhatsHap and LongPhase haplotagging, HP1/HP2 splitting and haplotype depth | 2.5.4.3 |
| `07_VACmap/` | Read extraction, VACmap realignment, and repeat of alignment metrics and haplotyping | 2.5.4.4 |
| `08_SV_calling/` | Sniffles2 and cuteSV on padded KIV-2 BAMs | 2.5.4.5 |
| `09_local_assembly/` | Read extraction by haplotype, Flye assembly, contig realignment, KIV-2 coverage script and Bandage graphs | 2.5.4.5.1, 2.5.4.5.2 |
| `10_methylation/` | modkit pileup of CpG 5mC calls, all reads and by haplotype | 2.5.4.6 |
| `11_variant_calling/` | Read group fix, GATK HaplotypeCaller, genotyping at known variants, filtering and preparation for VEP | 2.5.5.2, 2.5.5.3 |
| `12_KIV2_copy_number/` | KILDA, array depth, allele sizes from phasing SNPs, and validation against optical mapping | 2.5.6 |
| `trials/` | Methods tested but not used in the final analysis | see below |

Statistical analysis (Mann-Whitney U tests) was done in IBM SPSS Statistics
and is not included here.


## Trials not used in the final analysis

These scripts are kept so the full analysis process can be followed.

**VACmap small k-mer size (2.5.4.4).** Setting `-k 9` was tested to improve
resolution within the repeat. The `-k` parameter applies to the genome-wide
search stage, and the local refinement step already uses a fixed k-mer size
of 9 by default. Applying k=9 genome-wide made alignment computationally
infeasible, so default settings were used.

**VACmap internal sorted output (Appendix D).** Writing sorted output
directly with `-o output.sorted.bam` kept only a small fraction of input
reads. Piping VACmap output through `samtools sort` recovered all reads and
was used for all final alignments.

**Collapsed-reference copy number (2.5.6.2, Appendix G).** Reads were
realigned to a GRCh38 LPA region with the six KIV-2 copies collapsed into
one. This underestimated copy number and was replaced by array depth across
the six reference copies.


## Reference coordinates

KIV-2 is not annotated as a standalone feature in either assembly.
Coordinates were derived with BLASTn using the GRCh37 KIV-2 region
(chr6:161,033,785-161,066,618; Coassin et al., 2019) as the query, and the
KIV-2 interval was defined as the span from the first to the last hit. LPA
gene boundaries were cross-checked against the NCBI RefSeq annotation
(GeneID 4018).

### T2T-CHM13v2.0 (NC_060930.1)

| Region | Coordinates |
|---|---|
| LPA window | 161,500,000-162,300,000 |
| KIV-2 (unpadded) | 161,862,121-161,998,520 |
| KIV-2 (padded 50 kb) | 161,812,121-162,048,520 |
| Flanking, 5′ arm | 161,812,121-161,862,120 |
| Flanking, 3′ arm | 161,998,521-162,048,520 |
| Promoter (TSS ± 1 kb) | 162,010,763-162,012,762 |

### GRCh38 (chr6, NC_000006.12)

| Region | Coordinates |
|---|---|
| LPA window | 160,300,000-160,800,000 |
| KIV-2 (unpadded) | 160,612,753-160,645,586 |
| KIV-2 (padded 50 kb) | 160,562,753-160,695,586 |
| KIV-2 array, six reference copies (array depth) | 160,613,492-160,646,773 |
| Control, upstream | 160,320,001-160,330,000 |
| Control, LPA 3′ | 160,535,001-160,545,000 |
| Control, LPA 5′ | 160,680,001-160,690,000 |
| Control, downstream | 160,770,001-160,780,000 |

The array-depth interval is defined by KIV-2 unit boundaries (Behera et al.,
2025) and differs slightly from the BLAST-derived KIV-2 region.


## Script origin

Dr Aaron Stevens basecalled the raw reads, aligned them to both references,
merged them into one BAM per sample, filtered by qs≥10, and subset the merged
BAMs to the LPA locus. Dr Klaus Lehnert provided the short-read VCFs used to
check long-read genotypes. Some other scripts were adapted from scripts
written by Aaron Stevens. The origin of each script, including any AI
assistance, is recorded in Appendix B of the thesis.


## Citation

Yock, I. (2026). *Investigating LPA gene variants associated with premature
coronary artery disease in Fijian populations* [BBiomedSc(Hons) thesis,
University of Otago].


## Contact

Isabella Yock - isabellayock@gmail.com
