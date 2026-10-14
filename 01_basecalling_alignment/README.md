# 01_basecalling_alignment

Dorado basecalling and alignment (Section 2.5.1). Each script builds a list of
POD5 files, removes unassigned, failed and small (<200 kB) files, then submits
a GPU array job with one task per POD5 file.

| Script | What it does |
|---|---|
| `basecall_align_t2t.sh` | Basecalling with the 5mC/5hmC and 6mA models and alignment to T2T-CHM13v2.0. Written and run by Dr Aaron Stevens |
| `basecall_align_hg38.sh` | Canonical basecalling and alignment to GRCh38, adapted from the T2T script |
