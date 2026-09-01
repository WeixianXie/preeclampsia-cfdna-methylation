# Reproducible analysis code: cfDNA methylation signatures in pre-eclampsia

This repository contains the analysis pipeline for the manuscript
**"Maternal leukocyte composition, not placental pathology, underlies
pre-eclampsia-associated differences in plasma cell-free DNA methylation"**.

- **Corresponding author / PI**: Weixian Xie, MD (Master of Medicine, Guangzhou Medical University), Department of Anesthesiology, Affiliated Qingyuan Hospital, Guangzhou Medical University, Qingyuan, Guangdong, China. Email: 2026691021@gzhmu.edu.cn
- **ORCID**: 0009-0005-5285-8716

## How to cite this code

> Xie W. Analysis code and intermediate result tables for cfDNA methylation signatures in pre-eclampsia [Software]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX  *(DOI to be inserted after release)*

## Pipeline overview

The pipeline is organized in five phases. Each script is numbered to show the intended execution order.

| Phase | Scripts | Purpose |
|-------|---------|---------|
| 0. Data download and QC | `00_*` to `21_*` | Download GEO datasets, FinnGen, GoDMC, 1000 Genomes, GWAS summary statistics, build metadata and sample QC |
| 1. DMR discovery | `22_*` to `29_*` | Region-level cfRRBS methylation matrix, limma DMR calling, sensitivity analyses, final DMR set |
| 2. Phase 2: mQTL overlap | `30_*` to `39_*` | GoDMC cis-mQTL overlap, classifier, leukocyte deconvolution, immune score |
| 3. Phase 3: integration | `40_*` to `49_*` | Cross-cohort expression validation, enrichment, JASPAR motifs, placenta origin, MR preparation |
| 4. Manuscript and figures | `50_*` to `65_*` | Deconvolution plots, WGCNA, SMR/coloc, figures, tables, docx conversion |

A detailed script-to-output mapping is provided in `SCRIPTS.md`.

## System requirements

- **R** 4.6.1 or later with Bioconductor 3.20 or later
  - Core packages: `limma`, `WGCNA`, `EpiDISH`, `minfi`, `DMRcate`, `GenomicRanges`, `data.table`, `ggplot2`, `ggsignif`, `pROC`, `glmnet`, `TwoSampleMR` (or equivalent MR functions), `meta`, `dplyr`, `tidyr`, `readr`, `stringr`, `R.utils`, `httr`
- **Python** 3.13 or later
  - Core packages: `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `pillow`, `requests`, `pysam` / `pysam` substitute, `docx` (`python-docx`)
- **External tools**: `tabix` (optional; the project includes a pure-Python fallback `scripts/rtabix.py`)
- **OS**: Developed and tested on Windows 10/11 with Git Bash. Most scripts assume POSIX-style paths inside the project root.

## Data sources

All source data are public:

- cfRRBS discovery cohort: GEO **GSE282512**
- cfDNA WGBS longitudinal cohort: GEO **GSE154378**
- Maternal leukocyte HM27 cohort: GEO **GSE37722**
- Placental HM450 cohorts: **GSE57767**, **GSE73375**, **GSE75196**
- Early-pregnancy blood expression cohorts: **GSE85307**, **GSE86200**, **GSE98224**
- cis-meQTL summary statistics: **GoDMC**
- GWAS outcomes: FinnGen R10 (PE, gestational hypertension), Tyrmi 2023 PE meta-GWAS, Sakaue 2021 atlas, Keaton 2024 blood pressure

Exact accession numbers and GCST IDs are listed in the manuscript.

## How to run

1. Clone or download this repository and enter the project root.
2. Edit path variables at the top of the early scripts if your working directory differs.
3. Run scripts in numeric order. Most R scripts can be executed with `Rscript`; Python scripts with `python`.
4. Output figures are written to `figures/` and intermediate/result tables to `results/`.

Because raw GEO/FTP downloads are large, the scripts are split into download + analysis stages so that intermediate files can be reused.

## Repository contents

- `scripts/` — analysis code (R, Python, shell helper scripts)
- `results/` — selected intermediate result tables backing the manuscript tables and figures
- `figures/` — source PNGs and composed submission figures
- `docs/` — manuscript markdown/docx, cover letter, and EGA application materials

## License

This code is released under the MIT License. See `LICENSE` for details.

When using these scripts or intermediate data in a publication, please cite the
Zenodo DOI shown above and the associated manuscript.

## Contact

For questions about the analysis, please contact Weixian Xie at 2026691021@gzhmu.edu.cn.
