# Script-to-output mapping

This file maps each numbered analysis script to its main output(s) in `results/` or `figures/`.

## Phase 0 — Data download and QC

| Script | Output | Description |
|--------|--------|-------------|
| `01_download_geo.sh` / `01_download_geo.ps1` | `data/` (raw GEO tar files) | Downloads GEO series matrices and supplementary files |
| `02_download_gwas.sh` | `data/gwas/` | FinnGen / GWAS Catalog / GoDMC summary stats |
| `03_download_raw_tar.ps1` | `data/` | Raw tar retrieval helpers |
| `05_check_data_integrity.R` | `phase0_data_check.csv` | Sample and file integrity QC |
| `21_build_wgbs_metadata.R` | GSE154378 metadata | Metadata tables for external WGBS cohort |
| `22_select_subcohort.R` | `GSE282512_subcohort.csv` | GA-matched 32 PE / 32 CT subcohort |

## Phase 1 — DMR discovery

| Script | Output | Description |
|--------|--------|-------------|
| `23_build_region_matrix.R` | Region-level methylation matrix | Coverage-weighted autosomal regions |
| `24_dmr_discovery.R` | `GSE282512_dmr_final*.csv` | Primary limma DMR calls |
| `25_dmr_validate.R` | validation summaries | Permutation and fold-sensitivity checks |
| `27_dmr_sensitivity.R` | `GSE282512_dmr_sensitivity*.csv` | EDTA-only, early-onset, tube interaction |
| `28_dmr_final_call.R` / `28_dmr_final_v2.R` | final 166 DMR set | Final region list with annotations |
| `29_extract_mqtl_gwas.sh` | tabix-extracted GWAS | GoDMC / FinnGen variant extraction |

## Phase 2 — mQTL / classifier / deconvolution

| Script | Output | Description |
|--------|--------|-------------|
| `30_phase2_mqtl_overlap.R` | `GSE282512_phase2_*.csv` | GoDMC cis-mQTL overlap |
| `37_jaspar_motif.R` | `GSE282512_jaspar_*.csv` | JASPAR2024 motif enrichment |
| `38_leukocyte_deconv.R` / `39_cellscore_deconv.R` | `GSE282512_deconv_fractions*.csv`, `GSE282512_cellscore_*.csv` | Reference- and change-based deconvolution, immune score |
| `40_dmr_classifier.R` | `GSE282512_classifier_*.csv` | Elastic net classifier (LOOCV, permutation, SHAP, DCA) |

## Phase 3 — Integration and validation

| Script | Output | Description |
|--------|--------|-------------|
| `34_gse98224_expr_validate.R` / `35_cross_cohort_expr.R` | `GSE282512_expr_meta*.csv` | Placental and blood expression meta-analysis |
| `41_manuscript_docx.py` (predecessor) | `docs/manuscript_main_*.docx` | Manuscript drafting scripts |
| `42_dmr_enrichment.R` / `42b_kegg_enrich.R` | `GSE282512_dmr_enrichment_*.csv` | GO/KEGG enrichment |
| `43_model_deep.R` | `GSE282512_classifier_model_comparison.csv` | Deep model checks |
| `45_gse37722_leukocyte.R` | `GSE37722_probe_marker_map.csv`, `GSE37722_GSE154378_quasideconv.csv` | Leukocyte-cfDNA concordance |
| `46_smr_coloc.R` | `GSE282512_smr_*.csv`, `GSE282512_coloc_results.csv` | SMR / coloc sensitivity |
| `47_wgcna.R` | `GSE282512_wgcna_*.csv` | Weighted co-methylation modules |
| `48_immune_infiltration.R` | `GSE282512_immune_*.csv` | Immune infiltration scoring |
| `49d_placenta_discrimination.R` | `GSE282512_placenta_*.csv` | Placental origin testing |

## Phase 4 — MR and manuscript production

| Script | Output | Description |
|--------|--------|-------------|
| `49_gse154378_aggregate.py` / `50_gse154378_trajectory.R` | `GSE154378_*.csv` | External WGBS trajectory and subgroup analyses |
| `59_godmc_annotation.py` | `GSE282512_godmc_*.csv` | GoDMC mQTL annotation |
| `60_outcome_discovery.py` / `61_extract_outcome_effects.py` | outcome instrument sets | MR outcome preparation |
| `62_organ_mr.py` | `GSE282512_organ_mr_results.csv` | Organ-wide MR (112 tests) |
| `63_manuscript_docx.py` | `docs/manuscript_main_*.docx` | Markdown → docx converter (supports tables) |
| `64_figure1_design.py` | `figures/Figure1_study_design.png` | Study design diagram |
| `65_figure_compose.py` | `figures/Figure2_*` … `Figure8_*` | Multi-panel TIFF/PNG submission figures |

## Key helper modules

| Script | Purpose |
|--------|---------|
| `rtabix.py` | Pure-Python remote tabix reader for GWAS Catalog / EBI FTP harmonised summary statistics |
