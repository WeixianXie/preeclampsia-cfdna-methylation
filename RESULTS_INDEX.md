# Selected intermediate result tables

This directory contains a subset of the `results/` tables that directly support
the manuscript's tables and figures. Large raw/processed methylation matrices are
not included because of size; they can be regenerated from the public GEO data
using the numbered scripts in `scripts/`.

| File | Use |
|------|-----|
| `GSE282512_dmr_final.csv` | Main 166 DMR list (Table 2 basis) |
| `GSE282512_dmr_sensitivity_candidates.csv` | Sensitivity analyses |
| `GSE282512_dmr_site_level.csv` | Site-level methylation within DMRs |
| `GSE282512_dmr_enrichment_go.csv` | GO/KEGG enrichment results |
| `GSE282512_classifier_performance.csv` | Classifier LOOCV metrics |
| `GSE282512_classifier_panel.csv` | Selected model features |
| `GSE282512_classifier_oof_predictions.csv` | Out-of-fold predictions |
| `GSE282512_classifier_perm_null.csv` | Permutation null AUCs |
| `GSE282512_classifier_shap_importance.csv` | SHAP feature importance (Figure 3D) |
| `GSE282512_deconv_fractions*.csv` | Cell fraction estimates from deconvolution |
| `GSE282512_deconv_tests.csv` | Statistical tests from deconvolution |
| `GSE282512_cellscore_*.csv` | Cell-score / immune associations |
| `GSE282512_placenta_dmr_level.csv` | Placental tissue DMR testing |
| `GSE282512_expr_meta*.csv` | Cross-cohort expression meta-analysis |
| `GSE282512_godmc_*.csv` / `.txt` | GoDMC mQTL overlap and direction triangulation |
| `GSE282512_organ_mr_results.csv` | Organ-wide MR (Table S basis) |
| `GSE282512_wgcna_*.csv` | WGCNA module-trait associations |
| `GSE37722_GSE154378_quasideconv.csv` | Change-based deconvolution leukocyte-cfDNA concordance |
| `GSE154378_subgroup_deltaf*.csv` | Subgroup PreX/GDM/cHTN deconvolution contrasts (Figure 6) |
