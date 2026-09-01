# Maternal leukocyte composition, not placental pathology, underlies pre-eclampsia-associated differences in plasma cell-free DNA methylation: an integrative cell-free reduced-representation bisulfite sequencing study with tissue, trajectory, and genetic validation

Weixian Xie1, [Co-authors to be completed]

1 Department of Anesthesiology, Affiliated Qingyuan Hospital, Guangzhou Medical University, Qingyuan, Guangdong, China

*Correspondence: 2026691021@gzhmu.edu.cn

**Keywords:** pre-eclampsia; cell-free DNA; DNA methylation; leukocyte composition; reduced-representation bisulfite sequencing; Mendelian randomization

---

## Abstract

**Background.** Plasma cell-free DNA (cfDNA) methylation is an attractive non-invasive readout in pregnancy, but the cellular origin of pre-eclampsia (PE)-associated methylation differences is unresolved: placental trophoblast, maternal leukocytes, or both. The answer determines what a cfDNA methylation biomarker for PE would actually measure.

**Methods.** We analysed cell-free reduced-representation bisulfite sequencing (cfRRBS) of 369 plasma cfDNA samples from 279 pregnant women (GSE282512), detecting differentially methylated regions (DMRs) with coverage-weighted region-level methylation and limma, adjusted for blood collection tube and batch. Origin was tested in five independent layers: cross-cohort expression meta-analysis of DMR genes (8 cohorts, 6 placental and 2 maternal blood); direct testing of DMRs in placental tissue from three independent cohorts (105 samples); regulatory and transcription factor motif annotation; longitudinal trajectories in an external cfDNA cohort (GSE154378, 134 samples) and a maternal leukocyte cohort (GSE37722, 84 samples, five time points) using change-based deconvolution; and genetic analyses comprising cis-mQTL coverage (GoDMC), two-sample Mendelian randomization (FinnGen R10, Tyrmi 2023 meta-GWAS), and organ-wide MR across 16 blood cell and organ traits.

**Results.** We identified 166 candidate DMRs (105 hyper-, 61 hypomethylated), robust to collection tube (100% direction consistency; 63 replicated at FDR<0.05 in EDTA-only samples) and enriched in early-onset PE. An elastic net classifier reached a leave-one-out AUC of 0.932 (95% CI 0.866–0.998; permutation p=0.048); we report it as an internally validated candidate signature only: leave-one-out estimates in a 64-sample training subcohort are prone to optimism even with fold-internal tuning and a significant permutation null, and independent external validation is lacking. DMR genes showed no coordinated differential expression in placenta (0/66 genes at meta p<0.05). In placental tissue, DMR CpGs carried smaller PE effects than genomic background (mean |z| 0.346 vs 0.444; permutation p=1.0), agreed with cfDNA direction in only 46% of regions (chance level), and correlated with cfDNA effects at rho=−0.036, while a positive control set of the 500 strongest placental PE CpGs reached mean |z|=3.99. In pregnancy trajectories, placental-marker methylation rose with gestation (mean rho +0.542, p=6×10−51) while blood-cell markers fell; T-cell fractions declined in both cfDNA and leukocyte cohorts, with postpartum reversal. In the PE/hypertension subgroup, the neutrophil fraction was elevated at delivery (Δ +0.063, 95% CI 0.024–0.098) despite a lower placental fraction (Δ −0.126), and a shift in the same direction was already present in the first trimester (+0.070, CI 0.012–0.129); we treat this as hypothesis-generating because the subgroup is small (eight patients). DMR CpGs carried fewer cis-mQTLs than background (28.8% vs 37.7%, OR 0.669, p=2.1×10−9; hypermethylated DMRs 22.2%), indicating weak genetic regulation. MR identified cg19490609 in SLC17A1 (OR 0.898 per SD methylation, p=1.0×10−4) and cg15445000 (OR 1.089, p=2.2×10−4) as PE-associated loci in FinnGen R10; only cg15445000 was additionally significant in the partly overlapping Tyrmi 2023 meta-GWAS (OR 1.059, p=0.0024), with lineage-specific genetically instrumented effects on monocyte, lymphocyte, and neutrophil counts.

**Conclusions.** Across expression, tissue, trajectory, and genetic layers, PE-associated cfDNA methylation differences behave as dynamic fingerprints of maternal leukocyte composition rather than placental trophoblast lesions, with a small genetically regulated subset linking methylation to immune and metabolic traits. cfDNA methylation signatures for PE should therefore be interpreted as readouts of maternal inflammatory and haematological state, and validated accordingly. These conclusions rest on a small, predominantly European-ancestry discovery cohort (64 matched samples) with triangulated but not externally validated evidence; clinical application requires large, independent, and ancestrally diverse validation cohorts.

---

## Background

Pre-eclampsia affects 2–8% of pregnancies and remains a leading cause of maternal and perinatal morbidity and mortality. The only definitive treatment is delivery, so prediction and early detection carry most of the clinical weight. Protein-based markers (sFlt-1/PlGF ratio, PlGF alone) have entered practice, but their performance in first-trimester screening leaves room, and they say little about mechanism [1,2].

Circulating cell-free DNA offers a different entry point. In pregnancy, maternal plasma cfDNA is a mixture: the fetal fraction, mostly placental in origin, sits on top of a much larger maternal background contributed largely by haematopoietic cells. Because methylation patterns are cell-type specific, cfDNA methylation has been used to map the tissue origins of circulating DNA and to track placental contribution non-invasively [3,4]. Several groups, including a large Leuven programme, have reported PE-associated methylation differences in maternal plasma [5], and the working assumption in much of that literature is that disease-relevant signals come from the placenta, the organ central to PE pathology.

The assumption deserves scrutiny. Two facts complicate it. First, PE is a systemic maternal inflammatory state: neutrophil counts rise, lymphocyte proportions fall, and the maternal haematological profile shifts measurably, starting well before clinical onset [6,7]. Since the maternal leukocyte pool contributes the majority of cfDNA molecules, a change in its composition would move plasma methylation levels without any change in the placenta itself. Second, methylation differences observed in bulk cfDNA are diluted mixtures; without explicit decomposition, attributing them to the placenta is a guess.

This question is not academic. If PE-associated cfDNA methylation differences originate in maternal leukocytes, they are readouts of the maternal immune state, informative about pathophysiology and possibly about early detection, but they must be validated against leukocyte-relevant biology, and claims about placental dysfunction should not be built on them. If they originate in the placenta, they complement protein markers as direct indices of the diseased organ.

We therefore took an existing cfRRBS cfDNA dataset of pre-eclamptic and control pregnancies (GSE282512) through a deliberately origin-focused pipeline: (i) DMR discovery with sensitivity analyses for collection tube and disease subtype; (ii) a classifier evaluation kept deliberately conservative; (iii) origin testing in four independent evidence layers, each with its own positive control: cross-cohort expression meta-analysis, direct DMR testing in placental tissue, longitudinal methylation trajectories with cell-type-resolved deconvolution in two external cohorts, and genetic architecture analysis (cis-mQTL coverage, Mendelian randomization against PE and against blood cell and organ traits). We report the result of that triangulation.

## Methods

### Study data

The discovery dataset (GSE282512) comprises 369 plasma cfDNA profiles generated by cell-free reduced-representation bisulfite sequencing (cfRRBS; Ghent University, Belgium; nf-core/methylseq with the --rrbs flag, hg38), from 279 pregnant women sampled at various points during pregnancy; 357 profiles passed QC (≥1×10⁶ CpGs with depth ≥1 after standard filtering; median sequencing depth 10.8×, range 7.6–25.5×; median 4.2×10⁶ CpGs measured per sample; median 56.1% of measured CpGs at ≥5× coverage). For DMR discovery we used a gestational-age-matched subcohort of 32 PE cases (21 early-onset, 10 late-onset, 1 with missing onset; 28 severe, 4 non-severe) and 32 controls (BMI median 25.2 vs 23.7 kg/m²; sampling GA median 32 weeks, range 16–36; predominantly European-ancestry: 28/32 PE and 31/32 controls White, with 1 Asian, 1 Black, 2 Other among cases). Blood was drawn into EDTA or PAXgene tubes; tube type and sequencing batch were recorded and adjusted for. The study was conducted under the original dataset's ethics approval; secondary analysis used de-identified public data.

External data were all public: three placental DNA methylation cohorts (GSE57767, 31 PE/14 controls; GSE73375, 19/17; GSE75196, 8/16; Illumina HumanMethylation450), an independent cfDNA WGBS cohort (GSE154378: 134 samples from 4 pregnant groups plus 7 non-pregnant women, with first-trimester to delivery and postpartum sampling), a maternal leukocyte methylation cohort of healthy pregnancies (GSE37722, 84 samples across nulligravid, early, middle, delivery, and postpartum time points; HumanMethylation27), eight expression cohorts (six placental, two maternal blood), GoDMC cis-mQTL summary statistics, FinnGen R10 PE and gestational hypertension GWAS, a PE meta-GWAS (Tyrmi 2023, 16,743 cases/280,081 controls), and blood cell and organ trait GWAS (Sakaue 2021 atlas; Keaton 2024 blood pressure) [8–13]. The non-invasive fetal fraction markers and tissue references followed Sun et al. [3].

### DMR discovery and sensitivity analyses

Per-sample CpG methylation was aggregated into coverage-weighted region-level beta values over fixed 1-kb windows merged into larger candidate regions when adjacent windows co-differentiated (final median region length 4.0 kb; ~73,000 autosomal regions tested; region definition scripts included in the analysis code). Differential methylation was tested with limma [15], adjusting for collection tube and batch, and corroborated with Wilcoxon rank-sum tests. Regions with FDR<0.05 and |Δβ|>0.05 in the primary model were retained as candidates. Sensitivity analyses repeated the model in EDTA-only samples (20 PE/25 controls), in early-onset PE versus matched controls, and tested group×tube interaction genome-wide. A PAXgene-only exploratory analysis assessed direction consistency. Direction consistency rates and Δβ shrinkage were computed for all 166 candidates.

### Classifier evaluation

An elastic net logistic regression was trained on region-level methylation of the 166 candidates, with feature selection, imputation, and lambda tuning inside each leave-one-out fold to avoid information leakage. Performance was summarised by AUC with bootstrap 95% CI, sensitivity and specificity at the Youden index, and compared against top-k fold-selected panels, subtype-specific models (early- vs late-onset PE), and a tube-batch-only negative control. Label permutation (n=20) provided an empirical null. Calibration (slope, Hosmer–Lemeshow), decision-curve analysis, and SHAP importance were computed on out-of-fold predictions. As a cross-modal check in early pregnancy, a signed z-sum expression score of DMR genes was evaluated in two first-trimester blood transcriptome cohorts (GSE85307, GSE86200) against a transcriptome-wide positive control and a random-gene null. External methylation-level validation was not feasible in the present study because no suitable independent cfRRBS dataset for pre-eclampsia was publicly available at the time of analysis; access to EGAS00001007071 has been requested for future validation work.

### Placental origin testing

Three layers were used. (i) Expression: DMR genes were tested for differential expression in a random-effects meta-analysis across six placental cohorts and two maternal blood cohorts. (ii) Tissue methylome: 1,030 HM450 CpGs within 127 of the 166 DMRs (hg38 coordinates lifted to hg19) were tested for PE-vs-control differences in each placental cohort (Welch t per CpG), pooled by per-cohort z averaging, and compared with (a) genome-wide background CpGs and (b) a positive control of the 500 background CpGs with the strongest placental PE signal; significance of DMR-set effects was assessed by size-matched permutation (12,000 draws). Direction agreement between placental and cfDNA effects was tested binomially. (iii) Regulatory annotation: DMRs were overlapped with Ensembl regulatory builds [21] (enhancer, promoter, CTCF, open chromatin, TF binding) and scanned for JASPAR2024 CORE motif enrichment [20] against GC- and length-matched background regions.

### Longitudinal trajectories and change-based deconvolution

In GSE154378, per-read marker methylation was assigned to the Sun 2015 marker panel (5,820 markers; bin identity reconstructed from the sorted coordinate list, validated by Spearman correlation of non-pregnant plasma against the blood reference, rho=0.899). Marker-set beta values per time point were tracked across gestation. Cell fraction changes were estimated with a change-based variant of reference-based deconvolution [14] (Δβ of marker sets regressed on reference differences between cell types, least squares, bootstrap 400 for 95% CIs), which cancels platform offsets by construction. This approach focuses on temporal changes in methylation rather than absolute cell-type fractions; it cancels constant platform-specific offsets but cannot correct for cell-type-independent methylation shifts over gestation, so fraction-change estimates should be read as relative contrasts between groups rather than calibrated absolute fractions. The same marker panel, lifted to hg38, was applied to GSE37722 leukocyte-layer data (±5 kb flanking probes), and fraction changes were contrasted between the cfDNA and leukocyte cohorts at each time point. Subgroup analyses in GSE154378 compared PE/gestational hypertension (hereafter PreX; n=8 patients, a combined and very small subgroup whose results are hypothesis-generating only), gestational diabetes (GDM, n=7), and chronic hypertension (cHTN, n=2) against normal pregnancy (n=9) at each time point, with GDM serving as an internal positive control for placental fraction.

### Genetic analyses

cis-mQTL coverage: CpGs within DMRs were intersected with GoDMC significant cis-mQTLs (p<10−5) and compared with non-DMR HM450 background CpGs (probe annotation per Zhou et al. [22]; Fisher exact test). MR: cis-meQTL instruments (p<10−5) for seven DMR CpGs with shared mQTLs were LD-clumped (1000 Genomes EUR phase 3 [19], r²<0.1), allele-harmonised against outcomes (palindromic variants resolved by allele frequency, ambiguous variants dropped), and tested with IVW (fixed effect; DerSimonian–Laird random effect when heterogeneity Q p<0.05), single-instrument Wald ratios, MR-Egger intercept [18], weighted median [17], leave-one-out, and F statistics. Outcomes were FinnGen R10 PE and gestational hypertension, and the Tyrmi 2023 PE meta-GWAS; both are European-ancestry GWAS, whereas the discovery cohort is predominantly but not exclusively of European ancestry, so ancestry mismatch is a limitation of the genetic layer; we therefore interpret the MR results with caution for non-European populations. Direction triangulation compared MR effect direction with the observed cfDNA Δβ direction. Organ-wide MR applied the same instrument pool to 16 blood cell and organ traits (12 blood cell counts, eGFR, SBP/DBP/pulse pressure) extracted from GWAS Catalog harmonised summaries by remote tabix query, with Benjamini–Hochberg correction over the 112 tests; the pipeline reproduced the FinnGen PE estimates from local files as a direction check.

### Statistics and software

All analyses were performed in R 4.6.1 (limma, WGCNA [16], EpiDISH, ggplot2) and Python 3.13; scripts are archived at [repository to be completed]. Multiple testing used Benjamini–Hochberg FDR unless stated. Ethics approval and consent: original GEO deposits; this secondary analysis required no new consent.

## Results

### Scope of the evidence base

Two structural features of this analysis should be stated before the results, because they shape how the evidence should be weighed. First, DMR discovery and classifier evaluation rest on a gestational-age-matched subcohort of 64 samples (32 PE, 32 controls) drawn from the 357 QC-passing profiles of GSE282512; the full-cohort and subgroup analyses were used for sensitivity and replication checks, but the discovery statistics themselves are small-sample statistics and should be read as such. Second, the evidence is deliberately distributed across cohorts. Origin testing draws on placental tissue, expression, longitudinal trajectory, and genetic data from datasets other than the discovery cohort, because the cfRRBS design of GSE282512 measures only about 27% of the reference cell-type markers at these loci, too few for a reliable in-cohort cell-composition deconvolution (the corresponding in-cohort analysis is reported below and is underpowered). The cell-composition evidence therefore comes from two external cohorts (GSE154378 for cfDNA, GSE37722 for the maternal leukocyte layer), where it is consistent across cohort, platform, and modality; residual cross-cohort heterogeneity is a limitation that we address in the Discussion.

### DMR discovery and robustness

In the gestational-age-matched subcohort, limma identified 166 candidate DMRs: 105 hypermethylated and 61 hypomethylated in PE (median length 4.0 kb; distribution across all autosomes plus chrX, with clustering on chr1, chr17, and chr19). Effect sizes were moderate (median |Δβ| 0.116, range up to 0.225; the largest absolute Δβ in the set was 0.225).

The candidates were robust to technical and clinical stratification. In EDTA-only samples, all 166 regions kept their discovery direction (100%; no flips), 63 replicated at FDR<0.05 and 73 at nominal p<0.05, with median Δβ shrinkage of only 5%. PAXgene-only samples (12 PE/7 controls), though underpowered, agreed in direction for 99.4% of candidates. The genome-wide group×tube interaction test produced zero significant regions. In early-onset PE versus matched controls, direction consistency was again 100%, with 43 regions at FDR<0.05 and 122 at nominal p<0.05, indicating that the signal is stronger, not weaker, in the early-onset phenotype. A co-methylation network built from the most variable regions placed 99 candidates into modules; the PE-associated blue module (r=0.49, p=4.5×10−5) was enriched for DMRs (36/371 regions, p=7.3×10−17) and correlated with NK and naïve B-cell fractions, an early hint that cell composition was involved.

### Classifier performance

The elastic net classifier on 166 DMR features reached a leave-one-out AUC of 0.932 (95% CI 0.866–0.998), sensitivity 0.844 and specificity 0.969 at the out-of-fold Youden index, with calibration slope 0.923 and Hosmer–Lemeshow p=0.056; decision curves showed net benefit over treat-all/treat-none across clinically relevant thresholds. All 20 label permutations produced lower AUCs (maximum 0.791, empirical p=0.048). A model using only tube and batch identifiers reached AUC 0.574, ruling out a technical artefact of ascertainment. Subtype models gave 0.836 (early-onset) and 0.991 (late-onset). We treat this as an internally validated candidate signature only: leave-one-out cross-validation in a sample of this size is known to produce optimistic estimates even with fold-internal tuning and a significant permutation null, so the AUC of 0.932 should be read as an upper-bound reading of internal discriminative performance rather than an expected out-of-sample value; the 64-sample subcohort cannot support external claims, and independent methylation-level external validation remains the decisive test. As a cross-modal check in early pregnancy, a signed z-sum expression score of DMR genes reached AUC 0.616 (95% CI 0.524–0.706) in GSE85307, above a transcriptome-wide positive control (0.532) and a random-gene null band (0.479, 95% CI 0.346–0.585), while the small GSE86200 early-pregnancy set (6 PE) was uninformative. The direction information carried by DMR genes is real but modest at the expression level.

### Placental origin testing: three negative layers with positive controls

**Expression.** Of 66 DMR genes measurable across at least two placental cohorts, none reached meta p<0.05 for differential expression in placenta, and the direction-consistency rate across cohorts did not exceed chance (0/64 genes with >50% same-sign logFC in ≥2 cohorts). Maternal blood expression meta-analysis was similarly null.

**Placental tissue methylome.** 1,030 HM450 CpGs within 127 DMRs were measurable across the three placental cohorts (105 samples). DMR CpGs carried small PE effects (mean pooled |z| = 0.346), smaller than the genomic background (0.444; KS p=2.6×10−16, reflecting a deficit of large effects rather than an excess), and size-matched permutation gave p=1.0: no enrichment of placental PE signal whatsoever. Region-level placental effects did not differ from zero (Wilcoxon p=0.229), agreed with the cfDNA direction in only 46% of regions (binomial p=0.478), and were uncorrelated with cfDNA Δβ (rho=−0.036, p=0.69). The same pipeline applied to the 500 background CpGs with the strongest placental PE signal produced mean |z| = 3.99, i.e. the DMR effect is 8.7% of what a genuinely placental signal looks like in these data.

**Regulatory annotation.** DMRs were enriched in enhancers (OR 1.47, FDR 0.037) and depleted in promoters (OR 0.48, FDR 3.4×10−4), a pattern expected of cell-type identity regions rather than actively transcribed disease-response genes. Motif scanning found FOS/AP-1 motifs enriched 15-fold (FDR 0.022), consistent with myeloid inflammatory programmes, and CGGBP1 and ARNT::HIF1A motifs depleted.

### Longitudinal trajectories: the cfDNA signal follows blood cells, with a pre-symptomatic subgroup signal (hypothesis-generating)

In the external cfDNA cohort, marker-bin identity was first validated: non-pregnant plasma methylation matched the blood-cell reference at Spearman rho=0.899. Placenta-informative markers then rose monotonically across gestation (mean rho +0.542; 341/386 markers positive, p=6.0×10−51; β 0.221 at non-pregnant baseline to 0.327 at delivery), replicating the known increase in placental cfDNA fraction. Marker-level trajectories correlated with the reference tissue profile in the expected order: placenta +0.475, and negative values for neutrophils (−0.339), T cells (−0.289), and B cells (−0.251), i.e. the rising placental fraction dilutes blood-cell DNA, exactly as the mixture model predicts.

Change-based deconvolution made the cell composition explicit. In cfDNA, the T-cell fraction fell across gestation (delivery Δf −0.156, 95% CI −0.221 to −0.096) while the placental fraction rose to +0.258. The same markers in the maternal leukocyte layer (a different platform, a different cohort, and no plasma at all) reproduced the T-cell decline (middle-pregnancy Δf −0.034, CI −0.069 to −0.006), and both cohorts showed a postpartum reversal (leukocyte-layer T cells +0.041, CI 0.028–0.061), which excludes a platform artefact and confirms that the composition shift is physiological and reversible.

Against normal pregnancy, the PE/gestational-hypertension subgroup (PreX; eight patients with confirmed PE or gestational hypertension, 40 samples) showed a distinct pattern: at delivery, the neutrophil fraction was higher by +0.063 (CI 0.024–0.098) while the placental fraction was lower by 0.126 (CI −0.150 to −0.103). The neutrophil shift was already present in the first trimester (+0.070, CI 0.012–0.129), before clinical disease. GDM pregnancies showed the opposite placental pattern (delivery +0.067, CI 0.042–0.093), an internal positive control confirming that the placental fraction estimate responds when it should. The two chronic hypertension samples showed no neutrophil shift. We caution that the PreX subgroup comprises eight patients, so these subgroup contrasts are hypothesis-generating and require confirmation in larger independent cohorts before any clinical interpretation. Within the discovery cohort, region-level reference-based deconvolution was largely null (only NK and naïve CD4 T fractions reached significance in a 12-cell-type reference), which we attribute to the ~27% marker measurability of cfRRBS at these loci and the modest sample size; the composition evidence therefore rests on the external cohorts, where it is consistent across cohort, platform, and modality.

### Genetic architecture: weak cis-regulation overall, with two PE loci

Only 28.8% of DMR CpGs carried a significant cis-mQTL, compared with 37.7% of background CpGs (OR 0.669, p=2.1×10−9). Hypermethylated DMRs, which dominate the signature, were least genetically regulated (22.2% vs 50.0% for hypomethylated DMRs, OR 0.285, p=4.6×10−16). A dynamically regulated marker is what one expects of a composition fingerprint; a genetically hardwired one is not.

Where cis-mQTLs existed, they permitted Mendelian randomization. Seven DMR CpGs across three loci had usable instrument pools (4–58 SNPs per CpG after clumping; mean F 52–752). Two CpGs produced BH-significant IVW effects on PE in FinnGen R10: cg19490609 in an intron of SLC17A1 (OR 0.898 per SD increase in methylation, 95% CI 0.850–0.948, p=1.0×10−4), and cg15445000 at an intergenic chr17 locus (OR 1.089, CI 1.041–1.139, p=2.2×10−4). Of these two, only cg15445000 was additionally significant in the Tyrmi 2023 meta-GWAS (OR 1.059, CI 1.020–1.098, p=0.0024); cg19490609 did not replicate there (OR 0.980, CI 0.933–1.030, p=0.43). Because the Tyrmi meta-GWAS partly overlaps FinnGen (FinnGen R6 included), even the cg15445000 result should be read as triangulation rather than fully independent replication. Leave-one-out analyses were stable, and an approximate HEIDI test showed no heterogeneity, but MR-Egger intercepts were significant (indicating possible directional pleiotropy) and weighted median estimates did not reach 0.05; we therefore read these results as triangulating evidence rather than definitive causal effects, and as fully compatible with a non-causal interpretation of most DMRs. Direction comparison with the observed cfDNA Δβ was instructive: for cg19490609 the MR direction matched the observed hypomethylation in PE (consistent), while for cg15445000 the genetically driven effect ran opposite to the disease-associated change, the pattern expected of a marker that follows rather than causes disease.

### Organ-wide MR: methylation at these loci moves blood cell counts

Applying the same instruments to 16 blood cell and organ traits (112 tests) produced 14 BH-significant results with a coherent lineage structure. Methylation at cg19490609 (SLC17A1, hypomethylated in PE) increased monocyte (+0.049 SD), lymphocyte (+0.041 SD), eosinophil (+0.035 SD), and total WBC counts (+0.040 SD, all q<0.015). Hypermethylated-DMR CpGs moved counts downward: cg24308560 and cg27360567 each lowered WBC and neutrophil counts (q≤0.012), and cg24308560 raised eGFR (+0.017 SD, q=0.045). cg18129748 raised platelet count (+0.041 SD, q=0.014) and lowered basophil count. cg25753631 raised systolic blood pressure by 0.27 mmHg per SD methylation (q=0.05). The pipeline reproduced the FinnGen PE estimates from local summary files (cg19490609 OR 0.926 vs 0.898; cg15445000 OR 1.034 vs 1.089), confirming the direction of the remote extraction. In these analyses, genetically instrumented methylation changes at these DMR CpGs moved the same blood cell counts whose composition shifts in PE; given the pleiotropy caveats above, we interpret this as supporting, not proving, a haematological link for the genetically regulated subset.

## Discussion

### Principal finding

Four independent evidence layers converge on one answer to the origin question, for the 166 regions identified in this dataset. These DMRs are not differentially expressed in placenta, not differentially methylated in placental tissue (with an effect distribution shifted toward zero, and a positive control proving the placenta cohorts could detect real signals), enriched for enhancer and FOS/AP-1 regulatory grammar, tracked blood-cell rather than placental trajectories across gestation in an external cohort, showed a neutrophil-dominant composition shift in PE that is already present in the first trimester, are only weakly genetically regulated as a set, and, where instruments exist, are connected to blood cell counts through genetically instrumented effects rather than behaving as placental disease output. The most economical reading is that the majority of these 166 PE-associated DMRs derived from cfRRBS data behave as a fingerprint of maternal leukocyte composition. Two scope caveats follow. First, the conclusion applies to the DMRs we identified in maternal plasma; it does not claim that all PE-associated cfDNA methylation differences anywhere in the genome arise from leukocytes. Second, a minority of the DMR set is under genetic regulation and behaves differently, and those regions are discussed separately below.

### Comparison with the placental-origin assumption

This reading sits against the grain of the cfDNA methylation literature in PE, which has largely assumed placental origin, sometimes by construction (fetal-fraction or placenta-uniqueness filtering). Our placenta testing was designed to be able to detect placental signal and failed to find it in these regions, at 8.7% of positive-control strength. That does not mean placental methylation changes do not exist in PE; placenta clearly carries PE methylation differences elsewhere in the genome (our own background top-500 set), and fetal-fraction changes are real in PE. It means that the regions that differ in maternal plasma are, in the main, not those regions. For biomarker work, the distinction is consequential: a maternal-leukocyte fingerprint measured non-invasively is a readout of the systemic inflammatory state that precedes and accompanies PE, not of trophoblast injury, and its dynamics (for example the postpartum reversal we observed) will follow maternal haematology, not placental involution.

### The first-trimester neutrophil shift

The subgroup analysis carries the most clinically interesting result: in pregnancies that would develop PE or gestational hypertension, the neutrophil-derived cfDNA fraction was already elevated in the first trimester, before any clinical manifestation, and this occurred while the placental fraction estimate was flat rather than elevated. Neutrophil activation and neutrophil extracellular trap formation are documented early in PE pathophysiology [6,7]; our data suggest a detectable methylation correlate in circulation that predates disease. Three qualifications keep this in proportion. First, the PreX subgroup contains eight patients, so the estimate is imprecise and the first-trimester contrast in particular rests on a small number of samples. Second, the trajectory cohort and the discovery cohort are different studies, on different platforms, and the subgroup classification in GSE154378 combines PE with gestational hypertension. Third, no threshold or prediction claim follows from these numbers. We therefore present this as a hypothesis-generating observation: first-trimester cfDNA methylation can carry pre-symptomatic information through the maternal compartment, and a dedicated first-trimester cohort with prospective PE outcomes is required to test it.

### The genetic subset and its biology

The minority of DMRs under cis genetic regulation behave differently, and the MR results should be kept separate from the composition story. cg19490609 lies in SLC17A1, which encodes the renal urate transporter URAT1; hyperuricaemia is a classic PE accompaniment, and this locus ties cfDNA methylation to a urate-handling gene with a consistent causal direction toward PE risk in FinnGen R10, an effect that did not reach significance in the Tyrmi meta-GWAS. cg15445000, the only one of the two significant loci that was also significant in the partly overlapping Tyrmi 2023 meta-GWAS, shows the opposite MR-vs-observation pattern, compatible with reverse causation or with the disease-associated methylation change being a consequence of the leukocyte shift. The organ-wide MR then showed that instrumented methylation at these loci moves monocyte, lymphocyte, neutrophil, and platelet counts, and systolic pressure, in lineage-specific directions. The sensitivity caveats are real: MR-Egger intercepts flagged directional pleiotropy and weighted median estimates fell short of significance, so we treat the genetic layer as triangulation, not proof. What the genetic layer adds is an explanation of why some regions are methylation-quantitative-trait loci for haematological traits, and why the DMR set as a whole is not: a composition fingerprint need not be genetically fixed, but the machinery that sets methylation at these loci does interact with blood cell production.

### Strengths and limitations

The study's strength is the layered design with positive controls at each level: expression meta-analysis, tissue methylome testing with a demonstrated detection ceiling, trajectory mapping with an independent validation of the marker panel (rho=0.899), change-based deconvolution that cancels platform offsets, and MR with pipeline reproduction checks. Its limitations are equally concrete. The discovery subcohort is 64 samples, and the classifier is internally validated only; leave-one-out cross-validation is itself an optimistic estimator at this sample size, so the reported AUC of 0.932 is best read as an upper bound of internal performance, and neither the permutation test nor fold-internal tuning substitutes for validation in an independent sample. The early-pregnancy expression proxy confirmed only weak directional information (AUC 0.616) and the external early-pregnancy set had 6 cases. In-cohort deconvolution was underpowered due to cfRRBS marker coverage, so the cell-composition evidence rests on the external cohorts rather than on the discovery cohort itself; evidence for leukocyte composition shifts draws from independent external cohorts, and residual inter-cohort heterogeneity cannot be fully excluded. The trajectory and subgroup cohorts are small (8 PreX patients, a combined PE and gestational hypertension subgroup treated as hypothesis-generating; 2 cHTN patients). The genetic analyses rely on European-ancestry GWAS (FinnGen, Tyrmi, Sakaue EUR, Keaton EUR), while the discovery cohort is predominantly but not exclusively of European ancestry (28/32 PE and 31/32 controls White), and the two PE outcome datasets partly overlap. Liver enzymes, CRP, and urate itself could not be included in organ-wide MR for lack of indexed summary statistics. Finally, the definitive external validation, a first-trimester cfDNA methylome cohort with early-onset PE (EGAS00001007071), is under access application; until it reports, the classifier stays labelled a candidate signature. Should access not be granted, this study necessarily remains at the mechanism-exploration stage, with an internally validated candidate signature and no external validation; we state this contingency explicitly so that the scope of the present claims is not read as broader than the data allow.

Three directions follow from these limitations. First, completion of access to EGAS00001007071 would allow the cfDNA methylation signature to be tested in an independent, phenotypically different cohort (early-onset PE) and would move the classifier from internally validated candidate to externally validated. Second, the first-trimester neutrophil finding requires a dedicated, adequately powered first-trimester cohort with prospective PE outcomes; the present subgroup can only generate, not test, that hypothesis. Third, multi-ancestry GWAS and a discovery cohort with verified ancestry would be needed to determine whether the genetic layer transfers across populations; the current analysis is European-centric and its transferability is unknown.

### Conclusions

Pre-eclampsia-associated differences in maternal plasma cfDNA methylation, at least in this dataset and at these 166 regions, behave as a dynamic fingerprint of maternal leukocyte composition: invisible in placental tissue, tracking blood-cell trajectories, shifted toward neutrophils before clinical onset, weakly genetically regulated as a set, and linked to blood cell counts through genetically instrumented effects where instruments exist. The scope of the claim is deliberately narrow: it concerns the DMRs identified here, it does not exclude the existence of placental-origin cfDNA methylation signals in PE elsewhere in the genome, and it carries no clinical prediction claim until external validation in large, independent, and ancestrally diverse cohorts is available. A non-invasive readout of the maternal inflammatory state in pregnancy is valuable, for mechanism and possibly for early detection, but it should be designed, validated, and interpreted as what it is, not as a placental biopsy in a tube.

## Declarations

**Ethics approval and consent to participate.** This study is a secondary analysis of de-identified public data deposited in GEO, EGA, and GWAS Catalog. The original cohorts were collected under the ethics approvals reported in their primary publications (GSE282512: ethics approval of Ghent University Hospital; GSE154378, GSE37722, and the placental and expression cohorts: as reported in the corresponding GEO deposits). No new human data were collected, and no additional ethics approval was required for this analysis. The external validation cohort EGAS00001007071 will be analysed under its data access agreement upon approval.

**Consent for publication.** Not applicable.

**Availability of data and materials.** All source data are public: GSE282512, GSE154378, GSE37722, GSE57767, GSE73375, GSE75196, GSE85307, GSE86200, GSE98224 and additional expression cohorts via GEO; GoDMC via its data access portal; FinnGen R10 via the FinnGen portal; Tyrmi 2023 and blood cell/organ trait GWAS via GWAS Catalog (GCST90269903; GCST90002381, GCST90002393, GCST90018946–GCST90018978, GCST90019506; GCST90310294–96). The external cfDNA validation cohort EGAS00001007071 is available through EGA under a data access agreement (application in progress); if access is not granted, this study remains a mechanism-oriented analysis without external validation of the candidate signature, and no external validation claims are made. Analysis code and intermediate result tables: [repository/Zenodo DOI to be added at submission]; the complete analysis pipeline and all scripts generating the reported results will be made available to reviewers during peer review and publicly archived upon acceptance.

**Competing interests.** The authors declare none.

**Funding.** [To be completed at submission.]

**Authors' contributions.** [To be completed at submission: W.X. conceived the analysis, performed the bioinformatic and statistical analyses, and wrote the manuscript. Co-authors to confirm contributions.]

**Acknowledgements.** [To be completed at submission.]

## References

1. Levine RJ, Maynard SE, Qian C, et al. Circulating angiogenic factors and the risk of preeclampsia. N Engl J Med. 2004;350:672–83.
2. Zeisler H, Llurba E, Chantraine F, et al. Predictive value of the sFlt-1:PlGF ratio in women with suspected preeclampsia. N Engl J Med. 2016;374:13–22.
3. Sun K, Jiang P, Chan KCA, et al. Plasma DNA tissue mapping by genome-wide methylation sequencing for prenatal diagnosis, immunotherapy, and cancers. Proc Natl Acad Sci USA. 2015;112:E5503–12.
4. Lun FMF, Chiu RWK, Sun K, et al. Noninvasive prenatal methylomic analysis by genomewide bisulfite sequencing of maternal plasma DNA. Clin Chem. 2013;59:663–72.
5. De Borre M, Che H, Yu Q, et al. Cell-free DNA methylome analysis for early preeclampsia prediction. Nat Med. 2023;29(9):2206–15.
6. Redman CW, Sacks GP, Sargent IL. Preeclampsia: an excessive maternal inflammatory response to pregnancy. Am J Obstet Gynecol. 1999;180:499–506.
7. Gupta AK, Hasler P, Holzgreve W, Gebhardt S, Hahn S. Induction of neutrophil extracellular DNA lattices by placental microparticles and IL-8 and their presence in preeclampsia. Hum Immunol. 2005;66:1146–54.
8. Barrett T, Wilhite SE, Ledoux P, et al. NCBI GEO: archive for functional genomics data sets. Nucleic Acids Res. 2013;41:D991–5.
9. Min JL, Hemani G, Hannon E, et al. Genomic and phenotypic insights from an atlas of genetic effects on DNA methylation. Nat Genet. 2021;53:1311–21.
10. Kurki MI, Karjalainen J, Palta P, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature. 2023;613:508–18.
11. Tyrmi JS, Kaartokallio T, Lokki AI, et al. Genetic risk factors associated with preeclampsia and hypertensive disorders of pregnancy. JAMA Cardiol. 2023;8(7):674–83.
12. Sakaue S, Kanai M, Tanigawa Y, et al. A cross-population atlas of genetic associations for 220 human phenotypes. Nat Genet. 2021;53(10):1415–24.
13. Keaton JM, Kamali Z, Xie T, et al. Genome-wide analysis in over 1 million individuals of European ancestry yields improved polygenic risk scores for blood pressure traits. Nat Genet. 2024;56(5):778–91.
14. Teschendorff AE, Breeze CE, Zheng SC, Beck S. A comparison of reference-based methods for estimating cell-type proportions in DNA methylation data. PLoS One. 2017;12:e0189534.
15. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43:e47.
16. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008;9:559.
17. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent estimation in Mendelian randomization with some invalid instruments using a weighted median estimator. Genet Epidemiol. 2016;40:304–14.
18. Bowden J, Del Greco MF, Minelli C, et al. Assessing the suitability of summary data for Mendelian randomization analyses using MR-Egger regression: the role of the I2 statistic. Int J Epidemiol. 2016;45:1961–74.
19. 1000 Genomes Project Consortium, Auton A, Brooks LD, et al. A global reference for human genetic variation. Nature. 2015;526:68–74.
20. Fornes O, Castro-Mondragon JA, Khan A, et al. JASPAR 2024: the 9th release of the open-access database of transcription factor binding profiles. Nucleic Acids Res. 2024;52:D174–82.
21. Zerbino DR, Wilder SP, Johnson N, et al. The Ensembl regulatory build. Genome Biol. 2015;16:56.
22. Zhou W, Laird PW, Shen H. Comprehensive characterization, annotation and innovative analysis of Infinium methylation beadchip arrays. Nucleic Acids Res. 2017;45:e22.

## Figure legends (figures available in project `figures/` directory)

**Figure 1.** Study design and evidence layers. The discovery cohort (GSE282512) yielded 166 candidate DMRs and an internally validated classifier; the cellular origin of the signal was then tested in five independent evidence layers, each with its own positive control, triangulating between placental trophoblast and maternal leukocyte origins (`Figure1_study_design.png`).

**Figure 2.** DMR discovery and robustness. (A) Volcano plot of region-level methylation differences (PE vs controls; `dmr_volcano.png`). (B) Manhattan plot (`dmr_manhattan.png`). (C) EDTA-only sensitivity: Δβ concordance (`dmr_sensitivity_scatter.png`). (D) Top DMR site-level methylation (`dmr_site_top6.png`).

**Figure 3.** Classifier evaluation. (A) LOOCV ROC with permutation nulls (`roc_loocv.png`). (B) Calibration (`calibration.png`). (C) Decision curve analysis (`dca.png`). (D) SHAP feature importance (`shap_importance.png`).

**Figure 4.** Placental origin testing. DMR CpG effect distribution in placental tissue versus background and positive-control sets, permutation null, and cfDNA–placenta direction agreement (`GSE282512_placenta_direct.png`).

**Figure 5.** Longitudinal trajectories in the external cfDNA cohort. (A) Placenta-informative marker methylation across gestation (`GSE154378_typeI_trajectory.png`). (B) DMR-overlapping bins (`GSE154378_dmr_trajectory.png`). (C) Cell-fraction changes in cfDNA and maternal leukocyte cohorts with postpartum reversal (`GSE37722_GSE154378_cell_fraction.png`).

**Figure 6.** Subgroup change-based deconvolution (hypothesis-generating; the PreX subgroup comprises eight patients and combines pre-eclampsia with gestational hypertension). Neutrophil, T-cell, and placental fraction contrasts (PreX, GDM, cHTN vs normal pregnancy) across gestation (`GSE154378_subgroup_deltaf_forest.png`; `GSE154378_subgroup_deltaf_trajectory.png`). These subgroup contrasts require confirmation in larger independent cohorts and support no clinical prediction claim.

**Figure 7.** Mendelian randomization for PE and gestational hypertension. IVW forest plots with leave-one-out and sensitivity methods (`GSE282512_mr_forest.png`).

**Figure 8.** Organ-wide MR. Forest plot of the 25 strongest cis-meQTL-instrumented effects on blood cell, eGFR, and blood pressure traits (`GSE282512_organ_mr_forest.png`).

## Tables

**Table 1.** Datasets and cohorts.

| Dataset | Platform / data type | Sample composition | Role in this study |
|---|---|---|---|
| GSE282512 (discovery; Ghent University, Belgium) | Plasma cfDNA reduced-representation bisulfite sequencing (cfRRBS, hg38; nf-core/methylseq --rrbs) | 369 profiles from 279 women; 357 passed QC (median depth 10.8×, range 7.6–25.5×; median 4.2×10⁶ CpGs/sample); gestational-age-matched subcohort: 32 PE (21 early-onset, 10 late-onset, 1 onset unknown; 28 severe/4 non-severe) / 32 controls; EDTA or PAXgene tubes | DMR discovery, tube/subtype sensitivity, classifier |
| GSE57767 / GSE73375 / GSE75196 (placenta) | HM450 tissue methylome | PE/controls: 31/14, 19/17, 8/16 (105 samples total) | Direct testing of DMR CpGs in placental tissue |
| GSE154378 (external cfDNA) | WGBS with Sun 2015 marker deconvolution | 134 samples: normal pregnancy 9 patients (44 samples), PE/gestational hypertension 8 (40), GDM 7 (33), chronic hypertension 2 (10), non-pregnant 7 | Gestational trajectories, change-based deconvolution, subgroup contrasts |
| GSE37722 (maternal leukocytes) | HM27, leukocyte-layer DNA | 84 samples, healthy pregnancies at nulligravid, early, middle, delivery and postpartum time points | Leukocyte-layer concordance of cell-fraction changes |
| Expression cohorts (8) | RNA-seq / microarray | 6 placental, 2 maternal blood | DMR-gene differential expression meta-analysis |
| GSE85307 / GSE86200 | First-trimester maternal blood transcriptome | Early-pregnancy samples | Cross-modal expression proxy of DMR genes |
| GoDMC | cis-mQTL summary statistics (n=32,000) | — | mQTL coverage analysis, MR instruments |
| FinnGen R10 | GWAS | PE 7,965/211,852; gestational hypertension 9,535/211,957 | MR outcomes |
| Tyrmi 2023 meta-GWAS (GCST90269903) | GWAS | PE 16,743/280,081 | MR replication outcome |
| Sakaue 2021 atlas; Keaton 2024 | GWAS | 12 blood cell traits + eGFR (EUR ~350k); SBP/DBP/PP (EUR 1,028,980) | Organ-wide MR outcomes |

**Discovery subcohort clinical and technical characteristics (Table 1, continued).** PE cases n=32, controls n=32. BMI median (range): 25.2 (20.3–39.6) vs 23.7 (16.5–33.2) kg/m². Gestational age at sampling: median 32 weeks (16–36) in both groups (matched within ±2 weeks). Ancestry (as recorded in the GEO metadata): White 28/32 cases, 31/32 controls; among cases, with 1 Asian, 1 Black, and 2 Other; 1 control missing. Collection tube: EDTA 20/32 cases, 25/32 controls; PAXgene 12/32 cases, 7/32 controls. Sequencing batches 2–10 (case/control distributions reported in the analysis code). Fetal sex recorded in the metadata but not used in this analysis.

**Table 2.** Classifier performance (leave-one-out cross-validation; feature selection, imputation and lambda tuning inside each fold).

| Model | n | AUC | 95% CI | Sensitivity | Specificity | Accuracy |
|---|---|---|---|---|---|---|
| Elastic net, 166 DMR features | 64 | 0.932 | 0.866–0.998 | 0.844 | 0.969 | 0.906 |
| Top-10 fold-selected panel | 64 | 0.871 | 0.777–0.965 | 0.781 | 0.875 | 0.828 |
| Top-20 fold-selected panel | 64 | 0.891 | 0.799–0.982 | 0.781 | 1.000 | 0.891 |
| Elastic net, early-onset PE vs controls | 53 | 0.836 | 0.714–0.958 | — | — | — |
| Elastic net, late-onset PE vs controls | 42 | 0.991 | 0.970–1.000 | — | — | — |
| Tube and batch only (negative control) | 64 | 0.574 | — | — | — | — |

Label permutation (n=20): null AUC mean 0.632, maximum 0.791; observed AUC empirical p=0.048. Calibration slope 0.923, Hosmer–Lemeshow p=0.056. Sensitivity/specificity at the out-of-fold Youden index.

**Table 3.** Placental tissue testing of DMR CpGs (three cohorts, 105 samples; PE vs control; pooled per-cohort z averaging).

| CpG set | n CpGs | Mean \|pooled z\| | Comparison |
|---|---|---|---|
| DMR CpGs (127/166 DMRs covered) | 1,030 | 0.346 | 57.8th percentile of background |
| Genome-wide background | 484,547 | 0.444 | KS p=2.6×10⁻¹⁶; DMR set shifted toward smaller effects |
| Positive control: 500 background CpGs with strongest placental PE signal | 500 | 3.99 | 11.6-fold above DMR set; KS p≈0 vs background |

Size-matched permutation (12,000 draws): p=1.0 (no enrichment of placental PE signal in DMR CpGs). Region level (n=127): Wilcoxon vs 0, p=0.229; direction agreement with cfDNA 59/127 (46%), binomial p=0.478; cfDNA Δβ vs placental Δβ Spearman rho=−0.036 (p=0.69). Mean placental Δβ: hyper-DMR subset −0.0035 (n=84), hypo-DMR subset −0.0019 (n=43).

**Table 4.** Two-sample Mendelian randomization for pre-eclampsia (PE) and gestational hypertension (GH). OR per 1-SD increase in CpG methylation (IVW; DerSimonian–Laird random effects where heterogeneity Q p<0.05).

| CpG (locus) | Outcome | k IVs | OR (95% CI) | p | p (BH) | Q p | MR-Egger p | Weighted median p | Mean F |
|---|---|---|---|---|---|---|---|---|---|
| cg18129748 (chr3:49.90 Mb, promoter, hyper-DMR) | GH, FinnGen R10 | 5 | 0.859 (0.728–1.014) | 0.072 | 0.160 | 0.065 | 0.420 | 0.58 | 79 |
| cg24308560 (chr3:49.90 Mb, promoter, hyper-DMR) | GH, FinnGen R10 | 29 | 0.965 (0.919–1.012) | 0.145 | 0.170 | 2.8×10⁻⁷ | 0.240 | 0.79 | 577 |
| cg27360567 (chr3:49.90 Mb, promoter, hyper-DMR) | GH, FinnGen R10 | 20 | 0.924 (0.857–0.996) | 0.040 | 0.140 | 1.4×10⁻⁶ | 0.100 | 0.95 | 355 |
| cg18129748 | PE, FinnGen R10 | 5 | 0.864 (0.723–1.032) | 0.107 | 0.249 | 0.980 | 0.0037 | 1.5×10⁻⁶ | 79 |
| cg24308560 | PE, FinnGen R10 | 29 | 0.991 (0.964–1.018) | 0.501 | 0.501 | 0.420 | 0.450 | 0.25 | 577 |
| cg27360567 | PE, FinnGen R10 | 20 | 0.981 (0.943–1.021) | 0.352 | 0.501 | 0.170 | 0.470 | 0.27 | 355 |
| cg18129748 | PE, Tyrmi 2023 | 5 | 0.977 (0.842–1.134) | 0.761 | 0.990 | 0.970 | 0.550 | 8.3×10⁻⁴ | 79 |
| cg24308560 | PE, Tyrmi 2023 | 28 | 1.006 (0.983–1.029) | 0.642 | 0.990 | 0.290 | 0.680 | 0.83 | 589 |
| cg27360567 | PE, Tyrmi 2023 | 19 | 1.000 (0.965–1.035) | 0.990 | 0.990 | 0.740 | 0.990 | 0.82 | 367 |
| cg00387872 (chr6:25.7 Mb, SLC17A1 gene body, hypo-DMR) | GH, FinnGen R10 | 58 | 1.019 (0.996–1.043) | 0.098 | 0.160 | 0.050 | 0.110 | 0.48 | 733 |
| cg19490609 (chr6:25.8 Mb, SLC17A1 gene body, hypo-DMR) | GH, FinnGen R10 | 22 | 0.973 (0.902–1.050) | 0.486 | 0.486 | 0.013 | 0.260 | 0.49 | 179 |
| cg25753631 (chr6:25.7 Mb, SLC17A1 gene body, hypo-DMR) | GH, FinnGen R10 | 28 | 1.027 (0.993–1.063) | 0.115 | 0.160 | 0.068 | 0.045 | 0.69 | 341 |
| cg00387872 | PE, FinnGen R10 | 58 | 0.989 (0.962–1.017) | 0.433 | 0.501 | 0.0024 | 0.650 | 0.058 | 733 |
| cg19490609 | PE, FinnGen R10 | 22 | 0.898 (0.850–0.948) | 1.0×10⁻⁴ | 7.3×10⁻⁴ | 0.058 | 0.019 | 0.13 | 179 |
| cg25753631 | PE, FinnGen R10 | 28 | 1.020 (0.969–1.074) | 0.446 | 0.501 | 0.036 | 0.610 | 0.46 | 341 |
| cg00387872 | PE, Tyrmi 2023 | 56 | 1.001 (0.979–1.023) | 0.956 | 0.990 | 0.035 | 0.980 | 0.99 | 752 |
| cg19490609 | PE, Tyrmi 2023 | 22 | 0.980 (0.933–1.030) | 0.430 | 0.990 | 0.380 | 0.410 | 0.29 | 158 |
| cg25753631 | PE, Tyrmi 2023 | 28 | 1.021 (0.989–1.053) | 0.195 | 0.682 | 0.330 | 0.580 | 0.32 | 341 |
| cg15445000 (chr17:39.45 Mb, promoter, hypo-DMR) | GH, FinnGen R10 | 14 | 1.048 (1.006–1.093) | 0.026 | 0.140 | 0.160 | 0.093 | 0.23 | 436 |
| cg15445000 | PE, FinnGen R10 | 14 | 1.089 (1.041–1.139) | 2.2×10⁻⁴ | 7.6×10⁻⁴ | 0.140 | 0.010 | 0.15 | 436 |
| cg15445000 | PE, Tyrmi 2023 | 15 | 1.059 (1.020–1.098) | 0.0024 | 0.017 | 0.250 | 0.017 | 0.31 | 411 |

Instruments: GoDMC cis-meQTLs (p<10⁻⁵), LD-clumped at r²<0.1 (1000 Genomes EUR phase 3), allele-harmonised with palindrome resolution by allele frequency. Leave-one-out analyses stable throughout; approximate HEIDI heterogeneity test p=1 for all significant estimates.

**Note on outcome independence.** The Tyrmi 2023 meta-GWAS **partly overlaps FinnGen (FinnGen R6 is included)** and is therefore read as triangulation, **not fully independent replication**. Of the two BH-significant FinnGen R10 PE loci, **only cg15445000 reached significance in Tyrmi 2023** (OR 1.059, p=0.0024); **cg19490609 did not replicate** (OR 0.980, CI 0.933–1.030, p=0.43). All MR results should be interpreted as triangulating supportive evidence rather than definitive causal inference due to potential horizontal pleiotropy.

**Table 5.** Organ-wide MR: BH-significant (q<0.05) effects of 1-SD genetically instrumented methylation change on blood cell and organ traits (112 tests corrected). Blood cell and eGFR effects in SD units; blood pressure in mmHg.

| CpG (locus, DMR direction) | Trait | k IVs | Beta (95% CI) | p | q | Mean F |
|---|---|---|---|---|---|---|
| cg24308560 (chr3, hyper) | White blood cell count | 29 | −0.018 (−0.027, −0.009) | 8.1×10⁻⁵ | 0.009 | 577 |
| cg24308560 | Neutrophil count | 29 | −0.015 (−0.024, −0.007) | 4.7×10⁻⁴ | 0.012 | 577 |
| cg27360567 (chr3, hyper) | White blood cell count | 20 | −0.027 (−0.042, −0.012) | 4.9×10⁻⁴ | 0.012 | 355 |
| cg27360567 | Neutrophil count | 20 | −0.024 (−0.038, −0.011) | 5.0×10⁻⁴ | 0.012 | 355 |
| cg19490609 (SLC17A1, hypo) | Monocyte count | 24 | +0.049 (0.021, 0.077) | 5.3×10⁻⁴ | 0.012 | 193 |
| cg19490609 | Lymphocyte count | 24 | +0.041 (0.017, 0.065) | 6.6×10⁻⁴ | 0.012 | 193 |
| cg19490609 | White blood cell count | 24 | +0.040 (0.017, 0.064) | 7.3×10⁻⁴ | 0.012 | 193 |
| cg18129748 (chr3, hyper) | Platelet count | 4 | +0.041 (0.017, 0.066) | 9.8×10⁻⁴ | 0.014 | 94 |
| cg19490609 | Eosinophil count | 24 | +0.035 (0.014, 0.056) | 1.2×10⁻³ | 0.015 | 193 |
| cg24308560 | Lymphocyte count | 29 | −0.015 (−0.025, −0.005) | 2.2×10⁻³ | 0.023 | 577 |
| cg18129748 | Basophil count | 4 | −0.035 (−0.057, −0.012) | 2.3×10⁻³ | 0.023 | 94 |
| cg00387872 (SLC17A1, hypo) | Red blood cell count | 58 | −0.008 (−0.014, −0.003) | 4.4×10⁻³ | 0.041 | 737 |
| cg24308560 | eGFR | 29 | +0.017 (0.005, 0.030) | 5.3×10⁻³ | 0.045 | 577 |
| cg25753631 (SLC17A1, hypo) | Systolic blood pressure | 29 | +0.27 mmHg (0.08, 0.47) | 6.3×10⁻³ | 0.050 | 333 |

Traits: Sakaue 2021 atlas (GCST90002381, GCST90002393, GCST90018946–78, GCST90019506; EUR ~350k); Keaton 2024 (GCST90310294–96; EUR 1,028,980). The same pipeline applied to local FinnGen PE/GH files reproduced the estimates of Table 4 (cg19490609 PE OR 0.926 vs 0.898; cg15445000 PE 1.034 vs 1.089; GH 1.024 vs 1.048).
