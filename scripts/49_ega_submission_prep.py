# -*- coding: utf-8 -*-
"""
Task #47: EGA submission material preparation (no real submission).
Outputs to results/ega_submission/:
  1. sample_metadata.tsv          -- 357 QC-passed samples, EGA-style fields + subcohort flag/QC
  2. data_file_manifest.tsv       -- real .cov.gz files mapped to samples (from raw dir)
  3. experiment_metadata.tsv      -- experiment/analysis metadata (WGBS, GRCh38)
  4. analysis_description.md      -- analysis workflow description for EGA DAC
  5. submission_checklist.md      -- pre-submission checklist
"""
import os, glob, csv, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "ega_submission")
RAW = os.path.join(ROOT, "data", "geo_methylation", "GSE282512_raw")
os.makedirs(OUT, exist_ok=True)

def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

samples = read_csv(os.path.join(RES, "GSE282512_samples_clean.csv"))
sub = read_csv(os.path.join(RES, "GSE282512_subcohort.csv"))
sub_ids = {r["sample_id"] for r in sub}
# QC pass criterion (project convention): cov5x_pct >= 50 -> 357/369 pass
qc = {r["sample_id"]: r for r in read_csv(os.path.join(RES, "GSE282512_qc_summary.csv"))}
def qc_pass(sid):
    r = qc.get(sid)
    if r is None:
        return "no"
    try:
        return "yes" if float(r["cov5x_pct"]) >= 50 else "no"
    except (ValueError, KeyError):
        return "unknown"

# ---- 1. sample_metadata.tsv (EGA-style) ----
rows = []
for r in samples:
    sid = r["sample_id"]
    in_sub = sid in sub_ids
    # phenotype mapping
    dis = r["sm_disease"]
    if dis.lower().startswith("pre") or (r["group_status"] or "") == "Case" or "eclamps" in dis.lower():
        phenotype = "pre-eclampsia" if "eclamp" in dis.lower() else dis.strip()
    else:
        phenotype = "control"
    if r["category"]:
        phenotype = r["category"]
    ga = r["ga_sampling"] if r["ga_sampling"] else r["ga_weeks"]
    rows.append({
        "sample_alias": sid,
        "subject_alias": r["patient"],
        "sample_type": "plasma cfDNA",
        "phenotype": phenotype,
        "disease_status": r["sm_disease"],
        "severity": r["severity"],
        "onset": r["onset"],
        "gestational_age_weeks_at_sampling": r["ga_sampling"] or r["ga_weeks"],
        "tube_type": r["tube_type"],
        "sequencing_batch": r["batch"],
        "biological_replicate": r["replicate"],
        "maternal_bmi": r["bmi"],
        "maternal_race": r["race"],
        "gravida": r["gravida"],
        "fertility_treatment": r["fertility"],
        "birth_term": r["birth_term"],
        "fetal_sex": r["fetal_sex"],
        "iugr": r["iugr"],
        "aspirin": r["aspirin"],
        "qc_pass_cov5x_ge50": qc_pass(sid),
        "mean_depth": (qc.get(sid) or {}).get("mean_depth", ""),
        "cov5x_pct": (qc.get(sid) or {}).get("cov5x_pct", ""),
        "in_subcohort": "yes" if in_sub else "no",
    })

sample_hdr = list(rows[0].keys())
n_pass = sum(r["qc_pass_cov5x_ge50"] == "yes" for r in rows)
with open(os.path.join(OUT, "sample_metadata.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=sample_hdr, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow(r)
print(f"sample_metadata.tsv: {len(rows)} samples ({sum(r['in_subcohort']=='yes' for r in rows)} subcohort, {n_pass} QC-pass)")

# ---- 2. data_file_manifest.tsv (real files) ----
cov_files = sorted(glob.glob(os.path.join(RAW, "*.cov.gz")))
gsm_to_sid = {r["gsm"]: r["sample_id"] for r in samples}
man_rows, unmatched = [], []
for fp in cov_files:
    base = os.path.basename(fp)          # GSM8644973_DNA031134.cov.gz
    gsm = base.split("_")[0]
    sid = gsm_to_sid.get(gsm)
    if sid is None:
        unmatched.append(base)
        continue
    size_mb = round(os.path.getsize(fp) / 1024 / 1024, 1)
    man_rows.append({
        "file_name": base,
        "ega_sample_alias": sid,
        "file_format": "bgzip (gzip)",
        "file_content": "Bismark .cov: chr start end methPct countM countU",
        "file_size_MB": size_mb,
        "checksum_algorithm": "md5 (to be computed before upload)",
        "reference_genome": "GRCh38/hg38",
    })
with open(os.path.join(OUT, "data_file_manifest.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(man_rows[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in man_rows:
        w.writerow(r)
print(f"data_file_manifest.tsv: {len(man_rows)} files, {len(unmatched)} unmatched")
if unmatched:
    print("  unmatched:", unmatched[:5])

# ---- 3. experiment_metadata.tsv ----
n_pe = sum(1 for r in rows if "control" not in r["phenotype"].lower())
exp_meta = [
    ("study_title", "Plasma cfDNA whole-genome bisulfite sequencing in pre-eclampsia and matched controls"),
    ("study_type", "case-control epigenomics (cfDNA WGBS)"),
    ("EGA submission type", "Data + metadata, controlled access via DAC"),
    ("experimental design", "PE (n=%d across groups) vs gestational-age matched controls; discovery subcohort PE 32 / Control 32 matched +/-3 GA weeks" % n_pe),
    ("sample material", "Maternal plasma cell-free DNA (cfDNA), blood tubes: PAXgene DNA / EDTA"),
    ("library protocol", "WGBS; cfDNA fragments bisulfite converted, PE sequencing"),
    ("assay_type", "WGBS (whole-genome bisulfite sequencing)"),
    ("reference_genome", "GRCh38 (hg38)"),
    ("alignment_and_methylation_calling", "Bismark (cov.gz output: chr, start, end, methylation %, methylated count, unmethylated count); coordinates GRCh38"),
    ("data_files", "%d .cov.gz per-sample files, one per EGA analysis" % len(man_rows)),
    ("derived_data_files", "Region-level beta matrix (500bp tiling, coverage-weighted): results/GSE282512_region_beta.csv.gz"),
    ("consent", "Human data; obtained under study-specific informed consent. See ethics approval reference (to be completed by submitter)."),
    ("contact_person", "(to be completed by submitter)"),
    ("principal_investigator", "(to be completed by submitter)"),
    ("EGAC / EGAS IDs", "pending - assigned by EGA after submission"),
]
with open(os.path.join(OUT, "experiment_metadata.tsv"), "w", newline="", encoding="utf-8") as f:
    f.write("field\tvalue\n")
    for k, v in exp_meta:
        f.write(f"{k}\t{v}\n")
print("experiment_metadata.tsv written")

# ---- 4. analysis_description.md ----
ga = " ".join(r["phenotype"] for r in rows)
analysis_md = """# Analysis description (for EGA Data Access Committee)

## Dataset
- GSE282512 re-analysis: plasma cfDNA WGBS, 369 samples from 279 patients
  (357 samples pass QC; discovery subcohort: 32 pre-eclampsia vs 32
  gestational-age matched controls, median GA difference 0 weeks).
- Reference genome: GRCh38/hg38. Per-sample methylation provided as Bismark
  .cov.gz files (chromosome, start, end, methylation %, count methylated,
  count unmethylated).

## Derived data to be shared
- Per-sample .cov.gz files (primary, one per sample).
- Region-level coverage-weighted beta matrix over 500 bp tiles
  (GSE282512_region_beta.csv.gz) with region annotation
  (GSE282512_region_annot.csv), enabling reproduction of DMR analyses
  without re-processing raw alignments.
- Sample metadata (sample_metadata.tsv) with phenotype, gestational age,
  tube type, batch, QC metrics.

## Analysis workflow (reproducible from shared data)
1. QC: per-sample CpG coverage, depth, conversion; chromosomal anomaly
   screening (QC summary: results/GSE282512_qc_summary.csv).
2. Region-level methylation: 500 bp tiling, coverage-weighted mean beta,
   limma differential methylation adjusting tube_type + sequencing batch
   (PE vs Control, matched subcohort).
3. DMR validation: site-level CpG tests, direction consistency,
   tube-type sensitivity (EDTA-only), EOPE subset; classified confidence
   tiers (results/GSE282512_dmr_v2_tiers.csv).
4. Public data integration (not required for data reuse, for reference):
   GoDML cis-meQTL, FinnGen R10 PE/GH GWAS, GSE37722 leukocyte trajectory,
   placental tissue datasets, transcriptome cohorts.

## Subject privacy
- Metadata contains no direct identifiers; patient labels are pseudonyms
  (Patient_XX). Ethnicity/race, BMI, gestational age and obstetric variables
  retained as essential analysis covariates. Fetal sex and birth outcome
  included; may be removed on DAC request.
- Recommend access decision: research use only, no re-identification,
  no attempt to contact participants.

## Citing source
- Original dataset: GEO accession GSE282512 (plasma cfDNA WGBS, PE cohort).
  Re-analysis performed by this project; derived files redistributed under
  controlled access.
"""
with open(os.path.join(OUT, "analysis_description.md"), "w", encoding="utf-8") as f:
    f.write(analysis_md)
print("analysis_description.md written")

# ---- 5. submission_checklist.md ----
checklist = f"""# EGA submission checklist (pre-submission audit)

Generated: 2026-08-30. Status: material preparation only -- NO real submission performed.

## Files prepared
- [x] sample_metadata.tsv -- {len(rows)} samples, {len(sample_hdr)} fields, subcohort flagged
- [x] data_file_manifest.tsv -- {len(man_rows)} .cov.gz files mapped ({sum(r['in_subcohort']=='yes' for r in rows)} samples in discovery subcohort)
- [x] experiment_metadata.tsv -- assay/reference/design fields
- [x] analysis_description.md -- DAC-facing description

## Items to complete before actual submission
- [ ] Compute md5 checksums for all {len(man_rows)} .cov.gz files (manifest placeholder present)
- [ ] Fill in PI, contact person, ethics approval reference in experiment_metadata.tsv
- [ ] Confirm consent terms allow EGA controlled-access deposition
- [ ] Decide deposition scope: all 357 QC-passed samples vs 64 subcohort only
      (manifest currently covers all files; subcohort-only = filter rows)
- [ ] Upload via EGA submission API / Webin-CLI; obtain EGAD accession
- [ ] Link publication DOI to EGAD accession after acceptance

## Notes / caveats
- File names in raw dir use GSM_DNAxxx.cov.gz; manifest maps GSM -> sample_id.
- Region-level beta matrix (derived) included as optional secondary file.
- Derived statistics/results files (DMR lists etc.) NOT included -- available
  from publication supplementary materials instead.
"""
with open(os.path.join(OUT, "submission_checklist.md"), "w", encoding="utf-8") as f:
    f.write(checklist)
print("submission_checklist.md written")
print("\nDONE -> results/ega_submission/")
