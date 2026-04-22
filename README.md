# Mediated PRS Script — User Guide & Pre-Flight Checklist
`mediated_score_calc_01_2026.py`

---

## 1  Overview

This script computes a Mediated Polygenic Risk Score (Mediated PRS) that combines a standard GWAS-weighted PRS with a mediator-adjusted component [(Thom CS,...,Vujkovic M, Voight BF., 2020)](https://pubmed.ncbi.nlm.nih.gov/32833022/) . It accepts genotype data in PLINK binary format, a SNP weight sheet, and a sample phenotype/covariate sheet, then outputs per-sample PRS values and (optionally) classifier performance metrics and ROC plots. Prior to running this script you should have gather effect sizes from two GWASes conducted in datsets independent of your test cohort. One GWAS should be your Mediated Phenotype (MP) adjusted by Mediator (M).The other should be a GWAS of you Mediator (M)

- **Script file:** `mediated_score_calc_01_2026.py`
- **External dependency:** `genotype_filter_01_2026.sh` (must be in the working directory unless `--skip_raw` is used)

---

## 2  Required Input Files

### 2.1  SAMPLESHEET — `--samplesheet`
Tab-separated, any extension. One row per individual. Must have a header row. The column named `IID` is used as the row index and must match the sample IDs in the PLINK `.fam` file.

| Column | Required? | Notes |
|--------|-----------|-------|
| `IID` | Required | Must match `.fam` IDs exactly |
| `P` | Required | Phenotype column. Binary: 0/1 or 1/2 — must be consistent with `--type_of_pheno`. Quantitative: numeric. |
| `M` | Required | Mediator value per individual. Can be binary or quantitative |
| `PC*` | Required if `--use_pcs` | Principal components; must be numeric |

To include additonal covariates in samplesheet please pass comma seperated list to `--use-covars` argument

### 2.2  SNPSHEET — `--snpsheet`
Whitespace-separated, one row per SNP. Used to extract variants from PLINK files, build the effect-allele file, and weight genotypes.

| Column | Required? | Notes |
|--------|-----------|-------|
| `SNP_ID` | Required | Format: `chr<N>:<POS>` e.g. `chr1:12345`. Used to build PLINK extraction list and effect-allele file. |
| `BETA_MP` | Required | Effect size of the SNP on the mediated phenotype |
| `BETA_M` | Required | Effect size of the SNP on the mediator |
| `EA_MP` | Required | Effect allele for the mediated phenotype |
| `EA_M` | Required | Effect allele for the mediator. If `EA_M ≠ EA_MP`, `BETA_M` is sign-flipped automatically. |

### 2.3  PLINK Files — `--plink_prefix`
Provide the path prefix only — do not include the file extension (`.bed`/`.bim`/`.fam`).

- `SNP_ID` values in the SNPSHEET must resolve to variant IDs present in the `.bim` file after the script strips the `chr` prefix and reformats to `<CHR>:<POS>`.
- `IID` values in the SAMPLESHEET must match the second column of the `.fam` file.

---

## 3  Command-Line Arguments

| Argument | Description |
|----------|-------------|
| `--samplesheet PATH` | Path to the sample/phenotype sheet (tab-separated, `IID` index). |
| `--snpsheet PATH` | Path to the SNP weight sheet (whitespace-separated). |
| `--plink_prefix PATH` | Prefix for PLINK binary files (`.bed`/`.bim`/`.fam`). Do not include the extension. |
| `--skip_raw` | *(flag)* Skip re-running the PLINK filtering shell script and use an existing `*_filtered.raw` file. |
| `--calculate_auc` | *(flag)* After scoring, run a logistic or linear regression classifier and output performance metrics and a ROC plot. |
| `--type_of_pheno VALUE` | `binary` or `quantitative`. Required when `--calculate_auc` is used. |
| `--use_pcs` | *(flag)* Include `PC1 - PC*` from the SAMPLESHEET as covariates in the classifier model where * is an integer passed to this argument. |
| `--use-covars` |  *(flag)* Include additonal covariate columns from SAMPLESHEET provided as comma seperated list (e.g 'AGE,SEX,BMI') |

---

## 4  Pre-Flight Checklist

### 4.1 Envionment Setup
- Use the provided `environment.yml` file to build a conda environment that contains all paclages required to succesfully run this script

`conda create -f environment.yml`
### 4.2  Variant Alignment
- [ ] `SNP_ID` format in SNPSHEET matches the PLINK `.bim` IDs after `chr`-prefix stripping (script converts `chr1:12345` → `1:12345`).
- [ ] No duplicate `SNP_ID` values in SNPSHEET.
- [ ] `EA_MP` and `EA_M` columns are populated for every row — missing values cause silent score errors.
- [ ] `BETA_MP` and `BETA_M` values are on the same scale (e.g., both log-odds, or both linear effect sizes).
- [ ] Effect alleles are aligned to the same strand as the PLINK reference. The script does not perform strand-flip correction.

### 4.3  Sample Alignment
- [ ] `IID` values in SAMPLESHEET match the second column of the `.fam` file exactly (case-sensitive).
- [ ] No duplicate `IID` values in SAMPLESHEET (duplicates are detected and printed; both entries will be used unless removed).
- [ ] No missing `IID` values (NaN rows are dropped at merge).

### 4.4  Phenotype / Mediator Sanity
- [ ] Phenotype column `P` is coded consistently: binary 0/1 (or 1/2) — confirm the exact coding expected by your downstream analysis.
- [ ] Mediator column `M` is numeric with no unexpected string values (the script filters rows where `P == 'SUS'` before AUC calculation).
- [ ] **If `--use_pcs` is set:** `PC1`, `PC2`, `PC3`, and `PC4` columns are present in SAMPLESHEET and are numeric.

### 4.5  Runtime / File Dependencies
- [ ] `genotype_filter_01_2026.sh` is in the working directory (required unless `--skip_raw` is used).
- [ ] If `--skip_raw` is set, the file `<plink_prefix>_filtered.raw` must already exist in the working directory.

---

## 5  Running the Script

### 5.1  Basic Usage (score only)
```bash
python mediated_score_calc_01_2026.py \
  --samplesheet /path/to/samplesheet.txt \
  --snpsheet /path/to/snpsheet.txt \
  --plink_prefix /path/to/plink/prefix
```

### 5.2  With AUC / ROC Output (binary phenotype)
```bash
python mediated_score_calc_01_2026.py \
  --samplesheet /path/to/samplesheet.txt \
  --snpsheet /path/to/snpsheet.txt \
  --plink_prefix /path/to/plink/prefix \
  --calculate_auc \
  --type_of_pheno binary \
  --use_pcs 5
```

### 5.3  Skip PLINK Filtering (re-use existing `.raw` file)
```bash
python mediated_score_calc_01_2026.py \
  --samplesheet /path/to/samplesheet.txt \
  --snpsheet /path/to/snpsheet.txt \
  --plink_prefix /path/to/plink/prefix \
  --skip_raw
```

---

## 6  Outputs

| File | Contents |
|------|----------|
| `<snpsheet>_Mediated_PRS_output.txt` | Tab-separated. Per-sample scores: `PRS_P`, `PRS_M`, `PRS` (combined), plus `M` and `P` columns. |
| `Dropped_SNPs.txt` | `SNP_ID`s from SNPSHEET that had no match in the filtered genotype data. |
| `<plink_prefix>_filtered.raw` | Overwritten each run (unless `--skip_raw`). Dosage matrix with allele coding inverted (2 − dosage). |
| `<samplesheet>_metrics.txt` | Only with `--calculate_auc --type_of_pheno binary`. Contains AUROC, Accuracy, Sensitivity, Specificity, Brier score. |
| `<samplesheet>_AUC_ROC.png` | Only with `--calculate_auc --type_of_pheno binary`. ROC curve plot saved to the working directory. |

The console also prints the percentage of SNPSHEET variants matched and dropped.

**Score interpretation:** `PRS_P` is the sum of genotype × `BETA_MP` across matched SNPs. `PRS_M` is the individual's mediator value (`M`) multiplied by the sum of genotype × `BETA_M`. `PRS` = `PRS_P` + `PRS_M`.

---

## 7  Troubleshooting

### 7.1  Common Issues

| Symptom | Likely Cause & Fix |
|---------|--------------------|
| 0% of variants matched / all SNPs in `Dropped_SNPs.txt` | SNP_ID format mismatch. Confirm your `.bim` uses `<CHR>:<POS>` format. The script strips `chr` and reformats; cross-check a few IDs manually. |
| 0 samples in output / empty PRS table | IID mismatch between SAMPLESHEET and `.fam`. Check for leading zeros, case differences, or FID prefixes in the `.fam` file. |
| `KeyError` on `M` or `P` columns | Column name missing or mis-spelled in SAMPLESHEET. The script expects exact column names `M` and `P`. |
| `KeyError` on `PC1`–`PC4` | `--use_pcs` was set but PC columns are absent or mis-named in SAMPLESHEET. |
| `ValueError` in logistic regression / AUC calc | Ensure `P` is numeric and binary after filtering `SUS` rows. |
| Shell script not found | `genotype_filter_01_2026.sh` must be in the working directory. Use `--skip_raw` if you already have a `*_filtered.raw` file. |
| All-zero or constant PRS values | All SNP betas may be zero, or all genotypes are the same dosage. Check SNPSHEET BETA columns and the filtered `.raw` file. |

### 7.2  Diagnostic Commands

**Check SNP ID overlap:**
```python
import pandas as pd
snps = pd.read_table('snpsheet.txt', sep=r'\s+')['SNP_ID'].str.replace('chr','')
bim  = pd.read_table('plink.bim', header=None)[1]
print(snps.isin(bim).mean(), '% overlap')
```

**Check sample ID overlap:**
```python
ss  = pd.read_table('samplesheet.txt')['IID']
fam = pd.read_table('plink.fam', header=None, sep=r'\s+')[1]
print(ss.isin(fam).mean(), '% overlap')
```

**Inspect first lines of each input:**
```bash
head -3 samplesheet.txt snpsheet.txt plink.bim plink.fam
```

---

## 8  Score Methodology (Brief)

The script computes three scores per individual:

| Score | Formula |
|-------|---------|
| `PRS_P` | Σ (genotypeᵢ × BETA_MPᵢ) over matched SNPs |
| `PRS_M` | M (individual mediator) × Σ (genotypeᵢ × BETA_Mᵢ). `BETA_M` is sign-flipped when `EA_M ≠ EA_MP`. |
| `PRS` | `PRS_P` + `PRS_M` (combined score) |

**Dosage coding:** the script inverts raw PLINK dosages as `2 − dosage`, so the score is weighted relative to the non-reference (effect) allele count. SNPs with `BETA_MP = 0` are excluded from `PRS_P`; SNPs with `BETA_M = 0` are excluded from `PRS_M`.

**Classifier model (binary):** logistic regression on the mediator score (`M` column), optionally with PC covariates, using a 75/25 train/test split. Note: the classifier uses `M` as the predictor, not the combined `PRS`.

---

## 9  Contact / Maintenance

![Lab Logo](./assets/Setia-Verma-Logo_Horizontal-Color.png)

*From the minds in Setia-Verma Lab at the University of Pennsylania* 

For questions or issues contact [Lanna Caruth](mailto:lanna.caruth@pennmedicine.upenn.edu)

Documentation generated with help of [Claude](https://claude.ai)

---

