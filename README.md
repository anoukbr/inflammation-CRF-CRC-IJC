# Analysis Code: Inflammation and Fatigue in Colorectal Cancer Patients

This repository contains all R analysis scripts accompanying the manuscript:

> **Longitudinal Associations Between Inflammation and Multi-Dimensional Fatigue up to Two Years after Colorectal Cancer Diagnosis**  
> Bruijnzeels, A.E.C., Mols, F., Van Deun, K, Schoormans, D.  
> *International Journal of Cancer*, 2026. DOI: 

---

## Study Overview

This study examined the relationship between inflammatory biomarkers and fatigue in cancer patients (PROCORE cohort), compared to a normative reference group. Data were collected at three timepoints: baseline (T0), 12-month follow-up (T1), and 24-month follow-up (T2).

**Key measures:**
- Inflammatory biomarkers: CRP, IL-6, IL-8, IL-10, IL-1α, IL-1β, IL-17A, IL-22, TNF-RI, TNF-RII, IFN-γ
- Fatigue: Multidimensional Fatigue Inventory (MFI) — five subscales
- Additional outcomes: BMI, pain, sleep quality (PSQI) (smoking, alcohol use, age, sex, number of comorbidities and days since diagnosis)

---

## Repository Structure

```
├── 01_descriptives/
│   ├── table1_patient_characteristics.R               # Table 1: Demographics & clinical characteristics
│   ├── table1_supplement_norm_comparisons.R           # Table 1: Demographics & clinical characteristics
│   └── table1_descriptives_norm_vs_procore.R          # Table 1: Comparison: norm vs. PROCORE
│
├── 02_main_analyses/
│   ├── figure2_fatigue_trajectories
|   ├── table2_LMM_fatigue_over_time.R                 # Table 2: LMM fatigue change over time
│   ├── table3_regression_inflammation_fatigue_norm.R  # Table 3: Regression inflammation → fatigue
│   └── tableS3_figure3_hybrid_LMM.R                   # Table S3 & Figure 3: Hybrid LMM models
│
├── 03_sensitivity_analyses/
│   ├── sensitivity_complete_vs_incomplete_tableS2_figureS1.R      # Table S2 & Figure S1: Complete vs. incomplete data comparison
│   ├── sensitivity_blood_vs_no_blood_tableS1.R                    # Table S1: Baseline: blood vs. no blood sample
|   └── correlations_preanalytical_influences_describedintext.R    # Correlations between preanalytical influences and biomarkers: described in-text
│
└── README.md
```

---

## Data Availability

The data used in these analyses contain sensitive patient information and are **not publicly available** due to privacy regulations (GDPR). Researchers interested in access to the PROCORE dataset may contact Anouk Bruijnzeels at a.e.c.bruijnzeels@tilburguniversity.edu. All data is available upon reasonable request from www.profilesregistry.nl.

The scripts reference the following datasets:
| Object name in script | Description |
|---|---|
| `X1_KEY_FILE_Final_met_inflammatoire_markers_en_bloedvragenlijst` | Full dataset including inflammatory markers and blood questionnaire |
| `X3_Inflammatoire_markers_x_behandeling_long_format_bloedvragenlijstje_zonder_DCRA_andere_imputatiemethode_1000_iteraties` | Long-format data, PROCORE patients, 1000-iteration imputation |
| `X4_Norm_PROCORE_inflammatoire_markers_x_behandeling` | Combined norm + PROCORE dataset |
| `X5_Norm_PROCORE_inflammatoire_markers_x_behandeling` | Extended combined norm + PROCORE dataset |

---

## Software & Packages

All analyses were conducted in **R** (version X.X.X).

| Package | Version | Purpose |
|---|---|---|
| `lme4` | — | Linear mixed models |
| `lmerTest` | — | p-values for mixed models |
| `emmeans` | — | Estimated marginal means & pairwise contrasts |
| `gtsummary` | — | Summary tables |
| `gt` | — | Table formatting |
| `dplyr` | — | Data manipulation |
| `tidyr` | — | Data reshaping |
| `purrr` | — | Functional programming |
| `broom` / `broom.mixed` | — | Tidy model output |
| `ggplot2` | — | Figures |
| `openxlsx` | — | Export to Excel |
| `correlation` | — | Correlation analyses |
| `Hmisc` | — | Statistical utilities |

---

## Script Overview

| Script | Output | Description |
|---|---|---|
| `01.descriptives/table1_patient_characteristics.R` | Table 1 | Demographics, clinical characteristics, BMI by cohort and treatment group |
| `01.descriptives/table1_supplement_norm_comparisons.R` | Table 1 supplement | One-sample t-tests: MFI subscales, BMI, pain, sleep vs. normative values |
| `01.descriptives/table1_descriptives_norm_vs_procore.R` | Descriptive output | Descriptive statistics comparing norm and PROCORE groups |
| `02.main_analyses/table2_LMM_fatigue_over_time.R` | Table 2 | LMMs per MFI subscale, pairwise contrasts across timepoints, FDR correction |
| `02.main_analyses/table3_regression_inflammation_fatigue_norm.R` | Table 3 | Linear regression: each biomarker predicting each MFI subscale, per group/timepoint, FDR correction |
| `02.main_analyses/tableS3_figure3_hybrid_LMM.R` | Tables S3, Figure 3 | Hybrid LMM models: inflammation × time interaction, forest plots |
| `02.main_analyses/figure2_fatigue_trajectories.R` | Figure 2 | Line plots: MFI subscale means ± SD for PROCORE and norm group over time |
| `03.sensitivity_analyses/sensitivity_complete_vs_incomplete_tableS2_figureS1.R` | Sensitivity Table | Comparison of patients with complete vs. incomplete data across timepoints |
| `03.sensitivity_analyses/sensitivity_blood_vs_no_blood_tableS1.R` | Sensitivity Table | Comparison of patients with vs. without blood sample at baseline |
| `03.sensitivity_analyses/correlations_preanalytical_influences_describedintext.R` | Correlation output | Inter-correlations among biomarkers and pre-analytical influences |

---

## Reproducibility Notes

- All biomarker variables are log-transformed prior to regression and LMM analyses (e.g., `log_CRP`, `log_IL6`).
- Multiple testing is corrected using the **False Discovery Rate (FDR)** method (`p.adjust(..., method = "fdr")`).
- Missing inflammatory biomarkerdata were handled via **multiple imputation** (1000 iterations).
- Analyses in scripts referencing `CRP_c > 0 & responder == 1` are restricted to patients with valid CRP measurements who completed the blood assessment.

---

## Contact

For questions about the analyses, please contact:  
**Anouk Bruijnzeels** — a.e.c.bruijnzeels@tilburguniversity.edu
Tilburg University, Department of Medical and Clinical Psychology
