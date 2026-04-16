# =============================================================================
# Table 1: Patient Characteristics
# =============================================================================
# Paper: Longitudinal Associations Between Inflammation and Multi-Dimensional Fatigue up to Two Years after Colorectal Cancer Diagnosis
# Authors: Bruijnzeels, A.E.C., Mols, F., Van Deun, K., Schoormans, D.
# Output: Table 1 — Sociodemographic and clinical characteristics of CRC
#         patients and the normative sample, by treatment group.
#         Includes chi-square tests (categorical variables) and ANOVAs
#         (continuous variables).
#
# Data: Combined normative sample + PROCORE dataset (one row per patient,
#       wide format with BMI at each timepoint).
#       Load your dataset below before running this script.
# =============================================================================

# Load your dataset here:
# df <- read.csv("your_data.csv")   # or: load("your_data.RData")
# Expected object name used in original analysis:
# X5_Norm_PROCORE_inflammatoire_markers_x_behandeling

# -------------------------
# PACKAGES
# -------------------------
library(dplyr)
library(tidyr)
library(forcats)
library(gtsummary)
library(gt)

# -------------------------
# Settings
# -------------------------
id_var        <- "patID"
time_var      <- "Time"        # 1,2,3
cohort_var    <- "norm"        # 0/1 of "Norm"/"Patient"
treatment_var <- "treatment"

# Variables for the table
cat_vars_input <- c("geslacht", "opleiding_3_cat_W1", "BET_B", "burgst",
                    "aant_comorb", "tumorsoort", "stage", "type_chemo")

# Labels
var_labels <- list(
  leeft_W1           ~ "Mean age (SD), years",
  BMI_T1             ~ "Mean BMI (SD), kg/m²\nT1",
  BMI_T2             ~ "T2",
  BMI_T3             ~ "T3",
  geslacht           ~ "Sex, No. (%)",
  opleiding_3_cat_W1 ~ "Education, No. (%)",
  BET_B              ~ "Employment status, No. (%)",
  burgst             ~ "Marital status, No. (%)",
  aant_comorb        ~ "Number of comorbidities at baseline, No. (%)",
  tumorsoort         ~ "Type carcinoma, No. (%)",
  stage              ~ "Stage, No. (%)",
  type_chemo         ~ "Type of chemotherapy, No. (%)"
)

theme_gtsummary_journal("jama")
theme_gtsummary_compact()

# -------------------------
# LONG → ANALYSE-DATA
# -------------------------

# 1) Baseline-rij (Time==1) for stable variables
base_cols <- c(id_var, cohort_var, treatment_var,
               "geslacht", "opleiding_3_cat_W1", "BET_B", "burgst",
               "aant_comorb", "tumorsoort", "stage", "type_chemo",
               "leeft_W1")

df_base <- df %>%
  mutate(!!time_var := as.integer(.data[[time_var]])) %>%
  filter(.data[[time_var]] == 1) %>%
  arrange(.data[[id_var]]) %>%
  group_by(across(all_of(id_var))) %>%
  slice_head(n = 1) %>%                         # 1 rij per id
  ungroup() %>%
  select(any_of(base_cols))

# 2) BMI per Time wide: BMI_T1, BMI_T2, BMI_T3
bmi_wide <- df %>%
  mutate(.TimeTmp = as.integer(.data[[time_var]])) %>%
  select(all_of(id_var), .TimeTmp, BMI = .data[["BMI"]]) %>%
  group_by(across(all_of(c(id_var, ".TimeTmp")))) %>%
  summarise(BMI = dplyr::first(BMI), .groups = "drop") %>%
  mutate(name = paste0("BMI_T", .TimeTmp)) %>%
  select(all_of(id_var), name, BMI) %>%
  pivot_wider(names_from = name, values_from = BMI)

# 3) Combine to analysis dataset
data_clean <- df_base %>%
  left_join(bmi_wide, by = id_var) %>%
  mutate(
    # norm → "Norm"/"Patient"
    !!cohort_var := {
      v <- .data[[cohort_var]]
      if (is.numeric(v)) if_else(v == 1, "Norm", "Patient") else as.character(v)
    },
    !!cohort_var := factor(.data[[cohort_var]], levels = c("Norm","Patient")),
    !!treatment_var := as.factor(.data[[treatment_var]]),
    across(all_of(cat_vars_input), ~ as.factor(.x)),
    # comorbidities categories 0/1/≥2
    aant_comorb = case_when(
      is.na(aant_comorb) ~ NA_character_,
      as.numeric(as.character(aant_comorb)) == 0 ~ "0",
      as.numeric(as.character(aant_comorb)) == 1 ~ "1",
      as.numeric(as.character(aant_comorb)) >= 2 ~ "≥ 2"
    ) |> factor(levels = c("0","1","≥ 2"))
  )

num_vars <- c("leeft_W1", "BMI_T1", "BMI_T2", "BMI_T3")
cat_vars <- cat_vars_input

# -------------------------
# TABLES
# -------------------------

# A) Norm vs All patients
tbl_norm_vs_pat <- data_clean %>%
  select(all_of(c(num_vars, cat_vars, cohort_var))) %>%
  tbl_summary(
    by = all_of(cohort_var),
    type = list(all_of(num_vars) ~ "continuous2",
                all_of(cat_vars) ~ "categorical"),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    missing_text = "Missing",
    label = var_labels
  ) %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Norm vs. Patients**") %>%
  add_p(test = list(
    all_continuous()  ~ "t.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  bold_labels()



tbl_norm_vs_pat
