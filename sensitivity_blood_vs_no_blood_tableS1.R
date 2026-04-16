# =============================================================================
# Sensitivity Analysis: Blood Sample vs. No Blood Sample at Baseline (Table S1)
# =============================================================================
# Paper: Longitudinal Associations Between Inflammation and Multi-Dimensional Fatigue up to Two Years after Colorectal Cancer Diagnosis
# Authors: Bruijnzeels, A.E.C., Mols, F., Van Deun, K., Schoormans, D.
# Output: Table S1 — Comparison of sociodemographic and clinical
#         characteristics between CRC patients who provided a blood sample
#         at baseline (n = 411) and those who did not (n = 46).
#         Chi-square tests for categorical variables, independent samples
#         t-tests for continuous variables.
#
# Data: Full PROCORE key file including inflammatory markers and blood
#       questionnaire. Load your dataset below before running this script.
# =============================================================================

# Load your dataset here:
# df_filtered is derived inside the script via filtering
# data <- read.csv("your_data.csv")   # or: load("your_data.RData")
# Expected object name used in original analysis:
# X1_KEY_FILE_Final_met_inflammatoire_markers_en_bloedvragenlijst

library(dplyr)
df_filtered <- X1_KEY_FILE_Final_met_inflammatoire_markers_en_bloedvragenlijst %>%
  filter(is.na(inclusie) | inclusie <= 2) %>%
  filter(responder_W1 == 1)

#crosstabs with chi-kwadraat
crosstab_chi <- function(data, var) {
  tbl <- table(data$bloed_W1, data[[var]])
  print(paste("Crosstab: bloed_W1 x", var))
  print(tbl)
  print(prop.table(tbl, margin = 1))  
  print(prop.table(tbl, margin = 2))  
  print(chisq.test(tbl))
}

# Crosstabs
crosstab_chi(df_filtered, "gesl")
crosstab_chi(df_filtered, "aant_comorb_W1")
crosstab_chi(df_filtered, "treatment")
crosstab_chi(df_filtered, "stage")
crosstab_chi(df_filtered, "opleiding_3_cat_W1")
crosstab_chi(df_filtered, "burgst_W1")
crosstab_chi(df_filtered, "tumorsoort")

# T-tests for continuous variables
t.test(time_since_diagnosis_W1 ~ bloed_W1, data = df_filtered)
t.test(leeftijd_W1 ~ bloed_W1, data = df_filtered)

library(gtsummary)

df_filtered %>%
  select(bloed_W1, gesl, aant_comorb_W1, treatment, stage,
         opleiding_3_cat_W1, burgst_W1, tumorsoort,
         time_since_diagnosis_W1, leeftijd_W1) %>%
  tbl_summary(by = bloed_W1,
              missing = "no",
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = list(
                all_continuous() ~ 2,
                all_categorical() ~ 0
              )) %>%
  add_p() %>%
  modify_fmt_fun(p.value ~ function(x) style_pvalue(x, digits = 3))
