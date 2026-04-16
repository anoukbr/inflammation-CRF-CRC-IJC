# =============================================================================
# Correlations: Pre-Analytical Influences and Biomarker Levels
# =============================================================================
# Paper: Longitudinal Associations Between Inflammation and Multi-Dimensional Fatigue up to Two Years after Colorectal Cancer Diagnosis
# Authors: Bruijnzeels, A.E.C., Mols, F., Van Deun, K., Schoormans, D.
# Output: Correlations between pre-analytical factors (caffeine intake,
#         physical activity, influenza vaccination) and inflammatory biomarker
#         levels, reported in the Exploratory Analyses section (Section 3.1).
#         Sensitivity analyses excluding participants with high pre-analytical
#         exposure and elevated biomarker levels.
#
# Data: PROCORE long-format dataset with inflammatory markers,
#       multiple imputation (1000 iterations).
#       Load your dataset below before running this script.
# =============================================================================

# Load your dataset here:
# data <- read.csv("your_data.csv")   # or: load("your_data.RData")
# Expected object name used in original analysis:
# X3_Inflammatoire_markers_x_behandeling_long_format_bloedvragenlijstje_zonder_DCRA_andere_imputatiemethode_1000_iteraties


library(dplyr)
library(correlation)
library(openxlsx)
library(Hmisc)

data <- X3_Inflammatoire_markers_x_behandeling_long_format_bloedvragenlijstje_zonder_DCRA_andere_imputatiemethode_1000_iteraties

data <- data %>%
  filter(CRP_c > 0, responder == 1)

vars <- c("CRP_c", "IFNgamma", "IL_1_alpha", "IL_1_beta", "IL_10",
          "IL_17A", "IL_22", "IL_6", "IL_8", "TNFRI",
          "TNFRII", "bv_griep", "bv_prik", "bv_insp", "bv_koffie",
          "bv_alcohol", "bv_overgang", "bv_horm")

# Splits data per Time-niveau and calculate correlations
time_levels <- unique(data$Time)

wb <- createWorkbook()

for (t in time_levels) {
  
  # Subset per Time-niveau
  df_sub <- data %>%
    filter(Time == t) %>%
    select(all_of(vars)) %>%
    as.matrix()
  
  # Calculate correlations + p-values
  res <- rcorr(df_sub, type = "pearson")
  r_mat <- round(res$r, 3)
  p_mat <- round(res$P, 4)
  n_mat <- res$n
  
  # Combine r and p into one readable matrix r(p)
  combined <- matrix(
    paste0(r_mat, " (p=", p_mat, ")"),
    nrow = nrow(r_mat),
    dimnames = dimnames(r_mat)
  )
  diag(combined) <- "1"
  
  # Divert to dataframe with variablenames in the first column
  combined_df <- as.data.frame(combined)
  combined_df <- cbind(Variabele = rownames(combined_df), combined_df)
  rownames(combined_df) <- NULL
  
  # New tab
  sheet_name <- paste0("Time_", t)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, combined_df)
  
  header_style <- createStyle(
    fontColour = "#FFFFFF",
    fgFill = "#4472C4",
    halign = "CENTER",
    textDecoration = "Bold"
  )
  addStyle(wb, sheet_name, header_style,
           rows = 1, cols = 1:ncol(combined_df), gridExpand = TRUE)
  
  setColWidths(wb, sheet_name, cols = 1:ncol(combined_df), widths = "auto")
  
  # Mark significant correlations
  sig_style <- createStyle(fgFill = "#E2EFDA")
  
  for (row_i in 1:nrow(p_mat)) {
    for (col_j in 1:ncol(p_mat)) {
      if (!is.na(p_mat[row_i, col_j]) && p_mat[row_i, col_j] < 0.05 && row_i != col_j) {
        addStyle(wb, sheet_name, sig_style,
                 rows = row_i + 1,   # +1 vanwege header
                 cols = col_j + 1,   # +1 vanwege eerste kolom met namen
                 stack = TRUE)
      }
    }
  }
}

# Save
saveWorkbook(wb, "correlaties_per_time_ruw.xlsx", overwrite = TRUE)
cat("Bestand opgeslagen als: correlaties_per_time_ruw.xlsx\n")
