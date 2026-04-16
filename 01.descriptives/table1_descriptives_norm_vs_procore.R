# =============================================================================
# Descriptive Comparison: Normative Sample vs. PROCORE Patients
# =============================================================================
# Paper: Longitudinal Associations Between Inflammation and Multi-Dimensional Fatigue up to Two Years after Colorectal Cancer Diagnosis
# Authors: Bruijnzeels, A.E.C., Mols, F., Van Deun, K., Schoormans, D.
# Output: Descriptive statistics (means, SDs) for inflammatory biomarkers,
#         fatigue subscales, and covariates, compared between the normative
#         sample and CRC patients across timepoints.
#
# Data: Combined normative sample + PROCORE dataset.
#       Load your dataset below before running this script.
# =============================================================================

# Load your dataset here:
# data <- read.csv("your_data.csv")   # or: load("your_data.RData")
# Expected object name used in original analysis:
# X4_Norm_PROCORE_inflammatoire_markers_x_behandeling

data <- X4_Norm_PROCORE_inflammatoire_markers_x_behandeling

data$Time <- as.numeric(as.character(data$Time))
data$treatment <- as.factor(data$treatment)
data$geslacht <- as.factor(data$geslacht)
levels(data$treatment)
data$alcohol <- as.factor(data$alcohol)
data$ROOK <- as.factor(data$ROOK)
data$aant_comorb <- as.factor(data$aant_comorb)
data$norm <- factor(data$norm)
data$PAIN <- as.numeric(data$PAIN)
data$MFI_AM <- as.numeric(data$MFI_AM)
data$MFI_GV <- as.numeric(data$MFI_GV)
data$MFI_LV <- as.numeric(data$MFI_LV)
data$MFI_VA <- as.numeric(data$MFI_VA)
data$MFI_VM <- as.numeric(data$MFI_VM)

library(dplyr)
library(purrr)
library(broom)
library(tidyr)
library(writexl)

vars <- c("log_CRP", "log_IL6", "log_IL8", "log_IL10", "log_IL1a", "log_IL1b", "log_IL17A", "log_IL22", "log_TNFRI", "log_TNFRII", "log_IFN")

#List of biomarkers
biomarkers <- c("log_CRP", "log_IL6", "log_IL8", "log_IL10", "log_IL22", 
                "log_IL17A", "log_IL1b", "log_IL1a", "log_TNFRI", "log_TNFRII", "log_IFN")

#List of MFI subscales
mfi_subscales <- c("MFI_AM", "MFI_GV", "MFI_LV", "MFI_VA", "MFI_VM")

#only the normative sample 
norm_df <- subset(data, norm == 1)

#Covariates
covariates_norm <- c("geslacht","leeft_W1","BMI","aant_comorb","PAIN","alcohol","ROOK","PSQI_total")

#I want a lineair regression analysis per subscale of the MFI and per biomarker
run_norm_lm <- function(df, mfi_var, biomarker) {
  needed <- c(mfi_var, biomarker, covariates_norm)
  df_cc  <- df %>% dplyr::select(all_of(needed)) %>% na.omit()
  if (nrow(df_cc) < 10) return(NULL)
  
  fml <- as.formula(paste0(mfi_var, " ~ ", paste(c(biomarker, covariates_norm), collapse = " + ")))
  fit <- lm(fml, data = df_cc)
  
  tidy_fit   <- broom::tidy(fit, conf.int = TRUE) %>% mutate(n = nobs(fit)) 
  glance_fit <- broom::glance(fit) %>% dplyr::select(r.squared, adj.r.squared, AIC, BIC)
  
  tidy_fit <- tidy_fit %>%
    mutate(
      outcome   = mfi_var,
      predictor = biomarker
    ) %>%
    relocate(outcome, predictor, n) %>%
    left_join(glance_fit, by = character())
  
  return(tidy_fit)
}

#run all the models
grid <- expand.grid(mfi_var = mfi_subscales,
                    biomarker = biomarkers,
                    stringsAsFactors = FALSE)

norm_results <- pmap_dfr(grid, ~ run_norm_lm(norm_df, ..1, ..2))

norm_biomarker_slopes <- norm_results %>%
  filter(term == predictor) %>%                      # slope only of the biomarker itself
  transmute(
    MFI_Subscale = outcome,
    Biomarker    = predictor,
    beta         = estimate,
    SE           = std.error,
    t            = statistic,
    p            = p.value,
    n,
    R2           = r.squared,
    adj_R2       = adj.r.squared,
    AIC, BIC
  ) %>%
  arrange(Biomarker, MFI_Subscale)

#Export to Excel
write_xlsx(
  list(
    "All_terms_full_models" = norm_results,
    "Biomarker_slopes_only" = norm_biomarker_slopes
  ),
  path = "LR_normsample_per_biomarker_per_MFI.xlsx"
)

print(norm_biomarker_slopes, n = Inf)

