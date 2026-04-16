# =============================================================================
# Table 1 Supplement: One-Sample T-Tests Against Normative Values
# =============================================================================
# Paper: Longitudinal Associations Between Inflammation and Multi-Dimensional Fatigue up to Two Years after Colorectal Cancer Diagnosis
# Authors: Bruijnzeels, A.E.C., Mols, F., Van Deun, K., Schoormans, D.
# Output: Supplement to Table 1 — One-sample t-tests comparing MFI subscale
#         scores (General, Mental, Physical Fatigue, Reduced Activity,
#         Reduced Motivation), BMI, pain (EORTC), and sleep quality (PSQI)
#         in CRC patients against published normative means, per timepoint.
#
# Data: PROCORE long-format dataset with inflammatory markers,
#       multiple imputation (1000 iterations).
#       Load your dataset below before running this script.
# =============================================================================

# Load your dataset here:
# data <- read.csv("your_data.csv")   # or: load("your_data.RData")
# Expected object name used in original analysis:
# X3_Inflammatoire_markers_x_behandeling_long_format_bloedvragenlijstje_zonder_DCRA_andere_imputatiemethode_1000_iteraties

data <- X3_Inflammatoire_markers_x_behandeling_long_format_bloedvragenlijstje_zonder_DCRA_andere_imputatiemethode_1000_iteraties

# normvalues per subscale (based on values generated from table1_patient_characteristics)
norms <- list(
  MFI_AM = 8.7,  # General Fatigue
  MFI_GV = 7.2,  # Mental Fatigue
  MFI_LV = 8.5,  # Physical Fatigue
  MFI_VA = 8.8, # Reduced Activity
  MFI_VM = 8.0   # Reduced Motivation
)

labels <- list(
  MFI_AM = "General Fatigue",
  MFI_GV = "Mental Fatigue",
  MFI_LV = "Physical Fatigue",
  MFI_VA = "Reduced Activity",
  MFI_VM = "Reduced Motivation"
)

timepoints <- c(1, 2, 3)

for (var in names(norms)) {
  mu <- norms[[var]] 
  
  for (t in timepoints) {
    subset_data <- subset(data, Time == t & !is.na(data[[var]]))
    scores <- subset_data[[var]]
    
    if (length(scores) >= 2) {
      test <- t.test(scores, mu = mu)
      
      cat("==", labels[[var]], "– T", t, "==\n", sep = "")
      cat("Vergelijking met normwaarde:", mu, "\n")
      cat("M =", round(mean(scores), 2), 
          "| SD =", round(sd(scores), 2), 
          "| n =", length(scores), "\n")
      cat("t =", round(test$statistic, 2), 
          "| p =", format.pval(test$p.value, digits = 4), "\n\n")
    } else {
      cat("==", labels[[var]], "– T", t, "==\n", sep = "")
      cat("Niet genoeg data voor t-test (n =", length(scores), ")\n\n")
    }
  }
}

# the same for BMI and pain
extra_vars <- list(
  BMI = 26.6,
  PAIN = 13.8
)
# labels
extra_labels <- list(
  BMI = "Body Mass Index (BMI)",
  PAIN = "Pain"
)
 
timepoints <- c(1, 2, 3)

for (var in names(extra_vars)) {
  mu <- extra_vars[[var]]  
  
  for (t in timepoints) {
    subset_data <- subset(data, Time == t & !is.na(data[[var]]))
    scores <- subset_data[[var]]
    
    if (length(scores) >= 2) {
      test <- t.test(scores, mu = mu)
      
      cat("==", extra_labels[[var]], "– T", t, "==\n", sep = "")
      cat("Vergelijking met normwaarde:", mu, "\n")
      cat("M =", round(mean(scores), 2), 
          "| SD =", round(sd(scores), 2), 
          "| n =", length(scores), "\n")
      cat("t =", round(test$statistic, 2), 
          "| p =", format.pval(test$p.value, digits = 4), "\n\n")
    } else {
      cat("==", extra_labels[[var]], "– T", t, "==\n", sep = "")
      cat("Niet genoeg data voor t-test (n =", length(scores), ")\n\n")
    }
  }
}

#and for sleep quality
extra_sleepquality <- list(PSQI_total = 3.2905)

extra_labelsl <- list(PSQI_total = "Sleep quality (PSQI)")

timepoints <- c(1, 2, 3)

for (var in names(extra_sleepquality)) {
  mu <- extra_sleepquality[[var]] 
  
  for (t in timepoints) {
    subset_data <- subset(data, Time == t & !is.na(data[[var]]))
    scores <- subset_data[[var]]
    
    if (length(scores) >= 2) {
      test <- t.test(scores, mu = mu)
      
      cat("==", extra_labelsl[[var]], "– T", t, "==\n", sep = "")
      cat("Vergelijking met normwaarde:", mu, "\n")
      cat("M =", round(mean(scores), 2), 
          "| SD =", round(sd(scores), 2), 
          "| n =", length(scores), "\n")
      cat("t =", round(test$statistic, 2), 
          "| p =", format.pval(test$p.value, digits = 4), "\n\n")
    } else {
      cat("==", extra_labelsl[[var]], "– T", t, "==\n", sep = "")
      cat("Niet genoeg data voor t-test (n =", length(scores), ")\n\n")
    }
  }
}
