# =============================================================================
# Table 2 & Table 3: Linear Mixed Models — Fatigue Trajectories Over Time
# =============================================================================
# Paper: Longitudinal Associations Between Inflammation and Multi-Dimensional Fatigue up to Two Years after Colorectal Cancer Diagnosis
# Authors: Bruijnzeels, A.E.C., Mols. F., Van Deun, K., Schoormans, D.
# Output: Table 2 — LMM fixed effects (time) per MFI subscale.
#
# Model: LMM with time as fixed effect, random intercept per patient.
# Subset: CRC patients with valid CRP measurement (CRP_c > 0) who
#         completed the blood assessment (responder == 1).
#
# Data: PROCORE long-format dataset with inflammatory markers,
#       multiple imputation (1000 iterations, excl. DCRA).
#       Load your dataset below before running this script.
# =============================================================================

# Load your dataset here:
# data <- read.csv("your_data.csv")   # or: load("your_data.RData")
# Expected object name used in original analysis:
# X3_Inflammatoire_markers_x_behandeling_long_format_bloedvragenlijstje_zonder_DCRA_andere_imputatiemethode_1000_iteraties

# Note: the subsetting below (CRP_c > 0 & responder == 1) is applied within
# the script — adjust variable names to match your dataset if needed.

data <- subset(X3_Inflammatoire_markers_x_behandeling_long_format_bloedvragenlijstje_zonder_DCRA_andere_imputatiemethode_1000_iteraties, CRP_c > 0 & responder == 1)

library(lme4)
library(lmerTest)
library(emmeans)

data$Time <- factor(data$Time, levels = c(1, 2, 3), labels = c("T0", "T1", "T2"))

#General fatigue
model_AM <- lmer(MFI_AM ~ Time + (1 | patID), data = data)
emm_AM <- emmeans(model_AM, ~ Time)
pairwise_contrasts_AM <- contrast(emm_AM, method = "pairwise")

raw_pvals_AM <- summary(pairwise_contrasts_AM)$p.value

fdr_pvals_AM <- p.adjust(raw_pvals_AM, method = "fdr")
results_AM <- as.data.frame(summary(pairwise_contrasts_AM))
results_AM$fdr_p <- fdr_pvals_AM

print(results_AM)

#Physical fatigue
model_LV <- lmer(MFI_LV ~ Time + (1 | patID), data = data)
emm_LV <- emmeans(model_LV, ~ Time)
pairwise_contrasts_LV <- contrast(emm_LV, method = "pairwise")

raw_pvals_LV <- summary(pairwise_contrasts_LV)$p.value

fdr_pvals_LV <- p.adjust(raw_pvals_LV, method = "fdr")
results_LV <- as.data.frame(summary(pairwise_contrasts_LV))
results_LV$fdr_p <- fdr_pvals_LV

print(results_LV)

#Mental fatigue
model_GV <- lmer(MFI_GV ~ Time + (1 | patID), data = data)
emm_GV <- emmeans(model_GV, ~ Time)
pairwise_contrasts_GV <- contrast(emm_GV, method = "pairwise")

raw_pvals_GV <- summary(pairwise_contrasts_GV)$p.value

fdr_pvals_GV <- p.adjust(raw_pvals_GV, method = "fdr")
results_GV <- as.data.frame(summary(pairwise_contrasts_GV))
results_GV$fdr_p <- fdr_pvals_GV

print(results_GV)

#Reduced activity
model_VA <- lmer(MFI_VA ~ Time + (1 | patID), data = data)
emm_VA <- emmeans(model_VA, ~ Time)
pairwise_contrasts_VA <- contrast(emm_VA, method = "pairwise")

raw_pvals_VA <- summary(pairwise_contrasts_VA)$p.value

fdr_pvals_VA <- p.adjust(raw_pvals_VA, method = "fdr")
results_VA <- as.data.frame(summary(pairwise_contrasts_VA))
results_VA$fdr_p <- fdr_pvals_VA

print(results_VA)

#Reduced motivation
model_VM <- lmer(MFI_VM ~ Time + (1 | patID), data = data)
emm_VM <- emmeans(model_VM, ~ Time)
pairwise_contrasts_VM <- contrast(emm_VM, method = "pairwise")

raw_pvals_VM <- summary(pairwise_contrasts_VM)$p.value

fdr_pvals_VM <- p.adjust(raw_pvals_VM, method = "fdr")
results_VM <- as.data.frame(summary(pairwise_contrasts_VM))
results_VM$fdr_p <- fdr_pvals_VM

print(results_VM)
