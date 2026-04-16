# =============================================================================
# Table 3: Linear Regression — Inflammation and Fatigue in the Normative Sample
# =============================================================================
# Paper: Longitudinal Associations Between Inflammation and Multi-Dimensional Fatigue up to Two Years after Colorectal Cancer Diagnosis
# Authors: Bruijnzeels, A.E.C., Mols, F., Van Deun, K., Schoormans, D.
# Output: Table 3 — Associations between log-transformed inflammatory
#         biomarkers and MFI subscales in the normative sample, using
#         separate linear regression models per biomarker-subscale combination,
#         adjusted for covariates (age, sex, BMI, comorbidities, alcohol use,
#         smoking, pain, sleep quality). FDR correction applied per subscale.
#         Complete cases only (n = 147).
#
# Biomarkers: log_CRP, log_IL6, log_IL8, log_IL10, log_IL1a, log_IL1b,
#             log_IL17A, log_IL22, log_TNFRI, log_TNFRII, log_IFN
# Outcomes:   MFI_AM (General), MFI_GV (Mental), MFI_LV (Physical),
#             MFI_VA (Reduced Activity), MFI_VM (Reduced Motivation)
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

vars <- c("log_CRP", "log_IL6", "log_IL8", "log_IL10", "log_IL1a", "log_IL1b", "log_IL17A", "log_IL22", "log_TNFRI", "log_TNFRII", "log_IFN")

#List of biomarkers
biomarkers <- c("log_CRP", "log_IL6", "log_IL8", "log_IL10", "log_IL22", 
                "log_IL17A", "log_IL1b", "log_IL1a", "log_TNFRI", "log_TNFRII", "log_IFN")

#List of MFI subscales
mfi_subscales <- c("MFI_AM", "MFI_GV", "MFI_LV", "MFI_VA", "MFI_VM")

#Make a grouping variable
filtered_data <- data %>%
  mutate(group = case_when(
    norm == 1 ~ "norm_1",
    norm == 0 & Time == 1 ~ "norm0_time1",
    norm == 0 & Time == 2 ~ "norm0_time2",
    norm == 0 & Time == 3 ~ "norm0_time3",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))

run_lm_biomarker_predicts_fatigue <- function(data, mfi_var, biomarker, group_name) {
  
  base_formula <- paste0(mfi_var, " ~ ", biomarker)
  
  #add covariates based on group membership
  if (group_name == "norm_1") {
    covariates <- c("geslacht", "leeft_W1", "BMI", "aant_comorb", 
                    "PAIN", "alcohol", "ROOK", "PSQI_total")
  } else {
    covariates <- c("geslacht", "leeft_W1", "BMI", "aant_comorb", 
                    "treatment", "PAIN", "alcohol", "ROOK", "PSQI_total", "time_since_diagnosis_W1")
  }
  
  full_formula <- paste(base_formula, paste(covariates, collapse = " + "), sep = " + ")
  
  # Fit the model
  model <- lm(formula = as.formula(full_formula), data = data)
  
  # tidy output
  tidy_model <- broom::tidy(model)
  return(tidy_model)
}

#export to excel
library(dplyr)
library(purrr)
library(tidyr)
library(writexl)

# all groups together
groups <- unique(filtered_data$group)

model_input <- expand.grid(
  group = groups,
  mfi_var = mfi_subscales,
  biomarker = biomarkers,
  stringsAsFactors = FALSE
)

# all regression equations in one dataframe
model_results_lm <- model_input %>%
  mutate(
    result = purrr::pmap(list(group, mfi_var, biomarker), function(g, mfi, bio) {
      df_sub <- filtered_data %>% filter(group == g)
      tryCatch(
        run_lm_biomarker_predicts_fatigue(df_sub, mfi, bio, g),
        error = function(e) {
          message("Fout bij: ", g, ", ", mfi, ", ", bio, " - ", e$message)
          return(NULL)
        }
      )
    })
  ) %>%
  unnest(cols = c(result))

write_xlsx(model_results_lm, "modelresultaten_biomarkers_fatigue.xlsx")

#put this in a table with beta (estimate) and SD to calculate z-test
library(dplyr)
library(tidyr)

df <- model_tabel5_resultaten_biomarkers_fatigue

biomarkers <- c("log_CRP", "log_IL6", "log_IL8", "log_IL10", "log_IL22", 
                "log_IL17A", "log_IL1b", "log_IL1a", "log_TNFRI", "log_TNFRII", "log_IFN")

df <- df %>% filter(term %in% biomarkers)

df_wide <- df %>%
  pivot_wider(
    id_cols = c("mfi_var", "biomarker", "term"),
    names_from = group,
    values_from = c("estimate (beta)", std.error),
    names_glue = "{.value}_{group}"
  )

groups_to_compare <- c("norm0_time1", "norm0_time2", "norm0_time3")

for (grp in groups_to_compare) {
  diff_col <- paste0("diff_", grp)
  z_col <- paste0("z_", grp)
  p_col <- paste0("p_", grp)
  
  est_grp <- paste0("estimate_", grp)
  se_grp <- paste0("std.error_", grp)

  if(all(c(est_grp, se_grp, "estimate_norm_1", "std.error_norm_1") %in% colnames(df_wide))) {
    df_wide <- df_wide %>%
      mutate(
        !!diff_col := .data[[est_grp]] - .data[["estimate_norm_1"]],
        !!z_col := (!!sym(diff_col)) / sqrt(.data[[se_grp]]^2 + .data[["std.error_norm_1"]]^2),
        !!p_col := 2 * pnorm(-abs(.data[[z_col]]))
      )
  } else {
    warning(paste("Kolommen voor groep", grp, "ontbreken in data"))
  }
}

library(writexl)
write_xlsx(df_wide, "vergelijking_norm1_vs_andere_groups.xlsx")

print(df_wide)

#now i have all estimates (betas) and standard errors in a row per group, but i want the p values of the z-tests

library(readxl)
library(dplyr)
library(tidyr)
library(writexl)

datanorm <- tabel5_vergelijking_norm1_vs_andere_groups

biomarkers <- c("log_CRP", "log_IL6", "log_IL8", "log_IL10", "log_IL22", 
                "log_IL17A", "log_IL1b", "log_IL1a", "log_TNFRI", "log_TNFRII", "log_IFN")

#rename columns without (beta) in it, because that is hard with the code
names(datanorm) <- names(datanorm) %>%
  gsub("estimate \\(beta\\)", "estimate", .) %>%
  gsub(" ", "_", .) 

groups_to_compare <- c("norm0_time1", "norm0_time2", "norm0_time3")

for (grp in groups_to_compare) {
  est_grp <- paste0("estimate_", grp)
  se_grp <- paste0("std.error_", grp)
  est_ref <- "estimate_norm_1"
  se_ref <- "std.error_norm_1"
  
  diff_col <- paste0("diff_", grp)
  z_col <- paste0("z_", grp)
  p_col <- paste0("p_", grp)
  
  datanorm <- datanorm %>%
    mutate(
      !!diff_col := .data[[est_grp]] - .data[[est_ref]],
      !!z_col := .data[[diff_col]] / sqrt(.data[[se_grp]]^2 + .data[[se_ref]]^2),
      !!p_col := 2 * pnorm(-abs(.data[[z_col]]))
    )
}

#filter results on significance
df_significant <- datanorm %>%
  filter(p_norm0_time1 < 0.05 | p_norm0_time2 < 0.05 | p_norm0_time3 < 0.05)

write_xlsx(datanorm, "tabel3_ztesten_output.xlsx")
write_xlsx(df_significant, "tabel3_significante_resultaten.xlsx")

#table output
library(dplyr)
library(tidyr)
library(writexl)

apa_round_p <- function(p) {
  ifelse(is.na(p), NA,
         ifelse(p < .001, "< .001", sprintf("%.3f", p)))
}

final_wide_clean <- datanorm %>%
  transmute(
    Biomarker   = biomarker,
    MFI_Subscale = mfi_var,
    Norm_beta   = round(estimate_norm_1, 2),
    Norm_SE     = round(std.error_norm_1, 2),
    
    T1_beta     = round(estimate_norm0_time1, 2),
    T1_SE       = round(std.error_norm0_time1, 2),
    p_T1        = apa_round_p(p_norm0_time1),
    
    T2_beta     = round(estimate_norm0_time2, 2),
    T2_SE       = round(std.error_norm0_time2, 2),
    p_T2        = apa_round_p(p_norm0_time2),
    
    T3_beta     = round(estimate_norm0_time3, 2),
    T3_SE       = round(std.error_norm0_time3, 2),
    p_T3        = apa_round_p(p_norm0_time3)
  ) %>%
  arrange(Biomarker, MFI_Subscale)

write_xlsx(final_wide_clean, "tabel3_eindversie_schoon.xlsx")

#tabel in format to word
# ===================== PACKAGES =====================
library(dplyr)
library(tidyr)
library(stringr)
library(flextable)
library(officer)
library(writexl)

# ===================== HELPERS & LABELS =====================
apa_p <- function(p) ifelse(is.na(p), NA,
                            ifelse(p < .001, "< .001", sprintf("%.3f", p)))

time_labels <- c("norm0_time1"="T0", "norm0_time2"="T1", "norm0_time3"="T2")

nice_bio <- c(
  "z_log_CRP"   = "CRP (mg/L)",
  "z_log_IFN"   = "IFN-γ (pg/mL)",
  "z_log_IL1a"  = "IL-1α (pg/mL)",
  "z_log_IL1b"  = "IL-1β (pg/mL)",
  "z_log_IL6"   = "IL-6 (pg/mL)",
  "z_log_IL8"   = "IL-8 (pg/mL)",
  "z_log_IL10"  = "IL-10 (pg/mL)",
  "z_log_IL17A" = "IL-17A (pg/mL)",
  "z_log_IL22"  = "IL-22 (pg/mL)",
  "z_log_TNFRI" = "sTNFRI (ng/mL)",
  "z_log_TNFRII"= "sTNFRII (ng/mL)"
)
nice_mfi <- c(
  "MFI_AM" = "General fatigue",
  "MFI_GV" = "Physical fatigue",
  "MFI_LV" = "Mental fatigue",
  "MFI_VA" = "Reduced activity",
  "MFI_VM" = "Reduced motivation"
)

#add lables
dn_lab <- datanorm %>%
  mutate(
    biomarker_lab = dplyr::recode(biomarker, !!!nice_bio),
    mfi_lab       = dplyr::recode(mfi_var, !!!nice_mfi)
  )

#long table with per (biomarker x subscale x timepont) the beta, SE and p
dn_long <- dn_lab %>%
  transmute(
    biomarker, biomarker_lab, mfi_var, mfi_lab,
    T1_beta = as.numeric(estimate_norm0_time1),
    T1_SE   = as.numeric(std.error_norm0_time1),
    T1_p    = as.numeric(p_norm0_time1),
    T2_beta = as.numeric(estimate_norm0_time2),
    T2_SE   = as.numeric(std.error_norm0_time2),
    T2_p    = as.numeric(p_norm0_time2),
    T3_beta = as.numeric(estimate_norm0_time3),
    T3_SE   = as.numeric(std.error_norm0_time3),
    T3_p    = as.numeric(p_norm0_time3)
  ) %>%
  pivot_longer(starts_with("T"),
               names_to = c("tp","stat"), names_sep = "_", values_to = "val") %>%
  mutate(tp = dplyr::recode(tp,
                            "T1" = time_labels["norm0_time1"],
                            "T2" = time_labels["norm0_time2"],
                            "T3" = time_labels["norm0_time3"])) %>%
  pivot_wider(names_from = stat, values_from = val) %>%
  mutate(
    beta = round(as.numeric(beta), 2),
    SE   = round(as.numeric(SE), 2),
    p    = apa_p(as.numeric(p)),
    tp   = factor(tp, levels = unname(time_labels))  # T0->T1->T2
  ) %>%
  arrange(biomarker_lab, mfi_lab, tp) %>%
  rename(Time = tp)

make_biomarker_table <- function(df_bio, title){
  tab <- df_bio %>%
    select(`MFI subscale` = mfi_lab, Time, `β` = beta, SE, `pᵇ` = p)
  
  ft <- flextable(tab)
  ft <- theme_booktabs(ft)
  ft <- autofit(ft)
  ft <- align(ft, align = "center", part = "all")
  ft <- align(ft, j = 1, align = "left", part = "all")
  ft <- merge_v(ft, j = "MFI subscale")
  ft <- valign(ft, j = "MFI subscale", valign = "top")
  
  ft <- add_header_lines(ft, values = title)
  ft <- bold(ft, part = "header")
  ft
}

#generate the tables
bio_list <- unique(dn_long$biomarker_lab)

tables <- lapply(bio_list, function(bnm){
  df_bio <- dn_long %>% filter(biomarker_lab == bnm)
  make_biomarker_table(df_bio, title = bnm)
})
names(tables) <- paste0("Table: ", bio_list)

#export to word
doc <- read_docx()
for (nm in names(tables)) {
  doc <- body_add_par(doc, nm, style = "heading 2")
  doc <- body_add_flextable(doc, tables[[nm]])
  doc <- body_add_par(doc, "")  # lege regel
}
print(doc, target = "Tables_per_biomarker_PROCORE.docx")

excel_sheets <- lapply(bio_list, function(bnm){
  dn_long %>%
    filter(biomarker_lab == bnm) %>%
    select(`MFI subscale` = mfi_lab, Time, `β` = beta, SE, `p` = p)
})
names(excel_sheets) <- bio_list
write_xlsx(excel_sheets, path = "Tables_per_biomarker_PROCORE.xlsx")

#different table for overview
library(dplyr)
library(tidyr)
library(flextable)
library(officer)
library(writexl)

fmt_num2 <- function(x, digits = 2){
  ifelse(is.na(x), "",
         paste0(ifelse(x < 0, "-", ""),
                sub("^0", "", formatC(abs(x), format = "f", digits = digits))))
}
# p significant if < .034 (FDR correction)
fmt_p3 <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < .001, "< .001", sub("^0", "", sprintf("%.3f", p))))
}

time_labels <- c("norm0_time1"="T0", "norm0_time2"="T1", "norm0_time3"="T2")

# labels left
nice_mfi <- c(
  "MFI_AM"="General fatigue",
  "MFI_GV"="Physical fatigue",
  "MFI_LV"="Mental fatigue",
  "MFI_VA"="Reduced activity",
  "MFI_VM"="Reduced motivation"
)

# biomarker
nice_bio_head <- c(
  "z_log_CRP"   = "CRP (mg/L)ᵃ",
  "z_log_IFN"   = "IFN-γ (pg/mL)ᵃ",
  "z_log_IL1a"  = "IL-1α (pg/mL)ᵃ",
  "z_log_IL1b"  = "IL-1β (pg/mL)ᵃ",
  "z_log_IL6"   = "IL-6 (pg/mL)ᵃ",
  "z_log_IL8"   = "IL-8 (pg/mL)ᵃ",
  "z_log_IL10"  = "IL-10 (pg/mL)ᵃ",
  "z_log_IL17A" = "IL-17A (pg/mL)ᵃ",
  "z_log_IL22"  = "IL-22 (pg/mL)ᵃ",
  "z_log_TNFRI" = "sTNFRI (ng/mL)ᵃ",
  "z_log_TNFRII"= "sTNFRII (ng/mL)ᵃ"
)

dn_long <- datanorm %>%
  mutate(mfi_lab = dplyr::recode(mfi_var, !!!nice_mfi)) %>%
  transmute(
    biomarker, mfi_lab,
    T1_beta = as.numeric(estimate_norm0_time1),
    T1_SE   = as.numeric(std.error_norm0_time1),
    T1_p    = as.numeric(p_norm0_time1),
    T2_beta = as.numeric(estimate_norm0_time2),
    T2_SE   = as.numeric(std.error_norm0_time2),
    T2_p    = as.numeric(p_norm0_time2),
    T3_beta = as.numeric(estimate_norm0_time3),
    T3_SE   = as.numeric(std.error_norm0_time3),
    T3_p    = as.numeric(p_norm0_time3)
  ) %>%
  pivot_longer(starts_with("T"),
               names_to = c("tp","stat"), names_sep = "_", values_to = "val") %>%
  mutate(tp = dplyr::recode(tp,
                            "T1" = time_labels["norm0_time1"],
                            "T2" = time_labels["norm0_time2"],
                            "T3" = time_labels["norm0_time3"])) %>%
  pivot_wider(names_from = stat, values_from = val) %>%
  mutate(tp = factor(tp, levels = unname(time_labels)))  # ordening T0->T1->T2

#tablebuilder
make_panel_table <- function(dn_long, bio_order, title_top = "Inflammatory marker") {
  
  wide <- dn_long %>%
    filter(biomarker %in% bio_order) %>%
    select(`MFI subscale` = mfi_lab, Time = tp, biomarker, beta, SE, p) %>%
    pivot_longer(c(beta, SE, p), names_to = "stat", values_to = "val") %>%
    unite("bio_stat", biomarker, stat, sep = "__") %>%
    pivot_wider(names_from = bio_stat, values_from = val)
  
  col_seq <- c("MFI subscale", "Time",
               as.vector(t(outer(bio_order, c("β","SE","pᵇ"), paste, sep="__"))))
  # map: beta->β, p->pᵇ
  names(wide) <- sub("__beta$", "__β", names(wide))
  names(wide) <- sub("__p$",   "__pᵇ", names(wide))
  miss <- setdiff(col_seq, names(wide))
  if (length(miss)) for (m in miss) wide[[m]] <- NA_real_
  wide <- wide[, col_seq]
  
  beta_cols <- grep("__β$",  names(wide), value = TRUE)
  se_cols   <- grep("__SE$", names(wide), value = TRUE)
  p_cols    <- grep("__pᵇ$", names(wide), value = TRUE)
  
  wide[beta_cols] <- lapply(wide[beta_cols], fmt_num2, digits = 2)
  wide[se_cols]   <- lapply(wide[se_cols],   fmt_num2, digits = 2)
  wide[p_cols]    <- lapply(wide[p_cols],    fmt_p3)
  
  col_keys <- names(wide)
  h1 <- c("MFI subscale", "", rep(title_top, length(col_keys) - 2))
  h2 <- c("", "", unlist(lapply(bio_order, function(b) rep(nice_bio_head[[b]], 3))))
  h3 <- c("", "", rep(c("β","SE","pᵇ"), length(bio_order)))
  header_df <- data.frame(col_keys = col_keys, h1 = h1, h2 = h2, h3 = h3, check.names = FALSE)
  
  ft <- flextable(wide, col_keys = col_keys)
  ft <- set_header_df(ft, mapping = header_df, key = "col_keys")
  ft <- merge_h(ft, part = "header")  
  ft <- theme_booktabs(ft); ft <- autofit(ft)
  ft <- align(ft, part="all", align="center"); ft <- align(ft, j=1, part="all", align="left")
  ft <- merge_v(ft, j = "MFI subscale"); ft <- valign(ft, j="MFI subscale", valign="top")
  
  ft <- border_remove(ft)
  b <- fp_border(color = "#999999", width = 0.75)
  ft <- hline_top(ft, part="header", border=b)
  ft <- hline(ft, i=1, part="header", border=b)
  ft <- hline(ft, i=2, part="header", border=b)
  ft <- hline_bottom(ft, part="header", border=b)
  ft <- bold(ft, part = "header")
  ft
}

#panels
panel1_bios <- c("z_log_CRP","z_log_IFN","z_log_IL1a","z_log_IL1b","z_log_IL6","z_log_IL8")
panel2_bios <- c("z_log_IL10","z_log_IL17A","z_log_IL22","z_log_TNFRI","z_log_TNFRII")

ft_panel1 <- make_panel_table(dn_long, bio_order = panel1_bios)
ft_panel2 <- make_panel_table(dn_long, bio_order = panel2_bios)

#word export
doc <- read_docx()
doc <- body_add_par(doc, "Table. Comparison of biomarker–fatigue associations (PROCORE vs. norm)", style = "heading 2")
doc <- body_add_flextable(doc, ft_panel1)
doc <- body_add_par(doc, "")
doc <- body_add_flextable(doc, ft_panel2)
doc <- body_add_par(doc, "")
doc <- body_add_par(doc,
                    "Notes. ᵃ All biomarker concentrations are log- and z-transformed. ᵇ p-values are from z-tests comparing PROCORE slopes to normative slopes, adjusted for covariates.",
                    style = "Normal")
print(doc, target = "Table_biomarker_fatigue_panels_nozero.docx")

