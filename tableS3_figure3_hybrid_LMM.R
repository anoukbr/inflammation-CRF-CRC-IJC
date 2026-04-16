# =============================================================================
# Table S3 & Figure 3: Hybrid Linear Mixed Models — Inflammation and Fatigue
# =============================================================================
# Paper: Longitudinal Associations Between Inflammation and Multi-Dimensional Fatigue up to Two Years after Colorectal Cancer Diagnosis
# Authors: Bruijnzeels, A.E.C., Mols, F., Van Deun, K., Schoormans, D.
# Output: Table S3 — Full results of overall, between-subject, and
#         within-subject LMM associations between each log-transformed
#         inflammatory biomarker and each MFI subscale in CRC patients.
#         Figure 3 — Forest plots of standardized regression coefficients (β)
#         with 95% CIs, separated by overall, between-subject, and
#         within-subject effects.
#
# Model: Hybrid LMM approach:
#   - Overall model: biomarker as time-varying predictor
#   - Hybrid model: person-mean (between-subject) + person-mean-centered
#     deviation (within-subject) entered simultaneously
# Covariates: sex, age, BMI, alcohol use, smoking, comorbidities, treatment,
#             days since diagnosis, sleep quality (PSQI), pain (EORTC)
# FDR correction applied per MFI subscale across biomarkers.
#
# Data: PROCORE long-format dataset with inflammatory markers,
#       multiple imputation (1000 iterations, excl. DCRA).
#       Load your dataset below before running this script.
# =============================================================================

# Load your dataset here:
# filtered_data <- read.csv("your_data.csv")   # or: load("your_data.RData")
# Expected object name used in original analysis:
# X3_Inflammatoire_markers_x_behandeling_long_format_bloedvragenlijstje_zonder_DCRA_andere_imputatiemethode_1000_iteraties
# Subset applied: CRP_c > 0 & responder == 1

# install.packages(c("lme4","lmerTest","dplyr","purrr","tidyr","broom.mixed","ggplot2","forcats","openxlsx","stringr"))

library(lme4)
library(lmerTest)
library(dplyr)
library(purrr)
library(tidyr)
library(broom.mixed)
library(ggplot2)
library(forcats)
library(openxlsx)
library(stringr)

filtered_data <- subset(X3_Inflammatoire_markers_x_behandeling_long_format_bloedvragenlijstje_zonder_DCRA_andere_imputatiemethode_1000_iteraties, CRP_c > 0 & responder == 1)

filtered_data$log_CRP <- log(filtered_data$CRP_c)
filtered_data$log_IL6 <- log(filtered_data$IL_6)
filtered_data$log_IL8 <- log(filtered_data$IL_8)
filtered_data$log_IL10 <- log(filtered_data$IL_10)
filtered_data$log_IL22 <- log(filtered_data$IL_22)
filtered_data$log_IL17A <- log(filtered_data$IL_17A)
filtered_data$log_IL1b <- log(filtered_data$IL_1_beta)
filtered_data$log_IL1a <- log(filtered_data$IL_1_alpha)
filtered_data$log_TNFRI <- log(filtered_data$TNFRI)
filtered_data$log_TNFRII <- log(filtered_data$TNFRII)
filtered_data$log_IFNy <- log(filtered_data$IFNgamma)


set.seed(1)
dat <- filtered_data %>%
  mutate(
    Time = factor(Time, levels = c("1","2","3"), ordered = TRUE),
    treatment = factor(treatment),
    geslacht = factor(geslacht),
    alcohol   = factor(alcohol),
    ROOK      = factor(ROOK),
    aant_comorb = factor(aant_comorb)
  )

dat <- filtered_data %>%
  mutate(
    BMI_categorie = cut(BMI,
                        breaks = c(-Inf, 18.5, 24.9, 29.9, Inf),
                        labels = c("Underweight","Healthy weight","Overweight","Obese")) %>%
      relevel(ref = "Healthy weight")
  )

biomarkers <- c("log_CRP","log_IL6","log_IL8","log_IL10","log_IL22",
                "log_IL17A","log_IL1b","log_IL1a","log_TNFRI","log_TNFRII","log_IFNy")

mfi_subscales <- c("MFI_AM","MFI_GV","MFI_LV","MFI_VA","MFI_VM")

mfi_labels <- c(
  MFI_GV = "General Fatigue",
  MFI_LV = "Physical Fatigue",
  MFI_VA = "Reduced Activity",
  MFI_VM = "Reduced Motivation",
  MFI_AM = "Mental Fatigue"
)

biomarker_labels <- c(
  log_CRP   = "CRP (mg/L)",
  log_IL6   = "IL-6 (pg/mL)",
  log_IL8   = "IL-8 (pg/mL)",
  log_IL10  = "IL-10 (pg/mL)",
  log_IL22  = "IL-22 (pg/mL)",
  log_IL17A = "IL-17A (pg/mL)",
  log_IL1b  = "IL-1β (pg/mL)",
  log_IL1a  = "IL-1α (pg/mL)",
  log_TNFRI = "TNF-RI (pg/mL)",
  log_TNFRII= "TNF-RII (ng/mL)",
  log_IFNy  = "IFN-γ (pg/mL)"
)

# make hybrid variables (pm = between, wc = within)
make_within_between <- function(df, id_var = "patID", bio) {
  pm   <- paste0(bio, "_pm")
  wc   <- paste0(bio, "_wc")
  df %>%
    group_by(.data[[id_var]]) %>%
    mutate(
      !!pm := mean(.data[[bio]], na.rm = TRUE),
      !!wc := .data[[bio]] - .data[[pm]]
    ) %>%
    ungroup()
}

for (b in biomarkers) {
  dat <- make_within_between(dat, id_var = "patID", bio = b)
}

# ----------------
# MODELFUNCTIONS!!!
# ----------------

#Hybrid LMM
fit_hybrid_lmm <- function(data, mfi_var, bio) {
  bio_pm <- paste0(bio, "_pm")
  bio_wc <- paste0(bio, "_wc")
  fml <- paste0(
    mfi_var, " ~ (", bio_wc, " + ", bio_pm, ") + Time + ",
    "geslacht + leeft_W1 + BMI_categorie + aant_comorb + treatment + ",
    "PAIN + alcohol + ROOK + PSQI_total + time_since_diagnosis_W1 + ",
    "(1|patID)"
  )
  mod <- lmer(
    formula = as.formula(fml),
    data = data,
    REML = TRUE
  )
  out <- broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE)
  out$mfi_var  <- mfi_var
  out$biomarker <- bio
  out
}

#Overall LMM (without splitted effects, overall effects, 'normal' LMM)
fit_overall_lmm <- function(data, mfi_var, bio) {
  fml <- paste0(
    mfi_var, " ~ ", bio, " + Time + ",
    "geslacht + leeft_W1 + BMI_categorie + aant_comorb + treatment + ",
    "PAIN + alcohol + ROOK + PSQI_total + time_since_diagnosis_W1 + ",
    "(1|patID)"
  )
  mod <- lmer(as.formula(fml), data = data, REML = TRUE)
  out <- broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE)
  out$mfi_var <- mfi_var; out$biomarker <- bio
  out
}

#run models and combine results
hybrid_results <- purrr::cross_df(list(
  mfi_var = mfi_subscales,
  biomarker = biomarkers
)) %>%
  pmap_dfr(~ fit_hybrid_lmm(dat, ..1, ..2))

# Label: within / between 
hybrid_results <- hybrid_results %>%
  mutate(
    effect_type = case_when(
      str_detect(term, "_wc$") ~ "within",
      str_detect(term, "_pm$") ~ "between",
      TRUE ~ "other"
    ),
    pretty_term = case_when(
      effect_type == "within"  ~ "Within (deviation from person mean)",
      effect_type == "between" ~ "Between (person mean)",
      TRUE ~ term
    )
  )

#Overall results
overall_results_raw <- purrr::cross_df(list(mfi_var = mfi_subscales, biomarker = biomarkers)) %>%
  pmap_dfr(~ fit_overall_lmm(dat, ..1, ..2))
overall_results <- overall_results_raw %>%
  filter(term == biomarker) %>%
  mutate(effect_type = "overall") %>%
  select(mfi_var, biomarker, term, estimate, std.error, conf.low, conf.high, p.value, effect_type)

# --------
#PLOTS !!!
# --------

#combine overall and splitted results with separate plots for all and one combined forest plot (which will be figure 3 in the paper)
hyb_keep <- hybrid_results %>%
  filter(effect_type %in% c("between","within")) %>%
  select(mfi_var, biomarker, term, estimate, std.error, conf.low, conf.high, p.value, effect_type)

plot_df <- bind_rows(hyb_keep, overall_results) %>%
  mutate(
    effect_type = factor(effect_type,
                         levels = c("overall","between","within"),
                         labels = c("Overall association","Between (person mean)","Within (deviation)")),
    biomarker = factor(biomarker, levels = rev(biomarkers)),
    mfi_var_label = recode(mfi_var, !!!mfi_labels),
    sig_symbol = case_when(
      p.value < 0.034 ~ "\u2020",
      TRUE ~ ""
    ),
    significant = p.value < 0.034
  )

#separate plot per MFI subscale
make_one_mfi_plot <- function(df, one_mfi, file_tag = NULL, include_overall = TRUE) {
  use_df <- df %>% filter(mfi_var == one_mfi)
  if (!include_overall) {
    use_df <- use_df %>% filter(effect_type %in% c("Between (person mean)","Within (deviation)"))
  }
  p <- ggplot(use_df, aes(x = estimate, y = biomarker)) +
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
    geom_text(aes(label = sig_symbol),
              vjust = -0.7, size = 3, na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(x = "Effect (beta, 95% CI)", y = NULL,
         title = paste0("Associations — ", unique(use_df$mfi_var_label))) +
    facet_grid(rows = vars(mfi_var_label), cols = vars(effect_type), scales = "free_y") +
    scale_y_discrete(labels = biomarker_labels) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text.y   = element_text(angle = 0, hjust = 0.5),
      panel.spacing.x= unit(10, "pt"),
      panel.spacing.y= unit(8, "pt")
    )
  if (is.null(file_tag)) file_tag <- one_mfi
  width <- if (include_overall) 11 else 8
  ggsave(paste0("forest_", file_tag, ".png"), p, width = width, height = 4.6, dpi = 300)
  return(p)
}

plots_with_overall <- lapply(mfi_subscales, function(mfi) {
  make_one_mfi_plot(plot_df, one_mfi = mfi, file_tag = mfi, include_overall = TRUE)
})

# --------------------
#BIG OVERVIEW PLOT (FIGURE 3)
# --------------------

forest_three_panel <- ggplot(plot_df, aes(x = estimate, y = biomarker, color = significant)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
  geom_text(aes(label = sig_symbol),
            vjust = -0.7, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#c1327b"),
                     labels = c("Non-significant", "Significant (p < .034)"),
                     name = NULL) +
  labs(x = "Effect (\u03B2, 95% CI)", y = NULL,
       title = "Associations per MFI subscale: Overall vs Between vs Within") +
  facet_grid(rows = vars(mfi_var_label), cols = vars(effect_type), scales = "free_y") +
  scale_y_discrete(labels = biomarker_labels) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0.5),
    panel.spacing.x = unit(10, "pt"),
    panel.spacing.y = unit(16, "pt"),
    plot.title = element_text(size = 11),
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  )

print(forest_three_panel)
showtext_opts(dpi = 1000)
ggsave("forest_plot_ARIAL.png", forest_three_panel,
       width = 10, height = 14, dpi = 1000)

# --------------------------------
#SUPPLEMENT TABLES EXPORT!
# --------------------------------

# 1) Full fixed-effects from hybrid models
write.xlsx(hybrid_results, file = "Supplement_Full_FixedEffects_HYBRID.xlsx")

# 2) Compact tables: within / between
tab_between <- hybrid_results %>%
  filter(effect_type == "between") %>%
  select(mfi_var, biomarker, term, estimate, std.error, conf.low, conf.high, p.value)

tab_within <- hybrid_results %>%
  filter(effect_type == "within") %>%
  select(mfi_var, biomarker, term, estimate, std.error, conf.low, conf.high, p.value)

# 3) Compact table: overall
tab_overall <- overall_results %>%
  select(mfi_var, biomarker, term, estimate, std.error, conf.low, conf.high, p.value)

wb <- createWorkbook()
addWorksheet(wb, "Between_effects")
addWorksheet(wb, "Within_effects")
addWorksheet(wb, "Overall_effects")
writeData(wb, "Between_effects", tab_between)
writeData(wb, "Within_effects",  tab_within)
writeData(wb, "Overall_effects", tab_overall)
saveWorkbook(wb, "Supplement_WithinBetweenOverall_Tables.xlsx", overwrite = TRUE)

#install.packages(c("gt","scales"))
library(gt)
library(scales)

fmt_num <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

# Combined table wide with columns: Overall β (95%CI), p ; Between β (95%CI), p ; Within β (95%CI), p
make_apa_table_data <- function(plot_df, one_mfi, digits_beta = 2, digits_p = 3) {
  df <- plot_df %>% dplyr::filter(mfi_var == one_mfi)
  
  # Overall
  ov <- df %>%
    dplyr::filter(effect_type == "Overall association") %>%
    dplyr::select(biomarker, estimate, conf.low, conf.high, p.value) %>%
    dplyr::rename(
      ov_est = estimate, ov_low = conf.low, ov_high = conf.high, ov_p = p.value
    )
  
  # Between
  be <- df %>%
    dplyr::filter(effect_type == "Between (person mean)") %>%
    dplyr::select(biomarker, estimate, conf.low, conf.high, p.value) %>%
    dplyr::rename(
      be_est = estimate, be_low = conf.low, be_high = conf.high, be_p = p.value
    )
  
  # Within
  wi <- df %>%
    dplyr::filter(effect_type == "Within (deviation)") %>%
    dplyr::select(biomarker, estimate, conf.low, conf.high, p.value) %>%
    dplyr::rename(
      wi_est = estimate, wi_low = conf.low, wi_high = conf.high, wi_p = p.value
    )
  
  out <- ov %>%
    dplyr::full_join(be, by = "biomarker") %>%
    dplyr::full_join(wi, by = "biomarker") %>%
    # Zorg voor vaste biomarker-volgorde
    dplyr::mutate(biomarker = factor(biomarker, levels = rev(biomarkers))) %>%
    dplyr::arrange(biomarker) %>%
    # Maak APA-stijl strings
    dplyr::mutate(
      overall_beta_ci = ifelse(!is.na(ov_est),
                               paste0(fmt_num(ov_est, digits_beta), " [", fmt_num(ov_low, digits_beta), ", ", fmt_num(ov_high, digits_beta), "]"),
                               NA
      ),
      between_beta_ci = ifelse(!is.na(be_est),
                               paste0(fmt_num(be_est, digits_beta), " [", fmt_num(be_low, digits_beta), ", ", fmt_num(be_high, digits_beta), "]"),
                               NA
      ),
      within_beta_ci = ifelse(!is.na(wi_est),
                              paste0(fmt_num(wi_est, digits_beta), " [", fmt_num(wi_low, digits_beta), ", ", fmt_num(wi_high, digits_beta), "]"),
                              NA
      ),
      overall_p = ifelse(!is.na(ov_p), scales::pvalue(ov_p, accuracy = 10^(-digits_p)), NA),
      between_p = ifelse(!is.na(be_p), scales::pvalue(be_p, accuracy = 10^(-digits_p)), NA),
      within_p  = ifelse(!is.na(wi_p), scales::pvalue(wi_p, accuracy = 10^(-digits_p)), NA)
    ) %>%
    dplyr::select(
      biomarker,
      `Overall β [95% CI]` = overall_beta_ci, `Overall p` = overall_p,
      `Between β [95% CI]` = between_beta_ci, `Between p` = between_p,
      `Within β [95% CI]`  = within_beta_ci,  `Within p`  = within_p
    )
  
  return(out)
}

#save as png
save_apa_table_png <- function(tab_df, mfi_code, mfi_label,
                               file_prefix = "Table_MFI_", width_px = 1800) {
  
  gt_tab <- tab_df %>%
    gt(rowname_col = "biomarker") %>%
    tab_header(
      title = md(paste0("**", mfi_label, "** — Overall vs Between vs Within"))
    ) %>%
    cols_label(.list = list(
      `Overall β [95% CI]` = md("Overall β [95% CI]"),
      `Overall p` = md("Overall *p*"),
      `Between β [95% CI]` = md("Between β [95% CI]"),
      `Between p` = md("Between *p*"),
      `Within β [95% CI]` = md("Within β [95% CI]"),
      `Within p` = md("Within *p*")
    )) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels(everything())
    ) %>%
    tab_options(
      table.font.size = px(12),
      data_row.padding = px(6)
    )
  
  # save as png
  out_file <- paste0(file_prefix, mfi_code, ".png")
  gtsave(gt_tab, filename = out_file, vwidth = width_px)
  message("Saved: ", out_file)
  
  return(gt_tab)
}

#Separate tables for overall, within and between effects separately
library(gt)
library(scales)

biomarker_labels <- c(
  log_CRP   = "CRP (log)",
  log_IL6   = "IL-6 (log)",
  log_IL8   = "IL-8 (log)",
  log_IL10  = "IL-10 (log)",
  log_IL22  = "IL-22 (log)",
  log_IL17A = "IL-17A (log)",
  log_IL1b  = "IL-1β (log)",
  log_IL1a  = "IL-1α (log)",
  log_TNFRI = "TNF-RI (log)",
  log_TNFRII= "TNF-RII (log)",
  log_IFNy  = "IFN-γ (log)"
)
fmt_num <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

make_effect_table_long <- function(plot_df, effect_label,
                                   digits_beta = 2, digits_p = 3) {
  
  plot_df %>%
    dplyr::filter(effect_type == effect_label) %>%
    dplyr::mutate(
      biomarker_lbl = dplyr::recode(biomarker, !!!biomarker_labels),
      mfi_lbl = dplyr::recode(mfi_var, !!!mfi_labels)
    ) %>%
    dplyr::arrange(mfi_lbl, biomarker_lbl) %>%
    dplyr::mutate(
      beta_ci = paste0(fmt_num(estimate, digits_beta),
                       " [", fmt_num(conf.low, digits_beta),
                       ", ", fmt_num(conf.high, digits_beta), "]"),
      p_fmt = scales::pvalue(p.value, accuracy = 10^(-digits_p))
    ) %>%
    dplyr::select(
      `MFI subscale` = mfi_lbl,
      Biomarker = biomarker_lbl,
      `β [95% CI]` = beta_ci,
      `p` = p_fmt
    )
}

save_effect_table_png <- function(tab_df, title_text,
                                  file_name, width_px = 1800) {
  gt_tab <- tab_df %>%
    gt() %>%
    tab_header(title = md(paste0("**", title_text, "**"))) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels(everything())
    ) %>%
    cols_align(
      align = "left",
      columns = c(`MFI subscale`, Biomarker)
    ) %>%
    tab_options(
      table.font.size = px(12),
      data_row.padding = px(6)
    )
  
  gtsave(gt_tab, filename = file_name, vwidth = width_px)
  message("Saved: ", file_name)
  invisible(gt_tab)
}

# 1) OVERALL-tabel
overall_tab <- make_effect_table_long(plot_df, "Overall association")
save_effect_table_png(overall_tab,
                      title_text = "Overall associations (biomarker → MFI)",
                      file_name  = "Table_Overall_Associations.png")

# 2) BETWEEN-tabel
between_tab <- make_effect_table_long(plot_df, "Between (person mean)")
save_effect_table_png(between_tab,
                      title_text = "Between-person associations (person mean of biomarker → MFI)",
                      file_name  = "Table_Between_Associations.png")

# 3) WITHIN-tabel
within_tab <- make_effect_table_long(plot_df, "Within (deviation)")
save_effect_table_png(within_tab,
                      title_text = "Within-person associations (deviation from person mean of biomarker → MFI)",
                      file_name  = "Table_Within_Associations.png")

# All in one Excel file 
wb_out <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_out, "Overall")
openxlsx::addWorksheet(wb_out, "Between")
openxlsx::addWorksheet(wb_out, "Within")
openxlsx::writeData(wb_out, "Overall", overall_tab)
openxlsx::writeData(wb_out, "Between", between_tab)
openxlsx::writeData(wb_out, "Within", within_tab)
openxlsx::saveWorkbook(wb_out, "APA_Tables_Overall_Between_Within.xlsx", overwrite = TRUE)
# Session info for reproducibility
sessionInfo()
