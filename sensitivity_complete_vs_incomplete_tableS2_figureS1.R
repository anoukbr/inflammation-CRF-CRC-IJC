# =============================================================================
# Sensitivity Analysis: Complete vs. Incomplete Data (Table S2 & Figure S1)
# =============================================================================
# Paper: Longitudinal Associations Between Inflammation and Multi-Dimensional Fatigue up to Two Years after Colorectal Cancer Diagnosis
# Authors: Bruijnzeels, A.E.C., Mols, F., Van Deun, K., Schoormans, D.
# Output: Table S2 & Figure S1 — Comparison of MFI subscale scores between
#         patients with complete data at all timepoints (PROCORE 2, n = 209)
#         and patients with incomplete data (PROCORE 1, n = 202 at baseline).
#         Independent samples t-tests per subscale per timepoint.
#
# Data: Full PROCORE key file including inflammatory markers and blood
#       questionnaire or file with inflammatory markers and fatigue subscales. 
# Load your dataset below before running this script.
# =============================================================================

library(tidyverse)
library(dplyr)
library(ggplot2)
library(showtext)
library(openxlsx)
library(haven)

df <- X3_Inflammatoire_markers_x_behandeling_long_format_bloedvragenlijstje_zonder_DCRA_andere_imputatiemethode_1000_iteraties

df <- df %>%
  filter(CRP_c > 0, responder == 1)

summary_stats <- df %>%
  group_by(Time, compleet) %>%
  summarise(
    n = n(),
    
    MFI_AM_mean = mean(MFI_AM, na.rm = TRUE),
    MFI_AM_sd   = sd(MFI_AM, na.rm = TRUE),
    
    MFI_GV_mean = mean(MFI_GV, na.rm = TRUE),
    MFI_GV_sd   = sd(MFI_GV, na.rm = TRUE),
    
    MFI_LV_mean = mean(MFI_LV, na.rm = TRUE),
    MFI_LV_sd   = sd(MFI_LV, na.rm = TRUE),
    
    MFI_VA_mean = mean(MFI_VA, na.rm = TRUE),
    MFI_VA_sd   = sd(MFI_VA, na.rm = TRUE),
    
    MFI_VM_mean = mean(MFI_VM, na.rm = TRUE),
    MFI_VM_sd   = sd(MFI_VM, na.rm = TRUE),
    
    .groups = "drop"
  )

summary_stats

#Time needs to be ordered
df$Time <- factor(df$Time,
                  levels = c(1, 2, 3),
                  labels = c("Baseline", "12-month FU", "24-month FU"))

df_long <- df %>%
  select(Time, compleet, MFI_AM, MFI_GV, MFI_LV, MFI_VA, MFI_VM) %>%
  pivot_longer(
    cols = c(MFI_AM, MFI_GV, MFI_LV, MFI_VA, MFI_VM),
    names_to = "Scale",
    values_to = "Score"
  ) %>%
  mutate(
    Group = recode(as.character(compleet),
                   "0" = "Incomplete data",
                   "1" = "Complete data"),
    Scale = recode(Scale,
                   "MFI_AM" = "General fatigue",
                   "MFI_GV" = "Reduced activity",
                   "MFI_LV" = "Reduced motivation",
                   "MFI_VA" = "Physical fatigue",
                   "MFI_VM" = "Mental fatigue")
  )

df_summary <- df_long %>%
  group_by(Time, Scale, Group) %>%
  summarise(
    Mean = mean(Score, na.rm = TRUE),
    SD   = sd(Score, na.rm = TRUE),
    N    = sum(!is.na(Score)),
    .groups = "drop"
  )

df_stats <- df_summary %>%
  pivot_wider(
    names_from = Group,
    values_from = c(Mean, SD, N),
    names_glue = "{.value}_{Group}"
  ) %>%
  rename_with(~ gsub(" ", "_", .x)) %>%
  mutate(
    Mean_diff = Mean_Complete_data - Mean_Incomplete_data,
    SE_diff   = sqrt((SD_Complete_data^2 / N_Complete_data) +
                       (SD_Incomplete_data^2 / N_Incomplete_data)),
    t_value   = Mean_diff / SE_diff,
    df_welch  = ((SD_Complete_data^2 / N_Complete_data +
                    SD_Incomplete_data^2 / N_Incomplete_data)^2) /
      (((SD_Complete_data^2 / N_Complete_data)^2) / (N_Complete_data - 1) +
         ((SD_Incomplete_data^2 / N_Incomplete_data)^2) / (N_Incomplete_data - 1)),
    p_value   = 2 * pt(-abs(t_value), df_welch),
    tcrit     = qt(0.975, df_welch),
    CI_low    = Mean_diff - tcrit * SE_diff,
    CI_high   = Mean_diff + tcrit * SE_diff,
    Cohens_d  = Mean_diff / sqrt((SD_Complete_data^2 + SD_Incomplete_data^2) / 2)
  ) %>%
  select(Time, Scale,
         Mean_Complete_data, Mean_Incomplete_data,
         Mean_diff, SE_diff, t_value, df_welch, p_value, CI_low, CI_high, Cohens_d) %>%
  arrange(Scale, Time)

df_stats_print <- df_stats %>%
  mutate(across(c(Mean_Complete_data, Mean_Incomplete_data, Mean_diff, SE_diff,
                  t_value, df_welch, p_value, CI_low, CI_high, Cohens_d),
                ~ round(.x, 3)))

df_stats_print
write.xlsx(df_stats_print, "sensitivity_results_completeversusincomplete.xlsx")

df_plot <- df_summary %>%
  mutate(Group = recode(Group,
                        "PROCORE 1" = "Incomplete data",
                        "PROCORE 2" = "Complete data"))

# Line plot per subscale
ggplot(df_plot, aes(x = Time, y = Mean, color = Group, group = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_text(aes(label = paste0("n =", N), y = 5.1),
            position = position_dodge(width = 0.3),
            size = 2.8, show.legend = FALSE) +
  ylim(5, 11.5) +
  facet_wrap(~ Scale, scales = "free_y") +
  labs(
    title = "Fatigue scores by completeness of data",
    subtitle = "Incomplete data vs Complete data across timepoints",
    y = "Mean fatigue score",
    x = "Time",
    color = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

#Splitting the labels because they overlap
df_plot_n <- df_plot %>%
  mutate(y_label = ifelse(Group == "Complete data", 5.1, 5.5))

ggplot(df_plot, aes(x = Time, y = Mean, color = Group, group = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_text(data = df_plot_n,
            aes(label = paste0("n=", N), y = y_label),
            size = 3, show.legend = FALSE) +
  ylim(5, 11.5) +
  facet_wrap(~ Scale, scales = "free_y") +
  labs(
    title = "Fatigue scores by completeness of data",
    subtitle = "Incomplete data vs Complete data across timepoints",
    y = "Mean fatigue score",
    x = "Time",
    color = "Group"
  ) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

#Separate labels for each pane in the figure due to journal standards
panel_labels <- c(
  "General fatigue" = "(A)",
  "Mental fatigue" = "(B)",
  "Physical fatigue" = "(C)",
  "Reduced activity" = "(D)",
  "Reduced motivation" = "(E)"
)

p <- ggplot(df_plot, aes(x = Time, y = Mean, color = Group, group = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_text(data = df_plot_n,
            aes(label = paste0("n=", N), y = y_label),
            size = 3, show.legend = FALSE) +
  ylim(5, 11.5) +
  facet_wrap(~ Scale, scales = "free_y",
             labeller = labeller(Scale = panel_labels)) +
  labs(
    y = "Mean fatigue score",
    x = "Time",
    color = "Group"
  ) +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    legend.position = "bottom",
    strip.text = element_text(hjust = 0)
  )

ggsave("figureS1_correct.png", plot = p,
       width = 20, height = 14, units = "cm",
       dpi = 600) 
