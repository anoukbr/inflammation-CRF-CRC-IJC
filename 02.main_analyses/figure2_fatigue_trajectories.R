# =============================================================================
# Figure 2: Fatigue Trajectories Over Time — PROCORE vs. Normative Sample
# =============================================================================
# Paper: Longitudinal Associations Between Inflammation and Multi-Dimensional Fatigue up to Two Years after Colorectal Cancer Diagnosis
# Authors: Bruijnzeels, A.E.C., Mols, F., Van Deun, K., Schoormans, D.
# Output: Figure 2 — Five line plots showing mean MFI subscale scores (± SD)
#         for CRC patients (PROCORE) across three timepoints (Baseline,
#         12-month FU, 24-month FU), with normative reference values shown
#         as dashed lines.
#         Subscales: General Fatigue, Mental Fatigue, Physical Fatigue,
#         Reduced Activity, Reduced Motivation.
#
# Note: Means and SDs entered manually based on descriptive output from the
#       main dataset. No external dataset required for this script.
# =============================================================================

df <- tribble(
  ~Time, ~Scale, ~Group, ~Mean, ~SD,
  "Baseline", "General fatigue", "PROCORE", 8.8, 4.4,
  "12-month FU", "General fatigue", "PROCORE", 8.8, 4.2,
  "24-month FU", "General fatigue", "PROCORE", 8.7, 4.1,
  "Baseline", "General fatigue", "Norm", 8.7, NA,
  "12-month FU", "General fatigue", "Norm", 8.7, NA,
  "24-month FU", "General fatigue", "Norm", 8.7, NA,
  
  "Baseline", "Mental fatigue", "PROCORE", 8.4, 3.7,
  "12-month FU", "Mental fatigue", "PROCORE", 7.7, 3.7,
  "24-month FU", "Mental fatigue", "PROCORE", 7.9, 3.9,
  "Baseline", "Mental fatigue", "Norm", 7.2, NA,
  "12-month FU", "Mental fatigue", "Norm", 7.2, NA,
  "24-month FU", "Mental fatigue", "Norm", 7.2, NA,
  
  "Baseline", "Physical fatigue", "PROCORE", 8.8, 4.4,
  "12-month FU", "Physical fatigue", "PROCORE", 8.9, 4.3,
  "24-month FU", "Physical fatigue", "PROCORE", 8.7, 4.1,
  "Baseline", "Physical fatigue", "Norm", 8.5, NA,
  "12-month FU", "Physical fatigue", "Norm", 8.5, NA,
  "24-month FU", "Physical fatigue", "Norm", 8.5, NA,
  
  "Baseline", "Reduced activity", "PROCORE", 10.1, 4.3,
  "12-month FU", "Reduced activity", "PROCORE", 9.5, 4.3,
  "24-month FU", "Reduced activity", "PROCORE", 9.5, 4.1,
  "Baseline", "Reduced activity", "Norm", 8.8, NA,
  "12-month FU", "Reduced activity", "Norm", 8.8, NA,
  "24-month FU", "Reduced activity", "Norm", 8.8, NA,
  
  "Baseline", "Reduced motivation", "PROCORE", 9.0, 3.9,
  "12-month FU", "Reduced motivation", "PROCORE", 8.4, 3.5,
  "24-month FU", "Reduced motivation", "PROCORE", 8.8, 3.7,
  "Baseline", "Reduced motivation", "Norm", 8.0, NA,
  "12-month FU", "Reduced motivation", "Norm", 8.0, NA,
  "24-month FU", "Reduced motivation", "Norm", 8.0, NA
)


df$Time <- factor(df$Time, levels = c("Baseline", "12-month FU", "24-month FU"))

library(ggplot2)

plot_mfi <- function(data, scale_name) {
  ggplot(data %>% filter(Scale == scale_name),
         aes(x = Time, y = Mean, group = Group, color = Group, linetype = Group)) +
    geom_line(size = 1.2) +
    geom_point(size = 2.5) +
    geom_errorbar(
      data = ~ .x %>% filter(Group == "PROCORE"),
      aes(ymin = Mean - SD, ymax = Mean + SD),
      width = 0.1
    ) +
    scale_color_manual(values = c("PROCORE" = "#c1327b", "Norm" = "gray40")) +
    scale_linetype_manual(values = c("PROCORE" = "solid", "Norm" = "dashed")) +
    scale_x_discrete(labels = c(
      "Baseline\nn = 411",
      "12-month FU\nn = 304",
      "24-month FU\nn = 252"
    )) +
    labs(y = "Score", x = "Time") +
    coord_cartesian(ylim = c(4, 14), clip = "off") +
    theme_minimal(base_size = 12, base_family = "Arial") +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(t = 5, r = 5, b = 40, l = 5)
    )
}


plot_mfi(df, "General fatigue")
plot_mfi(df, "Mental fatigue")
plot_mfi(df, "Physical fatigue")
plot_mfi(df, "Reduced activity")
plot_mfi(df, "Reduced motivation")

library(patchwork)

p1 <- plot_mfi(df, "General fatigue")
p2 <- plot_mfi(df, "Mental fatigue")
p3 <- plot_mfi(df, "Physical fatigue")
p4 <- plot_mfi(df, "Reduced activity")
p5 <- plot_mfi(df, "Reduced motivation")

combined <- (p1 | p2 | p3) / (p4 | p5 | plot_spacer()) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = list(c("(A)", "(B)", "(C)", "(D)", "(E)")),
    tag_prefix = "",
    tag_suffix = ""
  ) &
  theme(legend.position = "bottom",
        plot.tag = element_text(size = 12, face = "bold"))

print(combined)

showtext_opts(dpi = 1000)
ggsave("mfi_plots_figure2_panes_ARIAL.png", combined, width = 12, height = 8, dpi = 1000)

