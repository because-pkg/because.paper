# Supplementary Publication Figures for Misspecification Validation

library(ggplot2)
library(dplyr)
library(patchwork)

core_scenarios <- c(
  "S0_True",
  "S1_WrongDAG_MissingEdges",
  "S3a_OmittedConfounder_BodyMass",
  "S4_WrongHierarchy_DropSurvey",
  "S5a_WrongStruct_None"
)

joint_scenarios <- c(
  "S0_True",
  "S4_WrongHierarchy_DropSurvey",
  "S5a_WrongStruct_None"
)

tests_df <- read.csv(file.path(in_dir, "misspecification_numpyro_tests.csv"), stringsAsFactors = FALSE) %>%
  filter(Scenario %in% core_scenarios)
joint_df <- read.csv(file.path(in_dir, "misspecification_numpyro_joint.csv"), stringsAsFactors = FALSE) %>%
  filter(Scenario %in% joint_scenarios)

# Figure 1: Complete Standardized d-Separation & ROPE Tests (Baseline + 1..4)

base_scen_titles <- c(
  "S0_True"                        = "Baseline: True Model",
  "S1_WrongDAG_MissingEdges"       = "(1) Wrong/Missing Edges",
  "S3a_OmittedConfounder_BodyMass" = "(2) Drop Body_Mass",
  "S4_WrongHierarchy_DropSurvey"   = "(3) Drop Survey Scale",
  "S5a_WrongStruct_None"           = "(4) structure = NULL"
)

scenario_labels <- base_scen_titles[core_scenarios]

joint_scen_labels <- c(
  "S0_True"                      = "Baseline",
  "S4_WrongHierarchy_DropSurvey" = "(3) Drop Survey",
  "S5a_WrongStruct_None"         = "(4) No Structure"
)

format_dsep_tests <- function(df) {
  df %>%
    mutate(
      Orig_Order = seq_len(n()),
      Scale_Tag = ifelse(
        grepl("\\[Scale:\\s*[a-z]+\\]", Test, perl = TRUE),
        sub(".*\\[Scale:\\s*([a-z]+)\\].*", "\\1", Test, perl = TRUE),
        "other"
      ),
      Scale_Order = factor(Scale_Tag, levels = c("site", "survey", "species", "obs", "other")),
      Clean_Test = gsub("\\s*\\[Scale:.*?\\]\\s*$", "", Test),
      Clean_Test = gsub("\\s*\\|\\s*\\{\\s*\\}\\s*$", "", Clean_Test),
      Clean_Test = gsub("\\s+", " ", trimws(Clean_Test)),
      Clean_Test = gsub("I\\(Temperature - Thermal_Tol\\)", "Mismatch", Clean_Test),
      ROPE_low = -0.10,
      ROPE_high = 0.10,
      Status = case_when(
        !Rejected_CI ~ "Supported (95% CI includes 0)",
        Rejected_CI & P_ROPE_num < 5 ~ "Rejected (CI != 0 & P_ROPE < 5%)",
        TRUE ~ "95% CI != 0, >= 5% in ROPE"
      ),
      Status = factor(
        Status,
        levels = c(
          "Supported (95% CI includes 0)",
          "95% CI != 0, >= 5% in ROPE",
          "Rejected (CI != 0 & P_ROPE < 5%)"
        )
      ),
      Scenario_Label = recode(Scenario, !!!scenario_labels),
      Scenario_Label = factor(Scenario_Label, levels = unname(scenario_labels)),
      ROPE_Tag = sprintf("P_ROPE = %s", P_ROPE),
      Row_ID = paste(Scenario, Parameter, sep = "::")
    ) %>%
    arrange(Scenario_Label, Scale_Order, Orig_Order) %>%
    mutate(Row_ID = factor(Row_ID, levels = rev(Row_ID)))
}

status_cols <- c(
  "Supported (95% CI includes 0)"     = "#1b7837",
  "95% CI != 0, >= 5% in ROPE"        = "#d95f02",
  "Rejected (CI != 0 & P_ROPE < 5%)"  = "#b2182b"
)

build_dsep_panel <- function(df_sub,
                             title,
                             subtitle,
                             x_limits = c(-0.62, 1.44),
                             x_breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0),
                             legend_nrow = 2,
                             axis_y_size = 7.2,
                             strip_size = 7.4,
                             tag_size = 2.5) {
  df_sub$Scenario_Label <- droplevels(df_sub$Scenario_Label)
  df_sub$Row_ID <- droplevels(df_sub$Row_ID)
  label_map <- setNames(df_sub$Clean_Test, as.character(df_sub$Row_ID))

  ggplot(df_sub, aes(x = Estimate, y = Row_ID, color = Status)) +
    geom_rect(
      aes(xmin = ROPE_low, xmax = ROPE_high, ymin = -Inf, ymax = Inf),
      fill = "grey88", alpha = 0.6, inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(xmin = LowerCI, xmax = UpperCI),
                  orientation = "y", width = 0.28, linewidth = 0.6, show.legend = TRUE) +
    geom_point(size = 2.1, show.legend = TRUE) +
    geom_text(
      aes(x = pmax(UpperCI, ROPE_high) + 0.03, label = ROPE_Tag),
      hjust = 0, size = tag_size, show.legend = FALSE
    ) +
    facet_grid(Scenario_Label ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_discrete(labels = label_map) +
    scale_x_continuous(breaks = x_breaks) +
    scale_color_manual(values = status_cols, limits = names(status_cols), drop = FALSE) +
    coord_cartesian(xlim = x_limits) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Standardized Posterior Slope (95% CI)",
      y = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 9.5) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 10.6),
      plot.subtitle = element_text(size = 8.1, color = "grey30"),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold", size = strip_size),
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      axis.text.y = element_text(size = axis_y_size),
      legend.position = "bottom",
      legend.text = element_text(size = 7.8),
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(0.24, "lines")
    ) +
    guides(color = guide_legend(nrow = legend_nrow, byrow = TRUE))
}

tests_formatted <- format_dsep_tests(tests_df)

left_scenarios <- c(
  "S0_True",
  "S1_WrongDAG_MissingEdges",
  "S3a_OmittedConfounder_BodyMass"
)
right_scenarios <- c(
  "S4_WrongHierarchy_DropSurvey",
  "S5a_WrongStruct_None"
)

left_df  <- filter(tests_formatted, Scenario %in% left_scenarios)
right_df <- filter(tests_formatted, Scenario %in% right_scenarios)

p_all_left <- build_dsep_panel(
  left_df,
  title = sprintf("Baseline & DAG/Confounder Misspecifications (1)-(2) (%d tests)", nrow(left_df)),
  subtitle = NULL,
  x_limits = c(-0.64, 1.44),
  x_breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0),
  legend_nrow = 1,
  axis_y_size = 7.0,
  strip_size = 7.2,
  tag_size = 2.4
)

p_all_right <- build_dsep_panel(
  right_df,
  title = sprintf("Hierarchy & Covariance Misspecifications (3)-(4) (%d tests)", nrow(right_df)),
  subtitle = NULL,
  x_limits = c(-0.64, 1.44),
  x_breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0),
  legend_nrow = 1,
  axis_y_size = 7.0,
  strip_size = 7.2,
  tag_size = 2.4
)

fig_all <- (p_all_left | p_all_right) +
  plot_layout(widths = c(1.0, 1.0), guides = "collect") +
  plot_annotation(
    theme = theme(
      legend.position = "bottom",
      plot.margin = margin(8, 10, 8, 8)
    )
  ) &
  theme(legend.position = "bottom")

all_png <- file.path(in_dir, "FigS_all_dsep_tests.png")
all_pdf <- file.path(in_dir, "FigS_all_dsep_tests.pdf")

ggsave(all_png, fig_all, width = 18.0, height = 11.4, dpi = 300, bg = "white")
ggsave(all_pdf, fig_all, width = 18.0, height = 11.4, bg = "white")


# Figure 2: Joint Poisson SEM Diagnostics for Baseline, (3), and (4)

pb_df <- joint_df %>%
  filter(!is.na(True_Val)) %>%
  mutate(
    Path_Label = recode(
      Parameter,
      "beta_Metabolic_Rate_Body_Mass"                 = "Metabolic_Rate ~ Body_Mass\n(True = 0.80)",
      "beta_Thermal_Tol_Body_Mass"                    = "Thermal_Tol ~ Body_Mass\n(True = 0.60)",
      "beta_Flower_Cover_NDVI"                        = "Flower_Cover ~ NDVI\n(True = 0.60)",
      "beta_Flower_Cover_Elevation"                   = "Flower_Cover ~ Elevation\n(True = 0.30)",
      "beta_Temperature_Elevation"                    = "Temperature ~ Elevation\n(True = -1.50)",
      "beta_Abundance_Flower_Cover"                   = "Abundance ~ Flower_Cover\n(True = 0.20)",
      "beta_Abundance_Wind_Speed"                     = "Abundance ~ Wind_Speed\n(True = -0.20)",
      "beta_Abundance_Temperature_minus_Thermal_Tol"  = "Abundance ~ Mismatch\n(True = -0.30)"
    ),
    Path_Label = factor(
      Path_Label,
      levels = c(
        "Metabolic_Rate ~ Body_Mass\n(True = 0.80)",
        "Thermal_Tol ~ Body_Mass\n(True = 0.60)",
        "Flower_Cover ~ NDVI\n(True = 0.60)",
        "Flower_Cover ~ Elevation\n(True = 0.30)",
        "Temperature ~ Elevation\n(True = -1.50)",
        "Abundance ~ Flower_Cover\n(True = 0.20)",
        "Abundance ~ Wind_Speed\n(True = -0.20)",
        "Abundance ~ Mismatch\n(True = -0.30)"
      )
    ),
    Scenario_Short = factor(
      recode(Scenario, !!!joint_scen_labels),
      levels = rev(unname(joint_scen_labels))
    ),
    Coverage = ifelse(Covered_95, "95% CI Covers Truth", "95% CI Misses Truth")
  )

pb_truth <- pb_df %>%
  distinct(Path_Label, True_Val)

scen_cols <- c(
  "Baseline"        = "#1f78b4",
  "(3) Drop Survey" = "#b15928",
  "(4) No Structure" = "#e7298a"
)

p_joint_a <- ggplot(pb_df, aes(x = Mean, y = Scenario_Short, color = Scenario_Short, shape = Coverage)) +
  geom_vline(
    data = pb_truth,
    aes(xintercept = True_Val),
    linetype = "dashed", color = "grey30", linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    aes(xmin = LowerCI, xmax = UpperCI),
    orientation = "y",
    width = 0.0, linewidth = 0.65
  ) +
  geom_point(size = 2.2, stroke = 1.0) +
  facet_wrap(~ Path_Label, ncol = 4, scales = "free_x") +
  scale_x_continuous(n.breaks = 4, expand = expansion(mult = 0.14)) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  scale_shape_manual(
    values = c("95% CI Covers Truth" = 16, "95% CI Misses Truth" = 1),
    guide = guide_legend(order = 1, nrow = 1)
  ) +
  scale_color_manual(
    values = scen_cols, drop = FALSE,
    guide = guide_legend(order = 2, reverse = TRUE, nrow = 1)
  ) +
  labs(
    title = "A. Joint SEM Path Coefficient Recovery (Mean & 95% CI)",
    subtitle = "Dashed lines = true generating values; hollow circles miss truth in (3) & (4)",
    x = "Estimated Path Coefficient (95% CI)",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 8.4, color = "grey30"),
    strip.text = element_text(face = "bold", size = 7.4, lineheight = 1.0, margin = margin(t = 3, b = 3)),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    axis.text.x = element_text(size = 7.4),
    axis.text.y = element_text(size = 7.8, color = "grey20"),
    panel.spacing.x = unit(0.75, "lines"),
    panel.spacing.y = unit(0.55, "lines"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 8.0),
    panel.grid.minor = element_blank()
  )

lambda_params <- c(
  "lambda_Body_Mass_phylo",
  "lambda_Metabolic_Rate_phylo",
  "lambda_Thermal_Tol_phylo",
  "lambda_Temperature_spatial",
  "lambda_Flower_Cover_spatial"
)

struct_sd_params <- c(
  "sigma_Body_Mass_phylo",
  "sigma_Metabolic_Rate_phylo",
  "sigma_Thermal_Tol_phylo",
  "sigma_Abundance_phylo",
  "sigma_Temperature_spatial",
  "sigma_Flower_Cover_spatial",
  "sigma_Abundance_spatial"
)

resid_sd_params <- c(
  "sigma_Metabolic_Rate",
  "sigma_Thermal_Tol",
  "sigma_Body_Mass",
  "sigma_Temperature",
  "sigma_Flower_Cover"
)

pb_lambda_df <- joint_df %>%
  filter(Parameter %in% lambda_params) %>%
  mutate(
    Param_Label = recode(
      Parameter,
      "lambda_Metabolic_Rate_phylo" = "Metabolic Rate (phylo)",
      "lambda_Thermal_Tol_phylo"    = "Thermal Tol (phylo)",
      "lambda_Body_Mass_phylo"      = "Body Mass (phylo)",
      "lambda_Temperature_spatial"  = "Temperature (spatial)",
      "lambda_Flower_Cover_spatial" = "Flower Cover (spatial)"
    ),
    Param_Label = factor(
      Param_Label,
      levels = rev(c(
        "Metabolic Rate (phylo)",
        "Thermal Tol (phylo)",
        "Body Mass (phylo)",
        "Temperature (spatial)",
        "Flower Cover (spatial)"
      ))
    ),
    Scenario_Short = factor(
      recode(Scenario, !!!joint_scen_labels),
      levels = rev(unname(joint_scen_labels))
    )
  )

pb_struct_df <- joint_df %>%
  filter(Parameter %in% struct_sd_params) %>%
  mutate(
    Param_Label = recode(
      Parameter,
      "sigma_Metabolic_Rate_phylo"  = "Metabolic Rate (phylo SD)",
      "sigma_Thermal_Tol_phylo"     = "Thermal Tol (phylo SD)",
      "sigma_Body_Mass_phylo"       = "Body Mass (phylo SD)",
      "sigma_Abundance_phylo"       = "Abundance (phylo SD)",
      "sigma_Temperature_spatial"   = "Temperature (spatial SD)",
      "sigma_Flower_Cover_spatial"  = "Flower Cover (spatial SD)",
      "sigma_Abundance_spatial"     = "Abundance (spatial SD)"
    ),
    Param_Label = factor(
      Param_Label,
      levels = rev(c(
        "Metabolic Rate (phylo SD)",
        "Thermal Tol (phylo SD)",
        "Body Mass (phylo SD)",
        "Abundance (phylo SD)",
        "Temperature (spatial SD)",
        "Flower Cover (spatial SD)",
        "Abundance (spatial SD)"
      ))
    ),
    Scenario_Short = factor(
      recode(Scenario, !!!joint_scen_labels),
      levels = rev(unname(joint_scen_labels))
    )
  )

pb_resid_df <- joint_df %>%
  filter(Parameter %in% resid_sd_params) %>%
  mutate(
    Param_Label = recode(
      Parameter,
      "sigma_Metabolic_Rate" = "Metabolic_Rate\n(species resid SD)",
      "sigma_Thermal_Tol"    = "Thermal_Tol\n(species resid SD)",
      "sigma_Body_Mass"      = "Body_Mass\n(species resid SD)",
      "sigma_Temperature"    = "Temperature\n(survey resid SD)",
      "sigma_Flower_Cover"   = "Flower_Cover\n(survey resid SD)"
    ),
    Param_Label = factor(
      Param_Label,
      levels = c(
        "Metabolic_Rate\n(species resid SD)",
        "Thermal_Tol\n(species resid SD)",
        "Body_Mass\n(species resid SD)",
        "Temperature\n(survey resid SD)",
        "Flower_Cover\n(survey resid SD)"
      )
    ),
    Scenario_Short = factor(
      recode(Scenario, !!!joint_scen_labels),
      levels = rev(unname(joint_scen_labels))
    )
  )

p_joint_b1 <- ggplot(pb_lambda_df, aes(x = Mean, y = Param_Label, color = Scenario_Short)) +
  geom_errorbar(
    aes(xmin = LowerCI, xmax = UpperCI),
    orientation = "y",
    position = position_dodge(width = 0.60),
    width = 0.0, linewidth = 0.65
  ) +
  geom_point(position = position_dodge(width = 0.60), size = 2.1) +
  scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.25, 0.50, 0.75, 1.0)) +
  scale_color_manual(values = scen_cols, drop = FALSE, guide = "none") +
  labs(
    title = "B. Spatial & Phylogenetic Signal and Variance-Component Distortion",
    subtitle = "B1. Signal strength (lambda; absent in (4))",
    x = "Posterior Signal lambda (95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 8.4, color = "grey25", face = "bold"),
    axis.text.x = element_text(size = 7.6),
    axis.text.y = element_text(size = 8.0),
    panel.grid.minor = element_blank()
  )

p_joint_b2 <- ggplot(pb_struct_df, aes(x = Mean, y = Param_Label, color = Scenario_Short)) +
  geom_errorbar(
    aes(xmin = LowerCI, xmax = UpperCI),
    orientation = "y",
    position = position_dodge(width = 0.62),
    width = 0.0, linewidth = 0.65
  ) +
  geom_point(position = position_dodge(width = 0.62), size = 2.1) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5)) +
  scale_color_manual(values = scen_cols, drop = FALSE, guide = "none") +
  labs(
    title = " ",
    subtitle = "B2. Structured random-effect SDs (absent in (4))",
    x = "Posterior Structured SD (95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 8.4, color = "grey25", face = "bold"),
    axis.text.x = element_text(size = 7.6),
    axis.text.y = element_text(size = 8.0),
    panel.grid.minor = element_blank()
  )

p_joint_b3 <- ggplot(pb_resid_df, aes(x = Mean, y = Scenario_Short, color = Scenario_Short)) +
  geom_errorbar(
    aes(xmin = LowerCI, xmax = UpperCI),
    orientation = "y",
    width = 0.0, linewidth = 0.65
  ) +
  geom_point(size = 2.1) +
  facet_wrap(~ Param_Label, ncol = 5, scales = "free_x") +
  scale_x_continuous(n.breaks = 4, expand = expansion(mult = 0.14)) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  scale_color_manual(values = scen_cols, drop = FALSE, guide = "none") +
  labs(
    subtitle = "B3. Residual SD distortion ((3) collapses survey CIs; (4) pools unmodeled spatial & phylo variance)",
    x = "Posterior Residual SD (95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    plot.title.position = "plot",
    plot.subtitle = element_text(size = 8.4, color = "grey25", face = "bold"),
    strip.text = element_text(face = "bold", size = 7.4, lineheight = 1.0, margin = margin(t = 3, b = 3)),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    axis.text.x = element_text(size = 7.2),
    axis.text.y = element_text(size = 7.8, color = "grey20"),
    panel.spacing.x = unit(0.9, "lines"),
    panel.grid.minor = element_blank()
  )

fig_joint <- (
  p_joint_a /
  ((p_joint_b1 | p_joint_b2) + plot_layout(widths = c(1.0, 1.08))) /
  p_joint_b3
) +
  plot_layout(heights = c(1.18, 0.88, 0.56), guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(8, 10, 8, 8))) &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "lines"),
    legend.text = element_text(size = 8.2)
  )

joint_png <- file.path(in_dir, "FigS_joint_diagnostics.png")
joint_pdf <- file.path(in_dir, "FigS_joint_diagnostics.pdf")

ggsave(joint_png, fig_joint, width = 12.2, height = 10.8, dpi = 300, bg = "white")
ggsave(joint_pdf, fig_joint, width = 12.2, height = 10.8, bg = "white")

cat("Saved:\n  ", all_png, "\n  ", all_pdf, "\n  ", joint_png, "\n  ", joint_pdf, "\n")
