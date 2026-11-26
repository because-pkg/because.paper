
# Misspecification Validation for Case Study 1 (Baseline + 4 Core Scenarios)
=========

library(because)
library(because.phybase)
library(ape)
library(MASS)

set.seed(49)

# Simulated data (from multiscale_simulation.R)

N_Species          <- 50
N_Site             <- 50
N_Surveys_per_Site <- 3
N_Surveys          <- N_Site * N_Surveys_per_Site

BETA_BM_MR          <-  0.80
BETA_MR_TT          <-  0.00
BETA_BM_TT          <-  0.60
BETA_ELEV_NDVI      <-  0.00
BETA_NDVI_FC        <-  0.60
BETA_ELEV_FC        <-  0.30
BETA_ELEV_TEMP      <- -1.50
BASELINE_INTERCEPT  <-  0.80
BETA_FC_ABUND       <-  0.20
BETA_WIND_ABUND     <- -0.20
BETA_MISMATCH_ABUND <- -0.30

# Phylogeny
species_tree <- rtree(N_Species)
species_tree$tip.label  <- paste0("Sp_", 1:N_Species)
species_tree$edge.length <- species_tree$edge.length / max(node.depth.edgelength(species_tree))
Vcv <- vcv.phylo(species_tree)

lambda <- 0.6
V_phylo <- lambda * Vcv + (1 - lambda) * diag(diag(Vcv))

# Species level
Body_Mass      <- rnorm(N_Species, 0, 1.0)
eps_MR         <- as.vector(mvrnorm(1, rep(0, N_Species), 0.3^2 * V_phylo))
Metabolic_Rate <- BETA_BM_MR * Body_Mass + eps_MR
eps_TT         <- as.vector(mvrnorm(1, rep(0, N_Species), 0.3^2 * V_phylo))
Thermal_Tol    <- BETA_BM_TT * Body_Mass + eps_TT

d_species <- data.frame(
  Species        = species_tree$tip.label,
  Body_Mass      = Body_Mass,
  Metabolic_Rate = Metabolic_Rate,
  Thermal_Tol    = Thermal_Tol
)

# Site level
site_coords <- cbind(runif(N_Site, 0, 10), runif(N_Site, 0, 10))
rownames(site_coords) <- paste0("Site_", 1:N_Site)
dist_mat  <- as.matrix(dist(site_coords))
V_spatial <- exp(-dist_mat / 5.0)

Elevation <- rnorm(N_Site, 0, 1.0)
NDVI      <- BETA_ELEV_NDVI * Elevation + rnorm(N_Site, 0, 1.0)

d_site <- data.frame(
  Site      = paste0("Site_", 1:N_Site),
  Elevation = Elevation,
  NDVI      = NDVI
)

# Survey level
site_of_survey <- rep(1:N_Site, each = N_Surveys_per_Site)
u_site_temp    <- as.vector(mvrnorm(1, rep(0, N_Site), 0.4^2 * V_spatial))
Temperature    <- BETA_ELEV_TEMP * Elevation[site_of_survey] + u_site_temp[site_of_survey] + rnorm(N_Surveys, 0, 1.0)
Wind_Speed     <- rnorm(N_Surveys, 0, 1.0)

u_site_FC    <- as.vector(mvrnorm(1, rep(0, N_Site), 0.4^2 * V_spatial))
Flower_Cover <- BETA_NDVI_FC * NDVI[site_of_survey] + BETA_ELEV_FC * Elevation[site_of_survey] + u_site_FC[site_of_survey] + rnorm(N_Surveys, 0, 1.0)

d_survey <- data.frame(
  Survey       = paste0("Survey_", 1:N_Surveys),
  Site         = paste0("Site_", site_of_survey),
  Temperature  = Temperature,
  Wind_Speed   = Wind_Speed,
  Flower_Cover = Flower_Cover
)

# Observation level (Poisson counts)
d_obs <- expand.grid(Survey = d_survey$Survey, Species = d_species$Species, stringsAsFactors = FALSE)
d_obs <- merge(d_obs, d_survey[, c("Survey", "Site", "Temperature", "Wind_Speed", "Flower_Cover")], by = "Survey")
d_obs <- merge(d_obs, d_species, by = "Species")
d_obs$obs_id <- seq_len(nrow(d_obs))

Thermal_Mismatch <- d_obs$Temperature - d_obs$Thermal_Tol

u_phylo_abund   <- as.vector(mvrnorm(1, rep(0, N_Species), 0.3^2 * V_phylo))
u_spatial_abund <- as.vector(mvrnorm(1, rep(0, N_Site), 0.3^2 * V_spatial))
names(u_phylo_abund)   <- paste0("Sp_", 1:N_Species)
names(u_spatial_abund) <- paste0("Site_", 1:N_Site)

log_abund <- BASELINE_INTERCEPT +
  BETA_FC_ABUND       * d_obs$Flower_Cover +
  BETA_WIND_ABUND     * d_obs$Wind_Speed   +
  BETA_MISMATCH_ABUND * Thermal_Mismatch   +
  u_phylo_abund[d_obs$Species]             +
  u_spatial_abund[d_obs$Site]

d_obs$Abundance <- rpois(nrow(d_obs), lambda = exp(log_abund))

# Base NumPyro inputs
data_true <- list(
  species = d_species,
  site    = d_site,
  survey  = d_survey[, c("Survey", "Site", "Temperature", "Wind_Speed", "Flower_Cover")],
  obs     = d_obs[, c("obs_id", "Survey", "Species", "Abundance")]
)

levels_true <- list(
  species = c("Body_Mass", "Metabolic_Rate", "Thermal_Tol", "Species"),
  site    = c("Elevation", "NDVI", "Site"),
  survey  = c("Temperature", "Wind_Speed", "Flower_Cover", "Survey"),
  obs     = c("Abundance")
)

eqs_true <- list(
  Flower_Cover   ~ NDVI + Elevation,
  Body_Mass      ~ 1,
  Metabolic_Rate ~ Body_Mass,
  Thermal_Tol    ~ Body_Mass,
  Temperature    ~ Elevation,
  Abundance      ~ Flower_Cover + Wind_Speed + I(Temperature - Thermal_Tol)
)

structures_true <- list(
  phylo   = species_tree,
  spatial = V_spatial
)

link_vars_true <- list(site = "Site", survey = "Survey", species = "Species")

# Output helpers (used after each explicit because() fit below)
true_betas <- c(
  beta_Metabolic_Rate_Body_Mass                = BETA_BM_MR,
  beta_Thermal_Tol_Body_Mass                   = BETA_BM_TT,
  beta_Flower_Cover_NDVI                       = BETA_NDVI_FC,
  beta_Flower_Cover_Elevation                  = BETA_ELEV_FC,
  beta_Temperature_Elevation                   = BETA_ELEV_TEMP,
  beta_Abundance_Flower_Cover                  = BETA_FC_ABUND,
  beta_Abundance_Wind_Speed                    = BETA_WIND_ABUND,
  beta_Abundance_Temperature_minus_Thermal_Tol = BETA_MISMATCH_ABUND
)

format_dsep_df <- function(sum_dsep, sc_name, sc_axis) {
  res_df <- sum_dsep$results
  res_df$P_ROPE_num        <- as.numeric(sub("%$", "", res_df$P_ROPE))
  res_df$Rejected_CI       <- (res_df$LowerCI > 0) | (res_df$UpperCI < 0)
  res_df$Rejected_ROPE_lt5 <- res_df$P_ROPE_num < 5
  res_df$Scenario          <- sc_name
  res_df$Axis              <- sc_axis
  res_df
}

format_joint_df <- function(sum_joint, sc_name, sc_axis) {
  j_res <- as.data.frame(sum_joint$results, stringsAsFactors = FALSE)
  j_df <- data.frame(
    Scenario  = sc_name,
    Axis      = sc_axis,
    Parameter = rownames(j_res),
    Mean      = round(as.numeric(j_res[["Mean"]]), 4),
    SD        = round(as.numeric(j_res[["SD"]]), 4),
    LowerCI   = round(as.numeric(j_res[["2.5%"]]), 4),
    UpperCI   = round(as.numeric(j_res[["97.5%"]]), 4),
    CI_Width  = round(as.numeric(j_res[["97.5%"]]) - as.numeric(j_res[["2.5%"]]), 4),
    Rhat      = if ("Rhat" %in% names(j_res)) round(as.numeric(j_res[["Rhat"]]), 3) else NA_real_,
    n.eff     = if ("n.eff" %in% names(j_res)) round(as.numeric(j_res[["n.eff"]]), 0) else NA_real_,
    True_Val  = unname(true_betas[rownames(j_res)]),
    stringsAsFactors = FALSE
  )
  j_df$Bias       <- round(j_df$Mean - j_df$True_Val, 4)
  j_df$Covered_95 <- !is.na(j_df$True_Val) & (j_df$LowerCI <= j_df$True_Val) & (j_df$True_Val <= j_df$UpperCI)
  j_df
}

format_summary_row <- function(res_df, sc_name, sc_axis, dt_dsep, dt_joint = NA_real_) {
  rej_labels <- if (any(res_df$Rejected_CI, na.rm = TRUE)) {
    paste(res_df$Test[res_df$Rejected_CI], collapse = "; ")
  } else {
    "NONE"
  }
  data.frame(
    Scenario         = sc_name,
    Axis             = sc_axis,
    N_Tests_Eval     = nrow(res_df),
    N_Rejected_CI    = sum(res_df$Rejected_CI, na.rm = TRUE),
    N_Rejected_ROPE5 = sum(res_df$Rejected_ROPE_lt5, na.rm = TRUE),
    DAG_Rejected     = any(res_df$Rejected_CI, na.rm = TRUE),
    Max_Rhat         = round(max(res_df$Rhat, na.rm = TRUE), 3),
    Time_dsep_sec    = round(dt_dsep, 1),
    Time_joint_sec   = round(dt_joint, 1),
    Rejected_Claims  = rej_labels,
    stringsAsFactors = FALSE
  )
}

# Baseline original model

t0 <- Sys.time()
fit_dsep_baseline <- because(
  equations  = eqs_true,
  data       = data_true,
  levels     = levels_true,
  hierarchy  = "site > survey > obs; species > obs",
  link_vars  = link_vars_true,
  id_col     = "Species",
  family     = c(Abundance = "poisson"),
  structure  = structures_true,
  dsep       = TRUE,
  parallel   = TRUE,
  n.cores    = 3,
  n.iter     = 2000,
  engine     = "numpyro",
  quiet      = FALSE
)
dt_dsep_baseline <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

fit_dsep_baseline_std <- standardize(fit_dsep_baseline)
sum_dsep_baseline     <- summary(fit_dsep_baseline_std, rope = TRUE)
tests_baseline        <- format_dsep_df(sum_dsep_baseline, "S0_True", "Baseline (Correct Model)")
print(tests_baseline)

# joint final SEM
t1 <- Sys.time()
fit_joint_baseline <- because(
  equations  = eqs_true,
  data       = data_true,
  levels     = levels_true,
  hierarchy  = "site > survey > obs; species > obs",
  link_vars  = link_vars_true,
  id_col     = "Species",
  family     = c(Abundance = "poisson"),
  structure  = structures_true,
  dsep       = FALSE,
  WAIC       = FALSE,
  parallel   = TRUE,
  n.cores    = 3,
  n.iter     = 2000,
  engine     = "numpyro",
  quiet      = TRUE
)
dt_joint_baseline <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

sum_joint_baseline <- summary(fit_joint_baseline)
joint_baseline     <- format_joint_df(sum_joint_baseline, "S0_True", "Baseline (Correct Model)")


#  SCENARIO (1): INCORRECT DAG - WRONG & MISSING EDGES

eqs_wrong_edges <- list(
  Flower_Cover   ~ NDVI,
  Body_Mass      ~ 1,
  Metabolic_Rate ~ Body_Mass,
  Thermal_Tol    ~ Body_Mass,
  Temperature    ~ Elevation,
  Abundance      ~ NDVI + Wind_Speed + I(Temperature - Thermal_Tol)
)

t0 <- Sys.time()
fit_dsep_s1 <- because(
  equations  = eqs_wrong_edges,
  data       = data_true,
  levels     = levels_true,
  hierarchy  = "site > survey > obs; species > obs",
  link_vars  = link_vars_true,
  id_col     = "Species",
  family     = c(Abundance = "poisson"),
  structure  = structures_true,
  dsep       = TRUE,
  parallel   = TRUE,
  n.cores    = 3,
  n.iter     = 2000,
  engine     = "numpyro",
  quiet      = FALSE
)
dt_dsep_s1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

fit_dsep_s1_std <- standardize(fit_dsep_s1)
sum_dsep_s1     <- summary(fit_dsep_s1_std, rope = TRUE)
tests_s1        <- format_dsep_df(sum_dsep_s1, "S1_WrongDAG_MissingEdges", "1. Incorrect DAG (Wrong/Missing Edges)")
print(tests_s1)

# SCENARIO (2): OMITTED CONFOUNDER (Drop Body_Mass)

data_no_bm         <- data_true
data_no_bm$species <- d_species[, c("Species", "Metabolic_Rate", "Thermal_Tol")]

levels_no_bm         <- levels_true
levels_no_bm$species <- c("Metabolic_Rate", "Thermal_Tol", "Species")

eqs_no_bm <- list(
  Flower_Cover   ~ NDVI + Elevation,
  Metabolic_Rate ~ 1,
  Thermal_Tol    ~ 1,
  Temperature    ~ Elevation,
  Abundance      ~ Flower_Cover + Wind_Speed + I(Temperature - Thermal_Tol)
)

t0 <- Sys.time()
fit_dsep_s2 <- because(
  equations  = eqs_no_bm,
  data       = data_no_bm,
  levels     = levels_no_bm,
  hierarchy  = "site > survey > obs; species > obs",
  link_vars  = link_vars_true,
  id_col     = "Species",
  family     = c(Abundance = "poisson"),
  structure  = structures_true,
  dsep       = TRUE,
  parallel   = TRUE,
  n.cores    = 3,
  n.iter     = 2000,
  engine     = "numpyro",
  quiet      = FALSE
)
dt_dsep_s2 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

fit_dsep_s2_std <- standardize(fit_dsep_s2)
sum_dsep_s2     <- summary(fit_dsep_s2_std, rope = TRUE)
tests_s2        <- format_dsep_df(sum_dsep_s2, "S3a_OmittedConfounder_BodyMass", "2. Omitted Confounder (Body_Mass Dropped)")
print(tests_s2)

# SCENARIO (3): MISSPECIFIED HIERARCHY (Drop Survey Level)

data_no_survey <- list(
  species = d_species,
  site    = d_site,
  obs     = d_obs[, c("obs_id", "Site", "Species", "Temperature", "Wind_Speed", "Flower_Cover", "Abundance")]
)

levels_no_survey <- list(
  species = c("Body_Mass", "Metabolic_Rate", "Thermal_Tol", "Species"),
  site    = c("Elevation", "NDVI", "Site"),
  obs     = c("Temperature", "Wind_Speed", "Flower_Cover", "Abundance")
)

t0 <- Sys.time()
fit_dsep_s3 <- because(
  equations  = eqs_true,
  data       = data_no_survey,
  levels     = levels_no_survey,
  hierarchy  = "site > obs; species > obs",
  link_vars  = list(site = "Site", species = "Species"),
  id_col     = "Species",
  family     = c(Abundance = "poisson"),
  structure  = structures_true,
  dsep       = TRUE,
  parallel   = TRUE,
  n.cores    = 3,
  n.iter     = 2000,
  engine     = "numpyro",
  quiet      = FALSE
)
dt_dsep_s3 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

fit_dsep_s3_std <- standardize(fit_dsep_s3)
sum_dsep_s3     <- summary(fit_dsep_s3_std, rope = TRUE)
tests_s3        <- format_dsep_df(sum_dsep_s3, "S4_WrongHierarchy_DropSurvey", "3. Misspecified Hierarchy (Drop Survey Level)")
print(tests_s3)

t1 <- Sys.time()
fit_joint_s3 <- because(
  equations  = eqs_true,
  data       = data_no_survey,
  levels     = levels_no_survey,
  hierarchy  = "site > obs; species > obs",
  link_vars  = list(site = "Site", species = "Species"),
  id_col     = "Species",
  family     = c(Abundance = "poisson"),
  structure  = structures_true,
  dsep       = FALSE,
  WAIC       = FALSE,
  parallel   = TRUE,
  n.cores    = 3,
  n.iter     = 2000,
  engine     = "numpyro",
  quiet      = TRUE
)
dt_joint_s3 <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

sum_joint_s3 <- summary(fit_joint_s3)
joint_s3     <- format_joint_df(sum_joint_s3, "S4_WrongHierarchy_DropSurvey", "3. Misspecified Hierarchy (Drop Survey Level)")

#  SCENARIO (4): MISSPECIFIED COVARIANCE (structure = NULL)


t0 <- Sys.time()
fit_dsep_s4 <- because(
  equations  = eqs_true,
  data       = data_true,
  levels     = levels_true,
  hierarchy  = "site > survey > obs; species > obs",
  link_vars  = link_vars_true,
  id_col     = "Species",
  family     = c(Abundance = "poisson"),
  structure  = NULL,
  dsep       = TRUE,
  parallel   = TRUE,
  n.cores    = 3,
  n.iter     = 2000,
  engine     = "numpyro",
  quiet      = FALSE
)
dt_dsep_s4 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

fit_dsep_s4_std <- standardize(fit_dsep_s4)
sum_dsep_s4     <- summary(fit_dsep_s4_std, rope = TRUE)
tests_s4        <- format_dsep_df(sum_dsep_s4, "S5a_WrongStruct_None", "4. Misspecified Covariance (structure = NULL)")
print(tests_s4)

t1 <- Sys.time()
fit_joint_s4 <- because(
  equations  = eqs_true,
  data       = data_true,
  levels     = levels_true,
  hierarchy  = "site > survey > obs; species > obs",
  link_vars  = link_vars_true,
  id_col     = "Species",
  family     = c(Abundance = "poisson"),
  structure  = NULL,
  dsep       = FALSE,
  WAIC       = FALSE,
  parallel   = TRUE,
  n.cores    = 3,
  n.iter     = 2000,
  engine     = "numpyro",
  quiet      = TRUE
)
dt_joint_s4 <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

sum_joint_s4 <- summary(fit_joint_s4)
joint_s4     <- format_joint_df(sum_joint_s4, "S5a_WrongStruct_None", "4. Misspecified Covariance (structure = NULL)")


final_tests <- rbind(
  tests_baseline,
  tests_s1,
  tests_s2,
  tests_s3,
  tests_s4
)

final_joint <- rbind(
  joint_baseline,
  joint_s3,
  joint_s4
)

final_summary <- rbind(
  format_summary_row(tests_baseline, "S0_True",                        "Baseline (Correct Model)",                     dt_dsep_baseline, dt_joint_baseline),
  format_summary_row(tests_s1,       "S1_WrongDAG_MissingEdges",       "1. Incorrect DAG (Wrong/Missing Edges)",     dt_dsep_s1),
  format_summary_row(tests_s2,       "S3a_OmittedConfounder_BodyMass", "2. Omitted Confounder (Body_Mass Dropped)",  dt_dsep_s2),
  format_summary_row(tests_s3,       "S4_WrongHierarchy_DropSurvey",   "3. Misspecified Hierarchy (Drop Survey Level)", dt_dsep_s3,    dt_joint_s3),
  format_summary_row(tests_s4,       "S5a_WrongStruct_None",           "4. Misspecified Covariance (structure = NULL)", dt_dsep_s4,   dt_joint_s4)
)

all_results <- list(
  S0_True = list(dsep_summary = tests_baseline, joint_summary = joint_baseline),
  S1_WrongDAG_MissingEdges       = list(dsep_summary = tests_s1, joint_summary = NULL),
  S3a_OmittedConfounder_BodyMass = list(dsep_summary = tests_s2, joint_summary = NULL),
  S4_WrongHierarchy_DropSurvey   = list(dsep_summary = tests_s3, joint_summary = joint_s3),
  S5a_WrongStruct_None           = list(dsep_summary = tests_s4, joint_summary = joint_s4)
)

saveRDS(all_results, rds_path)
write.csv(final_summary, sum_path,   row.names = FALSE)
write.csv(final_tests,   tests_path, row.names = FALSE)
write.csv(final_joint,   joint_path, row.names = FALSE)
