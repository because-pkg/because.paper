# multiscale simulation

library(because)
library(because.phybase)
library(ape)
library(MASS)

set.seed(1234)

# SIMULATION PARAMETERS
N_Species          <- 50
N_Site             <- 50
N_Surveys_per_Site <- 3
N_Surveys          <- N_Site * N_Surveys_per_Site

# Causal paths
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


# PHYLOGENY
species_tree <- rtree(N_Species)
species_tree$tip.label  <- paste0("Sp_", 1:N_Species)
species_tree$edge.length <- species_tree$edge.length / max(node.depth.edgelength(species_tree))
Vcv <- vcv.phylo(species_tree)

lambda <- 0.6
V_phylo <- lambda * Vcv + (1 - lambda) * diag(diag(Vcv))

# 3. SPECIES LEVEL
Body_Mass <- rnorm(N_Species, 0, 1.0)   # Simple random normal
eps_MR    <- as.vector(mvrnorm(1, rep(0, N_Species), 0.3^2 * V_phylo))
Metabolic_Rate <- BETA_BM_MR * Body_Mass + eps_MR
eps_TT    <- as.vector(mvrnorm(1, rep(0, N_Species), 0.3^2 * V_phylo))

Thermal_Tol <- BETA_BM_TT * Body_Mass + eps_TT


d_species <- data.frame(
  Species        = species_tree$tip.label,
  Body_Mass      = Body_Mass,
  Metabolic_Rate = Metabolic_Rate,
  Thermal_Tol    = Thermal_Tol
)

# SPATIAL SITE-LEVEL
site_coords <- cbind(runif(N_Site, 0, 10), runif(N_Site, 0, 10))
rownames(site_coords) <- paste0("Site_", 1:N_Site)
dist_mat <- as.matrix(dist(site_coords))

# Spatial gradient for residuals
V_spatial <- exp(-dist_mat / 5.0)

Elevation <- rnorm(N_Site, 0, 1.0)
NDVI      <- BETA_ELEV_NDVI * Elevation + rnorm(N_Site, 0, 1.0)

d_site <- data.frame(
  Site         = paste0("Site_", 1:N_Site),
  Elevation    = Elevation,
  NDVI         = NDVI
)


# SURVEY LEVEL

site_of_survey  <- rep(1:N_Site, each = N_Surveys_per_Site)
u_site_temp <- as.vector(mvrnorm(1, rep(0, N_Site), 0.4^2 * V_spatial))
Temperature <- BETA_ELEV_TEMP * Elevation[site_of_survey] + u_site_temp[site_of_survey] + rnorm(N_Surveys, 0, 1.0)
Wind_Speed  <- rnorm(N_Surveys, 0, 1.0)

u_site_FC <- as.vector(mvrnorm(1, rep(0, N_Site), 0.4^2 * V_spatial))
Flower_Cover <- BETA_NDVI_FC * NDVI[site_of_survey] + BETA_ELEV_FC * Elevation[site_of_survey] + u_site_FC[site_of_survey] + rnorm(N_Surveys, 0, 1.0)

d_survey <- data.frame(
  Survey      = paste0("Survey_", 1:N_Surveys),
  Site        = paste0("Site_", site_of_survey),
  Temperature = Temperature,
  Wind_Speed  = Wind_Speed,
  Flower_Cover = Flower_Cover
)

# OBSERVATION LEVEL

d_obs <- expand.grid(Survey = d_survey$Survey, Species = d_species$Species, stringsAsFactors = FALSE)
d_obs <- merge(d_obs, d_survey[, c("Survey", "Site", "Temperature", "Wind_Speed", "Flower_Cover")], by = "Survey")
d_obs <- merge(d_obs, d_species, by = "Species")
d_obs$obs_id <- seq_len(nrow(d_obs))

Thermal_Mismatch <- d_obs$Temperature - d_obs$Thermal_Tol

# Generate the residuals
u_phylo_abund   <- as.vector(mvrnorm(1, rep(0, N_Species), 0.3^2 * V_phylo))
u_spatial_abund <- as.vector(mvrnorm(1, rep(0, N_Site), 0.3^2 * V_spatial))


#Add the names
names(u_phylo_abund)   <- paste0("Sp_", 1:N_Species)
names(u_spatial_abund) <- paste0("Site_", 1:N_Site)


# Calculate log_abund
log_abund <- BASELINE_INTERCEPT +
  BETA_FC_ABUND       * d_obs$Flower_Cover   +
  BETA_WIND_ABUND     * d_obs$Wind_Speed      +
  BETA_MISMATCH_ABUND * Thermal_Mismatch      +
  u_phylo_abund[d_obs$Species]                +
  u_spatial_abund[d_obs$Site]


# Generate Poisson data
d_obs$Abundance <- rpois(nrow(d_obs), lambda = exp(log_abund))

#Aproximate log abundance for d-sep tests
d_obs$log_Abund <- log(d_obs$Abundance + 1)


data_list <- list(
    species = d_species,
    site    = d_site,
    survey  = d_survey[, c("Survey", "Site", "Temperature", "Wind_Speed", "Flower_Cover")],
    obs     = d_obs[, c("obs_id", "Survey", "Species", "log_Abund")]
)

structures <- list(
    phylo   = species_tree,
    spatial = V_spatial
)

eqs <- list(
    Flower_Cover ~ NDVI + Elevation,
    Body_Mass    ~ 1,
    Metabolic_Rate ~ Body_Mass,
    Thermal_Tol    ~ Body_Mass,
    Temperature  ~ Elevation,
    log_Abund ~ Flower_Cover + Wind_Speed +
                I(Temperature - Thermal_Tol)
)


start_time <- Sys.time()
fit.jg_dsep <- because(
  eqs,
  data      = data_list,
  levels    = list(
    species = c("Body_Mass", "Metabolic_Rate", "Thermal_Tol", "Species"),
    site    = c("Elevation", "NDVI", "Site"),
    survey  = c("Temperature", "Wind_Speed", "Flower_Cover", "Survey"),
    obs     = c("log_Abund")),
  hierarchy    = "site > survey > obs; species > obs",
  link_vars    = list(site = "Site", survey = "Survey", species = "Species"),
  id_col       = "Species",
  structure    = structures,
  dsep         = TRUE,
  parallel     = TRUE,
  n.cores      = 3,
  n.iter= 5000,
  engine       = "jags"
)
end_time <- Sys.time()
time_diff <- end_time - start_time

#FIT FINAL JOINT MODEL

data_list2 <- list(
  species = d_species,
  site    = d_site,
  survey  = d_survey[, c("Survey", "Site", "Temperature", "Wind_Speed", "Flower_Cover")],
  obs     = d_obs[, c("obs_id", "Survey", "Species", "Abundance")]
)


eqs2 <- list(
  Flower_Cover ~ NDVI + Elevation,
  Body_Mass    ~ 1,
  Metabolic_Rate ~ Body_Mass,
  Thermal_Tol    ~ Body_Mass,
  Temperature  ~ Elevation,
  Abundance ~ Flower_Cover + Wind_Speed +
    I(Temperature - Thermal_Tol)
)


start_time <- Sys.time()
fit.jg <- because(
  eqs2,
  data      = data_list,
  levels    = list(
    species = c("Body_Mass", "Metabolic_Rate", "Thermal_Tol", "Species"),
    site    = c("Elevation", "NDVI", "Site"),
    survey  = c("Temperature", "Wind_Speed", "Flower_Cover", "Survey"),
    obs     = c("Abundance")),
  hierarchy    = "site > survey > obs; species > obs",
  link_vars    = list(site = "Site", survey = "Survey", species = "Species"),
  id_col       = "Species",
  family      = c(Abundance = "poisson"),
  structure    = structures,
  dsep         = FALSE,
  parallel     = TRUE,
  n.cores      = 3,
  n.iter       = 50000,
  engine       = "jags",
  n.thin = 1
)
end_time <- Sys.time()
time_diff4 <- start_time - end_time
time_diff4

# PLOTS

#DAG plot with background boxes

library(ggplot2)
coords <- list(
  Body_Mass      = c(0, 4),
  Metabolic_Rate = c(-0.5, 3),
  Thermal_Tol    = c(0.5, 3),
  Elevation      = c(3, 3),
  NDVI           = c(4, 3),
  Flower_Cover   = c(3.5, 1.5),
  Temperature    = c(2.5, 1.5),
  Wind_Speed     = c(4.5, 1.5),
  Abundance      = c(2.5, -0.5),
  "I(Temperature - Thermal_Tol)" = c(1, 1)
)
p <- plot_dag(eqs, coords = coords)

background_boxes <- list(
  # Species Level Box (Red)
  annotate("rect", xmin = -1, xmax = 1, ymin = 2.6, ymax = 4.4,
           alpha = 0.1, fill = "red"),

  # Site Level Box (Green)
  annotate("rect", xmin = 2.6, xmax = 4.5, ymin = 2.6, ymax = 3.4,
           alpha = 0.1, fill = "forestgreen"),

  # Survey Level Box (Blue)
  annotate("rect", xmin = -0.5, xmax = 5.0, ymin = 0.5, ymax = 1.9,
           alpha = 0.1, fill = "blue"),

  # Observation Level Box (Purple)
  annotate("rect", xmin = 2, xmax = 3, ymin = -0.9, ymax = 0,
           alpha = 0.1, fill = "purple")
)

p$layers <- c(background_boxes, p$layers)

p +
  annotate("text", x = 0, y = 4.6, label = "Species Level",
           fontface = "bold", color = "darkred") +

  annotate("text", x = 3.5, y = 3.6, label = "Site Level",
           fontface = "bold", color = "forestgreen") +

  annotate("text", x = 5.0, y = 2.1, label = "Survey Level",
           fontface = "bold", color = "darkblue", hjust = 1) +


  annotate("text", x = 2.5, y = -1, label = "Observation Level",
           fontface = "bold", color = "purple") +

  coord_cartesian(xlim = c(-1.5, 5.5), ylim = c(-1, 4))


#Coefficients plot with true values

c_plot <- plot_coef(fit.jg2, type = "raw")

true_vals <- data.frame(
  Path = c(
    "Metabolic_Rate ~ Body_Mass",
    "Thermal_Tol ~ Body_Mass",
    "Flower_Cover ~ NDVI",
    "Flower_Cover ~ Elevation",
    "Temperature ~ Elevation",
    "Abundance ~ Flower_Cover",
    "Abundance ~ Wind_Speed",
    "Abundance ~ Temperature_minus_Thermal_Tol"
  ),
  Estimate = c(
    0.8,    # BETA_BM_MR
    0.6,   # BETA_BM_TT
    0.6,    # BETA_NDVI_FC
    0.3,    # BETA_ELEV_FC
  -1.5,     # BETA_ELEV_TEMP
    0.2,    # BETA_FC_ABUND
    -0.2,   # BETA_WIND_ABUND
    -0.3    # BETA_MISMATCH_ABUND
)

c_plot + geom_point(
  data = true_vals,
  aes(x = Path, y = Estimate),
  color = "red",
  size = 3,
  shape = 1,       # 1 = hollow circle (use 16 if you want a solid red circle)
  stroke = 1.2,    # Thickness of the circle's line
  inherit.aes = FALSE # <--- This stops it from looking for Lower/Upper bounds!
)   + labs(title = NULL, subtitle = NULL)


#PP_check plots

library(patchwork)

p_mr    <- pp_check(fit.jg2, resp = "Metabolic_Rate") + labs(title = "Metabolic Rate")
p_tt    <- pp_check(fit.jg2, resp = "Thermal_Tol")    + labs(title = "Thermal Tol")
p_temp  <- pp_check(fit.jg2, resp = "Temperature")    + labs(title = "Temperature")         + labs(title = "NDVI")
p_fc    <- pp_check(fit.jg2, resp = "Flower_Cover")   + labs(title = "Flower Cover")
p_abund <- pp_check(fit.jg2, resp = "Abundance")      + labs(title = "Abundance")

combined_plot <- (p_mr | p_tt) /
  (p_temp | p_fc) /
  (p_abund | plot_spacer())

