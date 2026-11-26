library(because)
library(ggplot2)


coords <- list(
  Size       = c( 0.0,  2.2),
  Mass       = c( 2.4,  3.8),
  Brain      = c( 2.4,  2.2),
  HWI        = c( 2.4,  0.6),
  Clutch     = c( 4.8,  3.6),
  Migrant   = c( 4.8,  0.8),
  Life_span   = c( 7.2,  2.2),
  Sample_size = c( 7.2,  0.4)
)

# Model 1
model1 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Clutch ~ Mass,
  Life_span ~ Mass + Sample_size,
  Migrant  ~ Mass + HWI
)
plot_dag(model1, coords = coords, latent = "Size")

# Model 2
model2 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Clutch ~ Mass,
  Life_span ~ Mass + Sample_size,
  Migrant  ~ Brain + HWI
)
plot_dag(model2, coords = coords, latent = "Size")

# Model 3
model3 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Clutch ~ Mass,
  Migrant  ~ Brain + HWI,
  Life_span ~ Clutch + Sample_size
)
plot_dag(model3, coords = coords, latent = "Size")

# Model 4
model4 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Migrant ~ Brain + HWI,
  Clutch ~ Mass,
  Life_span ~ Migrant + Clutch + Sample_size
)
plot_dag(model4, coords = coords, latent = "Size")

# Model 5
model5 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Clutch ~ Brain,
  Migrant ~ Brain + HWI,
  Life_span ~ Brain + Sample_size
)
plot_dag(model5, coords = coords, latent = "Size")

# Model 6
model6 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Clutch ~ Brain,
  Migrant ~ Brain + HWI,
  Life_span ~ Clutch + Sample_size
)
plot_dag(model6, coords = coords, latent = "Size")

# Model 7
model7 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Clutch ~ Brain,
  Migrant ~ Brain + HWI,
  Life_span ~ Clutch + Migrant + Sample_size
)
plot_dag(model7, coords = coords, latent = "Size")

# Model 8
model8 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Clutch ~ Mass,
  Migrant ~ Mass + Brain + HWI + Clutch,
  Life_span ~ Brain + Clutch + Sample_size
)
plot_dag(model8, coords = coords, latent = "Size")

# Model 9
model9 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Clutch ~ Mass,
  Migrant ~ Mass + Brain + HWI + Clutch,
  Life_span ~ Brain + Clutch + HWI + Sample_size
)
plot_dag(model9, coords = coords, latent = "Size")

# Model 10
model10 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Clutch ~ Brain + Migrant,
  Migrant ~ Brain + Mass + HWI,
  Life_span ~ Brain + Clutch + HWI + Sample_size
)
plot_dag(model10, coords = coords, latent = "Size")

# Model 11
model11 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Clutch ~ Brain + Migrant,
  Migrant ~ Brain + HWI,
  Life_span ~ Brain + Clutch + Migrant +  HWI + Sample_size
)
plot_dag(model11, coords = coords, latent = "Size")

# Model 12
model12 <- list(
  Brain ~ Size,
  Mass ~ Size,
  HWI ~ Size,
  Clutch ~ Brain + Migrant,
  Migrant ~ Brain + HWI,
  Life_span ~ Brain + Clutch + Migrant + Sample_size
)
plot_dag(model12, coords = coords, latent = "Size")


# Faceted Grid: All 12 Models Side-by-Side

all_models <- list(
  `Model 1`  = model1, `Model 2`  = model2, `Model 3`  = model3, `Model 4`  = model4,
  `Model 5`  = model5, `Model 6`  = model6, `Model 7`  = model7, `Model 8`  = model8,
  `Model 9`  = model9, `Model 10` = model10, `Model 11` = model11,`Model 12` = model12
)

p_grid <- plot_dag(
  all_models,
  coords    = coords,
  latent    = "Size",
  node_size = 10,
  node_stroke = 0.7,
  text_size = 1.6,
  edge_width_range =c(0.5,0.5)
) +
  facet_wrap(~model_label, ncol = 4) +
  theme(strip.text = element_text(face = "bold", size = 11))

print(p_grid)
ggsave(
  filename = "bird_models.pdf",
  plot = p_grid,
  width = 297,
  height = 210,
  units = "mm"
)
