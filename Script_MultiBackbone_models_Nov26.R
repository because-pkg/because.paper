rm(list=ls())
setwd("~/Documents/Research Projects/Collaborations and MS ideas/Phil Trans Special Issue/Empirical paper/Dataframes & Analysis")

library(ape)
library(phylolm)
library(because.phybase)

bird.dat<-read.csv("Bird.datMatched8NovScaled.csv")
names(bird.dat)[names(bird.dat) == "Samplesize"]<-"Sample_size"

# Create binary Migrant variable (0 = Sedentary, 1 = Migratory)
# Sedentary birds are all at the standardized minimum of -1.02175 (0 km)
bird.dat$Migrant <- as.numeric(bird.dat$Migration > -1.02)

#Load the three phylogenies for analyses incorporating phylogenetic uncertainty:
McTavish.tree <-read.tree("McTavishmulti.tre")
Hacket.tree <-read.tree("JetzHacketmulti.tre")
Ericson.tree <-read.tree("JetzEricsonmulti.tre")

#Quick check (should get character(0))
setdiff(McTavish.tree[[3]]$tip.label, bird.dat$matched_name)
setdiff(Hacket.tree[[8]]$tip.label, bird.dat$matched_name)

#The McTavish tree has a root edge of 1, we will eliminate it so it does not affect lambda 
#estimates and for consistency as the Jetz backbones don't have it.
McTavish.tree <- lapply(McTavish.tree, function(x) { x$root.edge <- NULL; x })
class(McTavish.tree) <- "multiPhylo"
sum(sapply(McTavish.tree, function(x) !is.null(x$root.edge)))

#Test - controlling for phylogeny - if there is a relationship between sample size
#and lifespan after controlling for body mass (We'll do so with a single tree)
bird.dat2<-bird.dat
row.names(bird.dat2)<-bird.dat2$matched_name
phy.pgls<-McTavish.tree[[1]]

t1<-phylolm(Lifespan~Mass+Sample_size, data=bird.dat2, phy=phy.pgls, model="lambda")
summary(t1)

#Run models in two batches, first we will test the effect of migration per se (independently
#of the distance migrants travel). Then, based on the results we will analyse the relationship between
#migration distance and lifespan, for migrant species only. 

#Run models
#Mass influences
model1<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Lifespan~Mass+Sample_size,
             Migrant~Mass+HWI
)

#To check the model 
plot_dag(model1, latent="Size")

#Mass influences (except migration)

model2<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Lifespan~Mass+Sample_size,
             Migrant~Brain+HWI
)

#To check the model 
plot_dag(model2, latent="Size")

model3<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Migrant~Brain+HWI,
             Lifespan~Clutch+Sample_size
)

plot_dag(model3, latent="Size")

model4<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Migrant~Brain+HWI,
             Clutch~Mass,
             Lifespan~Migrant+Clutch+Sample_size
)

plot_dag(model4, latent="Size")

#Brain influences:

model5<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Brain,
             Migrant~Brain+HWI,
             Lifespan~Brain+Sample_size
)

plot_dag(model5, latent="Size")

model6<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Brain,
             Migrant~Brain+HWI,
             Lifespan~Clutch+Sample_size
)

plot_dag(model6, latent="Size")

model7<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Brain,
             Migrant~Brain+HWI,
             Lifespan~Clutch+Migrant+Sample_size
)

plot_dag(model7, latent="Size")

model8<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Migrant~Mass+Brain+HWI+Clutch,
             Lifespan~Brain+Clutch+Sample_size)

plot_dag(model8, latent="Size")

model9<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Migrant~Mass+Brain+HWI+Clutch,
             Lifespan~Brain+Clutch+HWI+Sample_size)

plot_dag(model9, latent="Size")

model10<-list(Brain~Size,
              Mass~Size,
              HWI~Size,
              Clutch~Brain+Migrant,
              Migrant~Brain+Mass+HWI,
              Lifespan~Brain+Clutch+HWI+Sample_size)

plot_dag(model10, latent="Size")

model11<-list(Brain~Size,
              Mass~Size,
              HWI~Size,
              Clutch~Brain+Migrant,
              Migrant~Brain+HWI,
              Lifespan~Brain+Clutch+Migrant+HWI+Sample_size)

plot_dag(model11, latent="Size")

model12<-list(Brain~Size,
               Mass~Size,
               HWI~Size,
               Clutch~Brain+Migrant,
               Migrant~Brain+HWI,
               Lifespan~Brain+Clutch+Migrant+Sample_size)

plot_dag(model12, latent="Size")

#First test conditional independencies

#Model 1
model1.dsep<-because(
  equations = model1,
  latent = "Size",
  data = bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col = "matched_name",
  dsep = TRUE,
  parallel = TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model1.dsep, rope=c(-0.1,0.1))
plot_rope(model1.dsep)

#Model 2
model2.dsep<-because(
  equations=model2,
  latent="Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model2.dsep, rope=c(-0.1,0.1))
plot_rope(model2.dsep)

#Model 3
model3.dsep<-because(
  equations=model3,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model3.dsep, rope=c(-0.1,0.1))
plot_rope(model3.dsep)

#Model 4
model4.dsep<-because(
  equations=model4,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model4.dsep, rope=c(-0.1,0.1))
plot_rope(model4.dsep)

#Model 5
model5.dsep<-because(
  equations=model5,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model5.dsep, rope = c(-0.1,0.1))
plot_rope(model5.dsep)

#Model 6
model6.dsep<-because(
  equations=model6,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500, 
  n.thin = 1,
  engine = "numpyro"
)

summary(model6.dsep, rope = c(-0.1,0.1))
plot_rope(model6.dsep)

#Model 7
model7.dsep<-because(
  equations=model7,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model7.dsep, rope = c(-0.1,0.1))
plot_rope(model7.dsep)

#Model 8
model8.dsep<-because(
  equations=model8,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model8.dsep, rope = c(-0.1,0.1))
plot_rope(model8.dsep)

#Model 9
model9.dsep<-because(
  equations=model9,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model9.dsep, rope = c(-0.1,0.1))
plot_rope(model9.dsep)

#Model 10
model10.dsep<-because(
  equations=model10,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model10.dsep, rope = c(-0.1,0.1))
plot_rope(model10.dsep)

model11.dsep<-because(
  equations=model11,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model11.dsep, rope = c(-0.1,0.1))
plot_rope(model11.dsep)

model12.dsep<-because(
  equations=model12,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model12.dsep, rope = c(-0.1,0.1))
plot_rope(model12.dsep)

#Can now run full models. Will run only Models 11 and 12 as the others have CIs
#that fail. 

start.time11<-Sys.time()
model11.run<-because(
  equations=model11,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=FALSE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  WAIC = TRUE,
  engine = "numpyro"
)

time.diff11<-Sys.time()-start.time11
time.diff11

summary(model11.run, rope = c(-0.1,0.1))
plot_coef(model11.run)

model11.run.std<-standardize(model11.run)
plot_rope(model11.run.std)

#To export summary to table in Word for Supplementary materials:   

library(flextable)

tab <- as.data.frame(summary(model11.run, rope = c(-0.1,0.1))$results)
tab <- tibble::rownames_to_column(tab, "Parameter")
tab[sapply(tab, is.numeric)] <- lapply(tab[sapply(tab, is.numeric)], round, 2)

greek <- c(alpha = "α", beta = "β", gamma = "γ", delta = "δ", epsilon = "ε",
           zeta = "ζ", eta = "η", theta = "θ", lambda = "λ", mu = "μ",
           nu = "ν", rho = "ρ", sigma = "σ", tau = "τ", phi = "φ",
           psi = "ψ", omega = "ω")

to_greek <- function(x) {
  for (nm in names(greek)) {
    x <- gsub(paste0("(^|_)", nm, "(?=_|$|\\[)"),
              paste0("\\1", greek[[nm]]), x, perl = TRUE)
  }
  x
}

tab$Parameter <- to_greek(tab$Parameter)
save_as_docx(autofit(flextable(tab)), path = "Model11Migrant.docx")

start.time12<-Sys.time()
model12.run<-because(
  equations=model12,
  latent = "Size",
  data=bird.dat,
  structure = list(McTavish = McTavish.tree,Hacket = Hacket.tree,Ericson = Ericson.tree),
  family = c(Migrant = "binomial"),
  id_col="matched_name",
  dsep=FALSE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  WAIC = TRUE,
  engine = "numpyro"
)

time.diff12<-Sys.time()-start.time12
time.diff12

summary(model12.run, rope = c(-0.1,0.1))
plot_coef(model12.run)

model12.run.std<-standardize(model12.run)
plot_coef(model12.run.std)
plot_rope(model12.run.std)

because_compare(model12.run, model11.run)

tab <- as.data.frame(summary(model12.run, rope = c(-0.1,0.1))$results)
tab <- tibble::rownames_to_column(tab, "Parameter")
tab[sapply(tab, is.numeric)] <- lapply(tab[sapply(tab, is.numeric)], round, 2)

greek <- c(alpha = "α", beta = "β", gamma = "γ", delta = "δ", epsilon = "ε",
           zeta = "ζ", eta = "η", theta = "θ", lambda = "λ", mu = "μ",
           nu = "ν", rho = "ρ", sigma = "σ", tau = "τ", phi = "φ",
           psi = "ψ", omega = "ω")

to_greek <- function(x) {
  for (nm in names(greek)) {
    x <- gsub(paste0("(^|_)", nm, "(?=_|$|\\[)"),
              paste0("\\1", greek[[nm]]), x, perl = TRUE)
  }
  x
}

tab$Parameter <- to_greek(tab$Parameter)
save_as_docx(autofit(flextable(tab)), path = "Model12Migrant.docx")



#Create plot of posterior predictive checks for model 

# 1. Open the PDF device
library(gridExtra)

p1<-pp_check(model12.run, resp = "Mass", type = "dens_overlay", ndraws = 100)
p2<-pp_check(model12.run, resp = "Brain", type = "dens_overlay", ndraws = 100)
p3<-pp_check(model12.run, resp = "HWI", type = "dens_overlay", ndraws = 100)
p4<-pp_check(model12.run, resp = "Migrant", type = "bars", ndraws = 100)
p5<-pp_check(model12.run, resp = "Clutch", type = "dens_overlay", ndraws = 100)
p6<-pp_check(model12.run, resp = "Lifespan", type = "dens_overlay", ndraws = 100)

pdf("PostPredChecks.pdf", width = 9, height = 12)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
dev.off()


#Plot of the proposed dags for the empirical example:

library(ggplot2)


coords <- list(
  Size       = c( 0.0,  2.2),
  Mass       = c( 2.4,  3.8),
  Brain      = c( 2.4,  2.2),
  HWI        = c( 2.4,  0.6),
  Clutch     = c( 4.8,  3.6),
  Migrant   = c( 4.8,  0.8),
  Lifespan   = c( 7.2,  2.2),
  Sample_size = c( 7.2,  0.4)
)

# Model 1
plot_dag(model1, coords = coords, latent = "Size")

# Model 2
plot_dag(model2, coords = coords, latent = "Size")

# Model 3
plot_dag(model3, coords = coords, latent = "Size")

# Model 4
plot_dag(model4, coords = coords, latent = "Size")

# Model 5
plot_dag(model5, coords = coords, latent = "Size")

# Model 6
plot_dag(model6, coords = coords, latent = "Size")

# Model 7
plot_dag(model7, coords = coords, latent = "Size")

# Model 8
plot_dag(model8, coords = coords, latent = "Size")

# Model 9
plot_dag(model9, coords = coords, latent = "Size")

# Model 10
plot_dag(model10, coords = coords, latent = "Size")

# Model 11
plot_dag(model11, coords = coords, latent = "Size")

# Model 12
plot_dag(model12, coords = coords, latent = "Size")


# Faceted Grid: All 10 Models Side-by-Side

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
  filename = "ProposedDAGs.pdf",
  plot = p_grid,
  width = 297,
  height = 210,
  units = "mm"
)

#Create PDF of the posterior distributions of the d-seps for the 12 models:

library(ggplot2)
library(patchwork)

# plots: a list of 12 ggplot objects
m1.rope<-plot_rope(model1.dsep)
m2.rope<-plot_rope(model2.dsep)
m3.rope<-plot_rope(model3.dsep)
m4.rope<-plot_rope(model4.dsep)
m5.rope<-plot_rope(model5.dsep)
m6.rope<-plot_rope(model6.dsep)
m7.rope<-plot_rope(model7.dsep)
m8.rope<-plot_rope(model8.dsep)
m9.rope<-plot_rope(model9.dsep)
m10.rope<-plot_rope(model10.dsep)
m11.rope<-plot_rope(model11.dsep)
m12.rope<-plot_rope(model12.dsep)
plots <- list(m1.rope,m2.rope,m3.rope,m4.rope,m5.rope,m6.rope,m7.rope,m8.rope,m9.rope,m10.rope,m11.rope,m12.rope)

#Wrap the titles
plots <- lapply(plots, function(p) {
  p + labs(title = stringr::str_wrap(p$labels$title, width = 35))
})

# 6 plots per page (3 rows x 2 cols) -> 2 pages
per_page <- 6
pages <- split(plots, ceiling(seq_along(plots) / per_page))

pdf("d-sep_ROPEs.pdf", width = 8.27, height = 11.69)   # A4 portrait
for (pg in pages) {
  print(
    wrap_plots(pg, ncol = 2, nrow = 3) &
  theme(plot.title.position = "plot",
        plot.margin = margin(t = 10, r = 35, b = 10, l = 10),
        plot.title = element_text(size = 9),
        legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm")
        )
  )
}
dev.off()


#####################################################
#Now repeat analyses using only migrant species and with
#migration distance as a continuous variable. Will need to
#subset the dataframe and prune the three backbones to include
#only the migrants.

#Subset dataframe to keep only migrant species (non-migrants have value of -1.02 after log transformation and
#standardization)
bird.datM<-subset(bird.dat, bird.dat$Migration > -1.02)

#Prune the phylogenies
to.dropMcT<-setdiff(McTavish.tree[[3]]$tip.label, bird.datM$matched_name)
to.dropH<-setdiff(Hacket.tree[[8]]$tip.label, bird.datM$matched_name)
to.dropE<-setdiff(Ericson.tree[[6]]$tip.label, bird.datM$matched_name)

McTavish.treeM<-drop.tip(McTavish.tree, to.dropMcT)
Hacket.treeM<-drop.tip(Hacket.tree, to.dropH)
Ericson.treeM<-drop.tip(Ericson.tree, to.dropE)

#These are the same models but with Migration as we want the migration distance
#Mass influences
model1M<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Lifespan~Mass+Sample_size,
             Migration~Mass+HWI
)

#To check the model 
plot_dag(model1M, latent="Size")

#Mass influences (except migration)

model2M<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Lifespan~Mass+Sample_size,
             Migration~Brain+HWI
)

#To check the model 
plot_dag(model2M, latent="Size")

model3M<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Migration~Brain+HWI,
             Lifespan~Clutch+Sample_size
)

plot_dag(model3M, latent="Size")

model4M<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Migration~Brain+HWI,
             Clutch~Mass,
             Lifespan~Migration+Clutch+Sample_size
)

plot_dag(model4M, latent="Size")

#Brain influences:

model5M<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Brain,
             Migration~Brain+HWI,
             Lifespan~Brain+Sample_size
)

plot_dag(model5M, latent="Size")

model6M<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Brain,
             Migration~Brain+HWI,
             Lifespan~Clutch+Sample_size
)

plot_dag(model6M, latent="Size")

model7M<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Brain,
             Migration~Brain+HWI,
             Lifespan~Clutch+Migration+Sample_size
)

plot_dag(model7M, latent="Size")

model8M<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Migration~Mass+Brain+HWI+Clutch,
             Lifespan~Brain+Clutch+Sample_size)

plot_dag(model8M, latent="Size")

model9M<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Migration~Mass+Brain+HWI+Clutch,
             Lifespan~Brain+Clutch+HWI+Sample_size)

plot_dag(model9M, latent="Size")

model10M<-list(Brain~Size,
              Mass~Size,
              HWI~Size,
              Clutch~Brain+Migration,
              Migration~Brain+Mass+HWI,
              Lifespan~Brain+Clutch+HWI+Sample_size)

plot_dag(model10M, latent="Size")

model11M<-list(Brain~Size,
              Mass~Size,
              HWI~Size,
              Clutch~Brain+Migration,
              Migration~Brain+HWI,
              Lifespan~Brain+Clutch+Migration+HWI+Sample_size)

plot_dag(model11M, latent="Size")

model12M<-list(Brain~Size,
              Mass~Size,
              HWI~Size,
              Clutch~Brain+Migration,
              Migration~Brain+HWI,
              Lifespan~Brain+Clutch+Migration+Sample_size)

plot_dag(model12M, latent="Size")

#Test the conditional independencies of the same models as above

#Model 1
model1M.dsep<-because(
  equations = model1M,
  latent = "Size",
  data = bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col = "matched_name",
  dsep = TRUE,
  parallel = TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model1M.dsep, rope=c(-0.1,0.1))
plot_rope(model1M.dsep)

#Model 2
model2M.dsep<-because(
  equations=model2M,
  latent="Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model2M.dsep, rope=c(-0.1,0.1))
plot_rope(model2M.dsep)

#Model 3
model3M.dsep<-because(
  equations=model3M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model3M.dsep, rope=c(-0.1,0.1))
plot_rope(model3M.dsep)

#Model 4
model4M.dsep<-because(
  equations=model4M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model4M.dsep, rope=c(-0.1,0.1))
plot_rope(model4M.dsep)

#Model 5
model5M.dsep<-because(
  equations=model5M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model5M.dsep, rope = c(-0.1,0.1))
plot_rope(model5M.dsep)

#Model 6
model6M.dsep<-because(
  equations=model6M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500, 
  n.thin = 1,
  engine = "numpyro"
)

summary(model6M.dsep, rope = c(-0.1,0.1))
plot_rope(model6M.dsep)

#Model 7
model7M.dsep<-because(
  equations=model7M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model7M.dsep, rope = c(-0.1,0.1))
plot_rope(model7M.dsep)

#Model 8
model8M.dsep<-because(
  equations=model8M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model8M.dsep, rope = c(-0.1,0.1))
plot_rope(model8M.dsep)

#Model 9
model9M.dsep<-because(
  equations=model9M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model9M.dsep, rope = c(-0.1,0.1))
plot_rope(model9M.dsep)

#Model 10
model10M.dsep<-because(
  equations=model10M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model10M.dsep, rope = c(-0.1,0.1))
plot_rope(model10M.dsep)

model11M.dsep<-because(
  equations=model11M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model11M.dsep, rope = c(-0.1,0.1))
plot_rope(model11M.dsep)

model12M.dsep<-because(
  equations=model12M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  engine = "numpyro"
)

summary(model12M.dsep, rope = c(-0.1,0.1))
plot_rope(model12M.dsep)

#Now compare models 11 and 12 both of which meet the conditional independencies

model11M.run<-because(
  equations=model11M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=FALSE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  WAIC = TRUE,
  engine = "numpyro"
)

summary(model11M.run, rope = c(-0.1,0.1))
plot_coef(model11M.run)
plot_rope(model11M.run)

model12M.run<-because(
  equations=model12M,
  latent = "Size",
  data=bird.datM,
  structure = list(McTavish = McTavish.treeM,Hacket = Hacket.treeM,Ericson = Ericson.treeM),
  id_col="matched_name",
  dsep=FALSE,
  parallel=TRUE,
  n.chains = 2,
  n.cores = 8,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 1,
  WAIC = TRUE,
  engine = "numpyro"
)

summary(model12M.run, rope = c(-0.1,0.1))
plot_coef(model12M.run)
plot_rope(model12M.run)

because_compare(model11M.run, model12M.run)

#Create a table with the results of model 11, which fits the data better than model 12

tab <- as.data.frame(summary(model11M.run, rope = c(-0.1,0.1))$results)
tab <- tibble::rownames_to_column(tab, "Parameter")
tab[sapply(tab, is.numeric)] <- lapply(tab[sapply(tab, is.numeric)], round, 2)

greek <- c(alpha = "α", beta = "β", gamma = "γ", delta = "δ", epsilon = "ε",
           zeta = "ζ", eta = "η", theta = "θ", lambda = "λ", mu = "μ",
           nu = "ν", rho = "ρ", sigma = "σ", tau = "τ", phi = "φ",
           psi = "ψ", omega = "ω")

to_greek <- function(x) {
  for (nm in names(greek)) {
    x <- gsub(paste0("(^|_)", nm, "(?=_|$|\\[)"),
              paste0("\\1", greek[[nm]]), x, perl = TRUE)
  }
  x
}

tab$Parameter <- to_greek(tab$Parameter)
save_as_docx(autofit(flextable(tab)), path = "Model11Migration.docx")


