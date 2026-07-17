
library(ape)
library(because.phybase)

bird.dat<-read.csv("Bird.datMatchedScaledTarsus.csv")
bird.phy.cut<-read.nexus("MultiPhyJetzMcTav.nex")

to.drop<-setdiff(bird.phy.cut[[1]]$tip.label, bird.dat$matched_name)
bird.phy.cut<-drop.tip(bird.phy.cut, to.drop)


#Run models
#Mass influences
model1<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Lifespan~Mass,
             Migration~Mass+HWI
)

#To check the model 
plot_dag(model1)

#Mass influences (except migration)

model2<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Lifespan~Mass,
             Migration~Brain+HWI
)

#To check the model 
plot_dag(model2)

model3<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Migration~Brain+HWI,
             Lifespan~Clutch
)

plot_dag(model3)

model4<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Migration~Brain+HWI,
             Clutch~Mass,
             Lifespan~Migration+Clutch
)

plot_dag(model4)

#Brain influences:

model5<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Brain,
             Migration~Brain+HWI,
             Lifespan~Brain
)

plot_dag(model5)

model6<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Brain,
             Migration~Brain+HWI,
             Lifespan~Clutch
)

plot_dag(model6)

model7<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Brain,
             Migration~Brain+HWI,
             Lifespan~Clutch+Migration
)

plot_dag(model7)

model8<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Migration~Mass+Brain+HWI+Clutch,
             Lifespan~Brain+Clutch)

plot_dag(model8)

model9<-list(Brain~Size,
             Mass~Size,
             HWI~Size,
             Clutch~Mass,
             Migration~Mass+Brain+HWI+Clutch,
             Lifespan~Brain+Clutch+HWI)

plot_dag(model9)


model10<-list(Brain~Size,
              Mass~Size,
              HWI~Size,
              Clutch~Brain+Migration,
              Migration~Brain+Mass+HWI,
              Lifespan~Brain+Clutch+HWI)

plot_dag(model10)



#First test conditional independencies

#Model 1
model1.dsep<-because(
  equations=model1,
  latent="Size",
  data=bird.dat,
  structure = bird.phy.cut,
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine = "jags"
)

summary(model1.dsep, rope=c(-0.1,0.1))
plot_rope(model1.dsep)

#Model 2
model2.dsep<-because(
  equations=model2,
  latent="Size",
  data=bird.dat,
  structure = bird.phy.cut,
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine = "jags",
  reuse_models = list(model1.dsep)
)

summary(model2.dsep, rope=c(-0.1,0.1))
plot_rope(model2.dsep)

#Model 3
model3.dsep<-because(
  equations=model3,
  latent = "Size",
  data=bird.dat,
  structure = bird.phy.cut,
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine = "jags",
  reuse_models = list(model1.dsep,model2.dsep)
)

summary(model3.dsep, rope=c(-0.1,0.1))
plot_rope(model3.dsep)

#Model 4
model4.dsep<-because(
  equations=model4,
  latent = "Size",
  data=bird.dat,
  structure = bird.phy.cut,
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine = "jags",
  reuse_models = list(model1.dsep,model2.dsep,model3.dsep)
)

summary(model4.dsep, rope=c(-0.1,0.1))
plot_rope(model4.dsep)

#Model 5
model5.dsep<-because(
  equations=model5,
  latent = "Size",
  data=bird.dat,
  structure = bird.phy.cut,
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine = "jags",
  reuse_models = list(model1.dsep,model2.dsep,model3.dsep,
                      model4.dsep)
)

summary(model5.dsep, rope = c(-0.1,0.1))
plot_rope(model5.dsep)

#Model 6
model6.dsep<-because(
  equations=model6,
  latent = "Size",
  data=bird.dat,
  structure = bird.phy.cut,
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine = "jags",
  reuse_models = list(model1.dsep,model2.dsep,model3.dsep,
                      model4.dsep)
)

summary(model6.dsep, rope = c(-0.1,0.1))
plot_rope(model6.dsep)

#Model 7
model7.dsep<-because(
  equations=model7,
  latent = "Size",
  data=bird.dat,
  structure = bird.phy.cut,
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine = "jags",
  reuse_models = list(model1.dsep,model2.dsep,model3.dsep,
                      model4.dsep,model6.dsep)
)

summary(model7.dsep, rope = c(-0.1,0.1))
plot_rope(model7.dsep)

#Model 8
model8.dsep<-because(
  equations=model8,
  latent = "Size",
  data=bird.dat,
  structure = bird.phy.cut,
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine = "jags",
  reuse_models = list(model1.dsep,model2.dsep,model3.dsep,
                      model4.dsep,model6.dsep,
                      model7.dsep)
)

summary(model8.dsep, rope = c(-0.1,0.1))
plot_rope(model8.dsep)

#Model 9
model9.dsep<-because(
  equations=model9,
  latent = "Size",
  data=bird.dat,
  structure = bird.phy.cut,
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine = "jags",
  reuse_models = list(model1.dsep,model2.dsep,model3.dsep,
                      model4.dsep,model6.dsep,
                      model7.dsep,model8.dsep)
)

summary(model9.dsep, rope = c(-0.1,0.1))
plot_rope(model9.dsep)

#Model 10
model10.dsep<-because(
  equations=model10,
  latent = "Size",
  data=bird.dat,
  structure = bird.phy.cut,
  id_col="matched_name",
  dsep=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine = "jags",
  reuse_models = list(model1.dsep,model2.dsep,model3.dsep,
                      model4.dsep,model6.dsep,
                      model7.dsep,model8.dsep,model9.dsep)
)

summary(model10.dsep, rope = c(-0.1,0.1))
plot_rope(model10.dsep)


#Can now run full models. Will run only from Models 9, 10 and 11 as the others have CIs
#that fail. I include model 8 only to test whether the causal link between HWI and lifespan 
#is significant, which I think is interesting to check.

start.time8<-Sys.time()
model8.run<-because(
  equations=model8,
  data=bird.dat,
  latent="Size",
  structure = bird.phy.cut,
  id_col="matched_name",
  WAIC=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine="numpyro"
)
time_diff8<-Sys.time() - start.time8
time_diff8

summary(model8.run)
plot_coef(model8.run)

#Run models for longer to see if the Rhat and ESS of the lambdas and other evol parameters improve

start.time9<-Sys.time()
model9.run<-because(
  equations=model9,
  data=bird.dat,
  latent="Size",
  structure = bird.phy.cut,
  id_col="matched_name",
  n.thin = 1,
  WAIC=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine="jags"
)
time_diff9<-Sys.time() - start.time9
time_diff9

summary(model9.run)
plot_coef(model9.run)

start.time10<-Sys.time()
model10.run<-because(
  equations=model10,
  data=bird.dat,
  latent="Size",
  structure = bird.phy.cut,
  id_col="matched_name",
  n.thin= 1,
  WAIC=TRUE,
  parallel=TRUE,
  n.cores=3,
  engine="jags"
)

time_diff10<-Sys.time() - start.time10
time_diff10

summary(model10.run)
plot_coef(model10.run)


#Create plot of posterior predictive checks for model 10

# 1. Open the PDF device
library(gridExtra)

p1<-pp_check(model10.run, resp = "Mass", type = "dens_overlay", ndraws = 100)
p2<-pp_check(model10.run, resp = "Brain", type = "dens_overlay", ndraws = 100)
p3<-pp_check(model10.run, resp = "HWI", type = "dens_overlay", ndraws = 100)
p4<-pp_check(model10.run, resp = "Migration", type = "dens_overlay", ndraws = 100)
p5<-pp_check(model10.run, resp = "Clutch", type = "dens_overlay", ndraws = 100)
p6<-pp_check(model10.run, resp = "Lifespan", type = "dens_overlay", ndraws = 100)

pdf("PostPredChecks.pdf", width = 9, height = 12)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
dev.off()


#Plot of the proposed dags for the empirical example: 

library(ggplot2)
library(gridExtra)

tighten <- theme(
  plot.margin = margin(1, 1, 1, 1, unit = "pt")  # top, right, bottom, left
)

m1 <- plot_dag(model1) + tighten + ggtitle("Model1")
m2 <- plot_dag(model2) + tighten + ggtitle("Model2")
m3 <- plot_dag(model3) + tighten + ggtitle("Model3")
m4 <- plot_dag(model4) + tighten + ggtitle("Model4")
m5 <- plot_dag(model5) + tighten + ggtitle("Model5")
m6 <- plot_dag(model6) + tighten + ggtitle("Model6")
m7 <- plot_dag(model7) + tighten + ggtitle("Model7")
m8 <- plot_dag(model8) + tighten + ggtitle("Model8")
m9 <- plot_dag(model9) + tighten + ggtitle("Model9")
m10 <- plot_dag(model10) + tighten + ggtitle("Model10")

pdf("ProposedDAGs.pdf", width = 11, height = 16)
grid.arrange(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, ncol = 2)
dev.off()
