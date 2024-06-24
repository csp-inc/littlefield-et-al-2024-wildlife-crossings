library(amt)
library(sf)
library(tidyverse)
library(glmmTMB)
library(car)
library(AICcmodavg)

# Load SSF dataset
load(file = "data/south-ute-deer-elk-data/elk-migration-ssf-data.rda", verbose = TRUE) 

# Examine correlations between covariates
mig_ssf %>% dplyr::select("cos_ta", "log_sl","rd_cross", "crossing", "elevation","slope",
                          "aspect","chili","water","dev","ag","forest_pcov","well",
                          "tmax","tmin","cmd","pas","tr_scale") %>% 
  filter(!is.na(cos_ta)) %>% st_drop_geometry() %>% cor()

# Model fitting function using glmmTMB and Muff et al 2020 approach
modFit <- function(data, model){
  mod <- glmmTMB(model, data = data, family = poisson, doFit = FALSE,
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
  mod$parameters$theta[1] <- log(1e3)
  ntheta = length(mod$parameters$theta)
  mod$mapArg <- list(theta=factor(c(NA, rep(1,ntheta-1))))
  mod_fit <- glmmTMB:::fitTMB(mod)
  return(mod_fit)
}

# #---------------------------
# # First identify best road covariates (i.e., traffic volume or just road crossing 0-1, size of crossing structures)
# rd1 <- formula(case ~ -1 + rd_cross * cross100 + log_sl + cos_ta + (1 | step_id) + (0 + cross100 + rd_cross | Animal_id))
# rd2 <- formula(case ~ -1 + rd_cross * cross500 + log_sl + cos_ta + (1 | step_id) + (0 + cross500 + rd_cross | Animal_id))
# rd3 <- formula(case ~ -1 + tr_scale * cross100 + log_sl + cos_ta + (1 | step_id) + (0 + cross100 + tr_scale | Animal_id))
# rd4 <- formula(case ~ -1 + tr_scale * cross500 + log_sl + cos_ta + (1 | step_id) + (0 + cross500 + tr_scale | Animal_id))
# 
# rd_mods <- list()
# rd_mods[[1]] <- modFit(mig_ssf, rd1)
# rd_mods[[2]] <- modFit(mig_ssf, rd2)
# rd_mods[[3]] <- modFit(mig_ssf, rd3)
# rd_mods[[4]] <- modFit(mig_ssf, rd4)
# names(rd_mods) <- c('rd100','rd500','tr100','tr500')
# rd_tab <- aictab(rd_mods)

#---------------------------
# Define additional models, building on top road submodel

# Model formulae
anthro1 <- formula(case ~ -1 + tr_scale + log_sl + cos_ta + (1 | step_id) + (0 + tr_scale | Animal_id))
anthro2 <- formula(case ~ -1 + tr_scale + dev + dev_sq + ag + ag_sq + log_sl + cos_ta + 
                     (1 | step_id) + (0 + tr_scale + dev + ag | Animal_id))
tercli1 <- formula(case ~ -1 + elevation + slope + aspect + log_sl + cos_ta + 
                   (1 | step_id) + (0 + elevation + slope + aspect | Animal_id))
tercli2 <- formula(case ~ -1 + slope + aspect + tmin + pas + cmd + log_sl + cos_ta + 
                   (1 | step_id) + (0 + slope + aspect + tmin + pas + cmd | Animal_id))
lc1 <- formula(case ~ -1 + dev + dev_sq + ag + ag_sq + water + forest_pcov + log_sl + cos_ta + 
                     (1 | step_id) + (0 + dev + ag + water + forest_pcov | Animal_id))
combo1 <- formula(case ~ -1 + tr_scale + dev + dev_sq + ag + ag_sq + slope + 
                    dev:slope + aspect + water + forest_pcov + log_sl + cos_ta + 
                     (1 | step_id) + (0 + tr_scale + dev + ag + slope + aspect + water + forest_pcov | Animal_id))
combo2 <- formula(case ~ -1 + tr_scale + dev +  ag + slope + aspect + water + forest_pcov + log_sl + cos_ta + 
                    (1 | step_id) + (0 + tr_scale + dev +  ag + slope + aspect + water + forest_pcov | Animal_id))
combo3 <- formula(case ~ -1 + tr_scale + dev + dev_sq + slope + dev:slope + 
                    aspect + forest_pcov + pas + log_sl + cos_ta + 
                    (1 | step_id) + (0 + tr_scale + dev + slope + aspect + forest_pcov + pas | Animal_id))
combo4 <- formula(case ~ -1 + tr_scale + dev + dev_sq + ag + ag_sq + slope + 
                    dev:slope + forest_pcov + log_sl + cos_ta + 
                    (1 | step_id) + (0 + tr_scale + dev + ag + slope + forest_pcov | Animal_id))
combo5 <- formula(case ~ -1 + tr_scale + dev + dev_sq + ag + ag_sq + slope + forest_pcov + log_sl + cos_ta + 
                    (1 | step_id) + (0 + tr_scale + dev + ag + slope + forest_pcov | Animal_id))
move1 <- formula(case ~ -1 + log_sl + cos_ta + (1 | step_id) + (0 + log_sl + cos_ta | Animal_id))

# Fit models
mig_mods <- list()
mig_mods[[1]] <- modFit(mig_ssf, anthro1)
mig_mods[[2]] <- modFit(mig_ssf, anthro2)
mig_mods[[3]] <- modFit(mig_ssf, tercli1)
mig_mods[[4]] <- modFit(mig_ssf, tercli2)
mig_mods[[5]] <- modFit(mig_ssf, lc1)
mig_mods[[6]] <- modFit(mig_ssf, combo1)
mig_mods[[7]] <- modFit(mig_ssf, combo2)
mig_mods[[8]] <- modFit(mig_ssf, combo3)
mig_mods[[9]] <- modFit(mig_ssf, combo4)
mig_mods[[10]] <- modFit(mig_ssf, move1)
mig_mods[[11]] <- modFit(mig_ssf, combo5)
names(mig_mods) <- c('anthro1','anthro2','tercli1','tercli2','lc1','combo1','combo2','combo3','combo4','move1','combo5')

# Save everything
save(mig_mods, file = "output/migration-ssf-mods-v2.rda")
# load(file = "output/migration-ssf-mods.rda")

# Create mod selection table and save top model only
mig_tab <- aictab(mig_mods)
mig_top <- mig_mods$combo4
save(list = c('mig_tab', 'mig_top'), file = "output/migration-ssf-top-model-v2.rda")
load(file = "output/migration-ssf-top-model-v2.rda", verbose = T)



