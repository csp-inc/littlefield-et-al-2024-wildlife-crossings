library(amt)
library(tidyverse)
library(raster)
library(sf)
library(glmmTMB)
library(car)
library(DHARMa)
library(AICcmodavg)

# Load RSF datasets
load(file = 'data/south-ute-deer-elk-data/elk-rsf-datsets-buff10km.rda', verbose = TRUE)

#------------------------
# Check covariate correlations
w_cor_tab<-winter %>% dplyr::select(chili, slope, aspect, PAS, PPT, Tave, Tmin, water,
                                    dev, forest, ag, well, traffic, rd_dist) %>% cor()
diag(w_cor_tab) <- NA
max(w_cor_tab, na.rm = T)

s_cor_tab<-summer %>% dplyr::select(chili, slope, aspect, PAS, PPT, Tave, water,
                                    dev, forest, ag, well, traffic, rd_dist) %>% cor()
diag(s_cor_tab) <- NA
max(s_cor_tab, na.rm = T)

#------------------------
# Define model fitting function
modFit <- function(data, model){
  mod <- glmmTMB(model, data = data, family = 'binomial', doFit = FALSE)
  mod$parameters$theta[1] <- log(1e3)
  ntheta = length(mod$parameters$theta)
  mod$mapArg <- list(theta=factor(c(NA, rep(1,ntheta-1))))
  mod_fit <- glmmTMB:::fitTMB(mod)
  return(mod_fit)
}

#------------------------
# Specify candidate models

# Terrain only
terr1 <- formula(used ~ chili + slope + aspect + (1 | nid) + (0 + chili + slope + aspect | nid))

# Climate only
s_clim1 <- formula(used ~ Tave + Tave_sq + (1 | nid) + (0 + Tave | nid))
s_clim2 <- formula(used ~ PPT + PPT_sq + (1 | nid) + (0 + PPT | nid))
w_clim1 <- formula(used ~ PAS+ PAS_sq + Tmin + Tmin_sq + (1 | nid) + (0 + PAS + Tmin | nid))
w_clim2 <- formula(used ~ PAS+ PAS_sq + (1 | nid) + (0 + PAS | nid))

# LC only
lc1 <- formula(used ~ water + forest + dev + dev_sq + ag + ag_sq + well + (1 | nid) + 
                   (0 + water + forest + dev + ag + well | nid))
lc2 <- formula(used ~ forest + dev + dev_sq + ag + ag_sq + (1 | nid) + 
                 (0 + forest + dev + ag | nid))

# Combo models - winter
w_combo1 <- formula(used ~ PAS + PAS_sq + Tmin + Tmin_sq + slope + chili + 
                      dev + dev_sq + ag + ag_sq + forest + (1 | nid) + 
                      (0 + PAS + Tmin + slope + chili + dev + ag + forest | nid))
w_combo2 <- formula(used ~ PAS + PAS_sq + slope + chili + 
                      dev + dev_sq + ag + ag_sq + forest + (1 | nid) + 
                      (0 + PAS + slope + chili + dev + ag + forest | nid))
w_combo3 <- formula(used ~ PAS + PAS_sq + Tmin + Tmin_sq + slope + 
                      dev + dev_sq + ag + ag_sq + (1 | nid) + 
                      (0 + PAS + Tmin + slope + dev + ag | nid))

# Combo models - summer
s_combo1 <- formula(used ~ Tave + Tave_sq + chili + slope + 
                      dev + dev_sq + ag + ag_sq + water + forest + (1 | nid) + 
                      (0 + Tave + chili + slope + dev + ag + forest + water | nid))
s_combo2 <- formula(used ~ PPT + PPT_sq + chili + slope + 
                      dev + dev_sq + ag + ag_sq + water + well + forest + (1 | nid) + 
                      (0 + PPT + chili + slope + dev + ag + forest + water + well | nid))
s_combo3 <- formula(used ~ Tave + Tave_sq + slope + 
                      dev + dev_sq + ag + ag_sq + water + (1 | nid) + 
                      (0 + Tave + slope + dev + ag + water | nid))
s_combo4 <- formula(used ~ PPT + PPT_sq + slope + 
                      dev + dev_sq + ag + ag_sq + water + (1 | nid) + 
                      (0 + PPT + slope + dev + ag + water | nid))

#-------------------------------------
#-------------------------------------
# WINTER MODEL SELECTION

# Fit model candidate set
winter_mods <- list()
winter_mods[[1]] <- modFit(winter, terr1)
winter_mods[[2]] <- modFit(winter, w_clim1)
winter_mods[[3]] <- modFit(winter, w_clim2)
winter_mods[[4]] <- modFit(winter, lc1)
winter_mods[[5]] <- modFit(winter, lc2)
winter_mods[[6]] <- modFit(winter, w_combo1)
winter_mods[[7]] <- modFit(winter, w_combo2)
winter_mods[[8]] <- modFit(winter, w_combo3)

# set list names
names(winter_mods) <- c('terr1', 'clim1', 'clim2', 'lc1', 'lc2', 'combo1', 'combo2', 'combo3')

# Save models
save(winter_mods, file = "output/winter-rsf-mods-v5-wbuff.rda")
# load(file = "output/winter-rsf-mods-v5-wbuff.rda", verbose = TRUE)

# Generate mod selection table
winter_tab <- aictab(winter_mods)

#-------------------------------------
#-------------------------------------
# SUMMER MODEL SELECTION

# Fit model candidate set
summer_mods <- list()
summer_mods[[1]] <- modFit(summer, terr1)
summer_mods[[2]] <- modFit(summer, s_clim1)
summer_mods[[3]] <- modFit(summer, s_clim2)
summer_mods[[4]] <- modFit(summer, lc1)
summer_mods[[5]] <- modFit(summer, lc2)
summer_mods[[6]] <- modFit(summer, s_combo1)
summer_mods[[7]] <- modFit(summer, s_combo2)
summer_mods[[8]] <- modFit(summer, s_combo3)
summer_mods[[9]] <- modFit(summer, s_combo4)

# set list names
names(summer_mods) <- c('terr1', 'clim1', 'clim2', 'lc1', 'lc2', 'combo1', 'combo2', 'combo3', 'combo4')

# Save models
save(summer_mods, file = "output/summer-rsf-mods-v5-wbuff.rda")
# load(file = "output/summer-rsf-mods-v5-wbuff.rda", verbose = TRUE)

# Generate mod selection table
summer_tab <- aictab(summer_mods)

#-------------------------------------
#-------------------------------------
# SAVE TOP MODELS ONLY

summer_top <- summer_mods$combo2
winter_top <- winter_mods$combo1
save(list = c('summer_tab','winter_tab','summer_top','winter_top'), file = "output/winter-summer-rsf-top-models.rda")
