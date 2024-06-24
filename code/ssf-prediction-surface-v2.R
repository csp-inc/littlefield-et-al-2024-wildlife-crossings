library(amt)
library(tidyverse)
library(terra)
library(sf)
library(glmmTMB)

# Load top SSF mod and dataset
load(file = "output/migration-ssf-top-model-v2.rda", verbose = T)
load(file = "data/south-ute-deer-elk-data/elk-migration-ssf-data.rda", verbose = T)

# Define elk AOI and CRS
su_crs <- st_crs(26913)
aoi50 <- st_read("data/south-ute-deer-elk-data/elk-mcp/elk-mcp.shp") %>% st_buffer(dist = 50000)

#----------------------
# PREP PRESENT DAY COVS
covs90 <- rast('data/elk-covs-present-90m.tiff') 
covs_course <- rast('data/elk-covs-lc-road-present-450m.tiff')
# Change resolution of course covariates to match other covs
dev_90m <- disagg(covs_course$dev_dist, fact = 5) %>% crop(covs90$elevation_300m)
ag_90m <- disagg(covs_course$ag_dist, fact = 5) %>% crop(covs90$elevation_300m)
water_90m <- disagg(covs_course$water_dist, fact = 5) %>% crop(covs90$elevation_300m)
well_90m <- disagg(covs_course$well_dist, fact = 5) %>% crop(covs90$elevation_300m)
lc_covs <- c(dev_90m, ag_90m, water_90m, well_90m)
# Combine covs
mig_covs_pres <- covs90[[c("slope_300m","aspect_300m","CHILI_300m","forest_pcov_300m",
                      "Tmax_annual_300m", "Tmin_annual_300m", "CMD_300m", "PAS_300m")]]
mig_covs_pres <- c(mig_covs_pres, lc_covs)

#----------------------
# PREP FUTURE COVS
# Get Covariates
covs90f <- rast('data/elk-covs-future-90m.tiff')
covs_coursef <- rast('data/elk-covs-lc-road-future-450m.tiff')
# Change resolution of course covariates to match other covs
dev_90f <- disagg(covs_coursef$dev_dist, fact = 5) %>% crop(covs90f$PAS_1km)
ag_90f <- disagg(covs_coursef$ag_dist, fact = 5) %>% crop(covs90f$PAS_1km)
fc1 <- c(dev_90f, ag_90f)
fc2 <- covs90f[[c("forest_pcov_300m","Tmax_annual_300m","Tmin_annual_300m","CMD_300m", "PAS_300m")]]
pc1 <- mig_covs_pres[[c("slope_300m", "aspect_300m", "CHILI_300m", "water_dist", "well_dist")]]
mig_covs_fut <- c(pc1, fc1, fc2)

#----------------------
# PREP ROAD RELATED COVS
rd_width = 90
rds_union <- st_read("data/co-nm-roads-traffic/co-nm-addt-w-future.shp") %>% 
  st_buffer(rd_width) %>% st_simplify(dTolerance = 10) %>% st_union()
rds <- st_read("data/co-nm-roads-traffic/co-nm-addt-w-future.shp") %>% st_buffer(rd_width) %>% vect()
xing_exist <- st_read("data/crossing-structures/crossings-final/crossings-exist-recc-final-20230311.shp") %>% 
  filter(existRecc=='exist') %>% st_buffer(1500) %>% st_intersection(rds_union) %>% vect()
xing_all <- st_read("data/crossing-structures/crossings-final/crossings-exist-recc-final-20230311.shp") %>% 
  st_buffer(1500) %>% st_intersection(rds_union) %>% vect()
# Rasterize road covariates
traffic_pres <- terra::rasterize(rds, mig_covs_pres[[1]], field = "AADTpres", background = 0, touches = TRUE)
traffic_fut <- terra::rasterize(rds, mig_covs_pres[[1]], field = "AADT2050", background = 0, touches = TRUE)
no_traffic <- rast(traffic_pres, val = 0)

# Rename and add to cov raster stack
names(traffic_pres) <- "traffic_pres"
names(traffic_fut) <- "traffic_fut"
names(no_traffic) <- "traffic_none"
mig_covs_pres <- c(mig_covs_pres, traffic_pres, traffic_fut, no_traffic)
mig_covs_fut <- c(mig_covs_fut, traffic_pres, traffic_fut, no_traffic)

#-------------------------------------
#-------------------------------------
# PREDICTION AND RASTERIZATION FUNCTIONS
pConstrain <- function(x, bounds){
  x <- ifelse(x < bounds[1], bounds[1], x)
  x <- ifelse(x > bounds[2], bounds[2], x)
  return(x)
}

# Function to extract covariate values at all cell locations in study area
gridCovs <- function(covs, st_area, elk_dat, traffic = c("pres","fut","none")){
  
  # Get cell numbers for raster cells overlapping study area
  pred_cell_nums <- terra::cells(covs, vect(st_area))[,2]
  # Extract covariate values from across the covs raster for use in pred surface
  pred_cov_temp <- terra::extract(covs, pred_cell_nums) 
  pred_xy <- terra::xyFromCell(covs, pred_cell_nums) %>% 
    as.data.frame() %>% 
    rename(utm_x = x, utm_y = y) 
  pred_covs <- cbind(pred_cov_temp, pred_xy) %>% 
    drop_na()
  
  # Define traffic covariate layer to use
  if(traffic == 'pres') names(pred_covs)[which(names(pred_covs)=='traffic_pres')] <- 'traffic'  
  if(traffic == 'fut') names(pred_covs)[which(names(pred_covs)=='traffic_fut')] <- 'traffic'
  if(traffic == 'none') names(pred_covs)[which(names(pred_covs)=='traffic_none')] <- 'traffic'
  
  # Rescale prediction covs based on range of values in rsf dataset
  covScale <- function(pc = pred_covs, dat = as.data.frame(elk_dat), cname){
    return((pc[,cname] - mean(dat[,cname], na.rm = TRUE))/sd(dat[,cname], na.rm = TRUE))
  }
  
  pred_covs_scale <- data.frame(chili = covScale(cname = "CHILI_300m"),
                                slope = covScale(cname = "slope_300m"),
                                aspect = covScale(cname = "aspect_300m"),
                                pas = covScale(cname = "PAS_300m"),
                                tmax = covScale(cname = "Tmax_annual_300m"),
                                tmin = covScale(cname = "Tmin_annual_300m"),
                                water = covScale(cname = "water_dist"),
                                dev = covScale(cname = "dev_dist"),
                                forest_pcov = covScale(cname = "forest_pcov_300m"),
                                ag = covScale(cname = "ag_dist"),
                                well = covScale(cname = "well_dist"),
                                tr_scale = covScale(cname = "traffic"),
                                log_sl = mean(elk_dat$log_sl, na.rm = TRUE),
                                cos_ta = mean(elk_dat$cos_ta, na.rm = TRUE)) %>% 
    mutate(slope_sq = slope^2, 
           dev_sq = dev^2,
           ag_sq = ag^2,
           "dev:slope" = dev * slope,
           utm_x = pred_covs$utm_x,
           utm_y = pred_covs$utm_y)

  # Return scaled prediction covs
  return(pred_covs_scale)
}

# Function to predict from model equation at all cell locations in study area
gridPredict <- function(gcovs, mod, scale = "exp", constrain = NULL){
  coefs <- fixef(mod)$cond # get model coefficients
  gcovs_mat <- gcovs %>% dplyr::select(names(coefs)) %>% data.matrix() # grab applicable columns from grid covs
  # Create linear prediction (and transform, as specified in function call)
  pred <- (gcovs_mat %*% coefs)
  if(!is.null(constrain)) pred <- pConstrain(pred, constrain)
  if(scale == "logit") pred <- pred
  if(scale == "exp") pred <- exp(pred)
  if(scale == "response") pred <- exp(pred)/(1+exp(pred))
  pred_out <- data.frame(prediction = pred, utm_x = gcovs$utm_x, utm_y = gcovs$utm_y)
  return(pred_out)
}

# Function to rasterize predictions
gridRast <- function(dat, filename=NULL){
  # Convert prediction dataset to sf and then terra SpatVector
  pred_sp <- dat %>% 
    st_as_sf(coords = c("utm_x", "utm_y"), crs = su_crs) 
  
  # # Create template raster
  r_temp <- rast()
  terra::ext(r_temp) <- terra::ext(pred_sp)
  terra::res(r_temp) <- 90
  terra::crs(r_temp) <- su_crs$wkt
  
  # Rasterize
  pred_sp <- vect(pred_sp)
  
  if(is.null(filename)) {
    pred_rast <- terra::rasterize(pred_sp, mig_covs_pres[[1]], field = 'prediction') %>% 
    crop(r_temp)
  } else {
    pred_rast <- terra::rasterize(pred_sp, mig_covs_pres[[1]], field = 'prediction', 
                                  filename = filename, overwrite = TRUE) %>% crop(r_temp)
    } 
  
  return(pred_rast)
}

#-------------------
# Resistance surface functions
# Define resistance function (from Keeley et al. 2016 Landscape Ecol)
rscale <- function(h, c){
  out <- 100-99*((1-exp(-c*h))/(1-exp(-c)))
  return(out)
}

resist_depricated <- function(hs, c = NULL, cross = c('pres', 'fut'), cross_mod){

  # Rescale btwn 0 and 1
  mm <- minmax(hs)
  hs01 <- (hs - mm[1])/(mm[2] - mm[1])

  # Either rescale linearly or using neg exp
  if(is.null(c)) res_temp <- (1-hs01)*100
  if(!is.null(c)) res_temp <- rscale(hs01, c)

  # Specify present vs future crossings to use and modify resistance
  if(cross == "pres") ctemp <- xing_pres else ctemp <- xing_fut
  crast <- terra::rasterize(ctemp, res_temp, field = cross_mod, background = 1, touches = TRUE)
  res_out <- res_temp * crast

  return(res_out)
}

resist <- function(hs, c = NULL, lc = c('pres', 'fut'), cross = c('exist', 'all', 'none'), cross_mod){
  
  # Burn in crossing, setting them cross_mod more resistance than (roadless) background
  if(cross == "exist"){
    if(lc == "pres") ctemp <- cross_pres_exist * cross_mod else ctemp <- cross_fut_exist * cross_mod
    hs[which.lyr(ctemp)]<-ctemp
  } 
  if(cross == "all"){
    if(lc == "pres") ctemp <- cross_pres_all * cross_mod else ctemp <- cross_fut_all * cross_mod
    hs[which.lyr(ctemp)]<-ctemp
  }
  if(cross == "none") hs <- hs
  
  # Rescale btwn 0 and 1
  mm <- minmax(hs)
  hs01 <- (hs - mm[1])/(mm[2] - mm[1])
  
  # Either rescale linearly or using neg exp
  if(is.null(c)) {
    res_out <- (1-hs01)*100
    res_out <- res_out %>% clamp(lower = 1, values = TRUE)
  }  
  if(!is.null(c)) res_out <- rscale(hs01, c)
  
  return(res_out)
}

#-------------------------------------
#-------------------------------------
# CROSSING RESISTANCE MODIFIERS

# Create "no traffic" surface for use in assigning crossing resistance values
no_tr_pres_covs <- gridCovs(mig_covs_pres, aoi50, mig_ssf, traffic = 'none')
no_tr_pres_pred <- gridPredict(no_tr_pres_covs, mig_top)
no_tr_pres_rast <- gridRast(no_tr_pres_pred, filename = "output/mig-noT-pres-selection.tif")

no_tr_fut_covs <- gridCovs(mig_covs_fut, aoi50, mig_ssf, traffic = 'none')
no_tr_fut_pred <- gridPredict(no_tr_fut_covs, mig_top)
no_tr_fut_rast <- gridRast(no_tr_fut_pred, filename = "output/mig-noT-fut-selection.tif")

cross_pres_exist <- mask(no_tr_pres_rast, xing_exist) # Present day land cover, existing crossings
cross_pres_all <- mask(no_tr_pres_rast, xing_all) # Present day land cover, all crossings

cross_fut_exist <- mask(no_tr_fut_rast, xing_exist) # Future land cover, existing crossings
cross_fut_all <- mask(no_tr_fut_rast, xing_all) # Future land cover, all crossings

#-------------------------------------
#-------------------------------------
# MAKE PREDICTION RASTERS AND RESISTANCE SURFACES

# Present and future habitat suitability rasters
pres_covs <- gridCovs(mig_covs_pres, aoi50, mig_ssf, traffic = 'pres')
pres_pred <- gridPredict(pres_covs, mig_top)
pres_rast <- gridRast(pres_pred, filename = "output/mig_pres_selection.tif")

fut_covs <- gridCovs(mig_covs_fut, aoi50, mig_ssf, traffic = 'fut')
fut_pred <- gridPredict(fut_covs, mig_top)
fut_rast <- gridRast(fut_pred, filename = "output/mig_fut_selection.tif")

# Resistance surfaces for all scenarios w/ neg exp c = 4
res_p_nc_ne4 <- resist(pres_rast, c = 4, lc = "pres", cross = 'none', cross_mod = 0.75) # present, no crossings
res_p_ec_ne4 <- resist(pres_rast, c = 4, lc = "pres", cross = 'exist', cross_mod = 0.75) # present, existing crossings
res_p_ac_ne4 <- resist(pres_rast, c = 4, lc = "pres", cross = 'all', cross_mod = 0.75) # present, all crossings

res_f_nc_ne4 <- resist(fut_rast, c = 4, lc = "fut", cross = 'none', cross_mod = 0.75) # future, no crossings
res_f_ec_ne4 <- resist(fut_rast, c = 4, lc = "fut", cross = 'exist', cross_mod = 0.75) # future, existing crossings
res_f_ac_ne4 <- resist(fut_rast, c = 4, lc = "fut", cross = 'all', cross_mod = 0.75) # future, all crossings


# Resistance surfaces for all scenarios w/ neg exp c = 8
res_p_nc_ne8 <- resist(pres_rast, c = 8, lc = "pres", cross = 'none', cross_mod = 0.75) # present, no crossings
res_p_ec_ne8 <- resist(pres_rast, c = 8, lc = "pres", cross = 'exist', cross_mod = 0.75) # present, existing crossings
res_p_ac_ne8 <- resist(pres_rast, c = 8, lc = "pres", cross = 'all', cross_mod = 0.75) # present, all crossings

res_f_nc_ne8 <- resist(fut_rast, c = 8, lc = "fut", cross = 'none', cross_mod = 0.75) # future, no crossings
res_f_ec_ne8 <- resist(fut_rast, c = 8, lc = "fut", cross = 'exist', cross_mod = 0.75) # future, existing crossings
res_f_ac_ne8 <- resist(fut_rast, c = 8, lc = "fut", cross = 'all', cross_mod = 0.75) # future, all crossings

# Save 'em all
terra::writeRaster(res_p_nc_ne4,filename = "output/connectivity/resistance/pres-no-cross-ne4.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)
terra::writeRaster(res_p_ec_ne4,filename = "output/connectivity/resistance/pres-exist-cross-ne4.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)
terra::writeRaster(res_p_ac_ne4,filename = "output/connectivity/resistance/pres-all-cross-ne4.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)

terra::writeRaster(res_f_nc_ne4,filename = "output/connectivity/resistance/fut-no-cross-ne4.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)
terra::writeRaster(res_f_ec_ne4,filename = "output/connectivity/resistance/fut-exist-cross-ne4.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)
terra::writeRaster(res_f_ac_ne4,filename = "output/connectivity/resistance/fut-all-cross-ne4.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)

terra::writeRaster(res_p_nc_ne8,filename = "output/connectivity/resistance/pres-no-cross-ne8.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)
terra::writeRaster(res_p_ec_ne8,filename = "output/connectivity/resistance/pres-exist-cross-ne8.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)
terra::writeRaster(res_p_ac_ne8,filename = "output/connectivity/resistance/pres-all-cross-ne8.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)

terra::writeRaster(res_f_nc_ne8,filename = "output/connectivity/resistance/fut-no-cross-ne8.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)
terra::writeRaster(res_f_ec_ne8,filename = "output/connectivity/resistance/fut-exist-cross-ne8.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)
terra::writeRaster(res_f_ac_ne8,filename = "output/connectivity/resistance/fut-all-cross-ne8.asc", 
                   filetype = "AAIGrid", overwrite = TRUE)

