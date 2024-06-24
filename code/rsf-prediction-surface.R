library(amt)
library(tidyverse)
# library(raster)
library(terra)
library(sf)
library(glmmTMB)
library(car)
library(AICcmodavg)

#-------------------------------------
#-------------------------------------
# DATA PREP

# Define elk AOI and CRS
su_crs <- st_crs(26913)
aoi <- st_read("data/south-ute-deer-elk-data/elk-mcp/elk-mcp.shp") %>% st_buffer(dist = 120000)
aoi50 <- st_read("data/south-ute-deer-elk-data/elk-mcp/elk-mcp.shp") %>% st_buffer(dist = 50000)

#-------------
# GET MODELS AND DATASETS
# Get seasonal RSF models
load(file = "output/winter-summer-rsf-top-models.rda", verbose = TRUE)

# Get seasonal RSF datasets (for calculating cov means and SDs)
load(file = 'data/south-ute-deer-elk-data/elk-rsf-datsets-buff10km.rda', verbose = TRUE)

#-------------
# PREP PRESENT DAY COVS
# Get Covariates
covs90 <- rast('data/elk-covs-present-90m.tiff')
covs_course <- rast('data/elk-covs-lc-road-present-450m.tiff') # all "distance to" covs and traffic volume
# Change resolution of course covariates to match other covs
rd_90m <- disagg(covs_course$road_dist, fact = 5) %>% crop(covs90$elevation_1km)
dev_90m <- disagg(covs_course$dev_dist, fact = 5) %>% crop(covs90$elevation_1km)
ag_90m <- disagg(covs_course$ag_dist, fact = 5) %>% crop(covs90$elevation_1km)
water_90m <- disagg(covs_course$water_dist, fact = 5) %>% crop(covs90$elevation_1km)
well_90m <- disagg(covs_course$well_dist, fact = 5) %>% crop(covs90$elevation_1km)
lc_covs_90m <- c(rd_90m, dev_90m, ag_90m, water_90m, well_90m)

# Winter covariates
wt_covs <- covs90[[which(str_detect(names(covs90), 
                                  paste(c('_sm', '_sp', '_at'), collapse = "|"), 
                                  negate = TRUE) & 
                         str_detect(names(covs90), "_300m", negate = TRUE))]]
wt_covs <- c(wt_covs, lc_covs_90m) # attach lc covariates

# Summer covariates
sm_covs <- covs90[[which(str_detect(names(covs90), 
                                  paste(c('_wt', '_sp', '_at'), collapse = "|"), 
                                  negate = TRUE) & 
                         str_detect(names(covs90), "_300m", negate = TRUE))]]
sm_covs <- c(sm_covs, lc_covs_90m) # attach lc covariates

#-------------
# PREP FUTURE COVS
# Get Covariates
covs90f <- rast('data/elk-covs-future-90m.tiff')
covs_coursef <- rast('data/elk-covs-lc-road-future-450m.tiff') # all "distance to" covs and traffic volume
# Change resolution of course covariates to match other covs
dev_90f <- disagg(covs_coursef$dev_dist, fact = 5) %>% crop(covs90f$PAS_1km)
ag_90f <- disagg(covs_coursef$ag_dist, fact = 5) %>% crop(covs90f$PAS_1km)
# Compile future covs with some from present day that are unchanged in 2050
extra_covs_90f <- c(dev_90f, ag_90f, water_90m, well_90m, rd_90m, covs90$slope_1km, 
                    covs90$aspect_1km, covs90$CHILI_1km, covs90$forest_pcov_1km)

# Winter covariates
wt_covs_fut <- covs90f[[which(str_detect(names(covs90f), 
                                    paste(c('_sm', '_sp', '_at'), collapse = "|"), 
                                    negate = TRUE) & 
                           str_detect(names(covs90f), "_300m", negate = TRUE))]]
wt_covs_fut <- c(wt_covs_fut, extra_covs_90f) # attach extra covariates

# Summer covariates
sm_covs_fut <- covs90f[[which(str_detect(names(covs90f), 
                                    paste(c('_wt', '_sp', '_at'), collapse = "|"), 
                                    negate = TRUE) & 
                           str_detect(names(covs90f), "_300m", negate = TRUE))]]
sm_covs_fut <- c(sm_covs_fut, extra_covs_90f) # attach lc covariates


#-------------------------------------
#-------------------------------------
# PREDICTION AND RASTERIZATION FUNCTIONS
pConstrain <- function(x, bounds){
  x <- ifelse(x < bounds[1], bounds[1], x)
  x <- ifelse(x > bounds[2], bounds[2], x)
  return(x)
}

# Function to extract covariate values at all cell locations in study area
gridCovs <- function(covs, st_area, elk_dat, season){
  
  # Get cell numbers for raster cells overlapping study area
  pred_cell_nums <- terra::cells(covs, vect(st_area))[,2]
  # Extract covariate values from across the covs raster for use in pred surface
  pred_cov_temp <- terra::extract(covs, pred_cell_nums) 
  pred_xy <- terra::xyFromCell(covs, pred_cell_nums) %>% 
    as.data.frame() %>% 
    rename(utm_x = x, utm_y = y) 
  pred_covs <- cbind(pred_cov_temp, pred_xy) %>% 
    drop_na()
  
  # Rescale prediction covs based on range of values in rsf dataset
  covScale <- function(pc = pred_covs, dat = as.data.frame(elk_dat), cname){
    return((pc[,cname] - mean(dat[,cname], na.rm = TRUE))/sd(dat[,cname], na.rm = TRUE))
  }
  if(season == 'winter'){
    pred_covs_scale <- data.frame(chili = covScale(cname = "CHILI_1km"),
                                  slope = covScale(cname = "slope_1km"),
                                  aspect = covScale(cname = "aspect_1km"),
                                  PAS = covScale(cname = "PAS_1km"),
                                  PPT = covScale(cname = "PPT_wt_1km"),
                                  Tave = covScale(cname = "Tave_wt_1km"),
                                  Tmin = covScale(cname = "Tmin_annual_1km"),
                                  water = covScale(cname = "water_dist"),
                                  dev = covScale(cname = "dev_dist"),
                                  forest = covScale(cname = "forest_pcov_1km"),
                                  ag = covScale(cname = "ag_dist"),
                                  well = covScale(cname = "well_dist"),
                                  rd_dist = covScale(cname = "road_dist")) %>% 
      mutate(PPT_sq = PPT^2, 
             Tave_sq = Tave^2,
             Tmin_sq = Tmin^2,
             PAS_sq = PAS^2,
             dev_sq = dev^2,
             forest_sq = forest^2,
             ag_sq = ag^2,
             utm_x = pred_covs$utm_x,
             utm_y = pred_covs$utm_y)
  }
  
  if(season == 'summer') {
    pred_covs_scale <- data.frame(chili = covScale(cname = "CHILI_1km"),
                                  slope = covScale(cname = "slope_1km"),
                                  aspect = covScale(cname = "aspect_1km"),
                                  PAS = covScale(cname = "PAS_1km"),
                                  PPT = covScale(cname = "PPT_sm_1km"),
                                  Tave = covScale(cname = "Tave_sm_1km"),
                                  water = covScale(cname = "water_dist"),
                                  dev = covScale(cname = "dev_dist"),
                                  forest = covScale(cname = "forest_pcov_1km"),
                                  ag = covScale(cname = "ag_dist"),
                                  well = covScale(cname = "well_dist"),
                                  rd_dist = covScale(cname = "road_dist")) %>% 
      mutate(PPT_sq = PPT^2, 
             Tave_sq = Tave^2,
             dev_sq = dev^2,
             forest_sq = forest^2,
             ag_sq = ag^2,
             utm_x = pred_covs$utm_x,
             utm_y = pred_covs$utm_y)
  }
    
  # Return scaled prediction covs
  return(pred_covs_scale)
}

# Function to predict from model equation at all cell locations in study area
gridPredict <- function(gcovs, mod, scale = "exp", keepInt = FALSE, constrain = NULL){
  if(keepInt) coefs <- fixef(mod)$cond else coefs <- fixef(mod)$cond[-1] # get model coefficients, removing the intercept if desired
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
gridRast <- function(dat, filename){
  # Convert prediction dataset to sf and then terra SpatVector
  pred_sp <- dat %>% 
    st_as_sf(coords = c("utm_x", "utm_y"), crs = su_crs) 
  
  # Create template raster
  r_temp <- rast()
  terra::ext(r_temp) <- terra::ext(pred_sp)
  terra::res(r_temp) <- 90
  terra::crs(r_temp) <- su_crs$wkt
  
  # Rasterize
  pred_sp <- vect(pred_sp)
  pred_rast <- terra::rasterize(pred_sp, r_temp, field = 'prediction', 
                                filename = filename, overwrite = TRUE)
  return(pred_rast)
}

#-------------------------------------
#-------------------------------------
# MAKE PREDICTION RASTERS

# Winter-present day prediction surface
wt_pres_covs = gridCovs(wt_covs, aoi50, winter, season = 'winter')
wt_pres_pred <- gridPredict(wt_pres_covs, winter_top, scale = "exp")
wt_pres_rast <- gridRast(wt_pres_pred, filename = "output/wt-present-exp-aoi50.tif")

# Summer-present day prediction surface
sm_pres_covs = gridCovs(sm_covs, aoi50, summer, season = 'summer')
sm_pres_pred <- gridPredict(sm_pres_covs, summer_top, scale = "exp", constrain = c(-Inf,6))
# sm_pres_pred <- gridPredict(sm_pres_covs, summer_top, scale = "logit")
sm_pres_rast <- gridRast(sm_pres_pred, filename = "output/sm-present-exp-aoi50.tif")

#---------
# Winter-future prediction surface
wt_fut_covs = gridCovs(wt_covs_fut, aoi50, winter, season = 'winter')
wt_fut_pred <- gridPredict(wt_fut_covs, winter_top, scale = "exp")
wt_fut_rast <- gridRast(wt_fut_pred, filename = "output/wt-future-exp-aoi50.tif")

# Summer-future prediction surface
sm_fut_covs = gridCovs(sm_covs_fut, aoi50, summer, season = 'summer')
sm_fut_pred <- gridPredict(sm_fut_covs, summer_top, scale = "exp", constrain = c(-Inf,6))
sm_fut_rast <- gridRast(sm_fut_pred, filename = "output/sm-future-exp-aoi50.tif")

