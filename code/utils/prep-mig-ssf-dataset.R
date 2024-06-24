library(amt)
library(tidyverse)
library(raster)
library(sf)

# Load elk dataset
load("data/south-ute-deer-elk-data/su-elk-final-w-season-20220615.rda")
su_crs <- st_crs(26913)

# Get Migration dataset
mig <- elkFinal %>% filter(season == 'migration', mig_quality == 1)
ids <- unique(mig$Animal_id)

# Get raster covariates
covs90 <- stack('data/elk-covs-present-90m.tiff')
covs_course <- stack('data/elk-covs-lc-road-present-450m.tiff') # all "distance to" covs and traffic volume
# Subset to just covs of interest
covs_fine <- covs90[[c("elevation_300m","slope_300m","aspect_300m","CHILI_300m","water_pcov_300m",
                      "dev_pcov_300m","forest_pcov_300m","ag_pcov_300m","well_dens_300m",
                      "Tmax_annual_300m", "Tmin_annual_300m", "CMD_300m", "PAS_300m")]]
covs_course <- covs_course[[c("dev_dist","ag_dist","forest_dist","water_dist","well_dist")]]

# Get road stuff
rds <- st_read("data/co-nm-roads-traffic/co-nm-addt-w-future.shp")
xing <- st_read("data/crossing-structures/crossings_existRecc_aoi_tidy.shp") %>% 
  filter(type=='exist')
xing_100 <- xing %>% st_buffer(100)
xing_500 <- xing %>% st_buffer(500)
# Function for extracting ids of roads crossed by a migration step
rd_id <- function(x){
  t <- which(x)
  if(length(t)==0) t <- NA
  if(length(t)>1) t <- t[1]
  return(t)
}

#-------------------------------------------------
# CREATE SSF DATA SET
# Loop through each elk and create used and available steps,
# extract raster covs, convert to line dataset, and extract
# road covs

ssf_ua <- list()
for(i in 1:length(ids)){
  print(paste("running elk #", i, "of", length(ids)))
  
  # Subset to focal animal and add "previous step length" to remove short (likely stationary) steps
  temp <- mig %>% filter(Animal_id == ids[i])
  temp$pstep <- c(NA, temp$step)[1:nrow(temp)]
  temp <- temp[which(is.na(temp$pstep) | temp$pstep > 20),]
  
  # Create used and available steps and extract raster covs at end points of each
  ua_temp <- temp  %>% make_track(x, y, date, id = Animal_id, sex = Sex, season = season, crs = 26913) %>% 
    track_resample(rate = hours(2), tolerance = minutes(30)) %>% 
    steps_by_burst() %>% 
    random_steps() %>% 
    extract_covariates(covs_fine, where = 'end') %>% 
    extract_covariates(covs_course, where = 'end') %>% 
    mutate(cos_ta = cos(ta_),
           log_sl = log(sl_),
           Animal_id = ids[i],
           Sex = temp$Sex[1])
  
  # Convert to line dataset
  pts <- ua_temp %>% dplyr::select(x1_, y1_, x2_, y2_) # get all start and end points
  lns <- apply(pts, 1, function(x) # Convert to lines
  {
    v <- as.numeric(x[c(1,3,2,4)])
    m <- matrix(v, nrow = 2)
    return(st_sfc(st_linestring(m), crs = 26913))
  })
  lns <- Reduce(c, lns)
  # Compile other data columns
  ldat <- ua_temp %>% mutate(burst = paste(Animal_id, burst_, sep = "-"),
                             step_id = paste(Animal_id, step_id_, sep = "-"),
                             case = as.numeric(case_)) %>% 
    dplyr::select(Animal_id, Sex, burst, step_id, case, 6:10, 13:32) 
  # Combine into line dataframe
  ua <- st_sf(ldat, geometry = lns)
  
  # Get road crossing (and traffic volume if crossed) 
  rdx <- st_intersects(rds, ua, sparse = F) %>% apply(2, rd_id) %>% unlist()
  traffic <- rds$AADTpres[rdx]
  traffic[which(is.na(traffic))] <- 0
  ua$rd_cross <- ifelse(is.na(rdx), 0, 1)
  ua$traffic <- traffic
  
  # Get crossing structure based on range of possible sizes
  cs_100 <- st_intersects(xing_100, ua, sparse = F) %>% apply(2, rd_id) %>% unlist()
  cs_500 <- st_intersects(xing_500, ua, sparse = F) %>% apply(2, rd_id) %>% unlist()
  ua$cross100 <- ifelse(is.na(cs_100), 0, 1)
  ua$cross500 <- ifelse(is.na(cs_500), 0, 1)
  
  # Add to list
  ssf_ua[[i]] <- ua
}
# Bind em
elk_ssf_dat <- do.call(rbind, ssf_ua)

save(elk_ssf_dat, file = "data/south-ute-deer-elk-data/elk-ssf-temp.rda")
load(file = "data/south-ute-deer-elk-data/elk-ssf-temp.rda", verbose = TRUE)

#-------------------------------
# Scale/transform some variables and test correlations
mig_ssf <- elk_ssf_dat %>% 
  mutate(elevation = scale(elevation_300m),
         slope = scale(slope_300m),
         slope_sq = slope^2,
         aspect = scale(aspect_300m),
         chili = scale(CHILI_300m),
         water = scale(water_dist),
         dev = scale(dev_dist),
         dev_sq = dev^2,
         ag = scale(ag_dist),
         ag_sq = ag^2,
         forest_pcov = scale(forest_pcov_300m),
         well = scale(well_dist),
         tmax = scale(Tmax_annual_300m),
         tmin = scale(Tmin_annual_300m),
         cmd = scale(CMD_300m),
         pas = scale(PAS_300m),
         tr_scale = scale(traffic))

save(mig_ssf, file = "data/south-ute-deer-elk-data/elk-migration-ssf-data.rda")
