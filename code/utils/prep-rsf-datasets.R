library(amt)
library(tidyverse)
library(raster)
library(sf)

#-------------------------------------
#-------------------------------------
# DATA PREP

# Define elk AOI and CRS
su_crs <- st_crs(26913)
aoi <- st_read("data/south-ute-deer-elk-data/elk-mcp/elk-mcp.shp") %>% st_buffer(dist = 120000)

# Get Covariates
covs90 <- stack('data/elk-covs-present-90m.tiff')
covs_course <- stack('data/elk-covs-lc-road-present-450m.tiff') # all "distance to" covs and traffic volume
# Summer/winter covariates
covs_fine <- covs90[[which(str_detect(names(covs90), 
                                      paste(c('_sp', '_at'), collapse = "|"), 
                                      negate = TRUE) & 
                             str_detect(names(covs90), "_300m", negate = TRUE))]]

# Load elk dataset
load("data/south-ute-deer-elk-data/su-elk-final-w-season-20220615.rda")

# Elk IDs to loop through
ids = unique(elkFinal$Animal_id)

#-------------------------------------
#-------------------------------------
# PREP USED AVAILABLE DATASETS
# Standardized time between used locations, create available locations, 
# and extract covariates for each elk
elk_ua <- list()
for(i in 1:length(ids)){
  print(paste("running elk #", i, "of", length(ids)))
  
  # organize and extract covs for used points
  used <- elkFinal %>% filter(Animal_id %in% ids[i], mig_quality == 1, season != "migration") 
  if(nrow(used)==0) next # skip animals where mig quality !=1
  used <- used %>% make_track(x, y, date, id = Animal_id, sex = Sex, season = season, crs = 26913) %>%
    track_resample(rate = hours(8), tolerance = minutes(30)) %>% 
    extract_covariates(covs_course) %>% 
    extract_covariates(covs_fine) %>% 
    mutate(used = 1) %>% 
    dplyr::select(-"t_", -"burst_") 
  
  # Make buffered home range
  hr <- hr_mcp(used) %>% hr_isopleths() %>% st_buffer(10000)
  
  # create random available points
  avail <- random_points(hr, n = nrow(used)*10, presence = used) %>% 
    filter(case_ == FALSE) %>% 
    extract_covariates(covs_course) %>%
    extract_covariates(covs_fine) %>%
    mutate(id = ids[i],
           sex = used$sex[1],
           season = 'rand',
           used = 0) %>% 
    dplyr::select(names(used))
  
  temp_out <- rbind(used, avail)
  elk_ua[[i]] <- temp_out
}
# Remove empty elements
elk_ua<-elk_ua[lapply(elk_ua,length)>0]
# Bind em
elk_rsf_dat <- do.call(rbind, elk_ua)

#-------------------------------------
#-------------------------------------
# PREP SEASONAL DATASETS

#--------------
# Create 1:10 ratio for used:available
eids <- unique(elk_rsf_dat$id)
wlist <- list()
slist <- list()
for(i in 1:length(eids)){
  # Subset to just animal of interest
  temp <- elk_rsf_dat %>% filter(id %in% eids[i])
  # Get available points
  a_temp <- temp %>% filter(season == 'rand')
  
  # Get 1:10 used:available ratio for winter
  uw <- temp %>% filter(season == 'winter')
  aw <- a_temp[sample(1:nrow(a_temp), nrow(uw)*10),]
  w_out <- rbind(uw, aw)
  wlist[[i]] <- w_out
  
  # Get 1:10 used:available ratio for summer
  us <- temp %>% filter(season == 'summer')
  as <- a_temp[sample(1:nrow(a_temp), nrow(us)*10),]
  s_out <- rbind(us, as)
  slist[[i]] <- s_out
}
winter <- do.call(rbind, wlist)
summer <- do.call(rbind, slist)

#--------------
# Create scaled versions of subset of vars
winter <- winter %>% 
  mutate(chili = scale(CHILI_1km),
         slope = scale(slope_1km),
         aspect = scale(aspect_1km),
         PAS = scale(PAS_1km),
         PAS_sq = PAS^2,
         PPT = scale(PPT_wt_1km),
         PPT_sq = PPT^2,
         Tave = scale(Tave_wt_1km),
         Tave_sq = Tave^2,
         Tmin = scale(Tmin_annual_1km),
         Tmin_sq = Tmin^2,
         water = scale(water_dist),
         dev = scale(dev_dist),
         dev_sq = dev^2,
         forest = scale(forest_pcov_1km),
         forest_sq = forest^2,
         ag = scale(ag_dist),
         ag_sq = ag^2,
         well = scale(well_dist),
         traffic = scale(rd_traffic_smooth),
         rd_dist = scale(road_dist)) %>% 
  mutate(nid = as.numeric(factor(id)),
         weight = 100000 * (1-used))
  
summer <- summer %>%
  mutate(chili = scale(CHILI_1km),
         slope = scale(slope_1km),
         aspect = scale(aspect_1km),
         PAS = scale(PAS_1km),
         PPT = scale(PPT_sm_1km),
         PPT_sq = PPT^2,
         Tave = scale(Tave_sm_1km),
         Tave_sq = Tave^2,
         water = scale(water_dist),
         dev = scale(dev_dist),
         dev_sq = dev^2,
         forest = scale(forest_pcov_1km),
         forest_sq = forest^2,
         ag = scale(ag_dist),
         ag_sq = ag^2,
         well = scale(well_dist),
         traffic = scale(rd_traffic_smooth),
         rd_dist = scale(road_dist)) %>% 
  mutate(nid = as.numeric(factor(id)),
         weight = 100000 * (1-used))


# Save seasonal datasets
save(list = c("elk_rsf_dat", "winter", "summer"), file = 'data/south-ute-deer-elk-data/elk-rsf-datsets-buff10km.rda')
