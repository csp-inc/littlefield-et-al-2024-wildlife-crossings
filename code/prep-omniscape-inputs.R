library(tidyverse)
library(terra)
library(sf)

# Define which seasonal hab suit layers to use as source and ground
era = "present"
sm_name <- paste0("sm-", era)
wt_name <- paste0("wt-", era)

# Grab relevant layers
res <- rast("output/connectivity/resistance/fut-all-cross-ne8.asc") # only used to set extent
sm_raw <- rast(paste0("output/", sm_name, "-exp-aoi50.tif")) 
ext(sm_raw) <- ext(res)
wt_raw <- rast(paste0("output/", wt_name, "-exp-aoi50.tif")) 
ext(wt_raw) <- ext(res)

# Rescale function
rescale01 <- function(r){
  mm <- minmax(r)
  rout <- (r - mm[1])/(mm[2] - mm[1])
  return(rout)
}

#---------------------------
# Prep continuous src and cnd layers

# Source
src_rs1 <- rescale01(sm_raw)
src_rs2 <- rescale01(wt_raw)
src_rs1[src_rs1 < 0.1] <- 0
src_rs2[src_rs2 < 0.1] <- 0
src <- src_rs1 + src_rs2
src <- src * 10
src_fname <- paste0("output/connectivity/os-inputs/", era, "-src.asc")
src %>% terra::writeRaster(filename = src_fname,
                          filetype = "AAIGrid", overwrite = TRUE)

# Conditional file for winter
cnd_wt <- rescale01(wt_raw)
cnd_wt[cnd_wt < 0.1] <- -100
cnd_wt[cnd_wt >= 0.1] <- 1
cnd_wt_fname <- paste0("output/connectivity/os-inputs/", wt_name,"-cnd.asc")
cnd_wt %>% terra::writeRaster(filename = cnd_wt_fname,
                          filetype = "AAIGrid", overwrite = TRUE)

# Conditional file for summer
cnd_sm <- rescale01(sm_raw)
cnd_sm[cnd_sm < 0.1] <- -100
cnd_sm[cnd_sm >= 0.1] <- 1
cnd_sm_fname <- paste0("output/connectivity/os-inputs/", sm_name,"-cnd.asc")
cnd_sm %>% terra::writeRaster(filename = cnd_sm_fname,
                              filetype = "AAIGrid", overwrite = TRUE)
