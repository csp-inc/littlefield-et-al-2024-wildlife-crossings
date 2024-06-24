library(tidyverse)
library(terra)
library(sf)

# Define connectivity scenario to use
scenario = "fut-no-cross"
base_path = "output/connectivity/omniscape/"

# Get relevant current flow maps
s2w <- rast(paste0(base_path, scenario, "/", scenario, "-s2w/cum_currmap.tif"))
w2s <- rast(paste0(base_path, scenario, "/", scenario, "-w2s/cum_currmap.tif"))

# Combine seasonal maps
mcurr <- mean(s2w, w2s)
scurr <- sum(s2w, w2s)

# Export
terra::writeRaster(mcurr, filename = paste0(base_path, scenario, "/mean-current.tif"),overwrite = TRUE)
terra::writeRaster(scurr, filename = paste0(base_path, scenario, "/sum-current.tif"),overwrite = TRUE)
