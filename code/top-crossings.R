library(tidyverse)
library(terra)
library(sf)

# Define connectivity scenario to use
scenario = "-no-cross"
base_path = "output/connectivity/omniscape/"

# Define AOI and CRS
su_crs <- st_crs(26913)
aoi50 <- st_read("data/south-ute-deer-elk-data/elk-mcp/elk-mcp.shp") %>% st_buffer(dist = 50000)

# Get relevant current flow maps
pcurr <- rast(paste0(base_path, "pres", scenario, "/", "mean-current.tif"))
fcurr <- rast(paste0(base_path, "fut", scenario, "/", "mean-current.tif"))

# Get Crossings
xing <- st_read("data/crossing-structures/crossings-final/crossings-exist-recc-final-20230311.shp") %>% 
  st_buffer(1500) %>% st_filter(aoi50, .predicate = st_within) %>% 
  mutate(ER01 = ifelse(existRecc == 'exist', 1, 0))

# Create crossing 'networks' by merging all overlapping crossings into a single polygon
# Done separately for existing and recommended
xing_exist_merge <- xing %>% 
  filter(existRecc == 'exist') %>% 
  st_union() %>% 
  st_cast("POLYGON") %>% 
  st_as_sf(data.frame(existRecc = 'exist'))
xing_recc_merge <- xing %>% 
  filter(existRecc == 'recc') %>% 
  st_union() %>% 
  st_cast("POLYGON") %>% 
  st_as_sf(data.frame(existRecc = 'recc'))
xing_merge <- rbind(xing_exist_merge, xing_recc_merge)

# Summarize current flow within crossings/crossing networks
xing$pres_mean <- extract(pcurr, vect(xing), fun = mean)[,2]
xing$pres_sum <- extract(pcurr, vect(xing), fun = sum)[,2]
xing$fut_mean <- extract(fcurr, vect(xing), fun = mean)[,2]
xing$fut_sum <- extract(fcurr, vect(xing), fun = sum)[,2]

xing_merge$pres_mean <- extract(pcurr, vect(xing_merge), fun = mean)[,2]
xing_merge$pres_sum <- extract(pcurr, vect(xing_merge), fun = sum)[,2]
xing_merge$fut_mean <- extract(fcurr, vect(xing_merge), fun = mean)[,2]
xing_merge$fut_sum <- extract(fcurr, vect(xing_merge), fun = sum)[,2]

# Identify top 10 crossing for present and future current flow
# Done for all crossings and just recommended crossings
recc_pres_top10 <- xing %>% filter(existRecc=='recc') %>% 
  arrange(desc(pres_mean)) %>% slice(1:10) %>% 
  summarise(t10 = min(pres_mean)) %>% pull(t10)
recc_fut_top10 <- xing %>% filter(existRecc=='recc') %>% 
  arrange(desc(fut_mean)) %>% slice(1:10) %>% 
  summarise(t10 = min(fut_mean)) %>% pull(t10)

all_pres_top10 <- xing %>% 
  arrange(desc(pres_mean)) %>% slice(1:10) %>% 
  summarise(t10 = min(pres_mean)) %>% pull(t10)
all_fut_top10 <- xing %>% 
  arrange(desc(fut_mean)) %>% slice(1:10) %>% 
  summarise(t10 = min(fut_mean)) %>% pull(t10)

xing$recc_pTop10 <- ifelse(xing$pres_mean >= recc_pres_top10, 1, 0)
xing$recc_fTop10 <- ifelse(xing$fut_mean >= recc_fut_top10, 1, 0)
xing$all_pTop10 <- ifelse(xing$pres_mean >= all_pres_top10, 1, 0)
xing$all_fTop10 <- ifelse(xing$fut_mean >= all_fut_top10, 1, 0)

pfTop <- which(xing$recc_pTop10 == 1 & xing$recc_fTop10 == 1)
xing$recc_pTop10[pfTop] <- 2
xing$recc_fTop10[pfTop] <- 2

# Export as new shapefile
st_write(xing, dsn = "output/crossings-w-current-summary", 
         layer = "all-crossings-w-current-summary",
         driver = "ESRI Shapefile", delete_layer = TRUE)


plot(pcurr)
xing %>% 
  filter(existRecc=='recc', recc_pTop10 == 1) %>%
   st_geometry() %>% plot(add = T, col = rgb(1,0,0,0.5))
xing %>% 
  filter(existRecc=='recc', recc_fTop10 == 1) %>%
  st_geometry() %>% plot(add = T, col = rgb(0,0,1,0.5))

