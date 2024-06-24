library(sf)
library(tidyverse)
library(adehabitatLT)
library(migrateR)
library(lubridate)

# set crs
su_crs <- st_crs(26913)

# Load cleaned data
load(file = "data/south-ute-deer-elk-data/su-elk-cleaned-20220601.rda", verbose = T)

# Vetted elk bursts based on visual inspection of nsd plot
vet <- read_csv("data/south-ute-deer-elk-data/elk-bursts-vetted.csv") %>% 
  dplyr::select(yr_burst, mod, keep, status)
pburst <- vet$yr_burst[vet$status=='params' & vet$keep == "keep"]
gburst <- vet$yr_burst[vet$status %in% c('good', 'force_mig') & vet$keep == "keep"]

# Add yearly "burst" column to separate animal years and only keep those burst
# with at least min_days of data
elkClean <- elkClean %>% 
  mutate(yr_burst = paste(Animal_id, Year, sep = "-"))
# Find bursts with at least 300 days of data
min_days <- 300 
burst_days <- elkClean %>% 
  group_by(yr_burst) %>% 
  summarise(n_days = length(unique(yday(date))))
burst_keep <- burst_days$yr_burst[burst_days$n_days >= min_days]
# Filter shorter bursts
elkBurst <- elkClean %>% filter(yr_burst %in% burst_keep) %>% 
  arrange(yr_burst)

#---------------------------------------------------------
#---------------------------------------------------------
# Run movement modeling and point classification for
# "good" bursts (i.e., those that were well fit by 
# general model parameters - see vet dataset)

# Create adehabitat trajectory object
gdat <- elkBurst %>% 
  filter(yr_burst %in% gburst)
gtraj <- as.ltraj(xy = gdat[,c('x','y')], 
                  date = gdat$date,
                  id = gdat$Animal_id,
                  burst = gdat$yr_burst,
                  infolocs = data.frame(locID = gdat$locID))

# Fit NSD models and refine
# Define some new starting parameters
# s.d  = starting value for dist between ranges, l.r  = min days on second range
pe_new <- pEst(s.d = 1000, l.r = 65)
gnsd <- mvmtClass(gtraj)
gnsd <- refine(gnsd, p.est = pe_new)
length(which(!fullmvmt(gnsd))) # recheck fits
# extract top model name for each burst
tmods <- topmvmt(gnsd, mdelta = 1000)
top_mods_g <- data.frame(yr_burst = attributes(tmods)$burst, 
                       mod = attributes(tmods)$names)
# Reclassify some top mods to force "migrant" as top
top_mods_g$mod[which(top_mods_g$yr_burst %in% c("84217-2018","84913-2014","98517-2018"))] <- "migrant"

# Classify each location as migratory or not and
# produce plots for model checking
pdf(file = "data/south-ute-deer-elk-data/nsd-test-good.pdf", width = 20, height = 8)
par(mfrow = c(1,3))
# spatmig(gtraj,gnsd)
new_data_g <- data.frame()
for(i in 1:length(gburst)){
  # get burst ID and top model type
  bid = gburst[i]
  tm = top_mods_g$mod[top_mods_g$yr_burst==bid]
  print(bid)
  
  # Plot mvmt model
  plot.mvmt(gnsd[[bid]]) 
  
  # Classify points based on mvmt mod into migration, winter and summer
  mod_type = ifelse(tm %in% c('migrant','mixmig'), tm, 'migrant')
  st_end <- mvmt2dt(gnsd[[bid]], mod = mod_type)[[1]]
  td_burst <- gdat %>% 
    filter(yr_burst == bid) %>% 
    mutate(migrating = ifelse(date >= st_end$date[1] & 
                                date <= st_end$date[2] | 
                                date >= st_end$date[3] &
                                date <= st_end$date[4],1,0))
  td_burst$season <- "NA"
  td_burst$season[td_burst$migrating==1] <- "migration"
  td_burst$season[td_burst$date < st_end$date[1] | td_burst$date > st_end$date[4]] <- "winter"
  td_burst$season[td_burst$date > st_end$date[2] & td_burst$date < st_end$date[3]] <- "summer"
  # test plots
  plot(td_burst$date, td_burst$NSD, col = as.numeric(as.factor(td_burst$season)), pch = 16, cex = 0.5, 
       xlab = 'Date', ylab = 'NSD', main = tm)
  plot(td_burst$x, td_burst$y, col = as.numeric(as.factor(td_burst$season)), pch = 16, cex = 0.5, axes = F,
       xlab = "", ylab = "")
  
  # Add to data set
  new_data_g <- rbind(new_data_g, td_burst)
}
dev.off()


#---------------------------------------------------------
#---------------------------------------------------------
# Run movement modeling and point classification for
# "param" bursts (i.e., those that require tweaking of
# model parameters to achieve better differentiation
# between ranging and migrating - see vet dataset)

# Create adehabitat trajectory object
pdat <- elkBurst %>% 
  filter(yr_burst %in% pburst)
ptraj <- as.ltraj(xy = pdat[,c('x','y')], 
                  date = pdat$date,
                  id = pdat$Animal_id,
                  burst = pdat$yr_burst,
                  infolocs = data.frame(locID = pdat$locID))

# Fit NSD models and refine
# Define some new starting parameters
# s.d  = starting value for dist between ranges, l.r  = min days on second range
pe_new <- pEst(s.d = 1000, l.r = 90, s.r = 100)
pnsd <- mvmtClass(ptraj, p.est = pe_new)
length(which(!fullmvmt(pnsd))) # recheck fits
# extract top model name for each burst
tmods <- topmvmt(pnsd, mdelta = 1000)
top_mods_p <- data.frame(yr_burst = attributes(tmods)$burst, 
                         mod = attributes(tmods)$names)
# Reclassify some top mods to force "migrant" as top
top_mods_p$mod[which(top_mods_p$yr_burst %in% c("84217-2018","84913-2014","98517-2018"))] <- "migrant"

# Classify each location as migratory or not and
# produce plots for model checking
pdf(file = "data/south-ute-deer-elk-data/nsd-test-param.pdf", width = 20, height = 8)
par(mfrow = c(1,3))
new_data_p <- data.frame()
for(i in 1:length(pburst)){
  # get burst ID and top model type
  bid = pburst[i]
  tm = top_mods_p$mod[top_mods_p$yr_burst==bid]
  print(bid)
  
  # Plot mvmt model
  plot.mvmt(pnsd[[bid]]) 
  
  # Classify points based on mvmt mod
  mod_type = ifelse(tm %in% c('migrant','mixmig'), tm, 'migrant')
  st_end <- mvmt2dt(pnsd[[bid]], mod = mod_type)[[1]]
  td_burst <- pdat %>% 
    filter(yr_burst == bid) %>% 
    mutate(migrating = ifelse(date >= st_end$date[1] & 
                                date <= st_end$date[2] | 
                                date >= st_end$date[3] &
                                date <= st_end$date[4],1,0))
  td_burst$season <- "NA"
  td_burst$season[td_burst$migrating==1] <- "migration"
  td_burst$season[td_burst$date < st_end$date[1] | td_burst$date > st_end$date[4]] <- "winter"
  td_burst$season[td_burst$date > st_end$date[2] & td_burst$date < st_end$date[3]] <- "summer"
  # test plots
  plot(td_burst$date, td_burst$NSD, col = as.numeric(as.factor(td_burst$season)), pch = 16, cex = 0.5, 
       xlab = 'Date', ylab = 'NSD', main = tm)
  plot(td_burst$x, td_burst$y, col = as.numeric(as.factor(td_burst$season)), pch = 16, cex = 0.5, axes = F,
       xlab = "", ylab = "")
  
  # Add to data set
  new_data_p <- rbind(new_data_p, td_burst)
}
dev.off()

#----------------------------------------
#----------------------------------------
# Combine data sets, apply seasonal classification, and 
# apply "migration quality" classification

elkFinal <- rbind(new_data_g, new_data_p) %>% 
  arrange(Animal_id, yr_burst, Fix) %>% 
  # dplyr::select(!c(step, angle, NSD, dt, dt_hr, speed, fast)) %>% 
  mutate(mig_quality = ifelse(yr_burst %in% gburst, 1, 2)) 
  
save(elkFinal, file = "data/south-ute-deer-elk-data/su-elk-final-w-season-20220615.rda")

