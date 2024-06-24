library(tidyverse)
library(sf)

traffic <- st_read("data/co-nm-roads-traffic/aadt_priSecRds_aoi.shp") %>% st_zm() %>% 
  mutate(state = ifelse(AADTYear == 2021, "NM", "CO"),
         yearPres = ifelse(AADTYear == 2021, 2021, 2020),
         yearFut = 2050,
         routeName = ifelse(is.na(RouteID), ROUTE, RouteID), 
         ID = 1:nrow(.)) %>% 
  rename(AADTpres = AADT) %>% 
  select(ID, state, routeName, yearPres, yearFut, AADTpres, ProjectedA)

# Fit linear regression between present and future AADT data for CO
co <- traffic %>% filter(state=='CO') %>% 
  rename(AADT2050 = ProjectedA) %>% 
  mutate(tdiff = AADT2050-AADTpres)
tmod <- lm(tdiff ~ AADTpres + I(AADTpres^2), data = co)
summary(tmod)

# Apply model predictions to NM
nm <- traffic %>% filter(state == 'NM') %>% select(-ProjectedA) 
nm$tdiff<-round(predict(tmod, newdata = data.frame(AADTpres = nm$AADTpres), type = "response"))
nm <- nm %>% 
  mutate(AADT2050 = AADTpres+tdiff)
# Recombine state-level data and export
traffic_out <- rbind(nm, co)
st_write(traffic_out, dsn = "data/co-nm-roads-traffic", layer = 'co-nm-addt-w-future',
         driver = 'ESRI Shapefile')

# Check out mod fit
plot(co$AADTpres, co$tdiff, xlab = 'CO traffic volume 2020', 
     ylab = 'Traffic volume increase (2050-2020)')
lines(sort(co$AADTpres), sort(predict(tmod)),lwd = 2, col = 'indianred')
text(3000,14000, "R^2 = 0.66")

# Plot predicted pres vs future traffic
plot(nm$AADTpres, nm$AADT2050, type = 'n', xlab = 'Present', ylab = 'Future')
points(co$AADTpres, co$AADT2050, pch = 16, col = rgb(0,0,1,0.4))
points(nm$AADTpres, nm$AADT2050, pch = 16, col = rgb(1,0,0,0.5))

