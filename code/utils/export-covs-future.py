# ELK COVARIATES - FUTURE
# Collect and export continuous covariates for elk models
# Script only exports covs (climate and land cover) that change by 2050
import ee
import math

# Initialize ee
ee.Initialize()

# Get AOI (elk MCP with 120 km buffer)
aoi = ee.FeatureCollection('projects/GEE_CSP/Pew-Connectivity/elk-mcp').first().geometry().buffer(120000)

# export scale and projection
scale = 90
projection = ee.Projection('ESPG:26913') # NAD83 UTM Zone 13N

# Choose radii for summarizing covariates
rad_large = 1000 # For summer/winter range models. Based on 80% quantile of 8-hr movement steps
name_large = "_1km"
rad_small = 300 # For migration model. Based on average step during migration at 2-hr interval 
name_small = "_300m"

#----------------------------------------
# FUNCTIONS

# Focal mean
def focal_mean(image, radius, unit, name):
    names = image.bandNames().getInfo()
    new_names = [s + name for s in names]
    return image.reduceNeighborhood(kernel = ee.Kernel.circle(radius, unit),
                                    reducer = ee.Reducer.mean()).rename(new_names)
# Gaussian density estimate
def gaussian_mean(image, radius, sigma, unit, name):
    names = image.bandNames().getInfo()
    new_names = [s + name for s in names]
    return image.reduceNeighborhood(kernel = ee.Kernel.gaussian(radius, sigma, unit), 
                                    reducer = ee.Reducer.mean()).rename(new_names)
# Focal sum
def focal_sum(image, radius, unit):
    return image.reduceNeighborhood(kernel = ee.Kernel.circle(radius, unit, False),
                                    reducer = ee.Reducer.sum())
# Focal count
def focal_count(image, radius, unit):
    return image.reduceNeighborhood(kernel = ee.Kernel.circle(radius, unit, False),
                                    reducer = ee.Reducer.count())
# Percent cover
def percent_cov(image, radius, unit, name):
    names = image.bandNames().getInfo()
    new_names = [s + name for s in names]
    isum = focal_sum(image, radius, unit)
    icount = focal_count(image, radius, unit)
    return isum.divide(icount).rename(new_names)

#---------------------------------------------
# CLIMATE COVS
clim1 = ee.Image("projects/GEE_CSP/Pew-Connectivity/ClimateNA-future-annual-tmax-tmin-2041_2070")\
    .rename(['Tmax_annual', 'Tmin_annual']).clip(aoi)
clim2 = ee.Image("projects/GEE_CSP/Pew-Connectivity/ClimateNA-future-climate-2041_2070")\
    .rename(['CMD','MAR','PAS','PPT_at','PPT_sm','PPT_sp','PPT_wt','Tave_at','Tave_sm','Tave_sp','Tave_wt']).clip(aoi)

# Combine and smooth
clim_all = ee.Image([clim1, clim2])
clim_large = focal_mean(clim_all, rad_large, "meters", name_large)
clim_small = focal_mean(clim_all, rad_small, "meters", name_small)

#---------------------------------------------
# LAND COVER
development_path = 'projects/GEE_CSP/Pew-Connectivity/future-development'
forest_path = 'projects/GEE_CSP/Pew-Connectivity/future-forest'
ag_path = 'projects/GEE_CSP/Pew-Connectivity/future-ag'
dev = ee.Image(development_path).rename('dev')
forest = ee.Image(forest_path).rename('forest')
ag = ee.Image(ag_path).rename('ag')

lc_all = ee.Image([dev, forest, ag])
lc_large = percent_cov(lc_all, rad_large, 'meters', '_pcov' + name_large)
lc_small = percent_cov(lc_all, rad_small, 'meters', '_pcov' + name_small)

#---------------------------------------------
# LANDCOVER DISTANCE COVARIATES
# Create distance to each lc type
dev_dist = dev.distance(ee.Kernel.euclidean(100000,'meters')).rename('dev_dist')
ag_dist = ag.distance(ee.Kernel.euclidean(100000,'meters')).rename('ag_dist')
forest_dist = forest.distance(ee.Kernel.euclidean(100000,'meters')).rename('forest_dist')

#---------------------------------------------
# ROADS
# Get roads layer and rasterize
rd = ee.FeatureCollection("projects/GEE_CSP/Pew-Connectivity/co-nm-traffic-w-future")
rd_img = rd.reduceToImage(properties = ['AADT2050'], reducer = ee.Reducer.max()).unmask(0).clip(aoi)

# Create kernel density smooth of traffic volume
rd_smooth = rd_img.reduceNeighborhood(kernel = ee.Kernel.gaussian(10000, 1000, 'meters'), 
                                      reducer = ee.Reducer.mean()).rename('rd_traffic_smooth')
#---------------------------------------------
# COMBINE AND EXPORT

# Pull all non-road covariates together 
covs_all = ee.Image([clim_large, clim_small, lc_large, lc_small])
# Convert all to float-32
covs_all = covs_all.float().clip(aoi)

task1 = ee.batch.Export.image.toDrive(image = covs_all,
                                     folder = 'elk-covs',
                                     description = 'elk-covs-future-' + str(scale) + "m",
                                     scale = scale,
                                     region = aoi,
                                     maxPixels = 1e13,
                                     crs = "EPSG:26913")
task1.start()

#----------
# Pull all coarser scale covs together
# These need to be exported at coarser resolution to deal with 
# max kernel size issues and need to cover full study area

covs_lc_rd = ee.Image([dev_dist, ag_dist, forest_dist, rd_smooth]).float().clip(aoi)

task2 = ee.batch.Export.image.toDrive(image = covs_lc_rd,
                                     folder = 'elk-covs',
                                     description = 'elk-covs-lc-road-future-450m',
                                     scale = 450,
                                     region = aoi,
                                     maxPixels = 1e13,
                                     crs = "EPSG:26913")
task2.start()