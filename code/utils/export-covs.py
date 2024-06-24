# ELK COVARIATES - PRESENT
# Collect and export continuous covariates for elk models
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
# Vector ruggedness measure
def compute_vrm(slope_img, aspect_img, radius, units):
    slope_sine = slope_img.sin()
    x_sum_sq = focal_sum(slope_sine.multiply(aspect_img.sin()), radius, units).pow(2)
    y_sum_sq = focal_sum(slope_sine.multiply(aspect_img.cos()), radius, units).pow(2)
    z_sum_sq = focal_sum(slope_img.cos(), radius, units).pow(2)
    n = focal_sum(ee.Image(1), radius, units)
    r = x_sum_sq.add(y_sum_sq).add(z_sum_sq).sqrt()
    vrm_img = ee.Image(1).subtract(r.divide(n))
    return vrm_img

#---------------------------------------------
# TERRAIN VARIABLES
# Create some topographic viables
dsm = ee.ImageCollection("JAXA/ALOS/AW3D30/V3_2") # Digital surface model. Native res 30m
projElev = dsm.first().select(0).projection()
elev = dsm.select("DSM").mosaic().setDefaultProjection(projElev).rename('elevation') # Elevation
slope = ee.Terrain.slope(elev).rename('slope') # Slope
aspect = ee.Terrain.aspect(elev).rename('aspect') # Aspect
# Get vector ruggedness metric
window_radius = 1000
slopeRad = slope.multiply(ee.Number(math.pi).divide(180))
aspectRad = slope.multiply(ee.Number(math.pi).divide(180))
vrm = compute_vrm(slopeRad, aspectRad, window_radius, "meters").rename('vrm')
chili = ee.Image("CSP/ERGo/1_0/US/CHILI").rename('CHILI')

# Combine and smooth
terrain_all = ee.Image([elev, slope, aspect, vrm, chili])
terrain_large = focal_mean(terrain_all, rad_large, "meters", name_large)
terrain_small = focal_mean(terrain_all, rad_small, "meters", name_small)

#---------------------------------------------
# CLIMATE COVS
clim1 = ee.Image("projects/GEE_CSP/Pew-Connectivity/ClimateNA-present-annual-tmax-tmin-1991_2020")\
    .rename(['Tmax_annual', 'Tmin_annual']).clip(aoi)
clim2 = ee.Image("projects/GEE_CSP/Pew-Connectivity/ClimateNA-present-climate-1991_2020")\
    .rename(['CMD','PAS','PPT_at','PPT_sm','PPT_sp','PPT_wt','Tave_at','Tave_sm','Tave_sp','Tave_wt']).clip(aoi)
clim3 = ee.Image("projects/GEE_CSP/Pew-Connectivity/climateNA-MAR_1981_2010").rename('MAR_81-10').clip(aoi)

# Combine and smooth
clim_all = ee.Image([clim1, clim2, clim3])
clim_large = focal_mean(clim_all, rad_large, "meters", name_large)
clim_small = focal_mean(clim_all, rad_small, "meters", name_small)

#---------------------------------------------
# LAND COVER
# NLCD-Based land cover covs
nlcd = ee.ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD')\
    .filter(ee.Filter.eq('system:index', '2016')).first()
nlcd_lc = nlcd.select('landcover')
water = nlcd_lc.eq(11).rename('water')
dev = ee.Image(0).where(nlcd_lc.gte(21).And(nlcd_lc.lte(24)), 1).rename('dev')
forest = ee.Image(0).where(nlcd_lc.gte(41).And(nlcd_lc.lte(43)), 1).rename('forest')
ag = ee.Image(0).where(nlcd_lc.gte(81).And(nlcd_lc.lte(82)), 1).rename('ag')

lc_all = ee.Image([water, dev, forest, ag])
lc_large = percent_cov(lc_all, rad_large, 'meters', '_pcov' + name_large)
lc_small = percent_cov(lc_all, rad_small, 'meters', '_pcov' + name_small)

#---------------------------------------------
# O&G WELL DENSITY
# Read in well data for each state
co_wells = ee.FeatureCollection("projects/GEE_CSP/Pew-Connectivity/CO_Oil_and_Gas_Locations")
nm_wells = ee.FeatureCollection("projects/GEE_CSP/HM/sources/NM_wells_2018")
# Select only active wells
def isWell(ft):
    return ft.set({"is_well":1})
co_wells_active = co_wells.filter(ee.Filter.eq("fac_status","AC")).map(isWell)
nm_wells_active = nm_wells.filter(ee.Filter.eq('status','A')).map(isWell)
# Combine state-level data
active_wells = ee.FeatureCollection([co_wells_active, nm_wells_active]).flatten().select("is_well")
# Create image
well_img = active_wells.reduceToImage(properties = ['is_well'], 
                                      reducer = ee.Reducer.max()).unmask(0).rename('well_dens')
# Get well density w/in large and small radii (defined above)
well_large = gaussian_mean(well_img, rad_large, rad_large*0.25, 'meters', name_large)
well_small = gaussian_mean(well_img, rad_small, rad_small*0.25, 'meters', name_small)

#---------------------------------------------
# LANDCOVER & WELL DISTANCE COVARIATES
# Create distance to each lc type
dev_dist = dev.distance(ee.Kernel.euclidean(100000,'meters')).rename('dev_dist')
ag_dist = ag.distance(ee.Kernel.euclidean(100000,'meters')).rename('ag_dist')
forest_dist = forest.distance(ee.Kernel.euclidean(100000,'meters')).rename('forest_dist')
water_dist = water.distance(ee.Kernel.euclidean(100000,'meters')).rename('water_dist')
well_dist = well_img.distance(ee.Kernel.euclidean(100000,'meters')).rename('well_dist')

#---------------------------------------------
# ROADS
# Get roads layer and rasterize
rd = ee.FeatureCollection("projects/GEE_CSP/Pew-Connectivity/co-nm-roads-traffic")
rd_img = rd.reduceToImage(properties = ['AADT'], reducer = ee.Reducer.max()).unmask(0).clip(aoi)

# Create kernel density smooth of traffic volume
rd_smooth = rd_img.reduceNeighborhood(kernel = ee.Kernel.gaussian(10000, 1000, 'meters'), 
                                      reducer = ee.Reducer.mean()).rename('rd_traffic_smooth')
# Create distance to road raster
rd_dist = rd_img.distance(ee.Kernel.euclidean(50000,'meters')).rename('road_dist')

#---------------------------------------------
# COMBINE AND EXPORT

# Pull all fine scale covariates together 
covs_all = ee.Image([terrain_large, terrain_small, clim_large, clim_small,
                     lc_large, lc_small, well_large, well_small, rd_dist])
# Convert all to float-32
covs_all = covs_all.float().clip(aoi)

task1 = ee.batch.Export.image.toDrive(image = covs_all,
                                     folder = 'elk-covs',
                                     description = 'elk-covs-present-' + str(scale) + "m",
                                     scale = scale,
                                     region = aoi,
                                     maxPixels = 1e13,
                                     crs = "EPSG:26913")
task1.start()

#----------
# Pull all coarser scale covs together
# These need to be exported at coarser resolution to deal with 
# max kernel size issues and need to cover full study area

covs_lc_rd = ee.Image([dev_dist, ag_dist, forest_dist, water_dist, well_dist, rd_dist, rd_smooth]).float().clip(aoi)

task2 = ee.batch.Export.image.toDrive(image = covs_lc_rd,
                                     folder = 'elk-covs',
                                     description = 'elk-covs-lc-road-present-450m',
                                     scale = 450,
                                     region = aoi,
                                     maxPixels = 1e13,
                                     crs = "EPSG:26913")
task2.start()