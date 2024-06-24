import ee

# Initialize ee
ee.Initialize()

# Some parameters
p_year = 2016 # "present" year
f_year = 2050 # future year
tef_prop = 0.1 # Percent of trailing edge forest lost by 2050 (per Parks et al. 2019)
sqm_per_acre = 4046.86
acre_area = ee.Image.pixelArea().divide(sqm_per_acre)

# Get AOI (elk MCP with 120 km buffer)
aoi = ee.FeatureCollection('projects/GEE_CSP/Pew-Connectivity/elk-mcp').first().geometry().buffer(120000) # **
# ** For final analysis, use buffer size of 50 km (not 120 km)

# NLCD-Based land cover covs
nlcd = ee.ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD')\
    .filter(ee.Filter.eq('system:index', '2016')).first()
nlcd_lc = nlcd.select('landcover')
nlcd_imp = nlcd.select('impervious_descriptor')
developed_pres = ee.Image(0).where(nlcd_lc.gte(21).And(nlcd_lc.lte(24)), 1).clip(aoi)
impervious_pres = ee.Image(0).where(nlcd_imp.gte(24).And(nlcd_imp.lte(26)), 1).clip(aoi)
forest_pres = ee.Image(0).where(nlcd_lc.gte(41).And(nlcd_lc.lte(43)), 1).clip(aoi)
ag_pres = ee.Image(0).where(nlcd_lc.gte(81).And(nlcd_lc.lte(82)), 1).clip(aoi)

#---------------------------------
# FUTURE DEVELOPMENT
# Calculate rate of natural areas loss
aoi_area = aoi.area().divide(sqm_per_acre).getInfo()

# Get state data
states = ee.FeatureCollection('TIGER/2018/States')
co_nm = states.filter(ee.Filter.inList('NAME', ['Colorado','New Mexico'])).union().geometry()
co_nm_area = co_nm.area().divide(sqm_per_acre).getInfo()
# Total natural area in those states
natLC = nlcd_lc.expression(
    "(b(0) > 31 & b(0) != 81 & b(0) != 82) ? 1" +
      ": 0").clip(co_nm)
natLC_area = natLC.multiply(acre_area).rename('area')\
    .reduceRegion(reducer = ee.Reducer.sum(), 
                  geometry = aoi,
                  scale = 30,
                  maxPixels = 1e13).get('area').getInfo()
# Yearly rate of nat area loss in CO and NM btwn 2006 and 2011
# From CSP disappearing west report
loss_per_yr = (150000+73000)/5
# Area lost per year as a proportion of total natural area in those states 
# times number of years between pres and future
prop_loss = (loss_per_yr/co_nm_area) * (f_year-p_year)
area_loss = natLC_area*prop_loss

# Create development density layer with actual development zeroed out
# Will serve to determine where new pixels fall
dev_density = developed_pres.reduceNeighborhood(
    kernel= ee.Kernel.gaussian(radius = 7500, sigma = 500, units = 'meters'),
    reducer= ee.Reducer.mean()).multiply(ee.Number(1000)).clip(aoi)
# Restrict to natural areas
dev_density = dev_density.updateMask(natLC).rename('density')

# Set target proportion of natural area to lose
dev_target_prop = (area_loss/natLC_area)*100

# Threshold dev density layer to just the top target acres
dev_cut = dev_density.reduceRegion(reducer = ee.Reducer.percentile([100-dev_target_prop]),
                                  maxPixels = 1e13,
                                  geometry = aoi,
                                  scale = 250).get('density').getInfo()
new_dev = dev_density.gte(ee.Number(dev_cut))

# Get area of new dev to check how close we are to target
nd_area = new_dev.multiply(acre_area).rename('area')\
    .reduceRegion(reducer = ee.Reducer.sum(), 
                  geometry = aoi,
                  scale = 250,
                  maxPixels = 1e13).get('area').getInfo()
print('CHECK: number of top acres = ', nd_area)
print('Target acres:', area_loss)

# Add future dev to existing dev
dev_fut = developed_pres.where(new_dev.eq(1),1)

#---------------------------------
# FUTURE AGRICULTURE

# Remove future development from existing ag
ag_fut = ag_pres.where(new_dev.eq(1),0)


#---------------------------------
# FUTURE FOREST

# Get trailing edge forest from Parks et al. 2019
# DECISION POINT - Calculte 10% area before or after aligning with existing forest?
tef = ee.Image('projects/GEE_CSP/Pew-Connectivity/trailing-edge-forest')
tef = tef.eq(-1).unmask().clip(aoi)
tef_area = tef.multiply(acre_area).rename('area')\
    .reduceRegion(reducer = ee.Reducer.sum(), 
                  geometry = aoi,
                  scale = 90,
                  maxPixels = 1e13).get('area').getInfo()
tef_loss_area = tef_area*tef_prop

# Remove new (future) development
# dev_fut = ee.Image('projects/GEE_CSP/Pew-Connectivity/future-development') # DELETE ME
forest_fut = forest_pres.where(dev_fut.eq(1),0)

# Get TEF and future forest overlap
tef_for = tef.updateMask(forest_fut).unmask().clip(aoi)
tef_for_a = tef_for.multiply(acre_area).rename('area')
tef_for_area = tef_for_a.reduceRegion(reducer = ee.Reducer.sum(), 
                  geometry = aoi,
                  scale = 250,
                  maxPixels = 1e13).get('area').getInfo()

# create non-forest layer and distance to forest edge
# Use this to preferentially remove TEF along forest edges
non_forest = forest_fut.eq(0).unmask().clip(aoi)
edge_dist = non_forest.distance(kernel = ee.Kernel.euclidean(7500, 'meters'))
edge_dist = ee.Image(1000).divide(edge_dist).where(forest_pres.eq(0),0)
tef_edge = edge_dist.updateMask(tef_for).rename('edge')
# add a bit of random noise to obtain sufficient separation for thresholding
rand = ee.Image.random(seed=930).updateMask(tef_for) 
tef_edge = tef_edge.add(rand)

tef_target_prop = (tef_loss_area/tef_for_area) * 100

# Threshold dev density layer to just the top target acres
tef_cut = tef_edge.reduceRegion(reducer = ee.Reducer.percentile([100-tef_target_prop]),
                                  maxPixels = 1e13,
                                  geometry = aoi,
                                  scale = 90).get('edge').getInfo()
tef_loss = tef_edge.gte(ee.Number(tef_cut))

# # Get area of tef_loss to check against target
# tl_area = tef_loss.multiply(acre_area).rename('area')\
#     .reduceRegion(reducer = ee.Reducer.sum(), 
#                   geometry = aoi,
#                   scale = 90,
#                   maxPixels = 1e13).get('area').getInfo()
# print('CHECK: number of top acres = ', tl_area)
# print('Target acres:', tef_loss_area)


# Remove tef loss from future forest
forest_fut = forest_fut.where(tef_loss.eq(1), 0)


#---------------------------------
# EXPORT

development_path = 'projects/GEE_CSP/Pew-Connectivity/future-development'
if ee.data.getInfo(development_path):
    ee.data.deleteAsset(development_path)
    
dev_task = ee.batch.Export.image.toAsset(image = dev_fut,
                                         description = 'future-development',
                                         assetId = development_path,
                                         scale = 30,
                                         region = aoi,
                                         maxPixels = 1e13)
dev_task.start()

ag_path = 'projects/GEE_CSP/Pew-Connectivity/future-ag'
if ee.data.getInfo(ag_path):
    ee.data.deleteAsset(ag_path)
    
ag_task = ee.batch.Export.image.toAsset(image = ag_fut,
                                         description = 'future-ag',
                                         assetId = ag_path,
                                         scale = 30,
                                         region = aoi,
                                         maxPixels = 1e13)
ag_task.start()


forest_path = 'projects/GEE_CSP/Pew-Connectivity/future-forest'
if ee.data.getInfo(forest_path):
    ee.data.deleteAsset(forest_path)

forest_task = ee.batch.Export.image.toAsset(image = forest_fut,
                                            description = 'future-forest',
                                            assetId = forest_path,
                                            scale = 30,
                                            region = aoi,
                                            maxPixels = 1e13)
forest_task.start()

