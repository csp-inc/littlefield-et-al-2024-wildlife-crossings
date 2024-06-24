# Littlefield et al. 2024. Evaluating and Elevating the Role of Wildlife Road Crossings in Climate Adaptation. Frontiers in Ecology and the Environment
## Elk Climate-informed Wildlife Crossings Case Study Analysis

### Organization: Conservation Science Partners, Inc.
### Contact: Dr. Justin Suraci (justin@csp-inc.org)

This repo contains all code required to replication the analyses described in [ADD LINK TO PAPER]. This work examined the effects of climate change and projected future development on elk habitat suitability on, and migratory connectivity between, summer and winter ranges, and the role of proposed wildlife crossing structures in supporting migratory connectivity under present and potential future conditions.
<br>
<br>
## *Analysis Workflow* <br>

## Step 1 - Use net squared displacement (NSD) models to tag elk locations as being on winter range, summer range, or migration routes
NSD models were run using `code/elk-fit-nsd.R`. This step required some iteration and manual intervention to visually inspect plotted NSD outputs and to determine which elk-years should be retained for downstream analysis (i.e., only those for which 'migration' or 'mixed-migration' NSD models were best fit for elk-year track).
<br>
<br>

## Step 2 - Prep environmental covariates
With the exception of traffic layers, all covariates used in habitat selection models and resistance surfaces were prepared in Google Earth Engine. We created present (2016) and future (2050) land cover layers (i.e., predicted changes in development and forest cover) using `code/utils/pres-future-land-cover.py` and exported those and other covariates for use in downstream analyses using `code/utils/export-covs.py` and `code/utils/export-covs-future.py`. For traffic layers, we used estimates of present day and future traffic volume on Colorado highways provided by CDOT, and predicted future volumes from present day volumes for New Mexico highways using `code/utils/prep-pres-future-traffic-volume.R`
<br>
<br>

## Step 3 - Seasonal range RSFs
We ran winter range and summer range resources selection function models using `code/utils/prep-rsf-datasets.R` and `code/fit-summer-winter-rsf-V5.R`. We then created prediction surfaces from top model outputs using `code/rsf-prediction-surface.R`
<br>
<br>

## Step 4 - Migration iSSA and resistance surface
We ran integrated Step Selection Function models for elk migratory movements using `code/utils/prep-mig-ssf-datasets.R` and `code/fit-migration-ssf-V2.R`. We then created a landscape resistance surface from top model output using `code/ssf-prediction-surface-v2.R`. Note that the latter script produces multiple resistances surfaces, which were explored in preliminary analyses, but the only surfaces carried forward for the final connectivity models were "output/connectivity/resistance/pres-no-cross-ne8.asc" and "output/connectivity/resistance/fut-no-cross-ne8.asc"
<br>
<br>

## Step 5 - Connectivity analyses
Connectivity models were run in Omniscape.jl via docker (image created using `docker pull vlandau/omniscape:latest`; container launched using `run_ominscape.sh`). Input files for Omniscape runs were preped using `code/prep-omniscape-inputs.R` and .ini files for individual model runs are stored in `code/conn-mod-runs/omniscape`. Omniscape outputs (i.e, winter-to-summer model and summer-to-winter model for each present and future time point) were summarized using `code/prep-omniscpae-outputs.R`. We compared connectivity value of proposed crossing structure locations using `code/top-crossings.R`



