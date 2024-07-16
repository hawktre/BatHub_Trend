## ---------------------------
##
## Script name: 00b_BatGrid_Covariates
##
## Purpose of script: Collect and format NABat grid covariates
##
## Author: Trent VanHawkins
##
## Date Created: 2024-05-07
##


## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(here)
require(sf)
require(terra)

# Read in data ------------------------------------------------------------
## nabat covariates already aggregated for Udell et al., 2022
nabat_covars <- st_read(dsn = here("DataRaw.nosync/covariates/NABat_grid_covariates/"),
                    layer = "NABat_grid_covariates")
## conus_grts key
conus10k <- read_sf(here("DataRaw.nosync/batgrid/complete_conus_mastersample_10km_attributed.shp"))

##Landfire gap cover 
landfire_or <- terra::rast(here("DataRaw.nosync/covariates/LandFire/LF2022_OR/LC22_EVT_230.tif"))
landfire_wa <- rast(here("DataRaw.nosync/covariates/LandFire/LF2022_WA/LC22_EVT_230.tif"))
landfire_id <- rast(here("DataRaw.nosync/covariates/LandFire/LF2022_ID/LC22_EVT_230.tif"))

# subset and plot data ---------------------------------------------------------------
## Define pnw
pnw <- c("Oregon", "Washington", "Idaho")

## subset and plot
nabat_pnw <- nabat_covars %>% 
  filter(admin1 %in% pnw | admin2 %in% pnw) %>% 
  mutate(admin1 = case_when(!(admin1 %in% pnw) ~ admin2,
                               TRUE ~ admin1))


plot(nabat_pnw['admin1'])

pnw_shp <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE)) %>% 
  filter(ID %in% c("oregon", "washington", "idaho")) %>% 
  st_transform(st_crs(nabat_pnw))

nabat_pnw <- nabat_pnw[pnw_shp,]
# Join to get CONUS -------------------------------------------------------

conus_grts_key <- conus10k %>% 
  select(CONUS_10KM, GRTS_ID)

## Join with conus and do some formatting
conus_pnw_covars <- nabat_pnw %>% 
  left_join(as.data.frame(conus_grts_key) %>% select(-geometry), by = "GRTS_ID") %>% 
  select(CONUS_10KM, GRTS_ID, 
         "state" = admin1,
         long, lat, karst, p_forest, p_wetland, mean_temp, precip, DEM_max, 
         physio_div, dist_mines, riverlake, geometry)

## rename to make it easier to call
covars <- conus_pnw_covars

plot(covars["state"])
# Get Cliff_Canyon ---------------------------------------------------------------
## Create a single layer
landfire <- terra::merge(landfire_or, landfire_wa, landfire_id, first = T)

## subset all cliffs and canyons categories
lf_vals <- as.data.frame(unique(landfire)) %>% 
  filter(str_detect(EVT_NAME, "Cliff") | str_detect(EVT_NAME, "Canyon"))

## Keep cliffs and canyons only
cliff_canyon <- landfire %in% lf_vals$EVT_NAME

## project everything to equal area projection
covars_ea <- covars %>% 
  st_transform(st_crs("EPSG:5070"))

cliff_canyon_ea <- cliff_canyon %>% 
  terra::project("EPSG:5070", method = "near")

## Convert to terra spatial vector for next operation
covars_ea_sv <- vect(covars_ea)

## resample cliff_canyon (takes a long time)
cliffcanyon_st <- terra::extract(cliff_canyon_ea, covars_ea_sv, mean)

## Join results of extract with sf object
covars <- covars %>% 
  rowid_to_column("ID") %>% 
  left_join(cliffcanyon_st, by = "ID")


# Write out the results ---------------------------------------------------
write_sf(covars, here("DataProcessed.nosync/occurrence/batgrid_covars.shp"))
