## ---------------------------
##
## Script name: 00_format_spOcc.R
##
## Purpose of script: Format data specifically for fitting models with spOcc
##
## Author: Trent VanHawkins
##
## Date Created: 2024-07-16
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(sf)
library(spOccupancy)
library(bayesplot)

covars <- read_sf(here("DataProcessed.nosync/occurrence/batgrid_covars.shp"))
dets <- readRDS(here("DataProcessed.nosync/detections/nw_nights.rds"))
possible_bats <- c("laci",
                   "lano",
                   "myev",
                   "epfu",
                   "myyu",
                   "myth",
                   "myci",
                   "myvo",
                   "anpa",
                   "pahe",
                   "euma",
                   "myca",
                   "mylu",
                   "coto")

anpa_range <- read_sf(here("Background/nabat-code-tutorial.nosync/RawData_BatRanges/ANPA/mPABAx_CONUS_Range_2001v1.shp")) #using this for parastrellus too
euma_range <- read_sf(here("Background/nabat-code-tutorial.nosync/RawData_BatRanges/EUMA/EUMA_WBWG_FINAL.shp"))
myci_range <- read_sf(here("Background/nabat-code-tutorial.nosync/RawData_BatRanges/MYCI/mWSFMx_CONUS_Range_2001v1.shp"))

# Create Data Array -------------------------------------------------------

## Join dets and covars
covars_join <- covars %>% 
  rename("cell" = CONUS_10KM,
         "cliff_canyon" = EVT_NAME) %>% 
  st_transform(crs = st_crs(anpa_range)) %>% 
  select(-c(lat, long)) %>%
  select(-riverlake) %>% 
  mutate(p_forest = log(p_forest + 1),
         precip = log(precip + 1),
         cliff_canyon = log(cliff_canyon+1),
         anpa_range = as.factor(if_else(lengths(st_intersects(., anpa_range)) > 0, 1, 0)),
         euma_range = as.factor(if_else(lengths(st_intersects(., euma_range)) > 0, 1, 0)),
         myci_range = as.factor(if_else(lengths(st_intersects(., myci_range)) > 0, 1, 0)),
         pahe_range = euma_range,
         across(karst:cliff_canyon, ~scale(.x)[,1])) 

## Write out the scaled covariates
write_sf(covars_join, here("DataProcessed.nosync/occurrence/batgrid_covars_scaled.shp"))

covars_join <- st_drop_geometry(covars_join)

dets <- dets %>%
  left_join(covars_join, by = "cell") %>% 
  #mutate(clutter = factor(clutter, levels = c(-1, 0, 1, 2, 3), labels = c(-1, 0, 1, 2, 3))) %>% 
  mutate(clutter = as.integer(clutter) - 2) %>% 
  group_by(cell, year) %>%
  mutate(replicate_id = as.numeric(factor(replicate, levels = unique(replicate)))) %>%
  ungroup() %>% 
  drop_na(DEM_max, p_forest, precip, cliff_canyon) %>% 
  arrange(year, cell, replicate_id) 

dets.plt <- dets %>% 
  pivot_longer(cols = all_of(possible_bats), names_to = "spp", values_to = "det") %>% 
  st_as_sf(coords = c("lon", "lat"), crs = "WGS84") %>% 
  filter(det == 1)

ggplot()+
  geom_sf(data = dets.plt, aes(color = spp))+
  facet_wrap(~spp)

pairs(dets %>% select(DEM_max, p_forest, mean_temp, precip, cliff_canyon))
# Pivot long --------------------------------------------------------------

dets_long <- dets %>%
  pivot_longer(cols = all_of(possible_bats), names_to = "spp", values_to = "occ") %>% 
  select(spp, cell, year, replicate_id, occ) %>% 
  arrange(spp, year, cell, replicate_id) %>% 
  mutate(across(1:4, as.factor)) 

# Detection Data Array  -------------------------------------------------------------
## use tapply to create array with dim(spp,cell,year,replicate)
y <- tapply(dets_long$occ, select(dets_long, spp, cell, year, replicate_id), identity)

# Occurrence Covariates --------------------------------------------------
occ.covs <- dets %>% 
  select(cell, DEM_max, p_forest, cliff_canyon, precip) %>% 
  distinct() %>% #Select covariates as list
  mutate(site.effect = 1:length(unique(dets$cell))) %>% #create variable for random site effect
  as.list() #covert each column to list element

occ.covs$year <- matrix(rep(unique(dets$year), each = length(occ.covs$cell)),
                        nrow = length(occ.covs$cell), ncol = length(unique(dets$year)))



# Detection Covariates (each covariate is a 3d array (cell, year, replicate)) ----------------------------------------------------
design.matrix.surveyed <- model.matrix(~ clutter + scale(tmin) + scale(daylight) + water_ind, data = dets)

det.covs <- list()
# det.covs[['clutter1']] <- tapply(design.matrix.surveyed[,2], select(dets, cell, year, replicate_id), identity)
# det.covs[['clutter2']] <- tapply(design.matrix.surveyed[,3], select(dets, cell, year, replicate_id), identity)
# det.covs[['clutter3']] <- tapply(design.matrix.surveyed[,4], select(dets, cell, year, replicate_id), identity)
# det.covs[['clutter4']] <- tapply(design.matrix.surveyed[,5], select(dets, cell, year, replicate_id), identity)
det.covs[['clutter']] <- tapply(design.matrix.surveyed[,2], select(dets, cell, year, replicate_id), identity)
det.covs[['tmin']] <- tapply(design.matrix.surveyed[,3], select(dets, cell, year, replicate_id), identity)
det.covs[['dayl']] <- tapply(design.matrix.surveyed[,4], select(dets, cell, year, replicate_id), identity)
det.covs[['water']] <- tapply(design.matrix.surveyed[,5], select(dets, cell, year, replicate_id), identity)

# Site coordinates (projected) ---------------------------------------------
sites.sp <- covars %>% 
  st_transform(crs = 26911) %>% 
  st_centroid() 
  

coords <- data.frame(CONUS_10KM = occ.covs$cell) %>% 
  left_join(sites.sp) %>% 
  st_as_sf() %>% 
  st_coordinates()

 
# Replicate Coordinates ---------------------------------------------------
dets.sp <- dets %>% 
  st_as_sf(coords = c("lon", "lat"), crs = "WGS84") %>% 
  st_transform(crs = 26911) 

write_sf(dets.sp, here("DataProcessed.nosync/spOccupancy/spatial_dets.shp"), )

# Store in a single object ------------------------------------------------
bat.dat <- list('y' = y,
                'occ.covs' = occ.covs,
                'det.covs' = det.covs,
                'coords' = coords)

saveRDS(bat.dat, here("DataProcessed.nosync/spOccupancy/bat_dat.rds"))


