## ---------------------------
##
## Script name: 00_format_spOcc.R
##
## Purpose of script: Format data specifically for fitting models with spOcc
##
## Author: Trent VanHawkins
##
## Date Created: 2024-08-13
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

covars <- read_sf(here("DataProcessed/occurrence/batgrid_covars_spocc.shp"))
dets <- readRDS(here("DataProcessed/detections/nw_nights.rds"))
possible_bats <- c("laci", "lano", "mylu")

# Create Data Array -------------------------------------------------------

## Format occupancy covariates and filter to OR/WA Only
covars_join <- covars %>% 
  filter(precip > 0,
         state %in% c("Oregon", "Washington")) %>% 
  rename("cell" = CONUS_10KM) %>% 
  select(cell, GRTS_ID, state, forest, cliff, elev, precip, rough, physio) %>% 
  mutate(forest = log(forest + 1),
         across(cliff:physio, ~scale(.x)[,1])) #mean-center and scale covars

plot(covars[c('forest', 'cliff', 'elev', 'precip', 'rough', 'physio')])
plot(covars_join[c('forest', 'cliff', 'elev', 'precip', 'rough', 'physio')])


# Join covariates with detections -----------------------------------------
covars_join <- st_drop_geometry(covars_join)

dets <- dets %>%
  left_join(covars_join, by = "cell") %>%
  filter(year <= 2018) %>% 
  # mutate(clutter = factor(clutter, levels = c(-1, 0, 1, 2, 3), labels = c(0, 1, 2, 3, 4))) %>% 
  group_by(cell, year) %>%
  mutate(replicate_id = as.numeric(factor(replicate, levels = unique(replicate)))) %>%
  ungroup() %>% 
  drop_na(elev, forest, cliff, precip, rough) %>% 
  arrange(year, cell, replicate_id) 

dets.plt <- dets %>% 
  pivot_longer(cols = all_of(possible_bats), names_to = "spp", values_to = "det") %>% 
  st_as_sf(coords = c("lon", "lat"), crs = "WGS84") %>% 
  filter(det == 1)

ggplot()+
  geom_sf(data = dets.plt, aes(color = spp))+
  facet_wrap(~spp)

pairs(dets %>% select(elev, cliff, forest, precip))
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
  select(cell, elev, forest, precip, rough) %>% 
  distinct() %>% #Select covariates as list
  as.list() #covert each column to list element

occ.covs$year <- matrix(rep(unique(dets$year), each = length(occ.covs$cell)),
                        nrow = length(occ.covs$cell), ncol = length(unique(dets$year)))



# Detection Covariates (each covariate is a 3d array (cell, year, replicate)) ----------------------------------------------------
design.matrix.surveyed <- model.matrix(~ clutter + tmin + daylight + water_ind, data = dets)

det.covs <- list()
# det.covs[['date']] <- tapply(design.matrix.surveyed[,2], select(dets, cell, year, replicate_id), identity)
det.covs[['clutter']] <- tapply(design.matrix.surveyed[,2], select(dets, cell, year, replicate_id), identity)
det.covs[['tmin']] <- tapply(design.matrix.surveyed[,3], select(dets, cell, year, replicate_id), identity)
det.covs[['dayl']] <- tapply(design.matrix.surveyed[,4], select(dets, cell, year, replicate_id), identity)
det.covs[['water']] <- tapply(design.matrix.surveyed[,5], select(dets, cell, year, replicate_id), identity)
# Store in a single object ------------------------------------------------
bat.dat <- list('y' = y,
                'occ.covs' = occ.covs,
                'det.covs' = det.covs)

saveRDS(bat.dat, here("DataProcessed/results/spOccupancy/laci_redo_2019/bat_dat.rds"))
st_write(covars_join, here("DataProcessed/results/spOccupancy/laci_redo_2019/covars_scaled.shp"))

