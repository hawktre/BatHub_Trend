## ---------------------------
##
## Script name: 01b_occurrence_modprep.R
##
## Purpose of script: Format Occurrence-Level Data for Stan Modeling
##
## Author: Trent VanHawkins
##
## Date Created: 2024-06-13
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## ---------------------------


## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(sf)

# Read in the data --------------------------------------------------------
nw_grid_shp <- read_sf(here("DataProcessed/occurrence/batgrid_covars.shp"))

nw_nights <- readRDS(here("DataProcessed/detections/nw_nights.rds"))


## rename nw_grid_shape to have the name cell
nw_grid_shp <- nw_grid_shp %>% 
  rename("cell" = CONUS_10KM,
         "cliff_cover" = EVT_NAME)
# Figure out what su was surveyed each year -------------------------------
## create sample history
samp_hist <- nw_nights %>% 
  select(cell, year) %>%
  distinct() %>% 
  mutate(surveyed = 1) %>% 
  pivot_wider(id_cols = cell, names_from = year, values_from = surveyed, values_fill = 0, names_prefix = "samp_") %>% 
  mutate(samp_all = 1)

# left join to nw_grid_shp
nw_grid_shp <- left_join(nw_grid_shp, samp_hist, by = "cell") %>% 
  select(-ID)

#replace na with 0
nw_grid_shp[is.na(nw_grid_shp)] <- 0


# Format like wilson and save out -----------------------------------------
nw_grid_shp <- nw_grid_shp %>% 
  arrange(desc(samp_all), cell)
  
saveRDS(nw_grid_shp, here("DataProcessed/occurrence/nw_grid_shp.rds"))
