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

require(tidyverse)
require(here)
require(sf)

covars <- read_sf(here("DataProcessed.nosync/occurrence/batgrid_covars.shp"))
dets <- readRDS(here("DataProcessed.nosync/detections/nw_nights.rds"))
