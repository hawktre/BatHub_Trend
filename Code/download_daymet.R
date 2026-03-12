## ---------------------------
##
## Script name: download_daymet.r
##
## Purpose of script:
##
## Author: Trent VanHawkins
##
## Date Created: 2025-11-15
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
library(daymetr)

#Download the daymet data (optional) ------------------------------
df_batch <- download_daymet_batch(file_location = here("DataRaw/covariates/daymet/daymet_batch.csv"),
                                  start = 2016,
                                  end = 2024,
                                  internal = T,
                                  simplify = T, 
                                  silent = F)

saveRDS(df_batch, here("DataRaw/covariates/daymet/all_daymet.rds"))
