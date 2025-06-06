## ---------------------------
##
## Script name: table1.R
##
## Purpose of script: Create a table 1 for report
##
## Author: Trent VanHawkins
##
## Date Created: 2025-05-18
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
library(table1)
library(gtsummary)

covars <- readRDS(here("DataProcessed/occurrence/nw_grid_shp.rds"))
dets <- readRDS(here("DataProcessed/detections/nw_nights.rds"))


# Original summary
survey_summary <- covars %>%
  group_by(state) %>%
  summarise(
    sampled_2016 = sum(samp_2016),
    sampled_2017 = sum(samp_2017),
    sampled_2018 = sum(samp_2018),
    sampled_2019 = sum(samp_2019),
    sampled_2020 = sum(samp_2020),
    sampled_2021 = sum(samp_2021),
    sampled_2022 = sum(samp_2022),
    .groups = "drop"
  ) %>% 
  st_drop_geometry()

# Transpose to long, then back to wide (year as rows)
transposed_summary <- survey_summary %>%
  pivot_longer(
    cols = starts_with("sampled_"),
    names_to = "Year",
    names_prefix = "sampled_",
    values_to = "n_sampled"
  ) %>%
  pivot_wider(
    names_from = state,
    values_from = n_sampled
  ) %>%
  mutate(Total = rowSums(across(where(is.numeric))))

# View or output the final table
transposed_summary %>% 
  kableExtra::kable(format = "latex", booktabs = T, caption = "Number of grid cells surveyed by state and year",
                    label = "tbl-grids")
  
