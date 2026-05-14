## ---------------------------
##
## Script name: 01a_modelprep.R
##
## Purpose of script: Prepare data for modeling steps for single-species, single-season
##
## Author: Trent VanHawkins
##
## Date Created: 2024-05-08
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
library(janitor)
library(sf)


# Read in the data --------------------------------------------------------
detections <- readRDS(here("DataProcessed/detections/detection_histories.rds"))
daymet <- read.csv(here("DataRaw/covariates/daymet/daymet_output.csv"))


# Join them --------------------------------------------------------------
detections <- left_join(
  detections,
  daymet |> mutate(night = date(night)),
  by = c("location_name", "night")
) |>
  clean_names()

# Reformat --------------------------------------------------------------------

# Spatial Replicates ------------------------------------------------------

## Visualize how many spatial replicates in each cell
detections %>%
  group_by(sample_unit_id, year) %>%
  summarise(spatial_reps = n()) %>%
  ungroup() %>%
  ggplot(aes(x = spatial_reps)) +
  geom_bar() +
  facet_wrap(~year) +
  theme_classic()

## what is the range of spatial replicates?
spatial_reps_summary <- detections %>%
  group_by(sample_unit_id, year) %>%
  summarise(spatial_reps = n())

range(spatial_reps_summary$spatial_reps)
# Find covariate NAs ------------------------------------------------------
## Find rows with NA values
na_rows <- detections %>%
  filter(!complete.cases(.))

## only 4 few rows. drop them
detections <- detections %>%
  filter(complete.cases(.))


saveRDS(detections, here("DataProcessed/detections/nw_nights.rds"))
