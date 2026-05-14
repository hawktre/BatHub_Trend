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
library(janitor)
library(sf)

# Read in the data --------------------------------------------------------
nw_grid <- read_sf(here("DataProcessed/occurrence/batgrid_covars.shp")) |>
  clean_names()
nw_nights <- readRDS(here("DataProcessed/detections/nw_nights.rds"))


## rename nw_grid_shape to have the name cell
nw_grid <- nw_grid %>%
  st_set_agr("constant") %>%
  rename("sample_unit_id" = conus_10km, "cliff_canyon" = evt_name)
# Figure out what su was surveyed each year -------------------------------
## create sample history
samp_hist <- nw_nights %>%
  select(sample_unit_id, year) %>%
  distinct() %>%
  mutate(surveyed = 1) %>%
  pivot_wider(
    id_cols = sample_unit_id,
    names_from = year,
    values_from = surveyed,
    values_fill = 0,
    names_prefix = "sampled_"
  ) %>%
  mutate(samp_all = 1)

# left join to nw_grid_shp
nw_grid <- left_join(nw_grid, samp_hist, by = "sample_unit_id") |> select(-id)

#replace na with 0
nw_grid[is.na(nw_grid)] <- 0


# Format and save out -----------------------------------------
nw_grid <- nw_grid %>%
  arrange(desc(samp_all), sample_unit_id)

saveRDS(nw_grid, here("DataProcessed/occurrence/nw_grid.rds"))
