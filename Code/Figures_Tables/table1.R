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
## ---------------------------

options(scipen = 6, digits = 4)

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(sf)
library(table1)
library(gtsummary)

# Load Data ---------------------------------------------------------------
covars <- readRDS(here("DataProcessed/occurrence/nw_grid.rds"))
dets   <- readRDS(here("DataProcessed/detections/nw_nights.rds"))

# Survey Summary ----------------------------------------------------------

## Pivot long on any sampled_YYYY columns, drop samp_all, summarize by state and year
survey_summary <- covars %>%
  st_drop_geometry() %>%
  pivot_longer(
    cols         = starts_with("sampled_"),
    names_to     = "Year",
    names_prefix = "sampled_",
    values_to    = "sampled"
  ) %>%
  group_by(state, Year) %>%
  summarise(n_sampled = sum(sampled), .groups = "drop") %>%
  pivot_wider(
    names_from  = state,
    values_from = n_sampled
  ) %>%
  mutate(Total = rowSums(across(where(is.numeric))))

# Output Table ------------------------------------------------------------
survey_summary %>%
  kableExtra::kable(
    format   = "latex",
    booktabs = TRUE,
    caption  = "Number of grid cells surveyed by state and year",
    label    = "tbl-grids"
  )