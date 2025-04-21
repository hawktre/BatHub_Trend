## ---------------------------
##
## Script name: 02_TrendFigs.R
##
## Purpose of script: Creating Figures and plots for trend summaries
##
## Author: Trent VanHawkins
##
## Date Created: 2025-03-31
##
##
## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(jagsUI)
library(sf)


# Read in the data --------------------------------------------------------

## Full Results

full <- readRDS(here("DataProcessed/results/jags/full/res_summary.rds"))

## Sensitivity Results
sensitivity  <- readRDS(here("DataProcessed/results/jags/ORWA_Only/res_summary.rds"))



# Concatenate Results for Plotting ----------------------------------------

combine_species_summaries <- function(full) {
  psi_bar_overall_list <- list()
  trend_overall_list <- list()
  psi_bar_bystate_list <- list()
  trend_bystate_list <- list()
  params_list <- list()
  
  for (species in names(full)) {
    species_data <- full[[species]]
    
    # overall$psi_bar
    if (!is.null(species_data$psi$overall$psi_bar)) {
      df <- species_data$psi$overall$psi_bar
      df$species <- species
      psi_bar_overall_list[[species]] <- df
    }
    
    # overall$trend
    if (!is.null(species_data$psi$overall$trend)) {
      df <- species_data$psi$overall$trend
      df$species <- species
      trend_overall_list[[species]] <- df
    }
    
    # by_state$psi_bar
    if (!is.null(species_data$psi$by_state$psi_bar)) {
      df <- species_data$psi$by_state$psi_bar
      df$species <- species
      psi_bar_bystate_list[[species]] <- df
    }
    
    # by_state$trend
    if (!is.null(species_data$psi$by_state$trend)) {
      df <- species_data$psi$by_state$trend
      df$species <- species
      trend_bystate_list[[species]] <- df
    }
    
    # params
    if (!is.null(species_data$params)) {
      df <- species_data$params
      df$species <- species
      params_list[[species]] <- df
    }
  }
  
  list(
    psi_bar_overall = dplyr::bind_rows(psi_bar_overall_list),
    trend_overall   = dplyr::bind_rows(trend_overall_list),
    psi_bar_bystate = dplyr::bind_rows(psi_bar_bystate_list),
    trend_bystate   = dplyr::bind_rows(trend_bystate_list),
    params          = dplyr::bind_rows(params_list)
  )
}

# Extract Summaries -------------------------------------------------------
## Full
summary_full <- combine_species_summaries(full)

### Psi-bar
psi_full <- summary_full$psi_bar_overall %>% mutate("analysis" = "OR|WA|ID")
psi_full_bystate <- summary_full$psi_bar_bystate %>% mutate("analysis" = "OR|WA|ID")

### Trend
trend_full <-  summary_full$trend_overall %>% mutate("analysis" = "OR|WA|ID")
trend_full_bystate <-  summary_full$trend_bystate %>% mutate("analysis" = "OR|WA|ID")

### Params
params_full <- summary_full$params %>% mutate("analysis" = "OR|WA|ID")

## Sensitivity
summary_sens <- combine_species_summaries(sensitivity)

### Psi-bar
psi_sens <- summary_sens$psi_bar_overall %>% mutate("analysis" = "OR|WA Only")
psi_sens_bystate <- summary_sens$psi_bar_bystate %>% mutate("analysis" = "OR|WA Only")

### Trend
trend_sens <-  summary_sens$trend_overall %>% mutate("analysis" = "OR|WA Only")
trend_sens_bystate <-  summary_sens$trend_bystate %>% mutate("analysis" = "OR|WA Only")

### Params
params_sens <- summary_sens$params %>% mutate("analysis" = "OR|WA Only")


# Create mega results -----------------------------------------------------

all_res <- list(psi = bind_rows(psi_full, psi_sens) %>% mutate(year = as.integer(year) + 2015),
                psi_bystate = bind_rows(psi_full_bystate, psi_sens_bystate) %>% mutate(year = year + 2015),
                trend = bind_rows(trend_full, trend_sens),
                trend_bystate = bind_rows(trend_full_bystate, trend_sens_bystate),
                params = bind_rows(params_full, params_sens))

saveRDS(all_res, here("DataProcessed/results/jags/all_res.rds"))