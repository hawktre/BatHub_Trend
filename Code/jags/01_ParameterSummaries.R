## ---------------------------
##
## Script name: 01_ParameterSummaries.R
##
## Purpose of script: Summarise posterior estimates from jags model fits
##
## Author: Trent VanHawkins
##
## Date Created: 2025-03-29
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
library(jagsUI)
library(MCMCvis)
library(sf)

# Read in Data--------------------------------------------------------

## Read in data objects
### Jags formatted
occ_data <- readRDS(here("DataProcessed/results/jags/occ_data.rds"))

### Satial Info
nw_grid <- readRDS(here("DataProcessed/occurrence/nw_grid_shp.rds"))

## Cliff Bats
cliff_bats <- c("anpa", "euma", "myci", "pahe")

# Create a site/state key -------------------------------------------------
## Get the ordering of sites in the analysis
sites <- rownames(occ_data$laci$dets) %>% as.numeric()

## Create a key
site_state_key <- nw_grid %>% 
  select(cell, state) %>% 
  filter(cell %in% sites) %>% 
  st_drop_geometry() %>% 
  mutate(site_index = row_number())


# Functions to summarise psi and params ---------------------------------------------------------

## Function to compute psi-bar (entire area and by state)
compute_psi_bar <- function(fit, site_info) {
  psi_array <- fit$sims.list$psi  # Dimensions: [iter, site, year]
  
  # Get array dimensions
  n_iter <- dim(psi_array)[1]
  n_sites <- dim(psi_array)[2]
  n_years <- dim(psi_array)[3]
  
  # Reshape to long format: each row is an iter-site-year triplet
  psi_df <- as.data.frame.table(psi_array, responseName = "psi") %>%
    rename(iter = Var1, site = Var2, year = Var3) %>%
    mutate(iter = as.integer(iter),
           site = as.integer(site),
           year = as.integer(year)) %>%
    left_join(site_info, by = c("site" = "site_index"))  # Ensure site_info has site_index and state
  
  # --- Overall psi_bar ---
  psi_bar_all_iter <- psi_df %>%
    group_by(year, iter) %>%
    summarize(psi_bar = mean(psi), .groups = "drop")
  
  psi_bar_all <- psi_bar_all_iter %>%
    group_by(year) %>%
    summarize(mean = mean(psi_bar),
              lci = quantile(psi_bar, 0.025),
              q25 = quantile(psi_bar, 0.25),
              q75 = quantile(psi_bar, 0.75),
              uci = quantile(psi_bar, 0.975),
              .groups = "drop")
  
  # --- psi_bar by state ---
  psi_bar_state_iter <- psi_df %>%
    group_by(state, year, iter) %>%
    summarize(psi_bar = mean(psi), .groups = "drop")
  
  psi_bar_by_state <- psi_bar_state_iter %>% 
    group_by(state, year) %>%
    summarize(mean = mean(psi_bar),
              lci = quantile(psi_bar, 0.025),
              q25 = quantile(psi_bar, 0.25),
              q75 = quantile(psi_bar, 0.75),
              uci = quantile(psi_bar, 0.975), .groups = "drop")
  
  # --- lambda_tot (overall trend) ---
  psi_start <- psi_bar_all_iter %>% filter(year == 1)
  psi_end   <- psi_bar_all_iter %>% filter(year == n_years)
  
  lambda_df <- psi_start %>%
    inner_join(psi_end, by = "iter", suffix = c("_start", "_end")) %>%
    mutate(lambda_ratio = psi_bar_end / psi_bar_start)
  
  lambda_tot <- lambda_df %>%
    summarize(mean = mean(lambda_ratio),
              lci = quantile(lambda_ratio, 0.025),
              q25 = quantile(lambda_ratio, 0.25),
              q75 = quantile(lambda_ratio, 0.75),
              uci = quantile(lambda_ratio, 0.975))
  
  # --- lambda_tot (by state) ---
  psi_start <- psi_bar_state_iter %>% filter(year == 1)
  psi_end   <- psi_bar_state_iter %>% filter(year == n_years)
  
  lambda_df_state <- psi_start %>%
    inner_join(psi_end, by = c("state", "iter"), suffix = c("_start", "_end")) %>%
    mutate(lambda_ratio = psi_bar_end / psi_bar_start)
  
  lambda_tot_state <- lambda_df_state %>%
    group_by(state) %>% 
    summarize(mean = mean(lambda_ratio),
              lci = quantile(lambda_ratio, 0.025),
              q25 = quantile(lambda_ratio, 0.25),
              q75 = quantile(lambda_ratio, 0.75),
              uci = quantile(lambda_ratio, 0.975),
              .groups = "drop")
  
  list(overall = list(psi_bar = psi_bar_all,
                      trend = lambda_tot), 
       by_state = list(psi_bar = psi_bar_by_state,
                       trend = lambda_tot_state))
}

## Function to summarise other params

summarize_params <- function(fit, spp) {
  # Extract monitored alphas and betas
  params <- fit$sims.list[c("alpha01","alphas", "beta0", "beta1", "beta2",
                            "beta3", "beta4", "beta5", "beta6", "beta7")]
  
  # Convert each param to a data frame (named list of tibbles)
  param_list <- lapply(params, function(x) {
    if (is.null(dim(x))) {
      # Scalar: single column tibble
      tibble(value = x)
    } else {
      # Vector/matrix: convert to long format
      as_tibble(x, .name_repair = "unique") |> 
        pivot_longer(cols = everything(), names_to = "index", values_to = "value")
    }
  })
  
  # Combine with parameter name using list_rbind
  summary_df <- list_rbind(param_list, names_to = "param") %>% 
    group_by(param, index) %>% 
    summarize(mean = mean(value),
              lci = quantile(value, 0.025),
              q25 = quantile(value, 0.25),
              q75 = quantile(value, 0.75),
              uci = quantile(value, 0.975),
              .groups = "drop")
    
    
    if(spp %in% cliff_bats) {
      param_names <- c("alpha01", "log_fc", "precip", "elev", "cliff_cover",
                       "beta0", "tmin", "dayl", "clut0", "clut1", "clut2", "clut3",
                       "water")
      
    }else{
      param_names <- c("alpha01", "log_fc", "precip", "elev",
                       "beta0", "tmin", "dayl", "clut0", "clut1", "clut2", "clut3",
                       "water")
    }
    
  summary_df <- summary_df %>% mutate(param = param_names) %>% select(-index)
}


# Loop through the files --------------------------------------------------
## Full 
### Set up file paths
dir <- "DataProcessed/results/jags/full/fits/"
filenames <- list.files(here(dir))
fit_paths <- paste0(dir, filenames)
spps <- str_split_i(filenames, "_", 1)

full_res <- map(fit_paths, function(path) {
  fit <- readRDS(path)
  spp <- str_split_i(tools::file_path_sans_ext(basename(path)), "_", 1)
  psi_summary <- compute_psi_bar(fit, site_info = site_state_key)
  param_summary <- summarize_params(fit, spp = spp)
  list(psi = psi_summary,
       params = param_summary)
})

names(full_res) <- spps

saveRDS(full_res, here("DataProcessed/results/jags/full/res_summary.rds"))

## Sensitivity
### Set up file paths
dir <- "DataProcessed/results/jags/ORWA_Only/fits/"
filenames <- list.files(here(dir))
fit_paths <- paste0(dir, filenames)
spps <- str_split_i(filenames, "_", 1)

sens_res <- map(fit_paths, function(path) {
  fit <- readRDS(path)
  spp <- str_split_i(tools::file_path_sans_ext(basename(path)), "_", 1)
  print(spp)
  psi_summary <- compute_psi_bar(fit, site_info = site_state_key)
  param_summary <- summarize_params(fit, spp = spp)
  list(psi = psi_summary,
       params = param_summary)
})

names(sens_res) <- spps

saveRDS(sens_res, here("DataProcessed/results/jags/ORWA_Only/res_summary.rds"))
