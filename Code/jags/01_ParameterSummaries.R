## ---------------------------
##
## Script name: 01_ParameterSummaries.R
##
## Purpose of script: Summarise posterior estimates from JAGS model fits
##
## Author: Trent VanHawkins
##
## Date Created: 2025-03-29
##
## ---------------------------

options(scipen = 6, digits = 4)

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(jagsUI)
library(MCMCvis)
library(sf)

# Load Data ---------------------------------------------------------------

## JAGS formatted data
occ_data <- readRDS(here("DataProcessed/results/jags/occ_data.rds"))

## Design matrices built in 00_jagsPrep.R
design_mats       <- readRDS(here("DataProcessed/results/jags/design_matrices.rds"))
xmata             <- design_mats$xmata
xmatc             <- design_mats$xmatc
design_matrix_det <- design_mats$design_matrix_det

## Site and year index keys built in 00_jagsPrep.R
index_keys <- readRDS(here("DataProcessed/results/jags/index_keys.rds"))
site_ids   <- index_keys$site_ids
year_ids   <- index_keys$year_ids

## Spatial info
nw_grid <- readRDS(here("DataProcessed/occurrence/nw_grid.rds"))

## Cliff-associated species
cliff_spp <- c("anpa", "euma", "myci", "pahe")

# Create Site/State Key ---------------------------------------------------

## Match site ordering used in the model to state labels
site_state_key <- nw_grid %>%
  select(sample_unit_id, state) %>%
  filter(sample_unit_id %in% site_ids) %>%
  st_drop_geometry() %>%
  arrange(sample_unit_id) %>%          # must match order used in tapply()
  mutate(site_index = row_number())

# Functions ---------------------------------------------------------------

## Compute psi-bar (mean occupancy) overall and by state, plus lambda trend
compute_psi_bar <- function(fit, site_info) {
  psi_array <- fit$sims.list$psi       # dim(n_iter, n_sites, n_years)
  
  n_iter  <- dim(psi_array)[1]
  n_sites <- dim(psi_array)[2]
  n_years <- dim(psi_array)[3]
  
  # Reshape to long format: one row per iter-site-year
  psi_df <- as.data.frame.table(psi_array, responseName = "psi") %>%
    rename(iter = Var1, site = Var2, year = Var3) %>%
    mutate(iter = as.integer(iter),
           site = as.integer(site),
           year = as.integer(year)) %>%
    left_join(site_info, by = c("site" = "site_index"))
  
  # Overall psi_bar per year
  psi_bar_all_iter <- psi_df %>%
    group_by(year, iter) %>%
    summarize(psi_bar = mean(psi), .groups = "drop")
  
  psi_bar_all <- psi_bar_all_iter %>%
    group_by(year) %>%
    summarize(mean = mean(psi_bar),
              lci  = quantile(psi_bar, 0.025),
              q25  = quantile(psi_bar, 0.25),
              q75  = quantile(psi_bar, 0.75),
              uci  = quantile(psi_bar, 0.975),
              .groups = "drop")
  
  # psi_bar by state per year
  psi_bar_state_iter <- psi_df %>%
    group_by(state, year, iter) %>%
    summarize(psi_bar = mean(psi), .groups = "drop")
  
  psi_bar_by_state <- psi_bar_state_iter %>%
    group_by(state, year) %>%
    summarize(mean = mean(psi_bar),
              lci  = quantile(psi_bar, 0.025),
              q25  = quantile(psi_bar, 0.25),
              q75  = quantile(psi_bar, 0.75),
              uci  = quantile(psi_bar, 0.975),
              .groups = "drop")
  
  # Overall lambda: ratio of psi_bar in last year vs first year
  lambda_tot <- psi_bar_all_iter %>%
    filter(year %in% c(1, n_years)) %>%
    pivot_wider(names_from = year, values_from = psi_bar,
                names_prefix = "year_") %>%
    mutate(lambda = .data[[paste0("year_", n_years)]] / year_1) %>%
    summarize(mean = mean(lambda),
              lci  = quantile(lambda, 0.025),
              q25  = quantile(lambda, 0.25),
              q75  = quantile(lambda, 0.75),
              uci  = quantile(lambda, 0.975))
  
  # Lambda by state
  lambda_tot_state <- psi_bar_state_iter %>%
    filter(year %in% c(1, n_years)) %>%
    pivot_wider(names_from = year, values_from = psi_bar,
                names_prefix = "year_") %>%
    mutate(lambda = .data[[paste0("year_", n_years)]] / year_1) %>%
    group_by(state) %>%
    summarize(mean = mean(lambda),
              lci  = quantile(lambda, 0.025),
              q25  = quantile(lambda, 0.25),
              q75  = quantile(lambda, 0.75),
              uci  = quantile(lambda, 0.975),
              .groups = "drop")
  
  list(
    overall  = list(psi_bar = psi_bar_all,    trend = lambda_tot),
    by_state = list(psi_bar = psi_bar_by_state, trend = lambda_tot_state)
  )
}

## Summarize occupancy and detection coefficients with names from design matrices
summarize_params <- function(fit, spp, xmat_use, pmat_use) {
  # Parameter names from design matrices
  occ_names <- c("alpha01", colnames(xmat_use))
  det_names <- colnames(pmat_use)
  all_names <- c(occ_names, det_names)
  
  # Extract posterior samples
  params <- fit$sims.list[c("alpha01", "alphas", "betas")]
  
  # Convert to long format
  param_list <- lapply(params, function(x) {
    if (is.null(dim(x))) {
      tibble(value = x)
    } else {
      as_tibble(x, .name_repair = "unique") %>%
        pivot_longer(cols = everything(), names_to = "index", values_to = "value")
    }
  })
  
  # Combine, summarize, and attach parameter names
  summary_df <- list_rbind(param_list, names_to = "param") %>%
    group_by(param, index) %>%
    summarize(
      mean = mean(value),
      lci  = quantile(value, 0.025),
      q25  = quantile(value, 0.25),
      q75  = quantile(value, 0.75),
      uci  = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(param = all_names) %>%
    select(-index)
  
  # Pull ESS and R-hat directly from the fit object
  diag_df <- tibble(
    param = all_names,
    ess   = unlist(fit$n.eff[c("alpha01", "alphas", "betas")]),
    rhat  = unlist(fit$Rhat[c("alpha01", "alphas", "betas")])
  )
  
  left_join(summary_df, diag_df, by = "param")
}

# Helper to run both summary functions for one species --------------------
summarize_species <- function(path, site_info, cliff_spp, xmata, xmatc, design_matrix_det) {
  fit  <- readRDS(path)
  spp  <- str_split_i(tools::file_path_sans_ext(basename(path)), "_", 1)
  cat("Summarizing:", spp, "\n")
  
  xmat_use <- if (spp %in% cliff_spp) xmatc else xmata
  
  list(
    psi    = compute_psi_bar(fit, site_info = site_info),
    params = summarize_params(fit, spp = spp,
                              xmat_use = xmat_use,
                              pmat_use = design_matrix_det)
  )
}

# Summarize Full Model Results --------------------------------------------
full_paths  <- list.files(here("DataProcessed/results/jags/full/fits/"),
                          full.names = TRUE)
full_spps   <- str_split_i(tools::file_path_sans_ext(basename(full_paths)), "_", 1)

full_res <- map(full_paths, summarize_species,
                site_info         = site_state_key,
                cliff_spp         = cliff_spp,
                xmata             = xmata,
                xmatc             = xmatc,
                design_matrix_det = design_matrix_det)

names(full_res) <- full_spps
saveRDS(full_res, here("DataProcessed/results/jags/full/res_summary.rds"))

# Summarize Sensitivity (OR/WA Only) Results ------------------------------
sens_paths <- list.files(here("DataProcessed/results/jags/ORWA_Only/fits/"),
                         full.names = TRUE)
sens_spps  <- str_split_i(tools::file_path_sans_ext(basename(sens_paths)), "_", 1)

sens_res <- map(sens_paths, summarize_species,
                site_info         = site_state_key,
                cliff_spp         = cliff_spp,
                xmata             = xmata,
                xmatc             = xmatc,
                design_matrix_det = design_matrix_det)

names(sens_res) <- sens_spps
saveRDS(sens_res, here("DataProcessed/results/jags/ORWA_Only/res_summary.rds"))