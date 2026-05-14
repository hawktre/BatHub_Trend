## ---------------------------
##
## Script name: 04_MoransI.R
##
## Purpose of script: Check for residual spatial autocorrelation in
##                    site-level mean psi estimates using Moran's I.
##                    Run after 01_ParameterSummaries.R.
##
## Author: Trent VanHawkins
##
## Date Created: 2026-04-05
##
## ---------------------------

options(scipen = 6, digits = 4)

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(sf)
library(ape)

# Load Data ---------------------------------------------------------------
nw_grid    <- readRDS(here("DataProcessed/occurrence/nw_grid.rds"))
index_keys <- readRDS(here("DataProcessed/results/jags/index_keys.rds"))
site_ids   <- index_keys$site_ids

# Build Site Coordinate Matrix --------------------------------------------

## Centroids of sampled grid cells in the order used by the model
coords <- nw_grid %>%
  filter(sample_unit_id %in% site_ids) %>%
  arrange(sample_unit_id) %>%
  st_centroid() %>%
  st_coordinates()

## Inverse distance weight matrix for full analysis
dists     <- as.matrix(dist(coords))
inv_dists <- 1 / dists
diag(inv_dists) <- 0

# Moran's I Function ------------------------------------------------------

## Reads raw fit file, extracts site-level mean psi, runs Moran's I
run_morans <- function(path, inv_dists_mat) {
  fit <- readRDS(path)
  spp <- str_split_i(tools::file_path_sans_ext(basename(path)), "_", 1)
  cat("Running Moran's I for:", spp, "\n")
  
  psi_array <- fit$sims.list$psi    # dim(n_iter, n_sites, n_years)
  psi_site  <- apply(psi_array, 2, mean)
  
  result <- Moran.I(psi_site, inv_dists_mat)
  
  tibble(
    species  = spp,
    observed = result$observed,
    expected = result$expected,
    sd       = result$sd,
    p_value  = result$p.value
  )
}

# Run Moran's I for Full Analysis -----------------------------------------
full_paths <- list.files(here("DataProcessed/results/jags/full/fits/"),
                         full.names = TRUE)

morans_full <- map_dfr(full_paths, run_morans, inv_dists_mat = inv_dists) %>%
  mutate(analysis = "OR|WA|ID")

# Run Moran's I for Sensitivity Analysis ----------------------------------
sens_paths <- list.files(here("DataProcessed/results/jags/ORWA_Only/fits/"),
                         full.names = TRUE)

## Check if sensitivity analysis has a different number of sites
## and rebuild the weight matrix if so
sens_fit     <- readRDS(sens_paths[1])
n_sites_sens <- dim(sens_fit$sims.list$psi)[2]

if (n_sites_sens != nrow(coords)) {
  cat("Sensitivity analysis has different site set - rebuilding weight matrix\n")
  
  coords_sens <- nw_grid %>%
    filter(sample_unit_id %in% site_ids, state != "Idaho") %>%
    arrange(sample_unit_id) %>%
    st_centroid() %>%
    st_coordinates()
  
  dists_sens     <- as.matrix(dist(coords_sens))
  inv_dists_sens <- 1 / dists_sens
  diag(inv_dists_sens) <- 0
} else {
  inv_dists_sens <- inv_dists
}

morans_sens <- map_dfr(sens_paths, run_morans, inv_dists_mat = inv_dists_sens) %>%
  mutate(analysis = "OR|WA Only")

# Combine and Summarize ---------------------------------------------------
morans_all <- bind_rows(morans_full, morans_sens) %>%
  mutate(significant = p_value < 0.05) %>%
  arrange(analysis, p_value)

cat("\n--- Moran's I Results ---\n")
print(morans_all, n = Inf)

cat("\nSpecies with significant spatial autocorrelation (p < 0.05):\n")
morans_all %>%
  filter(significant) %>%
  select(species, analysis, observed, p_value) %>%
  print(n = Inf)

# Save Results ------------------------------------------------------------
saveRDS(morans_all, here("DataProcessed/results/jags/morans_results.rds"))