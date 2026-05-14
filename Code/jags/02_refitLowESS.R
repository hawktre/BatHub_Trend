## ---------------------------
##
## Script name: 02_refitLowESS.R
##
## Purpose of script: Identify and refit species with insufficient ESS
##                    from initial dynamic occupancy model fits, for both
##                    full and sensitivity analyses.
##                    Run after 02_TrendFigs.R.
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
library(jagsUI)

# Load Combined Results ---------------------------------------------------
all_res <- readRDS(here("DataProcessed/results/jags/all_res.rds"))

# Identify Low-ESS Species by Analysis ------------------------------------
low_ess <- all_res$params %>%
  filter(ess < 400) %>%
  distinct(species, analysis)

cat("Species to refit:\n")
print(low_ess)

# Refit Settings ----------------------------------------------------------
ni_refit <- 100000
nb_refit <- 20000
nt_refit <- 10
nc_refit <- 4

# Analysis Path Key -------------------------------------------------------
analysis_paths <- c(
  "OR|WA|ID"   = "DataProcessed/results/jags/full/fits/",
  "OR|WA Only" = "DataProcessed/results/jags/ORWA_Only/fits/"
)

# Refit -------------------------------------------------------------------
for (i in seq_len(nrow(low_ess))) {
  spp      <- low_ess$species[i]
  analysis <- low_ess$analysis[i]
  path     <- here(paste0(analysis_paths[analysis], spp, "_jagsfit.rds"))
  
  cat("\nRefitting:", spp, "(", analysis, ")\n")
  cat("Start time:", format(Sys.time()), "\n")
  
  fit     <- readRDS(path)
  fit_new <- update(fit,
                    n.iter   = ni_refit,
                    n.burnin = nb_refit,
                    n.thin   = nt_refit,
                    n.chains = nc_refit)
  
  saveRDS(fit_new, path)
  cat("Finished:", spp, "at", format(Sys.time()), "\n")
}

cat("\nAll species refit. Re-run 01_ParameterSummaries.R and 02_TrendFigs.R to update summaries.\n")