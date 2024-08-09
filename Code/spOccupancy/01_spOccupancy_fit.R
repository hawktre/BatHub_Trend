## ---------------------------
##
## Script name: 01_spOccupancy_fit.R
##
## Purpose of script: Fit models in spOccupancy
##
## Author: Trent VanHawkins
##
## Date Created: 2024-07-22
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

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(sf)
library(spOccupancy)


# Read in the data --------------------------------------------------------
bat.dat <- readRDS(here("DataProcessed.nosync/spOccupancy/bat_dat.rds"))
covars <- read_sf(here("DataProcessed.nosync/occurrence/batgrid_covars.shp"))


# Define formulas ---------------------------------------------------------

## Detection Formula (always the same)
det.formula <- ~ clutter1 + clutter2 + clutter3 + clutter4 + tmin + dayl + water

fit.spOcc <- function(dat, occ.formula, det.formula){
  # Specify inits -----------------------------------------------------------
  z.inits <- apply(dat$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
  bat.inits <- list(beta = 0,         # occurrence coefficients
                    alpha = 0,        # detection coefficients
                    sigma.sq.psi = 1, # occurrence random effect variances
                    z = z.inits)      # latent occurrence values

  # Specify Priors ----------------------------------------------------------
  bat.priors <- list(beta.normal = list(mean = 0, var = 2.72), 
                     alpha.normal = list(mean = 0, var = 2.72), 
                     sigma.sq.psi.ig = list(a = 0.1, b = 0.1))
  

  # specify mcmc settings ---------------------------------------------------
  n.chains <- 3
  n.batch <- 600
  batch.length <- 50
  n.samples <- n.batch * batch.length
  n.burn <- 2000
  n.thin <- 8

  # Approx. run time: ~ 1.3 min
  fit <- tPGOcc(occ.formula = occ.formula,
                     det.formula = det.formula,
                     data = dat,
                     n.batch = n.batch,
                     batch.length = batch.length,
                     inits = bat.inits,
                     priors = bat.priors,
                     ar1 = T,
                     n.burn = n.burn,
                     n.thin = n.thin,
                     n.chains = n.chains,
                     n.report = 200)
  return(fit)
}

for(s in 1:length(dimnames(bat.dat$y)[[1]])){
  print(paste0("You are currently running the model for", dimnames(bat.dat$y)[[1]][s]))
  #Subset the array for species s and change name
  spp.dat <- bat.dat
  spp.dat$y <- spp.dat$y[dimnames(bat.dat$y)[[1]][s],,,]

  if(dimnames(bat.dat$y)[[1]][s] %in% c("anpa", "euma", "myci", "pahe")){
    occ.formula <- ~ scale(year) + p_forest + precip + cliff_canyon + DEM_max
  }
  else{
    occ.formula <- ~ scale(year) + p_forest + precip + DEM_max
  }
  
  ## Detection Formula (always the same)
  det.formula <- ~ clutter1 + clutter2 + clutter3 + clutter4 + tmin + dayl + water
  
  #Fit the model and return the result
  fit <- fit.spOcc(dat = spp.dat, occ.formula = occ.formula, det.formula = det.formula)
  
  print(summary(fit))
  
  filename <- paste0(dimnames(bat.dat$y)[[1]][s], "_fit.rds")
  saveRDS(fit, here(paste0("DataProcessed.nosync/spOccupancy/fits/", filename)))
}


