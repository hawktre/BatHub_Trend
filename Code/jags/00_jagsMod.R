## ---------------------------
##
## Script name: 00_jagsPrep.R
##
## Purpose of script: Perpare data for JAGS model
##
## Author: Trent VanHawkins
##
## Date Created: 2024-09-29
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
library(rjags)
library(jagsUI)

covars <- read_sf(here("DataProcessed/occurrence/batgrid_covars.shp"))
dets <- readRDS(here("DataProcessed/detections/nw_nights.rds"))

possible_bats <- c("laci",
                   "lano",
                   "myev",
                   "epfu",
                   "myyu",
                   "myth",
                   "myci",
                   "myvo",
                   "anpa",
                   "pahe",
                   "euma",
                   "myca",
                   "mylu",
                   "coto")
# Create Data Array -------------------------------------------------------

## Join Detections and Covariates
covars_join <- covars %>% 
  rename("cell" = CONUS_10KM,
         "cliff_canyon" = EVT_NAME) %>% 
  select(-c(lat, long)) %>%
  select(-riverlake) %>% 
  mutate(p_forest = log(p_forest + 1),
         cliff_canyon = log(cliff_canyon*100+1),
         across(karst:cliff_canyon, ~scale(.x)[,1])) 

## Drop geometry to make a data frame
covars_join <- st_drop_geometry(covars_join)

## Create out detections data frame
dets <- dets %>%
  left_join(covars_join, by = "cell") %>% 
  group_by(cell, year) %>%
  mutate(replicate_id = as.numeric(factor(replicate, levels = unique(replicate))),
         clutter = relevel(clutter, ref = "3")) %>%
  ungroup() %>% 
  drop_na(DEM_max, p_forest, precip, cliff_canyon) %>% 
  arrange(year, cell, replicate_id) 

## Pivot Longer to prep for array
dets_long <- dets %>%
  pivot_longer(cols = all_of(possible_bats), names_to = "spp", values_to = "occ") %>% 
  select(spp, cell, year, replicate_id, occ) %>% 
  arrange(spp, year, cell, replicate_id) %>% 
  mutate(across(1:4, as.factor)) 

## use tapply to create array with dim(spp,cell,year,replicate)
y <- tapply(dets_long$occ, select(dets_long, spp, cell, replicate_id, year), identity)

# Detection Covariates (each covariate is a 3d array (cell, replicate, year)) ----------------------------------------------------
design.matrix.surveyed <- model.matrix(~ clutter + scale(tmin) + scale(daylight) + water_ind, data = dets)

det.covs <- list()
det.covs[['clutter0']] <- tapply(design.matrix.surveyed[,2], select(dets, cell, replicate_id, year), identity)
det.covs[['clutter1']] <- tapply(design.matrix.surveyed[,3], select(dets, cell, replicate_id, year), identity)
det.covs[['clutter2']] <- tapply(design.matrix.surveyed[,4], select(dets, cell, replicate_id, year), identity)
det.covs[['clutter3']] <- tapply(design.matrix.surveyed[,5], select(dets, cell, replicate_id, year), identity)
det.covs[['tmin']] <- tapply(design.matrix.surveyed[,6], select(dets, cell, replicate_id, year), identity)
det.covs[['dayl']] <- tapply(design.matrix.surveyed[,7], select(dets, cell, replicate_id, year), identity)
det.covs[['water']] <- tapply(design.matrix.surveyed[,8], select(dets, cell, replicate_id, year), identity)


# Occurrence Covariates ---------------------------------------------------
nw_grid_shp <- readRDS(here("DataProcessed/occurrence/nw_grid_shp.rds"))

#divide shapefile into sampled and unsampled grid cells.
#also arrange these by conus_id
nw_grida <- nw_grid_shp %>%
  filter(samp_all == 1) %>%
  arrange(cell)
nw_gridb <- nw_grid_shp %>%
  filter(samp_all == 0) %>%
  arrange(cell)

#combine these back again, reordered. calculate scaled versions of each
#covariate and save these in matrices that will be used for predictions
nw_grid_all <- rbind(nw_grida, nw_gridb)
xmat_all <- nw_grid_all %>%
  st_drop_geometry() %>%
  mutate(log_fc = log(p_forest + 1)) %>% 
  select(log_fc, precip, DEM_max) %>%
  scale()

xmata <- xmat_all[which(nw_grid_all$samp_all == 1), ]
xmatb <- xmat_all[which(nw_grid_all$samp_all == 0), ]

#Repeat for cliff species design matrix
xmat_cliff <- nw_grid_all %>%
  st_drop_geometry() %>%
  mutate(log_fc = log(p_forest + 1),
         log_cliff = log(cliff_cover*100 + 1)) %>% 
  select(log_fc, precip, DEM_max, log_cliff) %>%
  scale()
xmatc <- xmat_cliff[which(nw_grid_all$samp_all == 1), ]
xmatd <- xmat_cliff[which(nw_grid_all$samp_all == 0), ]
# Setup remaining Stan data. ----------------------------------------------
#Need to add a few things to separate the different years and
#keep track of which each index corresponds to.
n_sites_total <- nrow(xmata)
n_years <- length(unique(dets$year))

#Number of covariates for occupancy and detection
n_xcovs <- ncol(xmata)
n_xcovs_c <- ncol(xmatc)

#Number of Visits
#Arrange nights by conus_id, then year
nw_nights_all <- dets %>%
  arrange(cell, year) %>% 
  filter(!is.na(state)) # filter out california data

## Detection Histories ##
#Number of visits per grid cell each year
n_visits_dfa <- nw_nights_all %>%
  group_by(cell, year) %>%
  summarise(n_visits = length(date))

#Need to fill in the zeros for this
n_visits_dfb <- pivot_wider(n_visits_dfa,
                            id_cols = cell,
                            names_from = year,
                            values_from = n_visits,
                            values_fill = list(n_visits = 0))

#reorder columns
n_visits_dfb <- select(n_visits_dfb, cell, '2016', '2017', 
                       '2018', '2019', '2020',
                       '2021','2022')

n_visits_jags <- n_visits_dfb %>%
  ungroup() %>%
  select(-cell) %>%
  simplify2array()

#Get naive occupancy, appropriately fill in 0 for years without visits
## Initialize an empty list for all species data
occ_data <- list()

for (i in possible_bats) {
  if (i %in% c("anpa", "euma", "myci", "pahe")) {
    occ_data[[i]] <- list('dets' = y[i,,,],
                          'tmin' = det.covs$tmin,
                          'dayl' = det.covs$dayl,
                          'clut0' = det.covs$clutter0,
                          'clut1' = det.covs$clutter1,
                          'clut2' = det.covs$clutter2,
                          'clut3' = det.covs$clutter3,
                          'wind' = det.covs$water,
                          'n_sites' = n_sites_total,
                          'xmat' = xmatc,
                          'n_xcovs' = ncol(xmatc),
                          'n_visits' = n_visits_jags,
                          'n_years' = n_years)
  }
  else{occ_data[[i]] <- list('dets' = y[i,,,],
                             'tmin' = det.covs$tmin,
                             'dayl' = det.covs$dayl,
                             'clut0' = det.covs$clutter0,
                             'clut1' = det.covs$clutter1,
                             'clut2' = det.covs$clutter2,
                             'clut3' = det.covs$clutter3,
                             'wind' = det.covs$water,
                             'n_sites' = n_sites_total,
                             'xmat' = xmata,
                             'n_xcovs' = ncol(xmata),
                             'n_visits' = n_visits_jags,
                             'n_years' = n_years)}
  
}

saveRDS(occ_data, here("DataProcessed/results/jags/occ_data.rds"))

#Specify initial values for z
max2 <- function(x){
  if(sum(is.na(x)) == length(x)){
    return(0)
  } else {
    return(max(x, na.rm = TRUE))
  }
}

#Fit the model for each species and output the results
# Loop over species (only 1 here)
for (i in seq(1:length(occ_data))) {
  # Subset data
  tmp <- occ_data[[i]]
  
  # Print species name
  print(names(occ_data)[i])
  print(Sys.time())
  
  # Set initial values for z
  inits <- function() {
    list(z = apply(tmp$dets, c(1, 3), max2))
  }
  
  # Run model with jagsUI
  occ_jags <- jags(data = tmp,
                   inits = inits,
                   parameters.to.save = c('alpha01', 'alphas',
                                          'beta0', 'beta1', 'beta2',
                                          'beta3', 'beta4', 'beta5', 
                                          'beta6', 'beta7', 'psi'),
                   model.file = here("Code/jags/occ_model_royle_onlypsi.jags"),
                   n.chains = 4,
                   n.iter = 10000,
                   n.burnin = 2000,
                   n.thin = 1,
                   parallel = TRUE)  # optional for speed
  
  # Save model output
  saveRDS(occ_jags, here(paste0("DataProcessed/results/jags/full/fits/", names(occ_data)[i], "_jagsfit.rds")))
}
