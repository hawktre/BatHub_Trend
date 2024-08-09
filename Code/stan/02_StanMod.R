## ---------------------------
##
## Script name: 02_StanMod.R
##
## Purpose of script: Format model matrices and run stan model 
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

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(sf)
library(rstan)

# Read in the data --------------------------------------------------------

## Detections 
nw_nights <- readRDS(here("DataProcessed.nosync/detections/nw_nights.rds"))

## Occurrence Covariates
nw_grid_shp <- readRDS(here("DataProcessed.nosync/occurrence/nw_grid_shp.rds"))


possible_bats <- c("laci",
                   "lano",
                   "myev",
                   "epfu",
                   "myyu",
                   "myth",
                   "myci",
                   "myvo",
                   "tabr",
                   "anpa",
                   "pahe",
                   "euma",
                   "myca",
                   "mylu",
                   "coto")
# Detections and Nightly Covariates ---------------------------------------
## Nightly Covariates ##
#Join nw_nights to find state
state_key <- nw_grid_shp %>% 
  select(cell, state)

nw_nights <- left_join(nw_nights, state_key, by = "cell")

#Arrange nights by conus_id, then year
nw_nights_all <- nw_nights %>%
  arrange(cell, year) %>% 
  filter(!is.na(state)) # filter out california data

nw_nights_spat <- nw_nights_all %>% 
  st_as_sf(coords = c('lon','lat'), crs = st_crs(nw_grid_shp))

ggplot()+
  geom_sf(data = nw_grid_shp, aes(fill = factor(samp_2018)))+
  geom_sf(data = nw_nights_spat, aes(color = state))

#calculated scaled versions of the nightly continuous covariates
vmat_temp_a <- nw_nights_all %>%
  select(tmin, daylight) %>%
  scale()

#Reformat some of the categorical detection covariates and
#add those into the matrix for detection covariates.
#Also have a version with treating clutter as continuous in
#the model (values -2, 1, 0, 1, 2).
det_cat <- nw_nights_all %>%
  select(clutter, water_ind)
det_cat$clutter <- det_cat$clutter %>%
  as.character() %>%
  as.numeric()

vmat_temp_b <- model.matrix(~clutter + water_ind, data = det_cat)

vmat <- cbind(vmat_temp_a, vmat_temp_b[, -1])

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

#Convert to a vector
n_visits_vec <- n_visits_dfb %>%
  ungroup() %>%
  select(-cell) %>%
  as.matrix() %>%
  t() %>%
  as.vector()

#Additional pieces
n_site_years <- length(n_visits_vec)


# Occurrence level covariates ---------------------------------------------

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

xmat_cliff <- nw_grid_all %>%
  st_drop_geometry() %>%
  mutate(log_fc = log(p_forest + 1)) %>% 
  select(log_fc, precip, DEM_max, cliff_cover) %>%
  scale()
xmatc <- xmat_cliff[which(nw_grid_all$samp_all == 1), ]
xmatd <- xmat_cliff[which(nw_grid_all$samp_all == 0), ]
# Setup remaining Stan data. ----------------------------------------------

#Setup some of the Stan data.
#Need to add a few things to separate the different years and
#keep track of which each index corresponds to.
n_sites_total <- nrow(xmata)
n_years <- length(unique(nw_nights_all$year))

#Number of covariates for occupancy and detection
n_xcovs <- ncol(xmata)
n_vcovs <- ncol(vmat)

#site factor
site_f <- factor(rep(c(1:n_site_years), times = n_visits_vec))

#Get naive occupancy, appropriately fill in 0 for years without visits
## Initialize an empty list for all species data
occ_data <- list()

for (i in possible_bats) {
  dets <- nw_nights_all[[i]]
  n_obs <- length(dets)
  naive_occ <- rep(NA, n_site_years)
  naivea <- as.vector(tapply(dets, site_f, max))
  naiveb <- rep(NA, n_site_years)
  naiveb[which(n_visits_vec > 0)] <- naivea
  naiveb[which(n_visits_vec == 0)] <- 0
  naive_occ <- naiveb
  
  if (i %in% c("anpa", "euma", "myci", "pahe")) {
    occ_data[[i]] <- list('n_sites_total' = n_sites_total,
                          'n_site_years' = n_site_years,
                          'n_years' = n_years,
                          'n_obs' = n_obs,
                          'dets' = dets,
                          'n_visits' = n_visits_vec,
                          'naive_ind' = naive_occ,
                          'n_covs1' = n_xcovs,
                          'xmat' = xmatc,
                          'n_covs2' = n_vcovs,
                          'vmat' = vmat)
  }
  else{occ_data[[i]] <- list('n_sites_total' = n_sites_total,
                             'n_site_years' = n_site_years,
                             'n_years' = n_years,
                             'n_obs' = n_obs,
                             'dets' = dets,
                             'n_visits' = n_visits_vec,
                             'naive_ind' = naive_occ,
                             'n_covs1' = n_xcovs,
                             'xmat' = xmata,
                             'n_covs2' = n_vcovs,
                             'vmat' = vmat)}
  
}

occ_data$laci

# Fit the stan model ------------------------------------------------------
## Load the model
docc_model1 <- stan_model(here('Code/docc_model1.stan'))

for (i in 1:length(occ_data)) {
  print(names(occ_data)[i])
  
  # Fit occupancy model in Stan
  occ_stan <- rstan::sampling(object = docc_model1,
                       data = occ_data[[i]])
  path_name <- paste0("DataProcessed.nosync/results/stan/full/fits/",names(occ_data)[i], "_stan.rds")
  
  saveRDS(occ_stan, here(path_name))
}
# Write out needed files --------------------------------------------------

saveRDS(xmat_all, here("DataProcessed.nosync/results/stan/full/xmat_all.rds"))
saveRDS(occ_data, here("DataProcessed.nosync/results/stan/full/occ_data.rds"))
st_write(nw_grid_all, here("DataProcessed.nosync/occurrence/nw_grid_all.shp"), append = F)

print(Sys.time())
