library(here)
library(tidyverse)
library(sf)

nw_nights <- readRDS(here("Background/nabat-code-tutorial.nosync/DataFiles/nw_nights.rds"))
nw_grid_shp <- readRDS(here("Background/nabat-code-tutorial.nosync/DataFiles/nw_grid_shp.rds"))

#Arrange nights by conus_id, then year
nw_nights_all <- nw_nights %>%
  arrange(conus_id, year)

#divide shapefile into sampled and unsampled grid cells.
#also arrange these by conus_id
nw_grid3a <- nw_grid_shp %>%
  filter(samp_all == 1) %>%
  arrange(conus_id)
nw_grid3b <- nw_grid_shp %>%
  filter(samp_all == 0) %>%
  arrange(conus_id)

#combine these back again, reordered. calculate scaled versions of each
#covariate and save these in matrices that will be used for predictions
nw_grid_all3 <- rbind(nw_grid3a, nw_grid3b)
xmat_all3 <- nw_grid_all3 %>%
  st_drop_geometry() %>%
  mutate(log_fc = log(forest_cover*100+1)) %>% 
  select(log_fc, annual_prcp, cliff_cover) %>%
  scale()
xmat3a <- xmat_all3[which(nw_grid_all3$samp_all == 1), ]
xmat3b <- xmat_all3[which(nw_grid_all3$samp_all == 0), ]

#calculated scaled versions of the nightly continuous covariates
vmat_temp3a <- nw_nights_all %>%
  select(tmin, daylight) %>%
  scale()

#Reformat some of the categorical detection covariates and
#add those into the matrix for detection covariates.
#Also have a version with treating clutter as continuous in
#the model (values -2, 1, 0, 1, 2).
det_cat3 <- nw_nights_all %>%
  select(clutter2, water_ind2)
det_cat3$clutter2 <- det_cat3$clutter2 %>%
  as.character() %>%
  as.numeric()

vmat_temp3b <- model.matrix(~clutter2 + water_ind2, data = det_cat3)

vmat3 <- cbind(vmat_temp3a, vmat_temp3b[, -1])

#Setup some of the Stan data.
#Need to add a few things to separate the different years and
#keep track of which each index corresponds to.
n_sites_total3 <- nrow(xmat3a)
dets3 <- nw_nights_all$myci
n_obs3 <- length(dets3)
n_years3 <- length(unique(nw_nights_all$year))

#Number of visits per grid cell each year
n_visits_df3a <- nw_nights_all %>%
  group_by(conus_id, year) %>%
  summarise(n_visits = length(date))