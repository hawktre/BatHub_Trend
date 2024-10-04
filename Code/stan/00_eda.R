## ---------------------------
##
## Script name: 00_eda.r
##
## Purpose of script: Conducting EDA before trying to compare to tom script
##
## Author: Trent VanHawkins
##
## Date Created: 2024-08-13
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
library(FedData)
library(terra)
library(ggcorrplot)

# # Calculate Topographic Roughness-------------------------------------
# covars <- st_read(here("DataProcessed/occurrence/batgrid_covars.shp"))
# elev <- elevatr::get_elev_raster(covars, z = 8, clip = "bbox")
# covars_vec <- terra::vect(covars)
# elev_rast <- rast(elev)
# topo_rough <- terra::extract(elev_rast, covars_vec, sd, bind = T)
# 
# covars <- st_as_sf(topo_rough)
# 
# # Plot spatial covariates -------------------------------------------------
# covars <- covars %>% rename("cliff" = EVT_NAME,
#                             "rough" = filec64e1f1b6e0a,
#                             "physio" = physio_div,
#                             "elev" = DEM_max,
#                             "forest" = p_forest)
# 
# st_write(covars, here("DataProcessed/occurrence/batgrid_covars_spocc.shp"), append = F)
# Read in the data --------------------------------------------------------
covars <- st_read(here("DataProcessed/occurrence/batgrid_covars_spocc.shp"))
dets <- readRDS(here("DataProcessed/detections/nw_nights.rds"))

covars_orwa <- covars %>% filter(state %in% c("Oregon", "Washington"))
plot(covars_orwa[c('cliff', 'precip', 'physio', "elev", 'forest', 'rough')])

plot(covars['state'])
plot(covars['precip'])
# View histograms and transformations -------------------------------------
## Pivot Longer
covars_long <- covars_orwa %>% pivot_longer(cols = c(cliff, precip, physio, elev, forest, rough),
                                       names_to = "covar", values_to = "value")
## Histogram
covars_long %>% 
  ggplot()+
  geom_histogram(aes(x = value))+
  facet_wrap(~covar, scales = "free")

## Log-transformed histogram
covars_long %>% 
  mutate(log_val = log(value + 1)) %>% 
  ggplot()+
  geom_histogram(aes(x = log_val))+
  facet_wrap(~covar, scales = "free")

## We have some with negative precip
covars_orwa <- covars_orwa[precip > 0,]
covars_long <- covars_long[covars_long$value > 0, ]

## Check Correlation
cor <- cor(covars_orwa %>% st_drop_geometry() %>% select(cliff, elev, forest, physio, precip, rough, mean_temp))
ggcorrplot(cor, lab = T)


# View raw detection probabilities over time ------------------------------
dets <- left_join(dets, covars, by = c("cell" = "CONUS_10KM"))
dets %>%
  filter(state %in% c("Oregon", "Washington")) %>% 
  group_by(year) %>% 
  summarise(laci_mean = mean(laci),
            mylu_mean = mean(mylu),
            lano_mean = mean(lano)) %>% 
  pivot_longer(cols = laci_mean:lano_mean) %>% 
  ggplot(aes(x = year, y = value))+
  geom_point(aes(color = name))+
  ylim(c(0,1))

dets_cor <- cor(dets %>% mutate(clutter = as.numeric(clutter)) %>%  select(tmin, daylight, clutter))
ggcorrplot(dets_cor, lab = T)
