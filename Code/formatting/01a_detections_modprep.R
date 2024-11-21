## ---------------------------
##
## Script name: 01a_modelprep.R
##
## Purpose of script: Prepare data for modeling steps for single-species, single-season 
##
## Author: Trent VanHawkins
##
## Date Created: 2024-05-08
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

require(tidyverse)
require(here)
require(sf)


# Read in the data --------------------------------------------------------

detections <- read_csv(here("DataProcessed/detections/detections_formatted_2016-2022.csv"), na = c("<null>", "NA"))

# Reformat --------------------------------------------------------------------
## Make name of "cell" same in all datasets
detections <- detections %>% 
  mutate(cell = as.numeric(sub("\\_.*", "", LocationName))) %>% 
  select(cell, "replicate" = LocationName, everything())

# compute cells we would drop following single season vignette (i. --------
cell_total <- detections %>% 
  group_by(Year) %>% 
  summarise(n_cell = n_distinct(cell)) %>% 
  ungroup()

## Compute how many cells we would drop using Kathy's criteria
cell_drop <- detections %>% 
  group_by(Year, replicate) %>% 
  reframe(cell = cell,
          n_nights = n_distinct(Night)) %>% 
  ungroup() %>% 
  filter(n_nights > 1) %>% 
  group_by(Year) %>% 
  summarise(n_drop = n_distinct(cell)) %>% 
  ungroup() %>% 
  left_join(cell_total, by = "Year") %>% 
  mutate(p_drop = n_drop/n_cell*100)

require(kableExtra)

cell_drop %>% gt::gt() %>% gt::fmt_number(decimals = 1, columns = 4)



# Take the first night for any site surveyed multiple nights --------------
## Create the new filtered dataset
detections_fn <- detections %>% 
  group_by(replicate, Year) %>% 
  filter(Night == min(Night)) %>% 
  ungroup()

## Make sure we don't have any more temporal replicates
detections_fn %>% 
  group_by(replicate, Year) %>% 
  summarise(N = n()) %>% 
  filter(N > 1) %>% 
  sum(.$N)

## We don't, so we are good to go. 


# Spatial Replicates ------------------------------------------------------

## Visualize how many spatial replicates in each cell
detections_fn %>% 
  group_by(cell, Year) %>% 
  summarise(spatial_reps = n()) %>% 
  ungroup() %>% 
  ggplot(aes(x = spatial_reps)) + 
  geom_bar() +
  facet_wrap(~Year)+
  theme_classic()

## what is the range of spatial replicates?
spatial_reps_summary <- detections_fn %>% 
  group_by(cell, Year) %>% 
  summarise(spatial_reps = n())

range(spatial_reps_summary$spatial_reps)

## We have min 1 and max 7 spatial replicates in any given year.


# Find covariate NAs ------------------------------------------------------
## Find rows with NA values
na.rows <- detections_fn %>% 
  select(-c(ClutterType, WaterBodyType)) %>% 
  filter(!complete.cases(.))

detections_fn %>% 
  skimr::skim()

## We have 87 missing clutter percent and 4 missing daymet; so we drop those
detections_fn <- detections_fn %>% 
  select(-c(ClutterType, WaterBodyType)) %>% 
  drop_na() %>% 
  select(cell, 
         replicate, 
         "lat" = Latitude, 
         "lon" = Longitude, 
         "date" = Night,
         "year" = Year,
         "tmin" = tmin_deg_c,
         "prcp" = prcp_mm_day,
         "vp" = vp_pa,
         "daylight" = dayl_s,
         "clutter" = ClutterPercent,
         "water_ind" = water_ind,
         everything()) %>% 
  select(-tmax_deg_c)

# Reformatting occurrence data for model ----------------------------------
## center clutter 
detections_fn <- detections_fn %>%
  mutate(clutter = factor(clutter, labels = c(-1, 0, 1, 2, 3))) %>% #center clutter percent
  arrange(year, cell) 


saveRDS(detections_fn, here("DataProcessed/detections/nw_nights.rds"))
