## ---------------------------
##
## Script name: 00b_DetectionDataFormat.R
##
## Purpose of script:
##
## Author: Trent VanHawkins
##
## Date Created: 2026-01-17
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

# Read in corrected data --------------------------------------------------
deployment <- readRDS(here("DataProcessed/detections/detections_to2024.rds"))

# Read in detection data --------------------------------------------------
## Read in
acoustics_to_2024 <- data.table::fread(here("DataRaw/tables/calls_to_2024.csv"))
acoustics_from_2024 <- data.table::fread(here("DataRaw/tables/calls_from_2024.csv"))

all_raw_acoustics <- bind_rows(acoustics_to_2024, acoustics_from_2024)

## Join with deployments and clean
acoustics <- left_join(deployment, all_raw_acoustics, by = c("id" = "DeploymentID")) %>% 
  ## Select the Columns we want to keep
  select(id, location_name, Night, latitude, longitude, clutter_percent, water_ind, ManualIDSpp1) %>% 
  ## Drop Blanks from Manual SPP ID
  drop_na(ManualIDSpp1) %>% 
  ## Fix values and ensure all in same case
  mutate(ManualIDSpp1 = case_when(ManualIDSpp1 == 'LASCIN' ~ 'LACI',
                                  ManualIDSpp1 == 'LASNOC' ~ 'LANO',
                                  ManualIDSpp1 == 'MYOCIL' ~ 'MYCI',
                                  ManualIDSpp1 == 'MYOEVO' ~ 'MYEV',
                                  ManualIDSpp1 == 'MYOLUC' ~ 'MYLU',
                                  ManualIDSpp1 == 'MYOYUM' ~ 'MYYU',
                                  ManualIDSpp1 == 'MYOCAL' ~ 'MYCA',
                                  ManualIDSpp1 == 'EPTFUS' ~ 'EPFU',
                                  ManualIDSpp1 == 'MYOTHY' ~ 'MYTH',
                                  TRUE ~ ManualIDSpp1),
         ManualIDSpp1 = tolower(ManualIDSpp1),
         Night = mdy_hms(Night),
         Year = year(Night))

## Create a list of possible bat IDs
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


# Remove WA TABR ----------------------------------------------------------

##find the record
bad_tabr <- acoustics %>% filter(ManualIDSpp1 == "tabr") %>% dplyr::slice_max(latitude) %>% .$id

##remove bad record
acoustics <- acoustics %>% filter(id != bad_tabr)


# Remove non-bats ---------------------------------------------------------
acoustics <- acoustics %>% filter(ManualIDSpp1 %in% possible_bats)


# Pivot Wider to get Spp Richness -------------------------------------------------------------
## Pivot Wider
acoustics_wide <- acoustics %>%
  select(-id) %>% 
  distinct() %>% 
  pivot_wider(names_from = ManualIDSpp1, values_from = ManualIDSpp1,
              values_fill = 0,
              values_fn = ~if_else(is.na(.), 0, 1))

# Get Daymet min temp (using daymetr)-----------------------------------------------------
## Write out daymet batch file
daymet_batch <- acoustics_wide %>% select(location_name, latitude, longitude, Night) %>% 
  rename("site" = location_name, 
         "night" = Night) %>% 
  distinct()

write.csv(daymet_batch, here("DataRaw/covariates/daymet/daymet_batch.csv"), row.names = F)



daymet_all <- readRDS(here("DataRaw/covariates/daymet/all_daymet.rds"))

#Clean up
daymet_wide <- daymet_all %>%
  filter(!str_detect(pattern = "Error", string = yday)) |> 
  mutate(value = as.numeric(value)) |> 
  pivot_wider(names_from = measurement, values_from = value)

daymet_wide$date <- make_date(as.numeric(daymet_wide$year)) + days(as.numeric(daymet_wide$yday) - 1)

daymet_clean <- daymet_wide |> 
  select(site, date, dayl..s., prcp..mm.day., tmax..deg.c., tmin..deg.c.) |> 
  rename("daylight" = dayl..s.,
         "precipitation" = prcp..mm.day.,
         "tmax" = tmax..deg.c.,
         "tmin" = tmin..deg.c.)


#Join with detection data
acoustics_wide <- acoustics_wide %>%
  left_join(daymet_clean, by = c("location_name" = "site", "Night" = "date"))



states <- spData::us_states |> filter(NAME %in% c("Oregon", "Washington", "Idaho")) |> st_transform(crs = "WGS84")


#Write out
saveRDS(acoustics_wide, here("DataProcessed/detections/detections_formatted_2016-2024.rds"))
