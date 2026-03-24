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
library(janitor)
# Read in deloyment data and join with daymet  --------------------------------------------------
deployment <- readRDS(here("DataProcessed/detections/deployments_to2024.rds"))

# Read in detection data --------------------------------------------------
## Read in
acoustics_to_2024 <- data.table::fread(here("DataRaw/tables/calls_to_2024.csv"))
acoustics_from_2024 <- data.table::fread(here("DataRaw/tables/calls_from_2024.csv"))

all_raw_acoustics <- bind_rows(acoustics_to_2024, acoustics_from_2024)

## Clean up
acoustics <- all_raw_acoustics %>% 
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

# Remove non-bats ---------------------------------------------------------
acoustics <- acoustics %>% 
  filter(ManualIDSpp1 %in% possible_bats)


# Join with deployments ---------------------------------------------------

all_detections <- left_join(acoustics, deployment, by = c("DeploymentID" = "id"))

## Select just the columns we want
detections <- all_detections %>% 
  select(names(deployment)[-1], ManualIDSpp1, Night, DeploymentID) %>%
  select(sample_unit_id, location_name, Night, everything()) %>% 
  select(-c(deployment_date, recovery_date)) %>% 
  clean_names()

##Take the first night in the case of multiple nights
detections <- detections %>% 
  group_by(location_name, year) %>% 
  slice_min(night) %>% 
  ungroup() %>% 
  drop_na()

## Drop NA's

## Check that we don't have any more temporal replicates in a year
detections %>% 
  select(location_name, year, night) %>% 
  distinct() %>% 
  group_by(location_name, year) %>% 
  summarise(N = n(), .groups = "drop") %>% 
  filter(N > 1)

## Save out sites for daymet
daymet_sites <- detections %>% 
  select(location_name, latitude, longitude, night) %>% 
  distinct() %>% 
  mutate(night = date(night))

write.csv(daymet_sites, here("DataRaw/covariates/daymet/daymet_sites.csv"), row.names = F)
# Remove WA TABR ----------------------------------------------------------

##find the record
bad_tabr <- detections %>% filter(manual_id_spp1 == "tabr") %>% dplyr::slice_max(latitude) %>% .$deployment_id
##remove bad record
detections <- detections %>% filter(deployment_id != bad_tabr)

# Pivot Wider to get Spp Richness -------------------------------------------------------------
## Pivot Wider
detections_wide <- detections %>%
  drop_na() %>% 
  distinct() %>% 
  pivot_wider(names_from = manual_id_spp1, values_from = manual_id_spp1,
              values_fill = 0,
              values_fn = ~if_else(is.na(.), 0, 1)) %>% 
  select(-deployment_id)


#Write out
saveRDS(detections_wide, here("DataProcessed/detections/detection_histories.rds"))
