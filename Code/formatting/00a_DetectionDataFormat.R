## ---------------------------
##
## Script name: 00a_join_dbtables.R
##
## Purpose of script: Read-in and format data for trend analysis
##
## Author: Trent VanHawkins
##
## Date Created: 2024-04-29
##
##
## ---------------------------

## view outputs in non-scientific notation

options(digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(janitor)

# Read in Required Data ---------------------------------------------------
tblDeployment <- read_csv(here("DataRaw/tables/tblDeployment.csv"))
tblPointLocation <- read_csv(here("DataRaw/tables/tblPointLocation.csv"))
tblSite <- read_csv(here("DataRaw/tables/tblSite.csv"))
tluClutter <- read_csv(here("DataRaw/tables/tluClutterType.csv")) 
tluWaterBodyType <- read_csv(here("DataRaw/tables/tluWaterBodyType.csv"))

# Join all tables together ------------------------------------------------
all_join <- left_join(tblDeployment, tblPointLocation, by = join_by(PointLocationID == ID)) %>% 
  left_join(., tblSite, by = join_by(SiteID == ID)) %>% 
  left_join(., tluClutter, by = join_by(ClutterTypeID == ID)) %>% 
  left_join(., tluWaterBodyType, by = join_by(WaterBodyTypeID == ID))

# Select only the columns we need -----------------------------------------

deployment <- all_join %>% select(ID,
                    SampleUnitID,
                    LocationName,
                    Latitude,
                    Longitude,
                    DeploymentDate,
                    RecoveryDate,
                    Label.x,
                    Label.y,
                    ClutterPercent) %>% 
  rename("ClutterType" = "Label.x",
         "WaterBodyType" = "Label.y") 



# Make Dates ----------------------------------------
deployment$DeploymentDate <- as_date(as_datetime(deployment$DeploymentDate, format = "%m/%d/%y %H:%M:%S"))
deployment$RecoveryDate <- as_date(as_datetime(deployment$RecoveryDate, format = "%m/%d/%y %H:%M:%S"))
deployment$year <- year(deployment$DeploymentDate)


# Correct  NAs and wrong values -----------------------------
## Clutter Percent
unique(deployment$ClutterPercent)

deployment <- deployment %>% mutate(ClutterPercent = factor(case_when(ClutterPercent == "0% (no structural interference, e.g., open habitat)" ~ "0",
                                                               ClutterPercent == "1 to 25%" ~ "1", 
                                                               ClutterPercent == "26 to 50%" | ClutterPercent == "26-50" ~ "2",
                                                               ClutterPercent == "<null>" ~ NA, 
                                                               TRUE ~ ClutterPercent)))
sum(is.na(deployment$ClutterPercent))

## Clutter Type
unique(deployment$ClutterType)
deployment <- deployment %>% mutate(ClutterType = if_else(ClutterType == "<null>", NA, ClutterType))

## Water Bodies
unique(deployment$WaterBodyType)
### Create Water Indicator
deployment <- deployment %>% mutate(water_ind = if_else(WaterBodyType == "None", 0, 1))
### Correct Water Indicator NA's 
deployment <- deployment %>% mutate(water_ind = if_else(is.na(water_ind) & str_detect(ClutterType, "Water"), 1, water_ind))


# How many rows are we dropping?  -----------------------------------------
## Write CSV for missing locations
missing_sites <- deployment %>% filter(is.na(ClutterPercent) | is.na(water_ind))
write.csv(missing_sites, here("DataProcessed/detections/sites_missing_covars.csv"))

## How many sites? 
n_drop <- sum(is.na(deployment$ClutterPercent) | is.na(deployment$water_ind))
p_drop <- n_drop/nrow(deployment)

cat("Dropping", n_drop, "rows missing clutter percent or waterbody indicator")

## Drop the missing rows
deployment <- deployment %>% drop_na(ClutterPercent, water_ind) %>% 
  select(-c(ClutterType, WaterBodyType))

deployment <- clean_names(deployment)

# Save out for verification in excel --------------------------------------
saveRDS(deployment, here("DataProcessed/detections/detections_to2024.rds"))

