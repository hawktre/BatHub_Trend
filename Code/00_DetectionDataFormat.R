## ---------------------------
##
## Script name: 00_DataCuration.R
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

require(tidyverse)
require(here)
require(sf)
require(skimr)
require(daymetr)
require(future.apply)
require(parallel)

# Read in Required Data ---------------------------------------------------
tblDeployment <- read_csv("DataRaw.nosync/database/tblDeployment.csv")
tblPointLocation <- read_csv("DataRaw.nosync/database/tblPointLocation.csv")
tblSite <- read_csv("DataRaw.nosync/database/tblSite.csv")
tluClutter <- read_csv("DataRaw.nosync/database/tluClutter.csv")
tluWaterBodyType <- read_csv("DataRaw.nosync/database/tluWaterBodyType.csv")

acoustics_example <- read_csv("Code/vignette-bayesian-site-occupancy-model-bat-acoustic-data.nosync/data/survey_data.csv")

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

## Fix Date Columns
deployment$DeploymentDate <- as_datetime(deployment$DeploymentDate, format = "%m/%d/%y %T") %>% as_date()
deployment$RecoveryDate <- as_datetime(deployment$RecoveryDate, format = "%m/%d/%y %T") %>% as_date()
deployment$year <- year(deployment$DeploymentDate)


# Read in detection data --------------------------------------------------
## Read in
acoustics <- data.table::fread(here("DataRaw.nosync/database/tblDeploymentDetection7.csv"))

## Join with deployments and clean
acoustics <- left_join(acoustics, deployment, by = c("DeploymentID" = "ID")) %>% 
  ## Select the Columns we want to keep
  select(ID, LocationName, Night, Latitude, Longitude, ClutterType, WaterBodyType, ClutterPercent, ManualIDSpp1) %>% 
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
bad_tabr <-acoustics %>% filter(ManualIDSpp1 == "tabr") %>% dplyr::slice_max(Latitude) %>% .$ID

##remove bad record
acoustics <- acoustics %>% filter(ID != bad_tabr)


# Remove non-bats ---------------------------------------------------------
acoustics <- acoustics %>% filter(ManualIDSpp1 %in% possible_bats)


# Pivot Wider to get Spp Richness -------------------------------------------------------------
## Pivot Wider
acoustics_wide <- acoustics %>%
  select(-ID) %>% 
  distinct() %>% 
  pivot_wider(names_from = ManualIDSpp1, values_from = ManualIDSpp1, 
              id_cols = c("LocationName", "Night", "Latitude", "Longitude",
                          "ClutterType", "WaterBodyType", "ClutterPercent", "Year"),
              values_fill = 0,
              values_fn = ~if_else(is.na(.), 0, 1))

##Create indicator for water and change 
acoustics_wide <- acoustics_wide %>% 
  mutate(water_ind = if_else(WaterBodyType == "None", 0, 1))



# Get Daymet min temp -----------------------------------------------------
## Read in study extent
conus10k <- read_sf(here("DataRaw.nosync/conus_10km_full/complete_conus_mastersample_10km_attributed.shp"))

## Crop to contain only PNW
pnw <- c("Oregon", "Washington", "Idaho")

conus10k_pnw <- conus10k %>% 
  filter(state_n_1 %in% pnw | state_n_2 %in% pnw)

## Check our study area
plot(conus10k_pnw["state_n_1"])

#create a function to get daymet data for every pointi
get_daymet <- function(i, dat){
  tmp <- download_daymet(site = dat$LocationName[i],
                         lat = dat$Latitude[i],
                         lon = dat$Longitude[i],
                         start = 2016,
                         end = 2022) 
    daymet <- tmp %>% 
      .$data %>% 
    as_tibble() %>% 
    mutate(date = as.Date(paste(year,yday, sep = "-"),"%Y-%j"),
           site = tmp$site) %>% 
    janitor::clean_names() %>% 
    select(site, date, dayl_s, prcp_mm_day, tmax_deg_c, tmin_deg_c, vp_pa) 
  
  return(daymet)
}

#format data
daymet_get <- acoustics_wide %>% 
  select(LocationName, Latitude, Longitude) %>% 
  distinct() %>% 
  arrange(LocationName)

#download data for all sites 
daymet_all <- lapply(1:nrow(daymet_get), function(x){get_daymet(i = x, dat = daymet_get)}) 
write_rds(daymet_all, here("DataRaw.nosync/daymet/daymet_full.rds"))

#Clean up
daymet_all <- daymet_all %>%
  bind_rows() %>% 
  mutate(location_date = paste0(site,"_",as.character(date)))

#Join with detection data
acoustics_wide <- acoustics_wide %>% 
  mutate(location_date = paste0(LocationName, "_", Night)) %>% 
  left_join(.,daymet_all, by = "location_date")
  
acoustics_wide <- acoustics_wide %>% 
  select(-c(site, date, location_date))

#Write out
write_csv(acoustics_wide, here("DataProcessed.nosync/detections_formatted_2016-2022.csv"))

detections <- st_as_sf(acoustics_wide, coords = c("Longitude", "Latitude"))

plot(detections)
