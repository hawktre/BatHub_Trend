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
tblDeployment <- read_csv(here("DataRaw/tables/tblDeployment.csv"))
tblPointLocation <- read_csv(here("DataRaw/tables/tblPointLocation.csv"))
tblSite <- read_csv(here("DataRaw/tables/tblSite.csv"))
tluClutter <- read_csv(here("DataRaw/tables/tluClutterType.csv"))
tluWaterBodyType <- read_csv(here("DataRaw/tables/tluWaterBodyType.csv"))


# Set the years you want to analyze  --------------------------------------

first_year <- 2016
last_year <- 2024

# Fix clutter 
unique(tblDeployment$ClutterPercent)

tblDeployment <- tblDeployment |> 
  mutate(ClutterPercent = case_when(ClutterPercent == "0% (no structural interference, e.g., open habitat)" ~ 0, 
ClutterPercent == "1 to 25%" ~ 1,
ClutterPercent == "26-50" ~ 2,
ClutterPercent == "26 to 50%" ~ 2,
ClutterPercent == "<null>" ~ NA_integer_, 
TRUE ~ as.integer(ClutterPercent)))

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

## Create Waterbody Indicator & infer missing clutter percent
deployment <- deployment |> mutate(water_ind = if_else(WaterBodyType != "None" | ClutterType == "Water", 1, 0),
ClutterPercent = if_else(is.na(ClutterPercent) & ClutterType == "None", 0, ClutterPercent)) #Infer missing clutter percent from clutter type

## Fix Date Columns
deployment$DeploymentDate <- as_datetime(deployment$DeploymentDate, format = "%m/%d/%y %T") %>% as_date()
deployment$RecoveryDate <- as_datetime(deployment$RecoveryDate, format = "%m/%d/%y %T") %>% as_date()
deployment$year <- year(deployment$DeploymentDate)


# Read in detection data --------------------------------------------------
## Read in
acoustics_to_2024 <- data.table::fread(here("DataRaw/tables/calls_to_2024.csv"))
acoustics_from_2024 <- data.table::fread(here("DataRaw/tables/calls_from_2024.csv"))

all_raw_acoustics <- bind_rows(acoustics_to_2024, acoustics_from_2024)

## Join with deployments and clean
acoustics <- left_join(all_raw_acoustics, deployment, by = c("DeploymentID" = "ID")) %>% 
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

# Get Daymet min temp (using daymetr)-----------------------------------------------------
## Write out daymet batch file
daymet_batch <- acoustics_wide %>% select(LocationName, Latitude, Longitude) %>% rename("site" = LocationName,
                                                                                        "latitude" = Latitude,
                                                                                        "longitude" = Longitude) %>% 
  distinct() |> 
  drop_na()

write.csv(daymet_batch, here("DataRaw/covariates/daymet/daymet_batch.csv"), row.names = F)

#Download the daymet data (optional) ------------------------------
# df_batch <- download_daymet_batch(file_location = here("DataRaw/covariates/daymet/daymet_batch.csv"), 
#                       start = first_year, 
#                       end = last_year, 
#                       internal = T,
#                     simplify = T)

# saveRDS(df_batch, here("DataRaw/covariates/daymet/all_daymet.rds"))
daymet_all <- readRDS(here("DataRaw/covariates/daymet/all_daymet.rds"))
#Clean up
daymet_wide <- daymet_all %>%
  filter(!str_detect(pattern = "Error", string = yday)) |> 
  mutate(value = as.numeric(value)) |> 
  pivot_wider(names_from = measurement, values_from = value)

daymet_wide$date <- make_date(as.numeric(daymet_wide$year)) + days(as.numeric(daymet_wide$yday) - 1)

daymet_wide <- daymet_wide |> 
  select(site, date, dayl..s., prcp..mm.day., tmax..deg.c., tmin..deg.c.) |> 
  rename("daylight" = dayl..s.,
  "precipitation" = prcp..mm.day.,
"tmax" = tmax..deg.c.,
"tmin" = tmin..deg.c.)

#Join with detection data
acoustics_wide <- acoustics_wide %>% 
  left_join(daymet_wide, by = c("LocationName" = "site", "Night" = "date"))

test <- acoustics_wide |> filter(is.na(daylight)) |> drop_na(Latitude, Longitude)  |>  st_as_sf(coords = c("Longitude", "Latitude"), crs = "WGS84")

states <- spData::us_states |> filter(NAME %in% c("Oregon", "Washington", "Idaho")) |> st_transform(crs = "WGS84")

test |> 
  ggplot()+
  geom_sf()+
  geom_sf(data = states, fill = "transparent")
#Write out
saveRDS(acoustics_wide, here("DataProcessed/detections/detections_formatted_2016-2024.rds"))
