##Data Setup##

####Overview####
#Code used in the Data Setup section of the code tutorial.

#R Version 3.6.6 (2020-02-29)
#06/16/2020
#Wilson Wright

####Setup####
#Load packages
library(tidyverse)
library(sf)
library(readxl)

#Specify the file paths
file_path1 <- paste0('RawData_NWBatHub/',
                     'OR_WA_2016-2018-BatHub-AllData_03272020.xlsx')
file_path2 <- paste0('RawData_NWBatHub/',
                     'OR_WA_2019-2020_BatHub_AllData_01202022.xlsx')

#Read in both files
nw_nights_raw1 <- read_excel(file_path1)
nw_nights_raw2 <- read_excel(file_path2)

#warnings()

#Define file path for the covariate data
file_covs <- paste0('RawData_NWBatHub/',
                    'OR-WA-FullGrid-Landcover-Covariates.xlsx')

#Read in file
nw_grid_raw <- read_excel(file_covs)

#View file
nw_grid_raw

nw_grid <- nw_grid_raw %>%
  select('conus_id' = starts_with("Sample"),
         'elev' = starts_with("Elevation_M"),
         'elev_std' = starts_with("Elevation_S"),
         'snag' = starts_with("Snags"), #Note which column we are using - 2017 here 
         'prime_prod' = starts_with("NetPrimary"),
         'cliff_cover' = starts_with("Cliffs"),
         'forest_cover' = starts_with("Forest"),
         'annual_prcp' = starts_with("MeanAnnual"),
         'state' = State)

names(nw_grid)

#Section 2.2

#Specify the file paths
file_path1 <- paste0('RawData_NWBatHub/',
                     'OR_WA_2016-2018-BatHub-AllData_03272020.xlsx')
file_path2 <- paste0('RawData_NWBatHub/',
                     'OR_WA_2019-2020_BatHub_AllData_01202022.xlsx')

#Read in both files
nw_nights_raw1 <- read_excel(file_path1)
nw_nights_raw2 <- read_excel(file_path2)

#Examine data structure
nw_nights_raw1

nw_nights1 <- nw_nights_raw1 %>%
  select('conus_id' = CONUS10km,
         'quad_id' = Site,
         'lat' = lat,
         'lon' = long,
         'date' = Start,
         'year' = Year,
         'state' = State,
         'tmin' = tmin,
         'prcp' = prcp,
         'vp' = vp,
         'daylight' = dayl,
         'clutter' = 'Clutter %',
         'c_type' = 'Clutter Type',
         'hab_local' = LocalHabitat,
         'hab_broad' = BroadHabitat,
         'det_target' = DetectionTarget,
         'det_detail' = TargetDescriptor,
         'anpa' = ANPA,
         'coto' = COTO,
         'epfu' = EPFU,
         'euma' = EUMA,
         'laci' = LACI,
         'lano' = LANO,
         'myca' = MYCA,
         'myci' = MYCI,
         'myev' = MYEV,
         'mylu' = MYLU,
         'myth' = MYTH,
         'myvo' = MYVO,
         'myyu' = MYYU,
         'pahe' = PAHE,
         'tabr' = TABR)

#Reformat the state variable to match that from the nw_grid data
nw_nights1$state[which(nw_nights1$state == 'OR')] <- 'Oregon'
nw_nights1$state[which(nw_nights1$state == 'WA')] <- 'Washington'

#Reformat date to just be a normal date object, not include time
nw_nights1$date <- as.Date(nw_nights1$date)

#Identify which rows are problematic,
#note that in the actual file, the first row has the column names
which(nw_nights_raw2$Start == 0)

#####################################################################################################
#I don't think we need this with the new 2019-2020 dataset!

#Read in four separate chunks, excluding rows 96, 104, and 324.
#These are actually rows 97, 105, and 325 in the actual file
#(when the row for the column names is included).
#nw_nights_raw2a <- read_excel(file_path2,
#                              range = cell_rows(1:96))
#nw_nights_raw2b <- read_excel(file_path2,
 #                             range = cell_rows(98:104),
  #                            col_names = names(nw_nights_raw2a))
#nw_nights_raw2c <- read_excel(file_path2,
#                              range = cell_rows(106:324),
#                              col_names = names(nw_nights_raw2a))
#nw_nights_raw2d <- read_excel(file_path2,
#                              range = cell_rows(326:794),
#                              col_names = names(nw_nights_raw2a))

#Combine in a single file, saving over the previous version
#nw_nights_raw2 <- bind_rows(nw_nights_raw2a,
#                            nw_nights_raw2b,
#                            nw_nights_raw2c,
#                            nw_nights_raw2d)
############################################################################################
#Now we can check that the structure matches the first file
nw_nights_raw2

#Select the needed columns and give them new names,
#note that some of the column names have changed from the previous
#file so we need to fix that too.
nw_nights2 <- nw_nights_raw2 %>%
  select('conus_id' = CONUS10km,
         'quad_id' = Site,
         'lat' = lat,
         'lon' = long,
         'date' = Start,
         'year' = Year,
         'state' = state,
         'tmin' = tmin,
         'prcp' = prcp,
         'vp' = vp,
         'daylight' = dayl,
         'clutter' = 'Clutter%',
         'c_type' = ClutterType,
         'hab_local' = LocalHabitat,
         'hab_broad' = BroadHabitat,
         'det_target' = DetectionTarget,
         'det_detail' = TargetDescriptor,
         'anpa' = ANPA,
         'coto' = COTO,
         'epfu' = EPFU,
         'euma' = EUMA,
         'laci' = LACI,
         'lano' = LANO,
         'myca' = MYCA,
         'myci' = MYCI,
         'myev' = MYEV,
         'mylu' = MYLU,
         'myth' = MYTH,
         'myvo' = MYVO,
         'myyu' = MYYU,
         'pahe' = PAHE,
         'tabr' = TABR)

#Reformat the state variable to match that from the nw_grid data,
#For some reason there are some observations from California.
#We will exclude those later
nw_nights2$state[which(nw_nights2$state == 'OR')] <- 'Oregon'
nw_nights2$state[which(nw_nights2$state == 'WA')] <- 'Washington'
#nw_nights2$state[which(nw_nights2$state == 'CA')] <- 'California'

#Reformat date to just be a normal date object
nw_nights2$date <- as.Date(nw_nights2$date)
#################################################
#Section 2.3 Shapefiles

#load the NABat grid shapefiles for Oregon and Washington,
#Note that these folders need to be placed in your working directory
oregon_shp <- read_sf(dsn = 'RawData_NABatGrid/Oregon',
                      layer = 'NABat_Oregon_attributed')
washington_shp <- read_sf(dsn = 'RawData_NABatGrid/Washington',
                          layer = 'NABat_Washington_attributed')

###########################
#NOTE: Add in eastside versions here for desert spp

oregon_east_shp <- read_sf(dsn = 'RawData_NABatGrid/Oregon',
                      layer = 'NABat_Oregon_attributed_east')
washington_east_shp <- read_sf(dsn = 'RawData_NABatGrid/Washington',
                          layer = 'NABat_Washington_attributed_east')


#View these shapefiles
oregon_shp

washington_shp

oregon_east_shp

washington_east_shp


#Combine shapefiles
combined_shp <- rbind(oregon_shp, washington_shp)
combined_east_shp <- rbind(oregon_east_shp, washington_east_shp)
#Identify duplicates
combined_shp$grid_duplicates <- duplicated(combined_shp)
combined_east_shp$grid_duplicates <- duplicated(combined_east_shp)

####################
#For east side spp - follow this path 
#rename east side

combined_shp<-combined_east_shp


#Next we do a series of steps:
#Filter duplicated grid cells,
#Select just a few of the columns (include geometry, the spatial information),
#Join with the covariate data (this excludes cells without covariates),
#Remove grid cells that have 0 for elevation standard deviation,
#Project data to lat-lon spatial coordinates
combined_shp2 <- combined_shp %>%
  filter(grid_duplicates == FALSE) %>%
  select('conus_id' = CONUS_10KM,
         'grts_id' = GRTS_ID,
         geometry) %>%
  left_join(nw_grid, by = 'conus_id') %>%
  filter(elev_std != 0) %>%
  st_transform(crs = 4326)

##################
#testing here 
###################
#combined_shp2 <- combined_shp %>%
#  filter(grid_duplicates == TRUE)



#View new file
combined_shp2

#Next add indicators for which grid cells were sampled each year
combined_shp2$samp2016 <- ifelse(combined_shp2$conus_id %in%
                                   filter(nw_nights1, year == 2016)$conus_id,
                                 1, 0)
combined_shp2$samp2017 <- ifelse(combined_shp2$conus_id %in%
                                   filter(nw_nights1, year == 2017)$conus_id,
                                 1, 0)
combined_shp2$samp2018 <- ifelse(combined_shp2$conus_id %in%
                                   filter(nw_nights1, year == 2018)$conus_id,
                                 1, 0)
combined_shp2$samp2019 <- ifelse(combined_shp2$conus_id %in%
                                   filter(nw_nights1, year == 2019)$conus_id,
                                 1, 0)

combined_shp2$samp2020 <- ifelse(combined_shp2$conus_id %in%
                                   nw_nights2$conus_id,
                                 1, 0)

#Add an overall indicator for sampled in at least one year
combined_shp2$samp_all <- combined_shp2 %>%
  select(samp2016, samp2017, samp2018, samp2019,samp2020) %>%
  st_drop_geometry() %>%
  apply(1, max)

#divide shapefile into sampled and unsampled grid cells.
#also arrange these by conus_id
combined_shp2a <- combined_shp2 %>%
  filter(samp_all == 1) %>%
  arrange(conus_id)
combined_shp2b <- combined_shp2 %>%
  filter(samp_all == 0) %>%
  arrange(conus_id)

#Combine these back together again in the new order and then
#save the shapefile as an .rds file to more easily read into R later
combined_shp2_all <- rbind(combined_shp2a, combined_shp2b)
saveRDS(combined_shp2_all, file = 'DataFiles/nw_grid_shp.rds')

#Modify clutter covariate. Note that we will treat it as numeric in the model
#so here we are just defining values for the different levels observed.
#I specified 0 for the 0-25% category because that was most common and
#could be thought of as the 'mean' clutter value.
levels(factor(nw_nights1$clutter))

levels(factor(nw_nights2$clutter))

nw_nights1$clutter2 <- factor(nw_nights1$clutter,
                              labels = c(-1, 0, 1, 2, 3))
nw_nights2$clutter2 <- factor(nw_nights2$clutter,
                              labels = c(-1, 0, 1, 2, 3))

#Combine the 2016-2018 data with the 2019 and 2020 data. Remove original clutter.
nw_nights_all <- bind_rows(select(nw_nights1, -clutter),
                           select(nw_nights2, -clutter))

#Add detection target indicator for water feature
levels(factor(nw_nights_all$det_target))

##############################################################################################
#Note I changed this to 1:19 and one less '1'
nw_nights_all$water_ind <- factor(nw_nights_all$det_target,
                                  labels = c(1:19))

#Redefining levels based on whether they included a water feature or not
levels(nw_nights_all$water_ind) <- c( '1', '1', '1', '1', '1', '1',
                                     '1', '1', '1', '1', '1', '1', '1', '1',
                                     '15', '15', '15', '15', '15')
nw_nights_all$water_ind[which(is.na(nw_nights_all$water_ind) == TRUE)] <- '1'
nw_nights_all$water_ind2 <- ifelse(nw_nights_all$water_ind == '15', 1, 0)

#Remove the rows for grid cells that are not included in the shapefile,
#this will remove all of the observations from California because those
#were not included in the shapefiles.
inc_idx <- which(nw_nights_all$conus_id %in% combined_shp2_all$conus_id)
nw_nights_all2 <- nw_nights_all[inc_idx, ]

#Arrange by year and conus id, then remove the columns we will not need
nw_nights_all3 <- nw_nights_all2 %>%
  arrange(year, conus_id) %>%
  select(-c_type, -hab_local, -hab_broad, -det_target, -det_detail, -water_ind)

#Save file
saveRDS(nw_nights_all3, file = 'DataFiles/nw_nights.rds')

#Check it out...
#Load files
nw_grid_shp <- readRDS('DataFiles/nw_grid_shp.rds')
nw_nights <- readRDS('DataFiles/nw_nights.rds')

#Summary of these files
nw_grid_shp

nw_nights

###################
#3 Exploratory Graphics
ggplot(nw_grid_shp) +
  geom_sf(aes(fill = samp_all), size = 0.1) +
  theme_bw() +
  labs(fill = "Eastside Sample") +
  viridis::scale_fill_viridis(option = 'cividis')
#change fill variable "elev" "elev_std""snag" "prime_prod" "cliff_cover" "forest_cover" "annual_prcp" 
ggplot(nw_grid_shp) +
  geom_sf(aes(fill = annual_prcp), size = 0.1) +
  theme_bw() +
  labs(fill = "30-yr Avg. Precipitation (mm)") +
  viridis::scale_fill_viridis(option = 'cividis')
#Look at TABR detection histories - raw

#Select columns for grid cell, year, and tabr detections,
#then group by grid cell and year,
#then summarize detections in each grid cell by year
naive_tabr <- nw_nights %>%
  select(conus_id, year, tabr) %>%
  group_by(conus_id, year) %>%
  summarise(naive = max(tabr))

#Next we need to restructure this to a wide format and
#then associate it with the shapefile.
#Note that grid cells without surveys for a particular year will be NAs.
#We go to wide format here because we need to join with all grid cells.
naive_tabr_wide <- pivot_wider(naive_tabr,
                               id_cols = conus_id,
                               names_from = year,
                               values_from = naive)
naive_tabr_shp <- left_join(select(nw_grid_shp, conus_id, geometry),
                            naive_tabr_wide,
                            by = 'conus_id')

#Restructure one more time for plotting, this is back to long format so
#we can use facets in ggplot (plotting variable needs to be in a single column).
naive_tabr_shp2 <- naive_tabr_shp %>%
  select('2016', '2017', '2018', '2019','2020') %>%
  #select('2019') %>%
  #gather(year, naive, -geometry)
  gather(year, naive, -geometry)
#Make plot
ggplot(naive_tabr_shp2) +
  geom_sf(aes(fill = factor(naive)), size = 0.15) +
  theme_bw(base_size = 16) +
  facet_wrap(. ~ year, ncol = 3) + 
  viridis::scale_fill_viridis(discrete = TRUE,
                              option = 'viridis',
                              na.value = '#d9d9d9')

#Make plot
ggplot(naive_tabr_shp2) +
  geom_sf(aes(fill = factor(year)), size = 0.15) +
  theme_bw(base_size = 16) +
  facet_wrap(. ~ naive, ncol = 3) + 
  viridis::scale_fill_viridis(discrete = TRUE,
                              option = 'viridis',
                              na.value = '#d9d9d9')
###################################
#This is the single-season formatting chunk in front of Stan implementation
#Clear environment (I'll also be repeating some names used previously)
#rm(list = ls())

#Load files
#nw_grid_shp <- readRDS('DataFiles/nw_grid_shp.rds')
#nw_nights <- readRDS('DataFiles/nw_nights.rds')

################################################################################################
#Check forest cover or other covars scaling
##########################################################################################################
#Add log forest cover covariate
nw_grid_shp$log_fc <- log(nw_grid_shp$forest_cover*100 + 1)
##########################################################################################################


#Filter observations for 2019 only
#nw_nights19 <- filter(nw_nights, year == 2019)

########################################################################################
#Filter covariates matrices into sampled grid cells
#and unsampled ones (to be used for predictions)
#divide shapefile into sampled and unsampled grid cells for 2018 data.
#also arrange these by conus_id
#nw_grid1a <- nw_grid_shp %>%
#  filter(samp2019 == 1) %>%
#  arrange(conus_id)
#nw_grid1b <- nw_grid_shp %>%
#  filter(samp2019 == 0) %>%
#  arrange(conus_id)

#combine these back again, reordered. calculate scaled versions of each
#covariate and save these in matrices that will be used for predictions
#nw_grid_all1 <- rbind(nw_grid1a, nw_grid1b)

#Check which covariates and also scaling! 
##################################################################################################
#xmat_all1 <- nw_grid_all1 %>%
#  st_drop_geometry() %>%
#  select(log_fc, annual_prcp, cliff_cover) %>%
#  scale() #Scaling all
#xmat1a <- xmat_all1[which(nw_grid_all1$samp2019 == 1), ]
#xmat1b <- xmat_all1[which(nw_grid_all1$samp2019 == 0), ]
#####################################################################################################
#calculated scaled versions of the nightly continuous covariates
#vmat_temp1a <- nw_nights19 %>%
#  select(tmin, daylight) %>%
#  scale()

#Reformat some of the categorical detection covariates and
#add those into the matrix for detection covariates.
#Also have a version with treating clutter as continuous in
#the model (values -2, 1, 0, 1, 2).
#det_cat1 <- nw_nights19 %>%
#  select(clutter2, water_ind2)
#det_cat1$clutter2 <- det_cat1$clutter2 %>%
#  as.character() %>%
#  as.numeric()
#
#vmat_temp1b <- model.matrix(~clutter2 + water_ind2, data = det_cat1)

#vmat1 <- cbind(vmat_temp1a, vmat_temp1b[, -1])


###############################################
#2
################################################
#Proceed here for trend

#This is the multi-season formatting chunk in front of Stan implementation...

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
###########################################################################################################
#combine these back again, reordered. calculate scaled versions of each
#covariate and save these in matrices that will be used for predictions
#cliff_cover is added for MYTH, EUMA, MYCI, ANPA, PAHE as beta 5
#adding snag (2017) in here Aug 30 2022 and will run for 7 spp
nw_grid_all3 <- rbind(nw_grid3a, nw_grid3b)
xmat_all3 <- nw_grid_all3 %>%
  st_drop_geometry() %>%
  select(elev,elev_std,annual_prcp,snag) %>%
  scale()
xmat3a <- xmat_all3[which(nw_grid_all3$samp_all == 1), ]
xmat3b <- xmat_all3[which(nw_grid_all3$samp_all == 0), ]
###########################################################################################################
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

#Need to fill in the zeros for this
n_visits_df3b <- pivot_wider(n_visits_df3a,
                             id_cols = conus_id,
                             names_from = year,
                             values_from = n_visits,
                             values_fill = list(n_visits = 0))

#reorder columns
n_visits_df3b <- select(n_visits_df3b, conus_id, '2016', '2017', '2018', '2019','2020')

#Convert to a vector
n_visits_vec3 <- n_visits_df3b %>%
  ungroup() %>%
  select(-conus_id) %>%
  as.matrix() %>%
  t() %>%
  as.vector()

#Additional pieces
n_site_years3 <- length(n_visits_vec3)

#Number of covariates for occupancy and detection
n_xcovs3 <- ncol(xmat3a)
n_vcovs3 <- ncol(vmat3)

#site factor
site_f3 <- factor(rep(c(1:n_site_years3), times = n_visits_vec3))

#Get naive occupancy, appropriately fill in 0 for years without visits
naive_occ3 <- rep(NA, n_site_years3)
naive3a <- as.vector(tapply(dets3, site_f3, max))
naive3b <- rep(NA, n_site_years3)
naive3b[which(n_visits_vec3 > 0)] <- naive3a
naive3b[which(n_visits_vec3 == 0)] <- 0
naive_occ3 <- naive3b

occ_data_ex3 <- list('n_sites_total' = n_sites_total3,
                     'n_site_years' = n_site_years3,
                     'n_years' = n_years3,
                     'n_obs' = n_obs3,
                     'dets' = dets3,
                     'n_visits' = n_visits_vec3,
                     'naive_ind' = naive_occ3,
                     'n_covs1' = n_xcovs3,
                     'xmat' = xmat3a,
                     'n_covs2' = n_vcovs3,
                     'vmat' = vmat3)

########################################################################################
#JAGS multi-season single-species
#Note that 7 and 8 are max visits respectively
#Reformat detections

#NOTE here is where we change spp

dets3_jags_temp <- nw_nights_all %>%
  select(conus_id, year, myev) %>%
  group_by(conus_id, year) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = c(year, visit),
              values_from = myev)

dets3_jags <- array(NA, dim = c(236, 8, 5))
dets3_jags[1:236, 1:7, 1] <- select(dets3_jags_temp,
                                    starts_with('2016')) %>%
  simplify2array()
dets3_jags[1:236, 1:8, 2] <- select(dets3_jags_temp,
                                    starts_with('2017')) %>%
  simplify2array()
dets3_jags[1:236, 1:8, 3] <- select(dets3_jags_temp,
                                    starts_with('2018')) %>%
  simplify2array()
dets3_jags[1:236, 1:7, 4] <- select(dets3_jags_temp,
                                    starts_with('2019')) %>%
  simplify2array()
dets3_jags[1:236, 1:7, 5] <- select(dets3_jags_temp,
                                    starts_with('2020')) %>%
  simplify2array()
########################################################################
#For eastside mask

#write.table(dets3_jags_temp,"det3_jags_temp_east.csv",sep=",")


dets3_jags <- array(NA, dim = c(140, 8, 5))
dets3_jags[1:140, 1:6, 1] <- select(dets3_jags_temp,
                                    starts_with('2016')) %>%
  simplify2array()
dets3_jags[1:140, 1:8, 2] <- select(dets3_jags_temp,
                                    starts_with('2017')) %>%
  simplify2array()
dets3_jags[1:140, 1:7, 3] <- select(dets3_jags_temp,
                                    starts_with('2018')) %>%
  simplify2array()
dets3_jags[1:140, 1:7, 4] <- select(dets3_jags_temp,
                                    starts_with('2019')) %>%
  simplify2array()
dets3_jags[1:140, 1:7, 5] <- select(dets3_jags_temp,
                                    starts_with('2020')) %>%
  simplify2array()

####################################################################################
#Do the same thing for each covariate
#Minimum temperature
tmin3_jags_temp <- data.frame(conus_id = nw_nights_all$conus_id,
                              year = nw_nights_all$year,
                              tmin = vmat3[, 'tmin']) %>%
  group_by(conus_id, year) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = c(year, visit),
              values_from = tmin)

tmin3_jags <- array(NA, dim = c(236, 8, 5))
tmin3_jags[, 1:7, 1] <- select(tmin3_jags_temp,
                               starts_with('2016')) %>%
  simplify2array()
tmin3_jags[, 1:8, 2] <- select(tmin3_jags_temp,
                               starts_with('2017')) %>%
  simplify2array()
tmin3_jags[, 1:8, 3] <- select(tmin3_jags_temp,
                               starts_with('2018')) %>%
  simplify2array()
tmin3_jags[, 1:7, 4] <- select(tmin3_jags_temp,
                               starts_with('2019')) %>%
  simplify2array()
tmin3_jags[, 1:7, 5] <- select(tmin3_jags_temp,
                               starts_with('2020')) %>%
  simplify2array()

#Daylight
dayl3_jags_temp <- data.frame(conus_id = nw_nights_all$conus_id,
                              year = nw_nights_all$year,
                              dayl = vmat3[, 'daylight']) %>%
  group_by(conus_id, year) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = c(year, visit),
              values_from = dayl)

dayl3_jags <- array(NA, dim = c(236, 8, 5))
dayl3_jags[, 1:7, 1] <- select(dayl3_jags_temp,
                               starts_with('2016')) %>%
  simplify2array()
dayl3_jags[, 1:8, 2] <- select(dayl3_jags_temp,
                               starts_with('2017')) %>%
  simplify2array()
dayl3_jags[, 1:8, 3] <- select(dayl3_jags_temp,
                               starts_with('2018')) %>%
  simplify2array()
dayl3_jags[, 1:7, 4] <- select(dayl3_jags_temp,
                               starts_with('2019')) %>%
  simplify2array()
dayl3_jags[, 1:7, 5] <- select(dayl3_jags_temp,
                               starts_with('2020')) %>%
  simplify2array()

#Clutter
clut3_jags_temp <- data.frame(conus_id = nw_nights_all$conus_id,
                              year = nw_nights_all$year,
                              clut = vmat3[, 'clutter2']) %>%
  group_by(conus_id, year) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = c(year, visit),
              values_from = clut)

clut3_jags <- array(NA, dim = c(236, 8, 5))
clut3_jags[, 1:7, 1] <- select(clut3_jags_temp,
                               starts_with('2016')) %>%
  simplify2array()
clut3_jags[, 1:8, 2] <- select(clut3_jags_temp,
                               starts_with('2017')) %>%
  simplify2array()
clut3_jags[, 1:8, 3] <- select(clut3_jags_temp,
                               starts_with('2018')) %>%
  simplify2array()
clut3_jags[, 1:7, 4] <- select(clut3_jags_temp,
                               starts_with('2019')) %>%
  simplify2array()
clut3_jags[, 1:7, 5] <- select(clut3_jags_temp,
                               starts_with('2020')) %>%
  simplify2array()

#Water ind
wind3_jags_temp <- data.frame(conus_id = nw_nights_all$conus_id,
                              year = nw_nights_all$year,
                              wind = vmat3[, 'water_ind2']) %>%
  group_by(conus_id, year) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = c(year, visit),
              values_from = wind)

wind3_jags <- array(NA, dim = c(236, 8, 5))
wind3_jags[, 1:7, 1] <- select(wind3_jags_temp,
                               starts_with('2016')) %>%
  simplify2array()
wind3_jags[, 1:8, 2] <- select(wind3_jags_temp,
                               starts_with('2017')) %>%
  simplify2array()
wind3_jags[, 1:8, 3] <- select(wind3_jags_temp,
                               starts_with('2018')) %>%
  simplify2array()
wind3_jags[, 1:7, 4] <- select(wind3_jags_temp,
                               starts_with('2019')) %>%
  simplify2array()
wind3_jags[, 1:7, 5] <- select(wind3_jags_temp,
                               starts_with('2020')) %>%
  simplify2array()

n_visits3_jags <- n_visits_df3b %>%
  ungroup() %>%
  select(-conus_id) %>%
  simplify2array()

#########################################
#East side dims

tmin3_jags_temp <- data.frame(conus_id = nw_nights_all$conus_id,
                              year = nw_nights_all$year,
                              tmin = vmat3[, 'tmin']) %>%
  group_by(conus_id, year) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = c(year, visit),
              values_from = tmin)

tmin3_jags <- array(NA, dim = c(140, 8, 5))
tmin3_jags[, 1:6, 1] <- select(tmin3_jags_temp,
                               starts_with('2016')) %>%
  simplify2array()
tmin3_jags[, 1:8, 2] <- select(tmin3_jags_temp,
                               starts_with('2017')) %>%
  simplify2array()
tmin3_jags[, 1:7, 3] <- select(tmin3_jags_temp,
                               starts_with('2018')) %>%
  simplify2array()
tmin3_jags[, 1:7, 4] <- select(tmin3_jags_temp,
                               starts_with('2019')) %>%
  simplify2array()
tmin3_jags[, 1:7, 5] <- select(tmin3_jags_temp,
                               starts_with('2020')) %>%
  simplify2array()

#Daylight
dayl3_jags_temp <- data.frame(conus_id = nw_nights_all$conus_id,
                              year = nw_nights_all$year,
                              dayl = vmat3[, 'daylight']) %>%
  group_by(conus_id, year) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = c(year, visit),
              values_from = dayl)

dayl3_jags <- array(NA, dim = c(140, 8, 5))
dayl3_jags[, 1:6, 1] <- select(dayl3_jags_temp,
                               starts_with('2016')) %>%
  simplify2array()
dayl3_jags[, 1:8, 2] <- select(dayl3_jags_temp,
                               starts_with('2017')) %>%
  simplify2array()
dayl3_jags[, 1:7, 3] <- select(dayl3_jags_temp,
                               starts_with('2018')) %>%
  simplify2array()
dayl3_jags[, 1:7, 4] <- select(dayl3_jags_temp,
                               starts_with('2019')) %>%
  simplify2array()
dayl3_jags[, 1:7, 5] <- select(dayl3_jags_temp,
                               starts_with('2020')) %>%
  simplify2array()

#Clutter
clut3_jags_temp <- data.frame(conus_id = nw_nights_all$conus_id,
                              year = nw_nights_all$year,
                              clut = vmat3[, 'clutter2']) %>%
  group_by(conus_id, year) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = c(year, visit),
              values_from = clut)

clut3_jags <- array(NA, dim = c(140, 8, 5))
clut3_jags[, 1:6, 1] <- select(clut3_jags_temp,
                               starts_with('2016')) %>%
  simplify2array()
clut3_jags[, 1:8, 2] <- select(clut3_jags_temp,
                               starts_with('2017')) %>%
  simplify2array()
clut3_jags[, 1:7, 3] <- select(clut3_jags_temp,
                               starts_with('2018')) %>%
  simplify2array()
clut3_jags[, 1:7, 4] <- select(clut3_jags_temp,
                               starts_with('2019')) %>%
  simplify2array()
clut3_jags[, 1:7, 5] <- select(clut3_jags_temp,
                               starts_with('2020')) %>%
  simplify2array()

#Water ind
wind3_jags_temp <- data.frame(conus_id = nw_nights_all$conus_id,
                              year = nw_nights_all$year,
                              wind = vmat3[, 'water_ind2']) %>%
  group_by(conus_id, year) %>%
  mutate(visit = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = c(year, visit),
              values_from = wind)

wind3_jags <- array(NA, dim = c(140, 8, 5))
wind3_jags[, 1:6, 1] <- select(wind3_jags_temp,
                               starts_with('2016')) %>%
  simplify2array()
wind3_jags[, 1:8, 2] <- select(wind3_jags_temp,
                               starts_with('2017')) %>%
  simplify2array()
wind3_jags[, 1:7, 3] <- select(wind3_jags_temp,
                               starts_with('2018')) %>%
  simplify2array()
wind3_jags[, 1:7, 4] <- select(wind3_jags_temp,
                               starts_with('2019')) %>%
  simplify2array()
wind3_jags[, 1:7, 5] <- select(wind3_jags_temp,
                               starts_with('2020')) %>%
  simplify2array()

n_visits3_jags <- n_visits_df3b %>%
  ungroup() %>%
  select(-conus_id) %>%
  simplify2array()




###################################################################################
#Save everything in a list
occ_data_jags3 <- list('dets' = dets3_jags,
                       'tmin' = tmin3_jags,
                       'dayl' = dayl3_jags,
                       'clut' = clut3_jags,
                       'wind' = wind3_jags,
                       'n_sites' = n_sites_total3,
                       'xmat' = xmat3a,
                       'n_xcovs' = n_xcovs3,
                       'n_visits' = n_visits3_jags,
                       'n_years' = 5)
#####################################
#Remove dups manually for eastside masking

xmat3a<-xmat3a[-c(100,104),]
n_sites_total3<-nrow(xmat3a)

######################################

#version without n_xcovs for use w modified models I made for informed priors and altered trend structure
occ_data_jags3 <- list('dets' = dets3_jags,
                       'tmin' = tmin3_jags,
                       'dayl' = dayl3_jags,
                       'clut' = clut3_jags,
                       'wind' = wind3_jags,
                       'n_sites' = n_sites_total3,
                       'xmat' = xmat3a,
                       'n_visits' = n_visits3_jags,
                       'n_years' = 5)
#############
#Specify initial values for z
max2 <- function(x){
  if(sum(is.na(x)) == length(x)){
    return(0)
  } else {
    return(max(x, na.rm = TRUE))
  }
}

jags_inits3 <- list('z' = apply(dets3_jags, c(1, 3), max2))

############################
#JAGS model

#alphas are tmin,dayl,clutter,w-ind
#betas 1, 2, 3 4 are elev, elev sd, precip, log forest
#5th is cliffs for those spp
# then replace forest # 4 w snags for those models 

library(rjags)

#Fit model to the  data
occ_adapt_ex3 <- jags.model(file = 'Tutorial/occ_model3a_MYEVa.jags',
                            data = occ_data_jags3,
                            n.chains = 3,
                            n.adapt = 1000,
                            inits = jags_inits3)

occ_jags_ex3 <- coda.samples(model = occ_adapt_ex3,
                             variable.names = c('alpha0','gam','phi','alpha4',
                                                'alpha1','alpha2','alpha3',
                                                'beta0', 'betas',
                                                'psi.hat','p.hat','lam.tot','lam.avg','turnover',
                                                'surv','col','avg.psi','lam.tot.avg','lam.avg.avg'),
                             n.iter = 5000,thin=3)

summary(occ_jags_ex3)
###########################################
#some diagnostics
library(MCMCvis)

plot(occ_jags_ex3)
heidel.diag(occ_jags_ex3)
raftery.diag(occ_jags_ex3)
effectiveSize(occ_jags_ex3)

#Can compare priors and posteriors
prior<-c(rnorm(1000,0,32) #dnorm(0,0.001) in JAGS

prior.beta1<-rnorm(5000,-0.52,0.29)
prior.beta2<-rnorm(5000,-0.08,0.21)
prior.beta3<-rnorm(5000,-0.41,0.30)
prior.beta4<-rnorm(5000,0.62,0.26)  
priors<-cbind(prior.beta1,prior.beta2,prior.beta3,prior.beta4)   
#prior.u<-runif(1000,0,1)
MCMCtrace(occ_jags_ex3,params='betas',priors=priors,Rhat=TRUE,n.eff=TRUE,pdf=FALSE)

#MCMCtrace(out,params=c('p'),priors=prior.u,Rhat=TRUE,n.eff=TRUE)

params=c("col[2]","col[3]","col[4]","col[5]",
         "psi.hat[1]","psi.hat[2]","psi.hat[3]","psi.hat[4]","psi.hat[5]")

#Nice posterior plot, no topography (beta8) for plotting scale
MCMCplot(object=occ_jags_ex3,params=c('col','psi.hat'),
         offset=0.1,col="black",
         labels=c("Colonization 2017","Colonization 2018","Colonization 2019","Colonization 2020",
                  "Pr[Occupancy] 2016","Pr[Occupancy] 2017","Pr[Occupancy] 2018","Pr[Occupancy] 2019","Pr[Occupancy] 2020"),
         main="Hoary Bat trend 2016-2020 [2010 priors]
         posterior estimates and 95% intervals",ref_ovl=TRUE)