## ---------------------------
##
## Script name: 01b_spOccupancy_predict.R
##
## Purpose of script: Make predictions across unsamples cells
##
## Author: Trent VanHawkins
##
## Date Created: 2024-07-29
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
library(spOccupancy)


# Read in the data --------------------------------------------------------

covars <- read_sf(here("DataProcessed.nosync/occurrence/batgrid_covars_scaled.shp"))
bat.dat <- readRDS(here("DataProcessed.nosync/spOccupancy/bat_dat.rds"))

coords <- covars %>% st_centroid() %>% st_coordinates()
# format covariates for prediction ----------------------------------------
files <- list.files(path = here("DataProcessed.nosync/spOccupancy/fits/"))
str_split(files, "_")
for (i in 1:length(files)){
  spp <- str_split(files, "_")[[i]][[1]]
  
  fit <- readRDS(here(paste0("DataProcessed.nosync/spOccupancy/fits/", files[i])))
  
  # Number of prediction sites.
  J.pred <- nrow(covars)
  # Number of prediction years.
  n.years.pred <- 2
  # Number of predictors (including intercept)
  p.occ <- ncol(fit$beta.samples) 
  
  # Get covariates and standardize them using values used to fit the model
  year.pred <- matrix(rep((c(2016, 2022) - mean(bat.dat$occ.covs$year)) / 
                            sd(bat.dat$occ.covs$year), 
                          J.pred), J.pred, n.years.pred, byrow = TRUE)
  
  if(spp %in% c("anpa", "euma", "myci", "pahe")){
    # Create three-dimensional array
    X.0 <- array(1, dim = c(J.pred, n.years.pred, p.occ))
    # Fill in the array
    # Years
    X.0[, , 2] <- year.pred
    #Forest
    X.0[, , 3] <- covars$p_forst
    
    # Precip
    X.0[, , 4] <- covars$precip
    
    # Cliff Canyon
    X.0[, , 5] <- covars$clff_cn
    
    X.0[, , 6] <- covars$DEM_max
  }
  
  else{
  # Create three-dimensional array
  X.0 <- array(1, dim = c(J.pred, n.years.pred, p.occ))
  # Fill in the array
  # Years
  X.0[, , 2] <- year.pred
  #Forest
  X.0[, , 3] <- covars$p_forst
  
  # Precip
  X.0[, , 4] <- covars$precip
  
  X.0[, , 5] <- covars$DEM_max
  }
  # Check out the structure
  str(X.0)
  
  #Year indices we wish to predict
  t.cols <- c(1,7)
  
  # Make predictions
  preds <- predict(fit, X.0 = X.0, t.cols = t.cols, ignore.RE = F, type = "occupancy")
  
  #Summarise mean for each cell at years 2016 and 2022
  year.means <- apply(preds$psi.0.samples,c(2,3), mean) %>% as.data.frame()
  names(year.means) <- c("mean_2016", "mean_2022")
  
  #Summarise credible intervals 
  all.ci <- apply(preds$psi.0.samples, c(2,3), quantile, probs = c(0.025, 0.975))
  all.width <- all.ci[2,,] - all.ci[1,,]
  
  #Construct a single data frame with all info
  preds.sp <- covars %>% mutate(mean_2016 = year.means$mean_2016,
                                mean_2022 = year.means$mean_2022,
                                width_2016 = all.width[,1],
                                width_2022 = all.width[,2]) %>% 
    st_transform(crs = "WGS84")
  
  #Long-format for plotting
  preds.plt <- preds.sp %>% pivot_longer(cols = mean_2016:width_2022, names_sep = "_", 
                                              names_to = c("type", "year"), values_to = "value")
  # Create a map of the mean posterior for each cell for each year 
  means_map <- ggplot()+
    geom_sf(data = preds.plt %>% filter(type == "mean"), aes(fill = value), size = 0.1, lwd = 0.05)+
    facet_wrap(~year)+
    viridis::scale_fill_viridis(limits = c(0,1))+
    theme(axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.x = element_blank(),
         axis.ticks.y = element_blank())+
    ggtitle(paste0(' Posterior Mean Occupancy', "(", spp, ")"),
            'Multi-season, single-species model')
  
  ggsave(filename = paste0(spp, "_post_means_map_spOcc.png"), plot = means_map, path = here("DataProcessed.nosync/maps/spOccupancy/means/"), 
         width = 3024, height = 1964, units = "px")
  
  # Create a map of the 95% CI for each cell and each year
  width_map <- ggplot()+
    geom_sf(data = preds.plt %>% filter(type == "width"), aes(fill = value))+
    facet_wrap(~year)+
    viridis::scale_fill_viridis(limits = c(0,1))+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(paste0('95% CI Posterior Mean Occupancy', "(", spp, ")"),
            'Multi-season, single-species model')
  
  ggsave(filename = paste0(spp, "_post_width_map_spOcc.png"), plot = width_map, path = here("DataProcessed.nosync/maps/spOccupancy/width/"), 
         width = 3024, height = 1964, units = "px")

}

psi.hat <- readRDS(here("DataProcessed.nosync/spOccupancy/psi_hat.rds"))

psi.plot <- psi.hat %>% 
  ggplot(aes(x = year, y = mean))+
  geom_errorbar(aes(ymin = ci_2.5, ymax = ci_97.5), width = 0, linewidth = 0.5)+
  geom_errorbar(aes(ymin = ci_25, ymax = ci_75, colour = "50% CI"), width = 0, linewidth = 1)+
  geom_point(cex = 0.5)+
  facet_grid(spp~., scales = "free")

ggsave(filename = "psi_plot.png", plot = psi.plot, path = here("Reports/spOccupancy/figures/"), width = 1964,
       height = 3024, units = "px")
 




