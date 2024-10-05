## ---------------------------
##
## Script name: 02a_jagsMod_diagnostics
##
## Purpose of script: Derive desired parameters from JAGS output and check diagnostics
##
## Author: Trent VanHawkins
##
## Date Created: 2024-10-01
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
library(rjags)
library(bayesplot)
library(coda)
library(MCMCvis)
# Read in the desired files -----------------------------------------------

# Read in Data
## Full 
dir <- "DataProcessed/results/jags/full/fits/"

filenames <- list.files(here(dir))

#read in the files into a list object
occ_jags <- lapply(filenames,function (x) readRDS(here(paste0(dir, x))))

# rename the list elements to have the correct names
names(occ_jags) <- str_split_i(filenames, "_", 1)

## Sensitivity
dir.sens <- "DataProcessed/results/jags/ORWA_Only/fits/"

filenames.sens <- list.files(here(dir.sens))

#read in the files into a list object
occ_jags_sens <- lapply(filenames.sens,function (x) readRDS(here(paste0(dir.sens, x))))

# rename the list elements to have the correct names
names(occ_jags_sens) <- str_split_i(filenames.sens, "_", 1)

## define some species we will exclude
exclude <- c("tabr")

## Cliff Bats
cliff_bats <- c("anpa", "euma", "myci", "pahe")

## Read in data objects
occ_data <- readRDS(here("DataProcessed/results/jags/occ_data.rds"))
occ_data_sens <- readRDS(here("DataProcessed/results/jags/occ_data_sens.rds"))

n_sites <- occ_data$laci$n_sites
n_sites_sens <- occ_data_sens$laci$n_sites
# Extract mcmc draws and convert to dataframe for summarization and plotting (outputs all.jags which contains all mcmc draws of all params for all species)
mcmc_extract <- function(occ_jags){
  for(i in 1:length(occ_jags)){
    cur.fit <- occ_jags[[i]]
    spp <- names(occ_jags)[i]
    
    for(j in 1:length(cur.fit)){
      tmp <- cur.fit[[j]] %>% 
        as.data.frame() %>% 
        mutate(iter = seq(1:nrow(cur.fit[[j]])),
               chain = j,
               spp = spp)
      if(j == 1){dat <- tmp}
      else{dat <- bind_rows(dat, tmp)}
    }
    if(i == 1){all.jags <- dat}
    else{all.jags <- bind_rows(all.jags, dat)}
  }
  return(all.jags)
}
## Use our function to extract our data
all.jags <- bind_rows(mcmc_extract(occ_jags) %>% mutate(analysis = "OR|WA|ID"), mcmc_extract(occ_jags_sens) %>% mutate(analysis = "OR|WA Only"))
  
# Predicting and creating maps --------------------------------------------

# Occurrence Covariates ---------------------------------------------------
nw_grid_shp <- readRDS(here("DataProcessed/occurrence/nw_grid_shp.rds"))

library(sf)
#divide shapefile into sampled and unsampled grid cells.
#also arrange these by conus_id
nw_grida <- nw_grid_shp %>%
  st_drop_geometry() %>% 
  filter(samp_all == 1) %>%
  arrange(cell)
nw_gridb <- nw_grid_shp %>%
  st_drop_geometry() %>% 
  filter(samp_all == 0) %>%
  arrange(cell)

#combine these back again, reordered. calculate scaled versions of each
#covariate and save these in matrices that will be used for predictions
nw_grid_all <- rbind(nw_grida, nw_gridb)

xmat_all <- nw_grid_all %>%
  st_drop_geometry() %>%
  mutate(log_fc = log(p_forest + 1)) %>% 
  select(log_fc, precip, DEM_max) %>%
  scale() %>% 
  as.matrix()

xmat_cliff <- nw_grid_all %>%
  st_drop_geometry() %>%
  mutate(log_fc = log(p_forest + 1),
         log_cliff = log(cliff_cover*100 + 1)) %>% 
  select(log_fc, precip, DEM_max, log_cliff) %>%
  scale() %>% 
  as.matrix()

n_years <- occ_data$myci$n_years

get_occ_post <- function(jags, xmat, n_years, grid){

  # Extract Occurrence Model Params. (Posterior Distributions)-----
  for (i in 1:length(unique(jags$spp))) {
    alphas <- jags %>% 
      filter(spp == unique(jags$spp)[i], analysis == "OR|WA|ID") %>% 
      select(dplyr::contains("alpha")) 
    
    if (!unique(jags$spp)[i] %in% c("anpa", "euma", "myci", "pahe")) {
      alphas <- alphas %>% select(-`alphas[4]`)
    }
    
    alphas <- as.matrix(alphas)
    gamma <- jags %>% 
      filter(spp == unique(jags$spp)[i], analysis == "OR|WA|ID") %>% 
      select(dplyr::contains("gamma")) %>% 
      as.matrix()
    
    phi <- jags %>% 
      filter(spp == unique(jags$spp)[i], analysis == "OR|WA|ID") %>% 
      select(dplyr::contains("phi")) %>% 
      as.matrix()
    
    ## initialize list for psi for each year
    psi_post <- list()
    
    #Summarize for each grid cell and each year
    for (j in 1:n_years) {
      if (j == 1) {
        psi_post[[j]] <- plogis(cbind(1, xmat) %*% t(alphas))
      }
      else{
        
        #Combine colonization and alpha parameters to get posterior
        alpha_post_col <- cbind(gamma[,j-1], alphas[,-1])
        alpha_post_surv <- cbind(phi[,j-1], alphas[,-1])
        
        #Get the posteriors
        psi_post_col <- plogis(cbind(1, xmat) %*% t(alpha_post_col))
        psi_post_surv <- plogis(cbind(1, xmat) %*% t(alpha_post_surv))
        
        psi_post[[j]] <- psi_post[[j-1]] * psi_post_surv + (1 - psi_post[[j-1]]) * psi_post_col
      }
    }
    
    #Posterior summary of Psi
    psi_summ <- lapply(psi_post, function(x){data.frame(mean = apply(x, 1, mean),
                                                        q2.5 = apply(x, 1, quantile, probs = 0.025),
                                                        q25 = apply(x, 1, quantile, probs = 0.25),
                                                        q75 = apply(x, 1, quantile, probs = 0.75),
                                                        q97.5 = apply(x, 1, quantile, probs = 0.975),
                                                        spp = unique(jags$spp)[i],
                                                        cell = grid$cell)})
    #Add in year column
    year.match <- 2016
    for(k in 1:length(psi_summ)){
      psi_summ[[k]]$year <- year.match
      
      year.match <- year.match + 1
    }
    #Merge into one big data frame
    if(i == 1){
      all_psi_summ <- psi_summ
    }else{
      all_psi_summ <- bind_rows(all_psi_summ, psi_summ)
    }
  }
  return(all_psi_summ)
}


# Run the function --------------------------------------------------------

## All non-cliff bats
psi.post <- get_occ_post(all.jags %>% filter(!spp %in% cliff_bats), xmat = xmat_all, n_years = n_years, grid = nw_grid_all)

## Cliff bats
psi.post.cliff <- get_occ_post(all.jags %>% filter(spp %in% cliff_bats), xmat = xmat_cliff, n_years = n_years, grid = nw_grid_all)

## Merge them
psi.post.all <- bind_rows(psi.post, psi.post.cliff) %>% mutate(width = q97.5 - q2.5)

psi.post.map <- left_join(nw_grid_shp, psi.post.all, by = "cell")

saveRDS(psi.post.map, here("DataProcessed/results/jags/psi_post_map.rds"))
# Pivot Longer and Create Traceplots --------------------------------------
all.jags.long <- all.jags %>% 
  pivot_longer(cols = -c(iter, chain, spp))
## Alphas
for (i in 1:length(occ_jags)) {
  
  spp.tmp <- names(occ_jags)[i]
  
  alpha_traceplot <- all.jags.long %>%
    filter(spp == spp.tmp,
           str_detect(name, 'alpha')) %>% 
    ggplot(aes(y = value, x = iter))+
    geom_line(aes(color = chain))+
    facet_wrap(~name, scales = "free", ncol = 1)+
    labs(title = paste0(toupper(spp.tmp), " MCMC Traceplot (Occurrence Model)"),
         x = "Iteration",
         y = "Value",
         color = "Chain")
  
  ggsave(filename = paste0(spp.tmp, "_traceplot.png"), plot = alpha_traceplot, device = "png", path = here("DataProcessed/results/jags/full/plots/traceplots/occurrence/"), width = 3840, height = 2160, units = "px", dpi = "retina")
}

## Betas
for (i in 1:length(occ_jags)) {
  
  spp.tmp <- names(occ_jags)[i]
  
  beta_traceplot <- all.jags.long %>%
    filter(spp == spp.tmp,
           str_detect(name, 'beta')) %>% 
    ggplot(aes(y = value, x = iter))+
    geom_line(aes(color = chain))+
    facet_wrap(~name, scales = "free", ncol = 1)+
    labs(title = paste0(toupper(spp.tmp), " MCMC Traceplot (Detection Model)"),
         x = "Iteration",
         y = "Value",
         color = "Chain")
  
  ggsave(filename = paste0(spp.tmp, "_traceplot.png"), plot = beta_traceplot, device = "png", path = here("DataProcessed/results/jags/full/plots/traceplots/detections/"), width = 3840, height = 2160, units = "px", dpi = "retina")
}

