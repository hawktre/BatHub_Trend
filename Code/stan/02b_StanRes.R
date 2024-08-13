## ---------------------------
##
## Script name: 02a_StanRes
##
## Purpose of script: Explore restults of Stan modeling
##
## Author: Trent VanHawkins
##
## Date Created: 2024-06-18
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
library(rstan)
library(sf)
library(shinystan)

# Read in necessary data --------------------------------------------------
## Bat Grid Shapefile
nw_grid_all <- st_read(here("DataProcessed.nosync/occurrence/nw_grid_all.shp"))

## Covariates model matrix (Occupancy)
xmat_all <- readRDS(here("DataProcessed.nosync/results/stan/full/xmat_all.rds"))
xmat_cliff <- readRDS(here("DataProcessed.nosync/results/stan/full/xmat_cliff.rds"))

## cliff bats
cliff_bats <- c("anpa", "euma", "myci", "pahe")
## Original Model Object
occ_data <- readRDS(here("DataProcessed.nosync/results/stan/full/occ_data.rds"))

## Stan Model
filenames <- list.files(path = here("DataProcessed.nosync/results/stan/full/fits/"), pattern = "*.rds", 
                        full.names = T)

occ_stan <- lapply(filenames,function(x) readRDS(here(x)))

names(occ_stan) <- sort(names(occ_data))

## Reorder occ_data to match occ_stan (alphabetically)
occ_data <- occ_data[names(occ_stan)]

all(names(occ_data) == names(occ_stan))  

dim(rstan::extract(occ_stan$anpa, 'logit_psi')$logit_psi)
# Summarise Alphas ---------------------------------------------------------
## Extract all model parameters into a dataframe 
### Alphas ###
for (i in 1:length(occ_stan)) {
  ## Extract all model parameters into a dataframe 
  ## Intercepts
  alpha01_post <- rstan::extract(occ_stan[[i]], 'alpha01')$alpha01
  alpha02_post <- rstan::extract(occ_stan[[i]], 'alpha02')$alpha02
  ## Covars
  alphas_post <- rstan::extract(occ_stan[[i]], 'alphas')$alphas
  
  ## Autoregressive Param
  alpha_auto_post <- rstan::extract(occ_stan[[i]], 'alpha_auto')$alpha_auto
  
  tmp_alphas <- cbind(alpha01_post, alpha02_post, alphas_post, alpha_auto_post) %>% as.data.frame()
  
  if(names(occ_stan)[i] %in% c("anpa", "euma", "myci", "pahe")){
  names(tmp_alphas) <- c("int01", "int02", "forest_cover", "precip", "elevation", "cliff_cover", "alpha_auto")
  }
  else{
    names(tmp_alphas) <- c("int01", "int02", "forest_cover", "precip", "elevation", "alpha_auto")
  }
  
  tmp_alphas <- tmp_alphas %>% mutate(spp = names(occ_stan)[i])
  if(i == 1){all_alphas <- tmp_alphas}
  else{all_alphas <- bind_rows(all_alphas, tmp_alphas)}
}

# Plot alphas
alpha_summ <- all_alphas %>% 
  pivot_longer(cols = -spp, names_to = "alpha", values_to = "value") %>% 
  group_by(spp, alpha) %>% 
  summarise(mean = mean(value),
            q2.5 = quantile(value, probs = c(0.025), na.rm = T),
            q25 = quantile(value, probs = c(0.25), na.rm = T),
            q50 = quantile(value, probs = c(0.5), na.rm = T),
            q75 = quantile(value, probs = c(0.75), na.rm = T),
            q97.5 = quantile(value, probs = c(0.975), na.rm = T)) %>% 
  ungroup() %>% 
  mutate(alpha = factor(alpha, levels = c("int01", "int02", "forest_cover", "precip", "elevation", "cliff_cover", "alpha_auto")))

saveRDS(alpha_summ, here("DataProcessed.nosync/results/stan/full/alpha_summ.rds"))

# Summarise Betas ---------------------------------------------------------
for (i in 1:length(occ_stan)) {
  ## Extract all model parameters into a dataframe 
  ## Intercepts
  beta0_post <- rstan::extract(occ_stan[[i]], 'beta0')$beta0
  ## Covars
  betas_post <- rstan::extract(occ_stan[[i]], 'betas')$betas
  
  tmp_betas <- cbind(beta0_post,betas_post) %>% as.data.frame()
  names(tmp_betas) <- c("Intercept","tmin", "daylight", "clutter", "water")
  
  tmp_betas <- tmp_betas %>% mutate(spp = names(occ_stan)[i])
  if(i == 1){all_betas <- tmp_betas}
  else{all_betas <- bind_rows(all_betas, tmp_betas)}
}

# Plot betas
beta_summ <- all_betas %>% 
  pivot_longer(cols = -spp, names_to = "betas", values_to = "value") %>% 
  group_by(spp, betas) %>% 
  summarise(mean = mean(value),
            q2.5 = quantile(value, probs = c(0.025)),
            q25 = quantile(value, probs = c(0.25)),
            q50 = quantile(value, probs = c(0.5)),
            q75 = quantile(value, probs = c(0.75)),
            q97.5 = quantile(value, probs = c(0.975))) %>% 
  ungroup() %>% 
  mutate(betas = factor(betas, levels = c("Intercept","tmin", "daylight", "clutter", "water")))

saveRDS(beta_summ, here("DataProcessed.nosync/results/stan/full/beta_summ.rds"))

# Summarise Psi -----------------------------------------------------------
# Create a function to create summaries of the posterior predictio --------

get_occ_post <- function(occ_stan, xmat, n_years){
  # Extract Occurrence Model Params. (Posterior Distributions)-----
  ## Intercepts
  alpha01_post <- rstan::extract(occ_stan, 'alpha01')$alpha01
  alpha02_post <- rstan::extract(occ_stan, 'alpha02')$alpha02
  ## Covars
  alphas_post <- rstan::extract(occ_stan, 'alphas')$alphas
  
  ## Autoregressive Param
  alpha_auto_post <- rstan::extract(occ_stan, 'alpha_auto')$alpha_auto
  
  # Combine posteriors and calculate posterior predictions ------------------
  alpha_posta <- cbind(alpha01_post, alphas_post)
  alpha_postb <- cbind(alpha02_post, alphas_post)
  alpha_postc <- cbind(alpha02_post + alpha_auto_post,
                       alphas_post)

  
  ##Add intercept column to x-matrix and multiply, then transform
  psi_posta <- plogis(cbind(1, xmat) %*% t(alpha_posta))
  psi_postb <- plogis(cbind(1, xmat) %*% t(alpha_postb))
  psi_postc <- plogis(cbind(1, xmat) %*% t(alpha_postc))
  
  ## initialize list for psi for each year
  psi_post <- list()
  
  #Summarize for each grid cell and each year
  for (i in 1:n_years) {
    if (i == 1) {
      psi_post[[i]] <- psi_posta
    }
    else{
      psi_post[[i]] <- psi_post[[i-1]] * psi_postc + (1 - psi_post[[i-1]]) * psi_postb
    }
  }
  
  #Posterior summary of Psi
  psi_summ <- lapply(psi_post, function(x){data.frame(mean = apply(x, 1, mean),
                                                      q2.5 = apply(x, 1, quantile, probs = 0.025),
                                                      q25 = apply(x, 1, quantile, probs = 0.25),
                                                      q75 = apply(x, 1, quantile, probs = 0.75),
                                                      q97.5 = apply(x, 1, quantile, probs = 0.975))})
  
  return(list("psi_summary" = psi_summ,
              "psi_post" = psi_post, 
              "n_years" = n_years))
}

# Get predicted psi for each site.  ---------------------------------------

# Extract parameters, get posterior predictions, and summarise ------------
n_years <- occ_data[[1]]$n_years

occ_post <- lapply(occ_stan[!names(occ_stan) %in% cliff_bats], function(x) {get_occ_post(x, n_years = n_years, xmat = xmat_all)})
occ_post_cliff <- lapply(occ_stan[cliff_bats], function(x) {get_occ_post(x, n_years = n_years, xmat = xmat_cliff)})

occ_all <- c(occ_post, occ_post_cliff)
# Create vector of years to name columns in resulting sf object -----------
years <- paste0(2016:2022)


# Create function to put summaries in map ---------------------------------
post_map <- function(grid, post_summ, years){
  ## Create Shell for each map
  means_map <- select(grid, 
                      cell, 
                      sm_2016,
                      sm_2017,
                      sm_2018,
                      sm_2019,
                      sm_2020,
                      sm_2021,
                      sm_2022,
                      geometry)
  ci_map <- select(grid, 
                   cell, 
                   sm_2016,
                   sm_2017,
                   sm_2018,
                   sm_2019,
                   sm_2020,
                   sm_2021,
                   sm_2022,
                   geometry)
  
  ## Add Means and 95% BCI for each year
  for(i in 1:length(years)){
    means_map <- cbind(means_map, post_summ$psi_summary[[i]]$mean) 
    ci_map <- cbind(ci_map, (post_summ$psi_summary[[i]]$q97.5 - post_summ$psi_summary[[i]]$q2.5))
  }
   names(means_map) <- c('cell', 'sm_2016', 'sm_2017', 'sm_2018', 'sm_2019', 
                         'sm_2020', 'sm_2021', 'sm_2022', years, 'geometry')
   names(ci_map) <- c('cell', 'sm_2016', 'sm_2017', 'sm_2018', 'sm_2019', 'sm_2020', 
                      'sm_2021', 'sm_2022', years, 'geometry')
  
  ## Create Long Format also
  means_map_long <- means_map %>%
    as.data.frame() %>%
    pivot_longer(cols = any_of(years), names_to = "years", values_to = "mean") %>% 
    pivot_longer(cols = starts_with("sm"), names_to = "sm_year", values_to = "sampled") %>% 
    mutate(sm_year = str_split_i(sm_year, "_", 2)) %>% 
    filter(years == sm_year) %>% 
    select(cell, years, mean, sampled, geometry) %>%
    st_as_sf()
  
  ci_map_long <- ci_map %>%
    as.data.frame() %>%
    pivot_longer(cols = any_of(years), names_to = "years", values_to = "width") %>% 
    pivot_longer(cols = starts_with("sm"), names_to = "sm_year", values_to = "sampled") %>% 
    mutate(sm_year = str_split_i(sm_year, "_", 2)) %>% 
    filter(years == sm_year) %>% 
    select(cell, years, width, sampled, geometry) %>% 
    st_as_sf()

  return(list("means_map" = means_map,
              "ci_map" = ci_map,
              "means_map_long" = means_map_long,
              "ci_map_long" = ci_map_long))
}

# Summarise in map form ---------------------------------------------------
occ_map <- lapply(occ_all, function(x) post_map(grid = nw_grid_all, 
                                                 post_summ = x, 
                                                 years = years))
  
## Save out the map file 
saveRDS(occ_map, here("DataProcessed.nosync/results/stan/full/occ_map.rds"))

for (i in 1:length(occ_map)) {
  sp_name <- toupper(names(occ_map)[i])
  
  # mean map
  means_map <-occ_map[[i]]$means_map_long %>% 
    ggplot()+
    geom_sf(aes(fill = mean), size = 0.1, lwd = 0.05) +
    facet_wrap(. ~ years, ncol = 4) +
    viridis::scale_fill_viridis(limits = c(0, 1)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(paste0(sp_name, ' Posterior Mean Occupancy'),
            'Multi-season, single-species model')
  
  ggsave(plot = means_map ,filename = paste0(sp_name, "_post_means_map.png"), path = here("DataProcessed.nosync/maps/stan/full/means/"), width = 3840, height = 2160, units = "px", dpi = "retina")
  
  # ci map
  # mean map
  ci_map <-occ_map[[i]]$ci_map_long %>% 
    ggplot()+
    geom_sf(aes(fill = width), size = 0.1, lwd = 0.05) +
    facet_wrap(. ~ years, ncol = 4) +
    viridis::scale_fill_viridis(limits = c(0, 1)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(paste0(sp_name, ' Posterior 95% Credible Interval (Width)'),
            'Multi-season, single-species model')
  
  ggsave(plot = ci_map ,filename = paste0(sp_name, "_post_ci_map.png"), path = here("DataProcessed.nosync/maps/stan/full/width/"), width = 3840, height = 2160, units = "px", dpi = "retina")
}
