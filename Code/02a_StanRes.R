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
## Stan Model
occ_stan <- readRDS(here("DataProcessed.nosync/results/stan_res.rds"))

## Bat Grid Shapefile
nw_grid_all <- st_read(here("DataProcessed.nosync/occurrence/nw_grid_all.shp"))

## Covariates model matrix (Occupancy)
xmat_all <- readRDS(here("DataProcessed.nosync/occurrence/xmat_all.rds"))

## Original Model Object
occ_data <- readRDS(here("DataProcessed.nosync/results/occ_data.rds"))


# Model Summaries ---------------------------------------------------------

#Summarize occupancy coefficients
print(summary(occ_stan, pars = c('alpha01', 'alpha02',
                                 'alphas', 'alpha_auto'),
              use_cache = FALSE)$summary, digits = 2)
traceplot(occ_stan)
#Summarize detection coefficients
print(summary(occ_stan, pars = c('beta0', 'betas'),
              use_cache = FALSE)$summary, digits = 2)


# Create a function to create summaries of the posterior predictio --------

get_occ_post <- function(occ_data, occ_stan){
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
  psi_posta <- plogis(cbind(1, xmat_all) %*% t(alpha_posta))
  psi_postb <- plogis(cbind(1, xmat_all) %*% t(alpha_postb))
  psi_postc <- plogis(cbind(1, xmat_all) %*% t(alpha_postc))
  
  ## initialize list for psi for each year
  psi_post <- list()
  
  #Summarize for each grid cell and each year
  for (i in 1:occ_data$n_years) {
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
  
  return(list("psi_post" = psi_post,
              "psi_summary" = psi_summ,
              "n_years" = occ_data$n_years))
}

myci_post <- get_occ_post(occ_data, occ_stan)

years <- paste0(2016:2022)

post_map <- function(grid, post_summ, years){
  ## Create Shell for each map
  means_map <- select(grid, cell, geometry)
  ci_map <- select(grid, cell, geometry)
  
  ## Add Means and 95% BCI for each year
  for(i in 1:post_summ$n_years){
    means_map <- cbind(means_map, post_summ$psi_summary[[i]]$mean) 
    ci_map <- cbind(ci_map, (post_summ$psi_summary[[i]]$q97.5 - post_summ$psi_summary[[i]]$q2.5))
  }
  names(means_map) <- c("cell", years, "geometry")
  names(ci_map) <- c("cell", years, "geometry")
  
  # ## Create Long Format also 
  means_map_long <- means_map %>% 
    as.data.frame() %>% 
    pivot_longer(cols = any_of(years), names_to = "years", values_to = "mean") %>% 
    st_as_sf()
  ci_map_long <- ci_map %>% 
    as.data.frame() %>% 
    pivot_longer(cols = any_of(years), names_to = "years", values_to = "width") %>% 
    st_as_sf()

  return(list("means_map" = means_map,
              "ci_map" = ci_map,
              "means_map_long" = means_map_long,
              "ci_map_long" = ci_map_long))
}

myci_map <- post_map(grid = nw_grid_all, post_summ = myci_post, years = years)

myci_map$means_map_long %>% 
  ggplot()+
  geom_sf(aes(fill = mean), size = 0.1) +
  facet_wrap(. ~ years, ncol = 2) +
  viridis::scale_fill_viridis(limits = c(0, 1)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  ggtitle('MYCI Posterior Mean Occupancy',
          'Multi-season, single-species model')

ggsave(filename = "myci_post_means.png", path = here("Reports/"), width = 32, height = 24, units = "in")


myci_map$ci_map %>% 
  pivot_longer(cols = all_of(years), names_to = "year", values_to = "Width") %>% 
  ggplot()+
  geom_sf(aes(fill = Width), lwd = 0) +
  facet_wrap(. ~ year, ncol = 4) +
  theme_bw(base_size = 16) +
  viridis::scale_fill_viridis(limits = c(0, 1), option = "magma")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  ggtitle('MYCI Occupancy 95% Interval Width',
          'Multi-season, single-species model')

ggsave(filename = "myci_post_ci.png", path = here("Reports/"), width = 3840, height = 2160, units = "px")
