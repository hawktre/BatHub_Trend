## ---------------------------
##
## Script name: 03_ModelDiagnostics.R
##
## Purpose of script:Model Diagnostics for STAN models
##
## Author: Trent VanHawkins
##
## Date Created: 2024-09-28
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
library(bayesplot)
library(rstan)
library(rstanarm)

## Extract all model parameters into a dataframe 
### Alphas ###
alpha_trace <- function(occ_stan, spp){
  ## Extract all model parameters into a dataframe 
  ## Intercepts
  alpha01_post <- rstan::extract(occ_stan, 'alpha01')$alpha01
  alpha02_post <- rstan::extract(occ_stan, 'alpha02')$alpha02
  
  ## Covars
  alphas_post <- rstan::extract(occ_stan, 'alphas')$alphas
  
  ## Autoregressive Param
  alpha_auto_post <- rstan::extract(occ_stan, 'alpha_auto')$alpha_auto
  
  ## Bind into a single data frame
  tmp_alphas <- cbind(alpha01_post, alpha02_post, alphas_post, alpha_auto_post) %>% as.data.frame()
  
  ## Name the columns of the dataframe
  if(spp %in% c("anpa", "euma", "myci", "pahe")){
    names(tmp_alphas) <- c("int01", "int02", "forest_cover", "precip", "elevation", "cliff_cover", "alpha_auto")
  }
  else{
    names(tmp_alphas) <- c("int01", "int02", "forest_cover", "precip", "elevation", "alpha_auto")
  }
  
  #Get the number of chains
  n_chains <- occ_stan@sim$chains
  chain_length <- (occ_stan@sim$iter - occ_stan@sim$warmup)/occ_stan@sim$thin
  
  tmp_alphas <- tmp_alphas %>% mutate(chain = paste0("Chain ", rep(seq(1:4), each = 1000)),
                                      iter = rep(seq(1:1000), 4))
  
  alpha_traceplot <- tmp_alphas %>% 
    pivot_longer(cols = -c(chain, iter)) %>% 
    ggplot(aes(y = value, x = iter))+
    geom_line(aes(colour = chain))+
    facet_wrap(~name, scales = "free", nrow = 6)+
    labs(title = paste0(toupper(spp), " MCMC Traceplot (Occurrence Model)"),
         x = "Iteration",
         y = "Value",
         color = "Chain")
  
  return(alpha_traceplot)
}

beta_trace <- function(occ_stan, spp){
  ## Extract all model parameters into a dataframe 
  ## Extract all model parameters into a dataframe 
  ## Intercepts
  beta0_post <- rstan::extract(occ_stan, 'beta0')$beta0
  ## Covars
  betas_post <- rstan::extract(occ_stan, 'betas')$betas
  
  tmp_betas <- cbind(beta0_post,betas_post) %>% as.data.frame()
  names(tmp_betas) <- c("Intercept","tmin", "daylight", "clutter-1","clutter0", "clutter1", "clutter2", "water")
  
  #Get the number of chains
  n_chains <- occ_stan@sim$chains
  chain_length <- (occ_stan@sim$iter - occ_stan@sim$warmup)/occ_stan@sim$thin
  
  tmp_betas <- tmp_betas %>% mutate(chain = paste0("Chain ", rep(seq(1:4), each = 1000)),
                                      iter = rep(seq(1:1000), 4))
  
  betas_traceplot <- tmp_betas %>% 
    pivot_longer(cols = -c(chain, iter)) %>% 
    ggplot(aes(y = value, x = iter))+
    geom_line(aes(colour = chain))+
    facet_wrap(~name, scales = "free", nrow = 4)+
    labs(title = paste0(toupper(spp), " MCMC Traceplot (Detection Model)"),
         x = "Iteration",
         y = "Value",
         color = "Chain")
  
  return(betas_traceplot)
}

beta_trace(occ_stan, "laci")




filenames <- list.files(here("DataProcessed/results/stan/full/fits/"))
