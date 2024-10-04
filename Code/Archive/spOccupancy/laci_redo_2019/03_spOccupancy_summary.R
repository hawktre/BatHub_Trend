## ---------------------------
##
## Script name: 01a_spOccupancy_summary.R
##
## Purpose of script: Summarise Results and Diagnostics for models 
##
## Author: Trent VanHawkins
##
## Date Created: 2024-07-25
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
bat.dat <- readRDS(here("DataProcessed/results/spOccupancy/full/bat_dat.rds"))


# Trace Plots -------------------------------------------------------------
trc_plt <- function(fit, spp){
  n.chains <- fit$n.chains
  #Extract occurrence mod samples
  occ.samps <- as.data.frame(fit$beta.samples)
  det.samps <- as.data.frame(fit$alpha.samples)
  
  #Extract posterior sample index and chain
  ## Occurrence
  occ.samps$x <- rep(1:(nrow(occ.samps)/n.chains), n.chains)
  occ.samps$chain <- as.factor(rep(1:4, each =nrow(occ.samps)/n.chains))
  
  ## Detection
  det.samps$x <- rep(1:(nrow(det.samps)/n.chains), n.chains)
  det.samps$chain <- as.factor(rep(1:4, each =nrow(det.samps)/n.chains))
  
  # convert to long for plotting
  ## Occurrence
  names(occ.samps)[1] <- "Intercept" 
  occ.samps.long <- gather(
    occ.samps, parameter, value,
    Intercept:DEM_max, factor_key = T )
  
  ## Detection
  names(det.samps)[1] <- "Intercept"
  det.samps.long <- gather(det.samps, parameter, value, Intercept:water, factor_key = T )
  
  occ_plt <- ggplot(occ.samps.long) + geom_line(
    aes(x = x, y = value,
        col = chain, group = chain)) +
    labs(x ="Iteration") +
    facet_wrap(~ parameter, nrow = length(unique(occ.samps.long$parameter)), scales = "free_y")+
    theme_bw()
  
  
  det_plt <- ggplot(det.samps.long) + geom_line(
    aes(x = x, y = value,
        col = chain, group = chain)) +
    labs(x ="Iteration") +
    facet_wrap(~ parameter, nrow = length(unique(det.samps.long$parameter)), scales = "free_y")+
    theme_bw()
  
  ggsave(filename = paste0(spp, "_occ_traceplot.png"), plot = occ_plt, path = here("DataProcessed/results/spOccupancy/full/figures/traceplots/Occurrence/"), width = 3024, height = 1964, units = "px")
  
  ggsave(filename = paste0(spp, "_det_traceplot.png"), plot = det_plt, path = here("DataProcessed/results/spOccupancy/full/figures/traceplots/Detection/"), width = 3024, height = 1964, units = "px")
}


# Moran's I (within year) -------------------------------------------------
correl_fun <- function(dist.mat, resid.vec, dist.incr){
  # function to calculate moran's from Wright et al 2021
  n <- length(resid.vec)
  resid.bar <- mean(resid.vec)
  resid.var <- var(resid.vec)
  resid.scale <- resid.vec - resid.bar
  resid.pairs <- outer(resid.scale, resid.scale)  
  resid.pairs <- resid.pairs[lower.tri(resid.pairs)] 
  
  dgrp.mat <- ceiling(dist.mat / dist.incr) 
  dgrp.mat <- dgrp.mat[lower.tri(dgrp.mat)]
  
  dist.mat <- dist.mat[lower.tri(dist.mat)]

  moran <- tapply(resid.pairs, dgrp.mat, mean, na.rm = TRUE) * n / (n - 1) / resid.var 
  dist.means <- tapply(dist.mat, dgrp.mat, mean, na.rm = TRUE)
  return(list(dist = dist.means, moran = moran))
}


morans_occ <- function(
    out, #fitted model output to get residuals from
    data,
    niter){ #number of randomly selected MCMC samples to perform this over
  ###  build occurrence distance matrix
  occ_dist_m <- data$coords %>% 
    dist() %>%
    as.matrix
  
  ## convert to km
  occ_dist_km <- occ_dist_m / 1000
  
  # create storage
  out_ <- list()
  out_by_year <- list()
  for (yr in 1:dim(out$z.samples)[[3]]) {
    
    z.rep <- as.data.frame(out$z.samples[,,yr]) # samples of latent z state across all sites 
    psi <- as.data.frame(out$psi.samples[,,yr]) # samples of latent psi across all sites 
    
    # select iterations
    total_samples <- nrow(z.rep)
    iters <- sample(1:total_samples, size = niter, replace = F)
    
    # loop through iterations and produce occurrence residual vector at one mcmc iter: 
    for(ndx in 1:niter){    
      
      iter_ <- iters[ndx]
      
      # get residual at ONE mcmc sample:
      resid_raw_occ <- as.numeric(z.rep[iter_, ] - psi[iter_, ])
      
      # calculate moran's i with Wright's function
      moran <- correl_fun(occ_dist_km, resid_raw_occ, 15) #Use distance lag of 15 km
      
      out_[[ndx]] <- tibble(
        dist = moran$dist,
        moran = moran$moran
      )  %>%
        mutate(
          iter = iter_
        )
    }
    out_by_year[[yr]] <- do.call("rbind", out_) 
  }
  names(out_by_year) <- as.character(dimnames(data$y)[[3]])
  return(out_by_year)
  
}

# Plot Moran's I ----------------------------------------------------------

morans_plt <- function(morans, type, spp){
  
  morans <- bind_rows(morans, .id = "year")
  
  morans.plt <- morans %>% 
    mutate(type = type) %>% 
    filter(dist <= 150) %>% 
    ggplot()+
    geom_line(aes(x = dist, y = moran, group = iter))+
    theme_bw()+
    geom_hline(aes(yintercept = 0), linetype = "dotdash")+
    facet_wrap(~year, scales = "free")+
    labs(y = "Moran's I",
         x = "Distance (km)")
  
  if (type == "Occurrence") {
    ggsave(filename = paste0(spp, "_", type, "_morans.png"), plot = morans.plt, path = here("DataProcessed/results/spOccupancy/full/figures/correlation/Occurrence/"), width = 3024, height = 1964, units = "px")  
  }
  else{
    ggsave(filename = paste0(spp, "_", type, "_morans.png"), plot = morans.plt, path = here("DataProcessed/results/spOccupancy/full/figures/correlation/Detection"), width = 3024, height = 1964, units = "px")  
  }
}
#
#Summary of occurrence covars --------------------------------------------

covars_summary <- function(fit, spp){
  occ_quantiles <- fit$beta.samples %>% as.data.frame() %>% 
    pivot_longer(cols = everything(), names_to = "Covariate", values_to = "value") %>% 
    group_by(Covariate) %>% 
    summarise(mean = mean(value),
              q2.5 = quantile(value, probs = c(0.025)),
              q25 = quantile(value, probs = c(0.25)),
              q75 = quantile(value, probs = c(0.75)),
              q97.5 = quantile(value, probs = c(0.975))) %>% 
    ungroup() %>% 
    mutate(type = "Occurrence",
           spp = spp)
  

  det_quantiles <- fit$alpha.samples %>% as.data.frame() %>% 
    pivot_longer(cols = everything(), names_to = "Covariate", values_to = "value") %>% 
    group_by(Covariate) %>% 
    summarise(mean = mean(value),
              q2.5 = quantile(value, probs = c(0.025)),
              q25 = quantile(value, probs = c(0.25)),
              q75 = quantile(value, probs = c(0.75)),
              q97.5 = quantile(value, probs = c(0.975))) %>% 
    ungroup() %>% 
    mutate(type = "Detection",
           spp = spp)
  return(bind_rows(occ_quantiles, det_quantiles))
}

files <- list.files(path = here("DataProcessed/results/spOccupancy/full/fits/"))

covars <- list()

for (i in 1:length(files)){
  spp <- str_split(files, "_")[[i]][[1]]
  
  cur.fit <- readRDS(here(paste0("DataProcessed/results/spOccupancy/full/fits/", files[i])))
  
  
  #traceplots
  #trc_plt(cur.fit, spp)
  
  
  # # #Morans (occurrence)
  morans.occ <- morans_occ(cur.fit, bat.dat, niter = 100)
  morans_plt(morans.occ, type = "Occurrence", spp = spp)


  #covars[[i]] <- covars_summary(cur.fit, spp)
  
}

covars.full <- bind_rows(covars)

saveRDS(covars.full, here("DataProcessed/results/spOccupancy/laci_redo_2019/covars_summ.rds"))


# Plot Covariates ---------------------------------------------------------

occ_plt <- covars.full %>% 
  filter(type == "Occurrence") %>% 
  ggplot(aes(x = Covariate, y = mean))+
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, size = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, colour = Covariate), width = 0,
                size = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~spp, ncol = 1, scales = "free")+
  geom_hline(yintercept = 0, lty = 2) +
  xlab('Species') + 
  ylab('Posterior Distribution (Log-Odds)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  ggtitle('Comparing occupancy coefficients (Mean & 95% CI)',
          'Single-season, multi-species model')+
  labs(color = "Legend") 

ggsave(filename = "all_occ_covars.png", plot = occ_plt, path = here("DataProcessed/results/spOccupancy/laci_redo_2019/figures/covariates/Occurrence"), width = 3024, height = 1964, units = "px")

det_plt <- covars.full %>% 
  filter(type == "Detection") %>% 
  ggplot(aes(x = Covariate, y = mean))+
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, size = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, colour = Covariate), width = 0,
                size = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~spp, ncol = 1, scales = 'free') +
  geom_hline(yintercept = 0, lty = 2) +
  xlab('Species') + 
  ylab('Posterior Distribution (Log-Odds)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  ggtitle('Comparing detection coefficients (Mean & 95% CI)',
          'Single-season, multi-species model')+
  labs(color = "Legend") 

ggsave(filename = "all_det_covars.png", plot = det_plt, path = here("DataProcessed/results/spOccupancy/laci_redo_2019/figures/covariates/Detection/"), width = 3024, height = 1964, units = "px")

for (i in 1:length(files)){
  spp <- str_split(files, "_")[[i]][[1]]
  
  cur.fit <- readRDS(here(paste0("DataProcessed/results/spOccupancy/laci_redo_2019/fits/", files[i])))
  
  posterior_summary <- data.frame(mean = apply(cur.fit$psi.samples, c(2,3), mean) %>% colMeans(),
                                  q2.5 = apply(cur.fit$psi.samples, c(2,3), quantile, probs = c(0.025)) %>% colMeans(),
                                  q25 = apply(cur.fit$psi.samples, c(2,3), quantile, probs = c(0.25)) %>% colMeans(),
                                  q75 = apply(cur.fit$psi.samples, c(2,3), quantile, probs = c(0.75)) %>% colMeans(),
                                  q97.5 = apply(cur.fit$psi.samples, c(2,3), quantile, probs = c(0.975)) %>% colMeans(),
                                  spp = spp,
                                  year = rep(2016:2018)
  )
  
  if(!exists("psi.hat")){
    psi.hat <- posterior_summary
  }
  else{
    psi.hat <- bind_rows(psi.hat, posterior_summary)
  }
}

saveRDS(psi.hat, here("DataProcessed/results/spOccupancy/psi_hat.rds"))
library(latex2exp)
psi.plt <- psi.hat %>% 
  mutate(year = factor(year)) %>% 
  ggplot(aes(x = year, y = mean))+
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, size = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, colour = spp), width = 0,
                size = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_wrap(~spp, ncol = 1, scales = 'free_y') +
  theme_bw()+
  ylim(c(0.5,1))+
  xlab('Species') + 
  ylab(TeX("$\\hat{psi}$")) +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ggtitle('Comparing detection coefficients (Mean & 95% CI)',
          'Single-season, multi-species model')+
  labs(color = "Legend") 
  
ggsave(filename = "all_psi_hat.png", plot = psi.plt, path = here("DataProcessed/results/spOccupancy/laci_redo_2019/figures/"), width = 3024, height = 1964, units = "px")
