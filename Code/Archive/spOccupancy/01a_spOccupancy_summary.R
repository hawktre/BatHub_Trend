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
bat.dat <- readRDS(here("DataProcessed.nosync/spOccupancy/bat_dat.rds"))
dets.sp <- read_sf(here("DataProcessed.nosync/spOccupancy/spatial_dets.shp"))


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
    Intercept:precip, factor_key = T )
  
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
  
  ggsave(filename = paste0(spp, "_occ_traceplot.png"), plot = occ_plt, path = here("Reports/spOccupancy/figures/traceplots/occurrence/"), width = 3024, height = 1964, units = "px")
  
  ggsave(filename = paste0(spp, "_det_traceplot.png"), plot = det_plt, path = here("Reports/spOccupancy/figures/traceplots/detection/"), width = 3024, height = 1964, units = "px")
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


# Detection Residuals (Spatial Autocorrelation) ---------------------------

morans_det <- function(
    out, #fitted model output to get residuals from
    spp,
    data_spatial, #occurrence data with spatial information for each detection
    niter){ #number of randomly selected MCMC samples to perform this over

  # create storage
  out_ <- list()
  out_by_year <- list()

  for (yr in 1:dim(out$y)[[2]]) {
    years <- dimnames(out$y)[[2]]
    
    #Get the design matrix for that year
    design.matrix.surveyed <- model.matrix(~ clutter + scale(tmin) + scale(daylght) + factor(watr_nd), 
                                           data = data_spatial %>% 
                                             st_drop_geometry() %>% 
                                             arrange(year,cell,rplct_d) %>% 
                                             filter(year == years[yr]))
    # select iterations
    iters <- sample(1:nrow(out$alpha.samples), size = niter, replace = F)
    
    # loop through iterations
    pb <- txtProgressBar(min = 0, max = niter, style = 3,  width = 50, char = "=") 
    for(ndx in 1:niter){
      iter_ <- iters[ndx]
      
      # construct replicated z vector
      z.rep <- out$z.samples[iter_,,yr] # samples of latent z state across all sites 
      psi <- out$psi.samples[iter_,,yr] # samples of latent psi across all sites 
      
      # have to expand across the visits so each z is replicated for its corresponding visits:
      all_combinations <- expand.grid(cell = as.numeric(dimnames(out$y)[[1]]),
                                      year = as.numeric(dimnames(out$y)[[2]]))
      
      
      # Group by cell and year, count replicates
      visits_by_cell_year <- data_spatial %>% 
        st_drop_geometry() %>% 
        select(cell, year, rplct_d) %>% 
        group_by(cell, year) %>% 
        summarize(nvisit = n()) %>% 
        arrange(year,cell) %>% 
        ungroup()
      
      # Left join to include all combinations, fill missing with 0
      nvisits <- left_join(all_combinations, visits_by_cell_year, by = c("cell", "year")) %>%
        mutate(nvisit = replace_na(nvisit, 0)) %>% 
        filter(year == years[yr]) %>% 
        select(nvisit) %>% 
        unlist() %>% unname()
      
      z.rep.exp <- rep(z.rep, nvisits)
      
      y_ <- data_spatial %>% st_drop_geometry() %>% 
        arrange(cell, year, rplct_d) %>% 
        filter(year == years[yr]) %>% 
        .[[spp]]
      
      # filter to z == 1
      keep_ndx <- which(z.rep.exp == 1)
      
      # grab y's associated with filtered z's
      y_filtered <- y_[keep_ndx]
      
      # grab p's associated with filtered z's
      alpha <- as.vector(out$alpha.samples[iter_,])
      mult <- design.matrix.surveyed[keep_ndx,] %*% alpha 
      
      # convert to probability
      p <- c(exp(mult) / (1 + exp(mult)))
      resid_raw_det <- y_filtered - p 
  
      # create distance matrix
      det_dist_m <- data_spatial %>%
        arrange(cell,year,rplct_d) %>% 
        filter(year == years[yr]) %>% 
        slice(keep_ndx) %>%
        mutate(Cell = as.character(cell)) %>%
        st_as_sf %>%
        mutate(
          x_tmp = st_coordinates(.)[,1],
          y_tmp = st_coordinates(.)[,2]
        ) %>%
        as_tibble %>%
        dplyr::select(x_tmp, y_tmp) %>%
        as.matrix %>%
        dist() %>%
        as.matrix
      det_dist_km <- det_dist_m/1000
      
      # calculate moran's i with Wright's function
      moran <- correl_fun(det_dist_km, resid_raw_det, 15)
      out_[[ndx]] <- tibble(
        dist = moran$dist,
        moran = moran$moran
      )  %>%
        mutate(
          iter = iter_
        )
      setTxtProgressBar(pb, ndx)
    }
    close(pb)
    out_by_year[[yr]] <- do.call("rbind", out_) 
  }
  names(out_by_year) <- as.character(dimnames(out$y)[[2]])
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
    ggsave(filename = paste0(spp, "_", type, "_morans.png"), plot = morans.plt, path = here("Reports/spOccupancy/figures/correlation/occurrence/"), width = 3024, height = 1964, units = "px")  
  }
  else{
    ggsave(filename = paste0(spp, "_", type, "_morans.png"), plot = morans.plt, path = here("Reports/spOccupancy/figures/correlation/detection/"), width = 3024, height = 1964, units = "px")  
  }
}
#
#Summary of occurrence covars --------------------------------------------

covars_summary <- function(fit, spp){
  occ.samps <- fit$beta.samples #this holds the full outputted samples
  #create a dataframe to store quantiles in:
  occ_quantiles <- as.data.frame(matrix(NA, ncol(occ.samps),6)) 
  colnames(occ_quantiles) <- c("covariate", "mean", "q2.5", "q25", "q75", "q97.5") 
  occ_quantiles$covariate <- colnames(occ.samps)
  if(spp %in% c("anpa", "euma", "myci", "pahe")){
  occ_quantiles$covariate <- c(
    "(Intercept)" , "Year", "Forest Cover", "Precipitation","Cliff/Canyon Cover", "Elevation"
  )
  }
  else{
    occ_quantiles$covariate <- c(
      "(Intercept)" , "Year", "Forest Cover", "Precipitation", "Elevation"
    )
  }
  occ_quantiles$mean <- colMeans(occ.samps)
  occ_quantiles[,3:6] <- t(apply(occ.samps ,2,quantile,probs=c(0.025, 0.25, 0.75, 0.975)))
  occ_quantiles <- occ_quantiles %>% mutate(type = "Occurrence",
                                            spp = spp)

  det.samps <- fit$alpha.samples
  det_quantiles <- as.data.frame(matrix(NA, ncol(det.samps),6))
  colnames(det_quantiles) <- c("covariate", "mean", "q2.5", "q25", "q75", "q97.5") 
  det_quantiles$covariate <- colnames(det.samps)
  det_quantiles$mean <- colMeans(det.samps)
  det_quantiles[,3:6] <- t(apply(det.samps ,2,quantile,probs=c(0.025, 0.25, 0.75, 0.975)))
  det_quantiles <- det_quantiles %>% mutate(type = "Detection",
                                            spp = spp)
  return(bind_rows(occ_quantiles, det_quantiles))
}

# Occurrence Residuals ----------------------------------------------------

occ_binned <- function( out, #fitted model output to get residuals from
                        cov_vec, #covariates from the data (just surveyed site-nights)
                        numbins,
                        niter){ # the number of randomly selected MCMC samples
  # create storage
  out_ <- list()
  
  # select iterations
  iters <- sample(1:nrow(out$z.samples), size = niter, replace = F)
  
  # loop through iterations
  for(ndx in 1:niter){
    iter_ <- iters[ndx]
    
    # construct replicated z vector
    z.rep <- out$z.samples[iter_,,] %>% as.vector() # samples of latent z state across all sites
    psi <- out$psi.samples[iter_,,] %>% as.vector()# samples of latent psi across all sites 
    #remember that model output and data are ordered to match. 
    
    # get residual at one mcmc sample:
    resid_raw_occ <- as.numeric(z.rep - psi)    
    
    ### now do binned residuals at that iteration: 
    
    #in each bin we have:
    binsize <- floor(length(psi)/numbins) 
    bin <- rep(1:numbins, each = binsize)
    leftover <- length(psi)%%numbins # with some leftover (why they're approx equal sized 
    # groups) which will be added to the first group
    #(adding to first groups because anticipating having big zero groups)
    if (leftover > 0) {bin <- c(rep(1, leftover), bin)} #add leftover to first group 
    
    # put it in a table so that each covariate value is staying with it's residual:
    binned <- as.data.frame(cbind(psi, resid_raw_occ, cov_vec))
    
    # 1) order the covariates (or in this case the psi)
    binned <- binned %>% arrange(psi)
    # 2) sort into bins with approx equal number of points in each 
    binned$bin <- bin
    # 3) average within the bins 
    df <- aggregate(x = binned, by = list(binned$bin), FUN = mean)    
    df$sample <- iter_
    
    ## store: 
    out_[[ndx]] <- tibble(
      bin = df$bin,
      psi = df$psi,
      cov_vec = df$cov_vec,
      resid_raw = df$resid_raw_occ,
      Iteration = df$sample) 
  }
  
  return(do.call("rbind", out_))
}


occ.covnames <- names(bat.dat$occ.covs)[c(2:5, 7)]
for (i in 1:length(occ.covnames)) {
  cov_vec <- bat.dat$occ.covs[[occ.covnames[i]]]
  if(occ.covnames[i] == "year"){
    cov_vec <- cov_vec %>% as.vector()
  }
  else{
  cov_vec <- rep(cov_vec, 7)
  }
  print(ggplot( occ_binned(all.fits$anpa, cov_vec,100, 6) ) +
    geom_point(aes(x = cov_vec, y = resid_raw), color = "blue") + geom_hline(yintercept = 0, lty = 2) +
    theme_bw() +
    facet_wrap(~ Iteration, labeller = "label_both") +
    labs( y = "Binned occurrence residuals", x = occ.covnames[[i]],
          title = as.character(occ.covnames[i])))
}



files <- list.files(path = here("DataProcessed.nosync/spOccupancy/fits/"))

covars <- list()

for (i in 1:length(files)){
  spp <- str_split(files, "_")[[i]][[1]]
  
  cur.fit <- readRDS(here(paste0("DataProcessed.nosync/spOccupancy/fits/", files[i])))
  
  
  #traceplots
  trc_plt(cur.fit, spp)
  
  
  # #Morans (occurrence)
  # morans.occ <- morans_occ(cur.fit, bat.dat, niter = 100)
  # morans_plt(morans.occ, type = "Occurrence", spp = spp)
  # #Morans (detections)
  # morans.det <- morans_det(cur.fit, spp, dets.sp, 100)
  # morans_plt(morans.det, type = "Detection", spp = spp)
  
  
  covars[[i]] <- covars_summary(cur.fit, spp)
  
}

covars.full <- bind_rows(covars)

saveRDS(covars.full, here("DataProcessed.nosync/spOccupancy/covars_summary.rds"))


# Plot Covariates ---------------------------------------------------------

occ_plt <- covars.full %>% 
  filter(type == "Occurrence") %>% 
  ggplot(aes(x = spp, y = mean))+
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, size = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, colour = spp), width = 0,
                size = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_grid(covariate ~ ., scales = 'free') +
  geom_hline(yintercept = 0, lty = 2) +
  xlab('Species') + 
  ylab('Posterior Distribution (Log-Odds)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  ggtitle('Comparing occupancy coefficients (Mean & 95% CI)',
          'Single-season, multi-species model')+
  labs(color = "Legend") 

ggsave(filename = "all_occ_covars.png", plot = occ_plt, path = here("Reports/spOccupancy/figures/covariates/occurrence/"), width = 3024, height = 1964, units = "px")

det_plt <- covars.full %>% 
  filter(type == "Detection") %>% 
  ggplot(aes(x = spp, y = mean))+
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, size = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, colour = spp), width = 0,
                size = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_grid(covariate ~ ., scales = 'free') +
  geom_hline(yintercept = 0, lty = 2) +
  xlab('Species') + 
  ylab('Posterior Distribution (Log-Odds)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  ggtitle('Comparing detection coefficients (Mean & 95% CI)',
          'Single-season, multi-species model')+
  labs(color = "Legend") 

ggsave(filename = "all_det_covars.png", plot = det_plt, path = here("Reports/spOccupancy/figures/covariates/detection/"), width = 3024, height = 1964, units = "px")

for (i in 1:length(files)){
  spp <- str_split(files, "_")[[i]][[1]]
  
  cur.fit <- readRDS(here(paste0("DataProcessed.nosync/spOccupancy/fits/", files[i])))
  
  posterior_summary <- data.frame(mean = apply(cur.fit$psi.samples, c(2,3), mean) %>% colMeans(),
                                  ci_2.5 = apply(cur.fit$psi.samples, c(2,3), quantile, probs = c(0.025)) %>% colMeans(),
                                  ci_25 = apply(cur.fit$psi.samples, c(2,3), quantile, probs = c(0.25)) %>% colMeans(),
                                  ci_75 = apply(cur.fit$psi.samples, c(2,3), quantile, probs = c(0.75)) %>% colMeans(),
                                  ci_97.5 = apply(cur.fit$psi.samples, c(2,3), quantile, probs = c(0.975)) %>% colMeans(),
                                  spp = spp,
                                  year = rep(2016:2022)
  )
  
  if(!exists("psi.hat")){
    psi.hat <- posterior_summary
  }
  else{
    psi.hat <- bind_rows(psi.hat, posterior_summary)
  }
}

saveRDS(psi.hat, here("DataProcessed.nosync/spOccupancy/psi_hat.rds"))

for (i in unique(covars.full$spp)) {
  occ_plt <- covars.full %>% 
    filter(type == "Occurrence", spp == i) %>% 
    ggplot(aes(x = covariate, y = mean))+
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, linewidth = 0.6)+
    geom_errorbar(aes(ymin = q25, ymax = q75, color = covariate), width = 0, linewidth = 1.5, show.legend = F)+
    geom_point()+
    theme_bw()+
    geom_hline(yintercept = 0, lty = 2, linewidth = 0.3)+
    xlab('Covariate')+
    ylab('Posterior Distribution')+
    ggtitle(paste0('Occurrence Coefficients (', toupper(i),')'),
            'Multi-season, Single-species model')
  
  ggsave(filename = paste0(i, "_occ_covars.png"), plot = occ_plt, path = here("Reports/spOccupancy/figures/covariates/occurrence/"), width = 3024, height = 1964, units = "px")
  
  det_plt <- covars.full %>% 
    filter(type == "Detection", spp == i) %>% 
    ggplot(aes(x = covariate, y = mean))+
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, linewidth = 0.6)+
    geom_errorbar(aes(ymin = q25, ymax = q75, color = covariate), width = 0, linewidth = 1.5, show.legend = F)+
    geom_point()+
    theme_bw()+
    geom_hline(yintercept = 0, lty = 2, linewidth = 0.3)+
    xlab('Covariate')+
    ylab('Posterior Distribution')+
    ggtitle(paste0('Detection Coefficients (', toupper(i),')'),
            'Multi-season, Single-species model')
  
  ggsave(filename = paste0(i, "_det_covars.png"), plot = det_plt, path = here("Reports/spOccupancy/figures/covariates/detection/"), width = 3024, height = 1964, units = "px")
  
}

