library(tidyverse)
library(here)
library(spOccupancy)

fit <- readRDS(here("DataProcessed/results/spOccupancy/full/fits/epfu_fit.rds"))
bat.dat <- readRDS(here("DataProcessed/results/spOccupancy/full/bat_dat.rds"))

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
  print(ggplot( occ_binned(fit, cov_vec,50, 6) ) +
          geom_point(aes(x = cov_vec, y = resid_raw), color = "blue") + geom_hline(yintercept = 0, lty = 2) +
          theme_bw() +
          facet_wrap(~ Iteration, labeller = "label_both") +
          labs( y = "Binned occurrence residuals", x = occ.covnames[[i]],
                title = "Epfu Binned Occurrence Residuals"))
}


