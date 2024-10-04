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

dir <- "DataProcessed/results/jags/full/fits/"

filenames <- list.files(here(dir))

#read in the files into a list object
occ_jags <- lapply(filenames,function (x) readRDS(paste0(dir, x)))

# rename the list elements to have the correct names
names(occ_jags) <- str_split_i(filenames, "_", 1)


# Extract mcmc draws and convert to dataframe for summarization and plotting (outputs all.jags which contains all mcmc draws of all params for all species)
for(i in 1:length(occ_jags)){
  cur.fit <- occ_jags[[i]]
  spp <- names(occ_jags)[i]
  
  for(j in 1:length(cur.fit)){
    tmp <- cur.fit[[j]] %>% 
      as.data.frame() %>% 
      mutate(iter = seq(1:nrow(cur.fit[[j]])),
             chain = j,
             spp = spp)
    if(j == 1){
      dat <- tmp 
    }
    else{
      dat <- bind_rows(dat, tmp)
    }
  }
  if(i == 1){
    all.jags <- dat
  }
  else{
    all.jags <- bind_rows(all.jags, dat)
  }
}
  


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

all.jags <- all.jags %>% 
  rename("tmin" = beta1,
         "daylight" = beta2,
         "clutter0" = beta3,
         "clutter1" = beta4,
         "clutter2" = beta5,
         "clutter3" = beta6,
         "water" = beta7,
         "log_fc" = `alphas[1]`,
         "precip" = `alphas[2]`,
         "elev" = `alphas[3]`,
         "log_cliff" = `alphas[4]`)

# Posterior Summaries (Coefficients)
param_summ <- all.jags %>%
  select(tmin, daylight, clutter0, clutter1, clutter2, clutter3,
         water, log_fc, precip, elev, log_cliff, spp) %>% 
  pivot_longer(cols = -spp, names_to = "param") %>% 
  group_by(spp, param) %>% 
  summarise(mean = mean(value),
            q2.5 = quantile(value, probs = c(0.025), na.rm = T),
            q25 = quantile(value, probs = c(0.25), na.rm = T),
            q50 = quantile(value, probs = c(0.5), na.rm = T),
            q75 = quantile(value, probs = c(0.75), na.rm = T),
            q97.5 = quantile(value, probs = c(0.975), na.rm = T)) %>% 
  ungroup() %>% 
  mutate(type = if_else(param %in% c("log_fc", "precip", "elev", "log_cliff"), "Occurrence", "Detection"))

saveRDS(param_summ, here("DataProcessed/results/jags/full/param_summ.rds"))

## Plot of Occurrence Coefficients
param_summ %>% 
  filter(type == "Occurrence") %>% 
ggplot(aes(x = spp, y = mean))+
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, size = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0,
                size = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_grid(param ~ ., scales = 'free') +
  geom_hline(yintercept = 0, lty = 2) +
  xlab('Species') + 
  ylab('Posterior Distribution (Log-Odds)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  ggtitle('Comparing occupancy coefficients (Mean & 95% CI)',
          'Multi-season, single-species model')+
  labs(color = "Legend") 

## Plot of Detection Coefficients
param_summ %>% 
  filter(type == "Detection") %>% 
  ggplot(aes(x = spp, y = mean))+
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, size = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0,
                size = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_grid(param ~ ., scales = 'free') +
  geom_hline(yintercept = 0, lty = 2) +
  xlab('Species') + 
  ylab('Posterior Distribution (Log-Odds)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  ggtitle('Comparing Detection coefficients (Mean & 95% CI)',
          'Multi-season, single-species model')+
  labs(color = "Legend") 


# psi summary -------------------------------------------------------------
psi_summ <- all.jags %>% 
  select(contains("psi"), spp) %>% 
  pivot_longer(-spp) %>% 
  mutate(year = case_when(str_detect(name, '1')~2016,
                          str_detect(name, '2')~2017,
                          str_detect(name, '3')~2018,
                          str_detect(name, '4')~2019,
                          str_detect(name, '5')~2020,
                          str_detect(name, '6')~2021,
                          str_detect(name, '7')~2022),
         type = if_else(str_detect(name, 'avg'), 'Average', 'Hat')) %>% 
  group_by(spp, type, year) %>% 
  summarise(mean = mean(value),
            q2.5 = quantile(value, probs = c(0.025), na.rm = T),
            q25 = quantile(value, probs = c(0.25), na.rm = T),
            q50 = quantile(value, probs = c(0.5), na.rm = T),
            q75 = quantile(value, probs = c(0.75), na.rm = T),
            q97.5 = quantile(value, probs = c(0.975), na.rm = T)) %>% 
  ungroup()

saveRDS(psi_summ, here("DataProcessed/results/jags/full/psi_summ.rds"))

psi_summ %>% 
  filter(type == "Average") %>% 
  ggplot(aes(x = year, y = mean))+
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, size = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0,
                size = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_grid(spp ~ ., scales = 'free') +
  xlab('Species') + 
  ylab('Posterior Distribution (Log-Odds)') +
  ylim(c(0,1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  ggtitle('Occurrence Probabilities Over Time',
          'Multi-season, single-species model')+
  labs(color = "Legend") 
