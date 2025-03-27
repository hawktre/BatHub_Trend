## ---------------------------
##
## Script name: 02c_EstimateTrend.R
##
## Purpose of script: Estimate the derived trend parameter
##
## Author: Trent VanHawkins
##
## Date Created: 2025-03-25
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## ---------------------------

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


# Format the data for computation of the trend parameter ------------------
mean_post_long <- all.jags %>% 
  select(contains("avg.psi"), chain, iter, spp, analysis) %>% 
  pivot_longer(cols = -c(iter, chain, spp, analysis)) %>% 
  mutate(year = as.integer(str_extract(name, "\\d+")))


# Step 2: Compute slope
slopes <- mean_post_long %>%
  group_by(spp, analysis, chain, iter) %>%  
  summarise(
    t_bar = mean(year),
    psi_bar = mean(value),
    numerator = sum((year - t_bar) * (value - psi_bar)),
    denominator = sum((year - t_bar)^2),
    slope = numerator / denominator
  ) %>% 
  ungroup()

# Step 3: Summarise slope posterior
slopes_summary <- slopes %>% 
  group_by(spp, analysis) %>% 
  summarise(mean_slope = mean(slope),
            sd_slope = sd(slope),
            q2.5 = quantile(slope, probs = 0.025),
            q25 =  quantile(slope, probs = 0.25),
            q50 =  quantile(slope, probs = 0.5),
            q75 =  quantile(slope, probs = 0.75),
            q97.5 =  quantile(slope, probs = 0.975))


# Plot it -----------------------------------------------------------------

slopes_summary %>%
  ggplot(aes(x = spp, y = mean_slope, group = analysis)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, linewidth = 0.75, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = analysis), width = 0, linewidth = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_minimal() + 
  labs(x = "Species",
       y = expression("Posterior Distribution" ~ hat(beta)),
       color = "Analysis")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle("Estiamted Trend (Linear) 2016-2022")  
  