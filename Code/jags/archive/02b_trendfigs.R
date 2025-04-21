## ---------------------------
##
## Script name: 02b_trendfigs.R
##
## Purpose of script: Create mock-up figures for presentations
##
## Author: Trent VanHawkins
##
## Date Created: 2025-02-10
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
library(latex2exp)
library(sf)
library(rjags)
library(coda)
library(gt)
library(kableExtra)

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

## Summary of psi values
psi_summ <- all.jags %>% 
  select(contains("psi"), spp, analysis) %>% 
  pivot_longer(-c(spp, analysis)) %>% 
  mutate(year = case_when(str_detect(name, '1')~2016,
                          str_detect(name, '2')~2017,
                          str_detect(name, '3')~2018,
                          str_detect(name, '4')~2019,
                          str_detect(name, '5')~2020,
                          str_detect(name, '6')~2021,
                          str_detect(name, '7')~2022),
         type = if_else(str_detect(name, 'avg'), 'Average', 'Hat')) %>% 
  group_by(spp, analysis, type, year) %>% 
  summarise(mean = mean(value),
            q2.5 = quantile(value, probs = c(0.025), na.rm = T),
            q25 = quantile(value, probs = c(0.25), na.rm = T),
            q50 = quantile(value, probs = c(0.5), na.rm = T),
            q75 = quantile(value, probs = c(0.75), na.rm = T),
            q97.5 = quantile(value, probs = c(0.975), na.rm = T)) %>% 
  ungroup()

## Summary of lambda values
lambda_summ <- all.jags %>% 
  select(contains("lam.tot"), spp, analysis) %>% 
  pivot_longer(-c(spp, analysis)) %>% 
  group_by(spp, analysis) %>% 
  summarise(mean = mean(value),
            q2.5 = quantile(value, probs = c(0.025), na.rm = T),
            q25 = quantile(value, probs = c(0.25), na.rm = T),
            q50 = quantile(value, probs = c(0.5), na.rm = T),
            q75 = quantile(value, probs = c(0.75), na.rm = T),
            q97.5 = quantile(value, probs = c(0.975), na.rm = T)) %>% 
  ungroup()

#Pick out our colors
colors <- c("Very Unstable" = '#D7191C',  # Red
            "Unstable" = '#FE9000',       # Orange
            "Less Stable" = '#FFD700',    # Yellow
            "Stable" = '#1A9641',         # Green
            "Uncertain" = '#D3D3D3')      # Light Gray

#Make a matching key
status_key <- lambda_summ %>% 
  filter(analysis == "OR|WA|ID") %>% 
  mutate(status = case_when(spp %in% c("anpa", 'pahe') ~ "Uncertain",
                            q75 > 1 ~ "Stable",
                            q75 < 1 & q97.5 > 1 ~ "Less Stable",
                            q97.5 < 1 & mean > 0.6 ~ "Unstable",
                            TRUE ~ "Very Unstable")) %>% 
  mutate(status = factor(status, levels = c("Very Unstable", "Unstable", "Less Stable", "Stable", "Uncertain")))

#Report card figure
status_key %>%
  filter(!spp %in% c("coto", "euma")) %>% 
  ggplot(aes(x = spp, y = mean)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, size = 0.75) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = status), width = 0, size = 1.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_minimal() +  # Keep only one theme
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Species",
       y = expression("Posterior Distribution" ~ hat(lambda)),
       color = "Status")+
  scale_color_manual(values = setNames(colors, levels(status_key$status)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle("Relative Occupancy Probability in 2022 vs. 2016")

## Subset just the columns we need
status_key_only <- status_key %>% 
  select(spp, status)

## Append to psi values
psi_summ_full <- psi_summ %>% 
  filter(analysis == "OR|WA|ID" & type == "Average") %>% 
  select(-c(analysis, type))

psi_status <- left_join(psi_summ_full, status_key_only, by = "spp")


#First, let's create a function to generate different shades of a base color
generate_shades <- function(base_color, n_shades) {
  colorRampPalette(c(base_color, "white"))(n_shades + 1)[1:n_shades]
}

# Get counts of species per status group to determine how many shades we need
species_per_status <- psi_status %>%
  select(spp, status) %>%
  distinct() %>%
  group_by(status) %>%
  summarise(count = n())

# Create a list to store the color palettes for each status
status_shades <- list()

# Generate shades for each status
for(s in unique(species_per_status$status)) {
  n <- species_per_status$count[species_per_status$status == s]
  status_shades[[s]] <- generate_shades(colors[s], n)
}

# Create a mapping from species to colors
species_color_map <- psi_status %>%
  select(spp, status) %>%
  distinct() %>%
  arrange(status, spp) %>%
  group_by(status) %>%
  mutate(shade_index = row_number()) %>%
  ungroup() %>%
  mutate(color = map2_chr(status, shade_index, ~status_shades[[.x]][.y]))

# Convert to a named vector for ggplot
species_colors <- setNames(species_color_map$color, species_color_map$spp)

psi_status %>%
  mutate(spp_ordered = factor(spp, 
                              levels = species_color_map %>%
                                arrange(status, spp) %>%
                                pull(spp))) %>% 
  ggplot(aes(x = year, y = mean, color = spp_ordered, fill = spp_ordered, group = spp_ordered)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0,1)) +
  facet_wrap(~ status) +
  scale_color_manual(values = species_colors,
                     # Rename legend title to "Species"
                     name = "Species",
                     # Optional: group species by status in the legend
                     breaks = species_color_map %>%
                       arrange(status, spp) %>%
                       pull(spp)) +
  scale_fill_manual(values = species_colors,
                    name = "Species",
                    breaks = species_color_map %>%
                      arrange(status, spp) %>%
                      pull(spp)) +
  theme_bw() +
  labs(
    title = "Bat Occupancy Trends by Status Group (2016-2022)",
    x = "Year",
    y = "Occupancy Probability"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
