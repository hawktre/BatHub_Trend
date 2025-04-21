## ---------------------------
##
## Script name: general_figs_and_tables.R
##
## Purpose of script: Tables and figures summarizing NABat efforts
##
## Author: Trent VanHawkins
##
## Date Created: 2025-03-28
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
library(kableExtra)

# Read in the data --------------------------------------------------------

nw_grid <- readRDS(here("DataProcessed/occurrence/nw_grid_shp.rds"))

## Pivot Longer for Plotting and tables
nw_grid_long <- nw_grid %>% 
  pivot_longer(cols = samp_2016:samp_2022, names_sep = "_", 
               names_to = c("samp", "year"),
               values_to = "sampled") %>% 
  mutate(sampled_txt = if_else(sampled == 1, "Yes", "No"))

## Create a table
nw_grid_long %>% 
  group_by(year) %>% 
  reframe(sites = sum(sampled)) %>% 
  kableExtra::kable(col.names = c("Year", "Sites Sampled")) %>% 
  kable_styling()

## Plot Version (All Years All States)
nw_grid_long %>% 
  ggplot() +
  geom_sf(
    aes(color = sampled_txt, linewidth = sampled_txt),
    fill = NA
  ) +
  scale_color_manual(
    values = c("Yes" = "dodgerblue", "No" = "black")
  ) +
  scale_linewidth_manual(
    values = c("Yes" = 0.5, "No" = 0.05)
  ) +
  facet_wrap(~year)+
  theme_minimal()+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank())+
  labs(color = "Sampled")+
  guides(linewidth = "none")  # Optional: hides legend for size if redundant

ggsave(
  filename = "sampled_grids.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 10, height = 6, dpi = 300,  # adjust size and resolution as needed
  bg = "transparent"  # makes the background transparent
)

state_cols <- wesanderson::wes_palette("FantasticFox1")

## Plot Version (All Years OR/WA)
nw_grid_long %>% 
  filter(state %in% c("Oregon", "Washington")) %>% 
  ggplot() +
  geom_sf(
    aes(fill = state)
  ) +
  labs(fill = "State")+
  scale_fill_manual(values = c("Oregon" = state_cols[3], "Washington" = state_cols[4]))+
  theme_minimal()+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")

ggsave(
  filename = "OR_WA_grids.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 10, height = 6, dpi = 300  # adjust size and resolution as needed
)

## Plot Version (All Years OR/WA/ID)
nw_grid_long %>% 
  ggplot() +
  geom_sf(
    aes(fill = state)
  ) +
  labs(fill = "State")+
  scale_fill_manual(values = c("Oregon" = state_cols[3], "Washington" = state_cols[4], "Idaho" = state_cols[5]))+
  theme_minimal()+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")

ggsave(
  filename = "OR_WA_ID_grids.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 10, height = 6, dpi = 300  # adjust size and resolution as needed
)
