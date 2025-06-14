## ---------------------------
##
## Script name: 02_ResFigs.R
##
## Purpose of script: Create summary figures for results
##
## Author: Trent VanHawkins
##
## Date Created: 2025-03-31
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

## Specify Color Palatte
state_colors <- wesanderson::wes_palette("FantasticFox1")
analysis_cols <- wesanderson::wes_palette("Darjeeling1")

possible_bats <- sort(c("laci",
                   "lano",
                   "myev",
                   "epfu",
                   "myyu",
                   "myth",
                   "myci",
                   "myvo",
                   "anpa",
                   "pahe",
                   "euma",
                   "myca",
                   "mylu",
                   "coto"))
# Read in the data --------------------------------------------------------
results <- readRDS(here("DataProcessed/results/jags/all_res.rds"))

# Parameter Figurs --------------------------------------------------------
## Alpha
results$params %>%
  filter(param %in% c("alpha01", "log_fc", "precip", "elev", "cliff_cover")) %>% 
  mutate(param = case_when(param == "log_fc" ~ "Log(Forest Cover)",
            param == "precip" ~ "Precipitation",
            param == "elev" ~ "Elevation",
            param == "cliff_cover" ~ "Log(Cliff Cover)",
            T ~ "(Intercept)")) %>% 
  ggplot(aes(x = species, y = mean, group = analysis))+
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, linewidth = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = analysis), width = 0,
                linewidth = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_grid(param ~ ., scales = 'free') +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = analysis_cols[1:2])+
  xlab('Species') + 
  ylab('Posterior Distribution (Log-Odds)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle('Occupancy Coefficients (Mean & 95% CI)')+
  labs(color = "Analysis")

ggsave(
  filename = "alphas.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 6, height = 6, dpi = 300# adjust size and resolution as needed
)

## Betas
results$params %>%
  filter(!param %in% c("alpha01", "log_fc", "precip", "elev", "cliff_cover")) %>% 
  mutate(param = case_when(param == "tmin" ~ "Min. Temp",
                           param == "clut0" ~ "Clutter 0",
                           param == "clut1" ~ "Clutter 1",
                           param == "clut2" ~ "Clutter 2",
                           param == "clut3" ~ "Clutter 3",
                           param == "dayl" ~ "Day Length",
                           param == "water" ~ "Waterbody",
                           T ~ "(Intercept)")) %>% 
  ggplot(aes(x = species, y = mean, group = analysis))+
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, linewidth = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = analysis), width = 0,
                linewidth = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_grid(param ~ ., scales = 'free') +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = analysis_cols[1:2])+
  xlab('Species') + 
  ylab('Posterior Distribution (Log-Odds)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle('Detection Coefficients (Mean & 95% CI)')+
  labs(color = "Analysis")

ggsave(
  filename = "betas.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 6, height = 7, dpi = 300# adjust size and resolution as needed
)


# Trend (Overall) -------------------------------------------------------------
trend_cols <- c("Decreasing" = '#D7191C',  # Red  
                          "Increasing" = '#1A9641', # Green
                          "Flat" = '#D3D3D3')      # Light Gray
trend_key <- results$trend %>% 
  mutate(trend = case_when(lci > 1 ~ "Increasing",
                           uci < 1 ~ "Decreasing",
                           T ~ "Flat")) %>% 
  select(species, trend) %>% 
  distinct()

## Lambda Plot (coto euma removed)
results$trend %>% 
  filter(!species %in% c("coto", "euma")) %>% 
  ggplot(aes(x = species, y = mean, group = analysis))+
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, linewidth = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = analysis), width = 0,
                linewidth = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  geom_hline(yintercept = 1, lty = 2) +
  scale_color_manual(values = analysis_cols[1:2])+
  xlab('Species') + 
  ylab(expression("Posterior Distribution (" * lambda[tot] * ")")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right") +
  labs(color = "Analysis")
  
ggsave(
  filename = "lambda.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 6, height = 3, dpi = 300# adjust size and resolution as needed
)


# Lambda Table ------------------------------------------------------------
library(readxl)

sp_names <- read_xlsx(here("DataRaw/SpeciesNatHistMatrix.xlsx"), sheet = "USBats")

sp_names <- sp_names %>% 
  select(fourlettername, SPECIES, `COMMON NAME`)
#Read in species names

lambda_restab <- results$trend %>% 
  filter(analysis == "OR|WA|ID") %>% 
  select( species, mean, lci, uci) %>% 
  mutate("CI_95" = paste0("(", round(lci, 3), ", ", round(uci,3), ")"), 
         "p_decline" = 100*(mean - 1),
         "CI_95_decline" = paste0("(", round(100*(lci - 1), 3), ", ", round(100*(uci - 1),3), ")")) %>% 
  left_join(sp_names, by = c("species" = "fourlettername"))

decline_rows <- which(lambda_restab$uci < 1)
increase_rows <- which(lambda_restab$lci > 1)

lambda_restab %>% 
  select(`COMMON NAME`, SPECIES, species, mean, CI_95, p_decline, CI_95_decline) %>% 
  kableExtra::kable(
    col.names = c("Common Name", "Species Name", "Code", "$\\hat\\lambda$", "95\\% Credible Interval", "\\% Change", "Change 95\\% CI"),
    digits = 3, format = "latex", booktabs = TRUE, escape = FALSE,
    caption = "Trend estimates reported as by the posterior mean $\\hat\\lambda$ (95\\% Credible Interval). Species with a 95% CI below one (red) are interpreted to have a decreasing trend while species with a 95% CI above one have an increasing trend.",
    label = "tbl-lambda", 
  ) %>% 
  kableExtra::row_spec(decline_rows, background = "#F4CCCC") %>%  # red
  kableExtra::row_spec(increase_rows, background = "#D9EAD3")
  

## Lambda Plot (All Species)
results$trend %>% 
  filter(species %in% c("coto", "euma")) %>% 
  ggplot(aes(x = species, y = mean, group = analysis))+
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, linewidth = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = analysis), width = 0,
                linewidth = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  geom_hline(yintercept = 1, lty = 2) +
  scale_color_manual(values = analysis_cols[1:2])+
  xlab('Species') + 
  ylab(expression("Posterior Distribution (" * lambda[tot] * ")")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle('Trend Assessment (Mean & 95% CI)',
          'Multi-season, single-species model')+
  labs(color = "Legend")

ggsave(
  filename = "lambda_eumacoto.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 8, height = 8, dpi = 300# adjust size and resolution as needed
)

## Lambda by state
results$trend_bystate %>% 
  filter(!species %in% c("coto", "euma"), analysis == "OR|WA|ID") %>% 
  ggplot(aes(x = species, y = mean, group = state))+
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, linewidth = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = state), width = 0,
                linewidth = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  geom_hline(yintercept = 1, lty = 2) +
  scale_color_manual(values = state_colors[3:5])+
  xlab('Species') + 
  ylab(expression("Posterior Distribution (" * lambda[tot] * ")")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle('Trend Assessment By State (Mean & 95% CI)',
          'Multi-season, single-species model')+
  labs(color = "Legend")

ggsave(
  filename = "lambda_bystate.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 8, height = 8, dpi = 300# adjust size and resolution as needed
)

## Psi by lambda  

#First, let's create a function to generate different shades of a base color
generate_shades <- function(base_color, n_shades) {
  colorRampPalette(c(base_color, "white"))(n_shades + 1)[1:n_shades]
}

# Get counts of species per status group to determine how many shades we need
species_per_status <- trend_key %>%
  select(species, trend) %>%
  distinct() %>%
  group_by(trend) %>%
  summarise(count = n())

# Create a list to store the color palettes for each status
status_shades <- list()

# Generate shades for each status
for(s in unique(species_per_status$trend)) {
  n <- species_per_status$count[species_per_status$trend == s]
  status_shades[[s]] <- generate_shades(trend_cols[s], n)
}

# Create a mapping from species to colors
species_color_map <- trend_key %>%
  select(species, trend) %>%
  distinct() %>%
  arrange(trend, species) %>%
  group_by(trend) %>%
  mutate(shade_index = row_number()) %>%
  ungroup() %>%
  mutate(color = map2_chr(trend, shade_index, ~status_shades[[.x]][.y]))

# Convert to a named vector for ggplot
species_colors <- setNames(species_color_map$color, species_color_map$species)

## psi by lambda (overall)
results$psi %>%
  left_join(trend_key, by = "species") %>% 
  filter(analysis == "OR|WA|ID") %>% 
  mutate(species_ordered = factor(species, 
                              levels = species_color_map %>%
                                arrange(trend, species) %>%
                                pull(species))) %>% 
  ggplot(aes(x = year, y = mean, color = species_ordered, 
             fill = species_ordered, 
             group = species_ordered)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0,1)) +
  facet_wrap( ~ trend, nrow = 3) +
  scale_color_manual(values = species_colors,
                     # Rename legend title to "Species"
                     name = "Species",
                     # Optional: group species by trend in the legend
                     breaks = species_color_map %>%
                       arrange(trend, species) %>%
                       pull(species)) +
  scale_fill_manual(values = species_colors,
                    name = "Species",
                    breaks = species_color_map %>%
                      arrange(trend, species) %>%
                      pull(species)) +
  theme_bw() +
  labs(
    title = "Bat Occupancy Trends by Status Group (2016-2022)",
    x = "Year",
    y = expression("Posterior Distribution (" * bar(psi)[t] * ")")
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "overall_trend.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 10, height = 8, dpi = 300# adjust size and resolution as needed
)  
## psi by lambda (subset)
results$psi %>%
  left_join(trend_key, by = "species") %>% 
  filter(analysis == "OR|WA|ID") %>% 
  mutate(species_ordered = factor(species, 
                                  levels = species_color_map %>%
                                    arrange(trend, species) %>%
                                    pull(species))) %>% 
  filter(species_ordered %in% c("myth", "mylu")) %>% 
  ggplot(aes(x = year, y = mean, color = species_ordered, 
             fill = species_ordered, 
             group = species_ordered)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0,1)) +
  facet_wrap( ~ trend, nrow = 3) +
  scale_color_manual(values = species_colors,
                     # Rename legend title to "Species"
                     name = "Species",
                     # Optional: group species by trend in the legend
                     breaks = species_color_map %>%
                       arrange(trend, species) %>%
                       pull(species)) +
  scale_fill_manual(values = species_colors,
                    name = "Species",
                    breaks = species_color_map %>%
                      arrange(trend, species) %>%
                      pull(species)) +
  theme_bw() +
  labs(
    x = "Year",
    y = expression("Posterior Distribution (" * hat(psi)[t] * ")")
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "overall_trend_subset.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 5, height = 3, dpi = 300# adjust size and resolution as needed
)  

## psi by lambda (undertain spp)
results$psi %>%
  left_join(trend_key, by = "species") %>% 
  filter(analysis == "OR|WA|ID") %>% 
  mutate(species_ordered = factor(species, 
                                  levels = species_color_map %>%
                                    arrange(trend, species) %>%
                                    pull(species))) %>% 
  filter(species_ordered %in% c("euma", "coto")) %>% 
  ggplot(aes(x = year, y = mean, color = species_ordered, 
             fill = species_ordered, 
             group = species_ordered)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0,1)) +
  scale_color_manual(values = state_colors[4:5]) +
  scale_fill_manual(values = state_colors[4:5]) +
  theme_bw() +
  labs(
    title = "Bat Occupancy Trends by Status Group (2016-2022)",
    x = "Year",
    y = expression("Posterior Distribution (" * bar(psi)[t] * ")"),
    color = "Species",
    fill = "Species"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "psi_trend_uncertain.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 10, height = 8, dpi = 300# adjust size and resolution as needed
)  
results$psi_bystate %>%
  left_join(trend_key, by = "species") %>% 
  filter(analysis == "OR|WA|ID") %>% 
  mutate(species_ordered = factor(species, 
                                  levels = species_color_map %>%
                                    arrange(trend, species) %>%
                                    pull(species))) %>% 
  ggplot(aes(x = year, y = mean, color = species_ordered, 
             fill = species_ordered, 
             group = species_ordered)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0,1)) +
  facet_grid(trend ~ state) +
  scale_color_manual(values = species_colors,
                     # Rename legend title to "Species"
                     name = "Species",
                     # Optional: group species by trend in the legend
                     breaks = species_color_map %>%
                       arrange(trend, species) %>%
                       pull(species)) +
  scale_fill_manual(values = species_colors,
                    name = "Species",
                    breaks = species_color_map %>%
                      arrange(trend, species) %>%
                      pull(species)) +
  theme_bw() +
  labs(
    title = "Bat Occupancy Trends by Status Group (2016-2022)",
    x = "Year",
    y = expression("Posterior Distribution (" * bar(psi)[t] * ")")
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "bystate_trend.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 10, height = 8, dpi = 300# adjust size and resolution as needed
)  

results$psi_bystate %>%
  left_join(trend_key, by = "species") %>% 
  filter(analysis == "OR|WA|ID",
         species %in% c("myth", "mylu")) %>% 
  mutate(species_ordered = factor(species, 
                                  levels = species_color_map %>%
                                    arrange(trend, species) %>%
                                    pull(species))) %>% 
  ggplot(aes(x = year, y = mean, color = species_ordered, 
             fill = species_ordered, 
             group = species_ordered)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0,1)) +
  facet_grid(trend ~ state) +
  scale_color_manual(values = species_colors,
                     # Rename legend title to "Species"
                     name = "Species",
                     # Optional: group species by trend in the legend
                     breaks = species_color_map %>%
                       arrange(trend, species) %>%
                       pull(species)) +
  scale_fill_manual(values = species_colors,
                    name = "Species",
                    breaks = species_color_map %>%
                      arrange(trend, species) %>%
                      pull(species)) +
  theme_bw() +
  labs(
    x = "Year",
    y = expression("Posterior Distribution (" * hat(psi)[t] * ")")
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "bystate_trend_subset.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 9, height = 3, dpi = 300# adjust size and resolution as needed
)  


# Table of overall growth -------------------------------------------------
results$trend %>% 
  filter(analysis == "OR|WA|ID") %>% 
  mutate(mean_change = (mean - 1) * 100,
         lci_change = (lci - 1) * 100,
         uci_change = (uci - 1) * 100) %>% 
  arrange(species)
  


# psi by lambda (remaining species) ---------------------------------------
## psi by lambda (subset)
results$psi %>%
  left_join(trend_key, by = "species") %>% 
  filter(analysis == "OR|WA|ID",
         species %in% possible_bats[1:6]) %>% 
  mutate(species_ordered = factor(species, 
                                  levels = species_color_map %>%
                                    arrange(trend, species) %>%
                                    pull(species))) %>% 
  ggplot(aes(x = year, y = mean, color = trend, 
             fill = trend, 
             group = trend)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0,1)) +
  facet_wrap(~species, ncol = 1) +
  scale_color_manual(values = trend_cols,
                     # Rename legend title to "Species"
                     name = "Trend") +
  scale_fill_manual(values = trend_cols,
                    name = "Trend") +
  theme_bw() +
  labs(
    x = "Year",
    y = expression("Posterior Distribution (" * hat(psi)[t] * ")")
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "overall_trend_allspp_part1.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 6, height = 8, dpi = 300# adjust size and resolution as needed
)  

## psi by lambda (subset)
results$psi %>%
  left_join(trend_key, by = "species") %>% 
  filter(analysis == "OR|WA|ID") %>% 
  mutate(species_ordered = factor(species, 
                                  levels = species_color_map %>%
                                    arrange(trend, species) %>%
                                    pull(species))) %>% 
  ggplot(aes(x = year, y = mean, color = trend, 
             fill = trend, 
             group = trend)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0,1)) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(values = trend_cols,
                     # Rename legend title to "Species"
                     name = "Trend") +
  scale_fill_manual(values = trend_cols,
                    name = "Trend") +
  theme_bw() +
  labs(
    x = "Year",
    y = expression("Posterior Distribution (" * hat(psi)[t] * ")")
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "overall_trend_allspp.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 9, height = 9, dpi = 300# adjust size and resolution as needed
)  


# psi_bystate_washington  -------------------------------------------------
results$psi_bystate %>%
  left_join(trend_key, by = "species") %>% 
  filter(analysis == "OR|WA|ID",
         state == "Washington") %>% 
  ggplot(aes(x = year, y = mean, color = trend, 
             fill = trend, 
             group = trend)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0,1)) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(values = trend_cols,
                     # Rename legend title to "Species"
                     name = "Trend") +
  scale_fill_manual(values = trend_cols,
                    name = "Trend") +
  theme_bw() +
  labs(
    x = "Year",
    y = expression("Posterior Distribution (" * hat(psi)[t] * ")")
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "bystate_trend_washington.png",
  path = here("Background/presentation_figs/"),
  plot = last_plot(),  # or replace with the name of your plot object
  width = 9, height = 9, dpi = 300# adjust size and resolution as needed
)  



