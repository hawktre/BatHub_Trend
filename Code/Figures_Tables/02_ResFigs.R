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
## ---------------------------

options(scipen = 6, digits = 4)

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(readxl)

# Color Palettes ----------------------------------------------------------
state_colors    <- wesanderson::wes_palette("FantasticFox1")
analysis_cols   <- wesanderson::wes_palette("Darjeeling1")

trend_cols <- c("Decreasing" = '#D7191C',
                "Increasing" = '#1A9641',
                "Flat"       = '#D3D3D3')

possible_bats <- sort(c("laci", "lano", "myev", "epfu", "myyu", "myth",
                        "myci", "myvo", "anpa", "pahe", "euma", "myca",
                        "mylu", "coto"))

# Read in Data ------------------------------------------------------------
results <- readRDS(here("DataProcessed/results/jags/all_res.rds"))

# Dynamic Year Range ------------------------------------------------------
## Derived from data so plots update automatically when years are added
year_range <- results$psi %>%
  summarise(min_year = min(year), max_year = max(year))

year_label <- paste0(year_range$min_year, "\u2013", year_range$max_year)

# Trend Key ---------------------------------------------------------------
trend_key <- results$trend %>%
  mutate(trend = case_when(lci > 1 ~ "Increasing",
                           uci < 1 ~ "Decreasing",
                           TRUE    ~ "Flat")) %>%
  select(species, trend) %>%
  distinct()

# Species Color Map -------------------------------------------------------
species_per_status <- trend_key %>%
  group_by(trend) %>%
  summarise(count = n())

generate_shades <- function(base_color, n_shades) {
  colorRampPalette(c(base_color, "white"))(n_shades + 1)[1:n_shades]
}

status_shades <- list()
for (s in unique(species_per_status$trend)) {
  n <- species_per_status$count[species_per_status$trend == s]
  status_shades[[s]] <- generate_shades(trend_cols[s], n)
}

species_color_map <- trend_key %>%
  arrange(trend, species) %>%
  group_by(trend) %>%
  mutate(shade_index = row_number()) %>%
  ungroup() %>%
  mutate(color = map2_chr(trend, shade_index, ~ status_shades[[.x]][.y]))

species_colors <- setNames(species_color_map$color, species_color_map$species)

# Parameter Figures -------------------------------------------------------

## Occupancy coefficients (alphas)
results$params %>%
  filter(param %in% c("alpha01", "log_fc", "precip", "dem_max", "log_cliff")) %>%
  mutate(param = case_when(param == "log_fc"      ~ "Log(Forest Cover)",
                           param == "precip"       ~ "Precipitation",
                           param == "dem_max"         ~ "Elevation",
                           param == "log_cliff"  ~ "Log(Cliff Cover)",
                           TRUE                    ~ "(Intercept)")) %>%
  ggplot(aes(x = species, y = mean, group = analysis)) +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, linewidth = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = analysis), width = 0,
                linewidth = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(param ~ ., scales = "free") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = analysis_cols[1:2]) +
  labs(title = "Occupancy Coefficients (Mean & 95% CI)",
       x = "Species", y = "Posterior Distribution (Log-Odds)", color = "Analysis") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom")

ggsave("alphas.png", path = here("Background/presentation_figs/"),
       width = 6, height = 6, dpi = 300)

## Detection coefficients (betas)
## Detection coefficients (betas)
results$params %>%
  filter(param %in% c("(Intercept)", "clutter_percent1", "clutter_percent2",
                       "clutter_percent3", "clutter_percent4",
                       "scale(tmin)", "scale(dayl)", "water_ind")) %>%
  mutate(param = case_when(param == "scale(tmin)"       ~ "Min. Temp",
                           param == "clutter_percent1"  ~ "Clutter 1",
                           param == "clutter_percent2"  ~ "Clutter 2",
                           param == "clutter_percent3"  ~ "Clutter 3",
                           param == "clutter_percent4"  ~ "Clutter 4",
                           param == "scale(dayl)"       ~ "Day Length",
                           param == "water_ind"         ~ "Waterbody",
                           TRUE                         ~ "(Intercept)")) %>%
  ggplot(aes(x = species, y = mean, group = analysis)) +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, linewidth = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = analysis), width = 0,
                linewidth = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(param ~ ., scales = "free") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = analysis_cols[1:2]) +
  labs(title = "Detection Coefficients (Mean & 95% CI)",
       x = "Species", y = "Posterior Distribution (Log-Odds)", color = "Analysis") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom")

ggsave("betas.png", path = here("Background/presentation_figs/"),
       width = 6, height = 7, dpi = 300)

# Trend Figures (Overall) -------------------------------------------------

## Lambda plot (coto and euma excluded)
results$trend %>%
  filter(!species %in% c("coto", "euma")) %>%
  ggplot(aes(x = species, y = mean, group = analysis)) +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, linewidth = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = analysis), width = 0,
                linewidth = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, lty = 2) +
  scale_color_manual(values = analysis_cols[1:2]) +
  labs(x = "Species",
       y = expression("Posterior Distribution (" * lambda[tot] * ")"),
       color = "Analysis") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right")

ggsave("lambda.png", path = here("Background/presentation_figs/"),
       width = 6, height = 3, dpi = 300)

## Lambda plot (coto and euma only)
results$trend %>%
  filter(species %in% c("coto", "euma")) %>%
  ggplot(aes(x = species, y = mean, group = analysis)) +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, linewidth = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = analysis), width = 0,
                linewidth = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, lty = 2) +
  scale_color_manual(values = analysis_cols[1:2]) +
  labs(title = glue::glue("Trend Assessment (Mean & 95% CI)"),
       subtitle = "Multi-season, single-species model",
       x = "Species",
       y = expression("Posterior Distribution (" * lambda[tot] * ")"),
       color = "Legend") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom")

ggsave("lambda_eumacoto.png", path = here("Background/presentation_figs/"),
       width = 8, height = 8, dpi = 300)

## Lambda by state
results$trend_bystate %>%
  filter(!species %in% c("coto", "euma"), analysis == "OR|WA|ID") %>%
  ggplot(aes(x = species, y = mean, group = state)) +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, linewidth = 0.75,
                position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q25, ymax = q75, color = state), width = 0,
                linewidth = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, lty = 2) +
  scale_color_manual(values = state_colors[3:5]) +
  labs(title = paste("Trend Assessment By State (Mean & 95% CI)", year_label, sep = " \u2013 "),
       subtitle = "Multi-season, single-species model",
       x = "Species",
       y = expression("Posterior Distribution (" * lambda[tot] * ")"),
       color = "Legend") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom")

ggsave("lambda_bystate.png", path = here("Background/presentation_figs/"),
       width = 8, height = 8, dpi = 300)

# Lambda Table ------------------------------------------------------------
sp_names <- read_xlsx(here("DataRaw/SpeciesNatHistMatrix.xlsx"), sheet = "USBats") %>%
  select(fourlettername, SPECIES, `COMMON NAME`)

lambda_restab <- results$trend %>%
  filter(analysis == "OR|WA|ID") %>%
  select(species, mean, lci, uci) %>%
  mutate(CI_95        = paste0("(", round(lci, 3), ", ", round(uci, 3), ")"),
         p_decline    = 100 * (mean - 1),
         CI_95_decline = paste0("(", round(100 * (lci - 1), 3), ", ",
                                round(100 * (uci - 1), 3), ")")) %>%
  left_join(sp_names, by = c("species" = "fourlettername"))

decline_rows  <- which(lambda_restab$uci < 1)
increase_rows <- which(lambda_restab$lci > 1)

lambda_restab %>%
  select(`COMMON NAME`, SPECIES, species, mean, CI_95, p_decline, CI_95_decline) %>%
  kableExtra::kable(
    col.names = c("Common Name", "Species Name", "Code",
                  "$\\hat\\lambda$", "95\\% Credible Interval",
                  "\\% Change", "Change 95\\% CI"),
    digits = 3, format = "latex", booktabs = TRUE, escape = FALSE,
    caption = paste0("Trend estimates (", year_label, ") reported as posterior mean $\\hat\\lambda$ ",
                     "(95\\% Credible Interval). Species with a 95\\% CI below one (red) are ",
                     "interpreted as declining; above one (green) as increasing."),
    label = "tbl-lambda"
  ) %>%
  kableExtra::row_spec(decline_rows,  background = "#F4CCCC") %>%
  kableExtra::row_spec(increase_rows, background = "#D9EAD3")

# Psi Trend Figures -------------------------------------------------------

## Helper to build the ordered species factor (used repeatedly)
ordered_species <- function(df) {
  df %>%
    mutate(species_ordered = factor(species,
                                    levels = species_color_map %>%
                                      arrange(trend, species) %>%
                                      pull(species)))
}

## Overall psi by trend group (all species)
results$psi %>%
  left_join(trend_key, by = "species") %>%
  filter(analysis == "OR|WA|ID") %>%
  ordered_species() %>%
  ggplot(aes(x = year, y = mean, color = species_ordered,
             fill = species_ordered, group = species_ordered)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0, 1)) +
  facet_wrap(~trend, nrow = 3) +
  scale_color_manual(values = species_colors, name = "Species",
                     breaks = species_color_map %>% arrange(trend, species) %>% pull(species)) +
  scale_fill_manual(values = species_colors, name = "Species",
                    breaks = species_color_map %>% arrange(trend, species) %>% pull(species)) +
  labs(title = paste("Bat Occupancy Trends by Status Group", year_label),
       x = "Year",
       y = expression("Posterior Distribution (" * bar(psi)[t] * ")")) +
  theme_bw() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave("overall_trend.png", path = here("Background/presentation_figs/"),
       width = 10, height = 8, dpi = 300)

## Overall psi subset (myth, mylu)
results$psi %>%
  left_join(trend_key, by = "species") %>%
  filter(analysis == "OR|WA|ID", species %in% c("myth", "mylu")) %>%
  ordered_species() %>%
  ggplot(aes(x = year, y = mean, color = species_ordered,
             fill = species_ordered, group = species_ordered)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0, 1)) +
  facet_wrap(~trend, nrow = 3) +
  scale_color_manual(values = species_colors, name = "Species",
                     breaks = species_color_map %>% arrange(trend, species) %>% pull(species)) +
  scale_fill_manual(values = species_colors, name = "Species",
                    breaks = species_color_map %>% arrange(trend, species) %>% pull(species)) +
  labs(x = "Year",
       y = expression("Posterior Distribution (" * hat(psi)[t] * ")")) +
  theme_bw() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave("overall_trend_subset.png", path = here("Background/presentation_figs/"),
       width = 5, height = 3, dpi = 300)

## Overall psi (uncertain species: euma, coto)
results$psi %>%
  left_join(trend_key, by = "species") %>%
  filter(analysis == "OR|WA|ID", species %in% c("euma", "coto")) %>%
  ordered_species() %>%
  ggplot(aes(x = year, y = mean, color = species_ordered,
             fill = species_ordered, group = species_ordered)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0, 1)) +
  scale_color_manual(values = state_colors[4:5]) +
  scale_fill_manual(values = state_colors[4:5]) +
  labs(title = paste("Bat Occupancy Trends by Status Group", year_label),
       x = "Year",
       y = expression("Posterior Distribution (" * bar(psi)[t] * ")"),
       color = "Species", fill = "Species") +
  theme_bw() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave("psi_trend_uncertain.png", path = here("Background/presentation_figs/"),
       width = 10, height = 8, dpi = 300)

## All species faceted by trend (all species, ncol = 3)
results$psi %>%
  left_join(trend_key, by = "species") %>%
  filter(analysis == "OR|WA|ID") %>%
  ordered_species() %>%
  ggplot(aes(x = year, y = mean, color = trend, fill = trend, group = trend)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0, 1)) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(values = trend_cols, name = "Trend") +
  scale_fill_manual(values = trend_cols, name = "Trend") +
  labs(x = "Year",
       y = expression("Posterior Distribution (" * hat(psi)[t] * ")")) +
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("overall_trend_allspp.png", path = here("Background/presentation_figs/"),
       width = 9, height = 9, dpi = 300)

## First 6 species only (ncol = 1)
results$psi %>%
  left_join(trend_key, by = "species") %>%
  filter(analysis == "OR|WA|ID", species %in% possible_bats[1:6]) %>%
  ordered_species() %>%
  ggplot(aes(x = year, y = mean, color = trend, fill = trend, group = trend)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0, 1)) +
  facet_wrap(~species, ncol = 1) +
  scale_color_manual(values = trend_cols, name = "Trend") +
  scale_fill_manual(values = trend_cols, name = "Trend") +
  labs(x = "Year",
       y = expression("Posterior Distribution (" * hat(psi)[t] * ")")) +
  theme_bw() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave("overall_trend_allspp_part1.png", path = here("Background/presentation_figs/"),
       width = 6, height = 8, dpi = 300)

# By-State Psi Figures ----------------------------------------------------

## All species by state and trend
results$psi_bystate %>%
  left_join(trend_key, by = "species") %>%
  filter(analysis == "OR|WA|ID") %>%
  ordered_species() %>%
  ggplot(aes(x = year, y = mean, color = species_ordered,
             fill = species_ordered, group = species_ordered)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0, 1)) +
  facet_grid(trend ~ state) +
  scale_color_manual(values = species_colors, name = "Species",
                     breaks = species_color_map %>% arrange(trend, species) %>% pull(species)) +
  scale_fill_manual(values = species_colors, name = "Species",
                    breaks = species_color_map %>% arrange(trend, species) %>% pull(species)) +
  labs(title = paste("Bat Occupancy Trends by Status Group", year_label),
       x = "Year",
       y = expression("Posterior Distribution (" * bar(psi)[t] * ")")) +
  theme_bw() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave("bystate_trend.png", path = here("Background/presentation_figs/"),
       width = 10, height = 8, dpi = 300)

## myth and mylu by state
results$psi_bystate %>%
  left_join(trend_key, by = "species") %>%
  filter(analysis == "OR|WA|ID", species %in% c("myth", "mylu")) %>%
  ordered_species() %>%
  ggplot(aes(x = year, y = mean, color = species_ordered,
             fill = species_ordered, group = species_ordered)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0, 1)) +
  facet_grid(trend ~ state) +
  scale_color_manual(values = species_colors, name = "Species",
                     breaks = species_color_map %>% arrange(trend, species) %>% pull(species)) +
  scale_fill_manual(values = species_colors, name = "Species",
                    breaks = species_color_map %>% arrange(trend, species) %>% pull(species)) +
  labs(x = "Year",
       y = expression("Posterior Distribution (" * hat(psi)[t] * ")")) +
  theme_bw() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave("bystate_trend_subset.png", path = here("Background/presentation_figs/"),
       width = 9, height = 3, dpi = 300)

## Washington only
results$psi_bystate %>%
  left_join(trend_key, by = "species") %>%
  filter(analysis == "OR|WA|ID", state == "Washington") %>%
  ggplot(aes(x = year, y = mean, color = trend, fill = trend, group = trend)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0, 1)) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(values = trend_cols, name = "Trend") +
  scale_fill_manual(values = trend_cols, name = "Trend") +
  labs(x = "Year",
       y = expression("Posterior Distribution (" * hat(psi)[t] * ")")) +
  theme_bw() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave("bystate_trend_washington.png", path = here("Background/presentation_figs/"),
       width = 9, height = 9, dpi = 300)

# Table of Overall Growth -------------------------------------------------
results$trend %>%
  filter(analysis == "OR|WA|ID") %>%
  mutate(mean_change = (mean - 1) * 100,
         lci_change  = (lci  - 1) * 100,
         uci_change  = (uci  - 1) * 100) %>%
  arrange(species)