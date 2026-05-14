## ---------------------------
##
## Script name: 09_spOcc_PPC.R
##
## Purpose of script: Posterior predictive checks for tPGOcc and stPGOcc
##                    model fits using spOccupancy's built-in ppcOcc()
##                    function. Bayesian p-values near 0.5 indicate good
##                    fit; values near 0 or 1 indicate misfit.
##                    Run after 05_stPGOcc.R and 06_tPGOcc.R.
##
## Author: Trent VanHawkins
##
## Date Created: 2026-04-05
##
## ---------------------------

options(scipen = 6, digits = 4)

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(spOccupancy)

# Load Data ---------------------------------------------------------------
index_keys <- readRDS(here("DataProcessed/results/jags/index_keys.rds"))
year_ids   <- index_keys$year_ids

# Functions ---------------------------------------------------------------

## Run PPC for all fits in a directory and return combined summary
run_ppcs <- function(fit_dir, out_dir, model_label) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  paths <- list.files(fit_dir, full.names = TRUE)
  spps  <- str_split_i(tools::file_path_sans_ext(basename(paths)), "_", 1)
  
  ppc_results <- list()
  
  for (i in seq_along(paths)) {
    spp <- spps[i]
    cat("\nRunning PPC for:", spp, "(", model_label, ")\n")
    
    fit <- readRDS(paths[i])
    
    ## Freeman-Tukey statistic grouped across sites (group = 1)
    ## fit.y and fit.y.rep are matrices: dim(n_samples, n_years)
    ## for both tPGOcc and stPGOcc
    ppc <- ppcOcc(fit, fit.stat = "freeman-tukey", group = 1)
    
    ## Save raw ppc object
    saveRDS(ppc, file.path(out_dir, paste0(spp, "_ppc.rds")))
    
    n_years <- ncol(ppc$fit.y)
    
    ## Overall p-value — collapse across all time periods
    overall_pval <- mean(ppc$fit.y.rep > ppc$fit.y)
    
    ## Per time period p-values — compare column by column
    year_pvals <- sapply(seq_len(n_years), function(t) {
      mean(ppc$fit.y.rep[, t] > ppc$fit.y[, t])
    })
    
    ppc_results[[spp]] <- tibble(
      species     = spp,
      model       = model_label,
      time_period = c("Overall", as.character(year_ids)),
      bayes_pval  = c(overall_pval, year_pvals)
    )
    
    cat("  Overall Bayesian p-value:", round(overall_pval, 3), "\n")
  }
  
  bind_rows(ppc_results)
}

# Run PPCs ----------------------------------------------------------------
tpg_ppc  <- run_ppcs(
  fit_dir     = here("DataProcessed/results/tPGOcc/fits/"),
  out_dir     = here("DataProcessed/results/tPGOcc/ppc/"),
  model_label = "tPGOcc"
)

stpg_ppc <- run_ppcs(
  fit_dir     = here("DataProcessed/results/stPGOcc/fits/"),
  out_dir     = here("DataProcessed/results/stPGOcc/ppc/"),
  model_label = "stPGOcc"
)

# Combine and Save --------------------------------------------------------
ppc_summary <- bind_rows(tpg_ppc, stpg_ppc) %>%
  mutate(
    fit = case_when(
      bayes_pval < 0.1 | bayes_pval > 0.9 ~ "Poor",
      bayes_pval < 0.2 | bayes_pval > 0.8 ~ "Moderate",
      TRUE                                  ~ "Good"
    )
  )

saveRDS(ppc_summary,
        here("DataProcessed/results/ppc_summary.rds"))

# Print Summary -----------------------------------------------------------
cat("\n--- Overall Bayesian P-values by Species and Model ---\n")
ppc_summary %>%
  filter(time_period == "Overall") %>%
  select(species, model, bayes_pval, fit) %>%
  pivot_wider(names_from  = model,
              values_from = c(bayes_pval, fit)) %>%
  arrange(species) %>%
  print(n = Inf)

cat("\nSpecies with poor overall fit (p < 0.1 or p > 0.9):\n")
ppc_summary %>%
  filter(time_period == "Overall", fit == "Poor") %>%
  arrange(model, species) %>%
  print(n = Inf)

cat("\nTime periods with poor fit (p < 0.1 or p > 0.9):\n")
ppc_summary %>%
  filter(time_period != "Overall", fit == "Poor") %>%
  arrange(model, species, time_period) %>%
  print(n = Inf)

# Visualize ---------------------------------------------------------------

## PPC by model — faceted by species, side by side for tPGOcc vs stPGOcc
ppc_summary %>%
  mutate(
    time_period = factor(time_period,
                         levels = c("Overall", as.character(year_ids))),
    model = factor(model, levels = c("tPGOcc", "stPGOcc"))
  ) %>%
  ggplot(aes(x = time_period, y = bayes_pval,
             group = interaction(species, model),
             color = model,
             linetype = model)) +
  geom_hline(yintercept = c(0.1, 0.9), lty = 2, color = "gray50") +
  geom_hline(yintercept = 0.5, lty = 1, color = "gray80") +
  geom_line(alpha = 0.6) +
  geom_point(aes(shape = fit), size = 2.5) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(values = c("tPGOcc"  = "#E69F00",
                                "stPGOcc" = "#0072B2")) +
  scale_linetype_manual(values = c("tPGOcc"  = "dashed",
                                   "stPGOcc" = "solid")) +
  scale_shape_manual(values = c("Good"     = 16,
                                "Moderate" = 17,
                                "Poor"     = 4)) +
  lims(y = c(0, 1)) +
  labs(title    = "Posterior Predictive Check — Bayesian P-values",
       subtitle = "Dashed lines at 0.1 and 0.9 — values outside indicate misfit",
       x        = "Time Period",
       y        = "Bayesian P-value",
       color    = "Model",
       linetype = "Model",
       shape    = "Fit") +
  theme_bw() +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave("spOcc_ppc.png",
       path  = here("Background/presentation_figs/"),
       width = 10, height = 10, dpi = 300)

## Overall p-value comparison: tPGOcc vs stPGOcc scatter plot
ppc_summary %>%
  filter(time_period == "Overall") %>%
  select(species, model, bayes_pval) %>%
  pivot_wider(names_from = model, values_from = bayes_pval) %>%
  ggplot(aes(x = tPGOcc, y = stPGOcc, label = species)) +
  geom_abline(intercept = 0, slope = 1, lty = 2, color = "gray50") +
  geom_vline(xintercept = c(0.1, 0.9), lty = 3, color = "gray70") +
  geom_hline(yintercept = c(0.1, 0.9), lty = 3, color = "gray70") +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3) +
  lims(x = c(0, 1), y = c(0, 1)) +
  labs(title    = "Overall PPC: tPGOcc vs stPGOcc",
       subtitle = "Points above diagonal = stPGOcc closer to 0.5 (better fit)",
       x        = "tPGOcc Bayesian P-value",
       y        = "stPGOcc Bayesian P-value") +
  theme_bw()

ggsave("spOcc_ppc_scatter.png",
       path  = here("Background/presentation_figs/"),
       width = 6, height = 6, dpi = 300)
