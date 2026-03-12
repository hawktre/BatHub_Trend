# Code/Figures_Tables

Scripts for generating all publication-quality figures and summary tables from model results. Run these after the full JAGS pipeline (`Code/jags/`) is complete.

## Files

| File | Input | Output | Description |
|------|-------|--------|-------------|
| `02_ResFigs.R` | `DataProcessed/results/jags/all_res.rds` | `Reports/` figures (`alphas.png`, `betas.png`, `lambda.png`, `overall_trend.png`, etc.) | Generates occupancy trend plots (ψ̄ₜ by year and species), trend ratio (λ) figures, and occupancy/detection coefficient plots. Covers full OR/WA/ID analysis and OR/WA sensitivity analysis. |
| `general_figs_and_tables.R` | `DataProcessed/occurrence/nw_grid_shp_to2024.rds` | `Reports/sampled_grids.png`, `OR_WA_grids.png`, etc. | Sampling effort maps showing which 10 km NABat grid cells were surveyed each year. |
| `table1.R` | `DataProcessed/occurrence/nw_grid_shp_to2024.rds` | LaTeX/HTML table | Table 1 — survey effort by state and year (grid cells surveyed). |

## Color Palettes Used
- States: `wesanderson::wes_palette("FantasticFox1")`
- Analysis comparison: `wesanderson::wes_palette("Darjeeling1")`
- Trend classification: Red (Decreasing λ < 1), Green (Increasing λ > 1), Gray (Flat)

## Trend Classification
- **Decreasing**: 95% CI upper bound < 1.0
- **Increasing**: 95% CI lower bound > 1.0
- **Flat**: CI spans 1.0
