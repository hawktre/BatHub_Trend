# JAGS code

This folder contains all the code relevant to model-fitting and summarizing results.

Details about the files in this folder:

File | Description
---|---------------------------------------------------------------------
`00_jagsMod.R` | Prepares data arrays for JAGS and fits the dynamic occupancy model for each species (full OR/WA/ID analysis). **WARNING: Takes ~2 hours per species on a standard laptop.** Outputs per-species fit objects to `DataProcessed/results/jags/full/fits/`.
`00a_jagsMod_sensitivity.R` | Same as `00_jagsMod.R` but restricted to Oregon and Washington only. Outputs to `DataProcessed/results/jags/ORWA_Only/fits/`.
`01_ParameterSummaries.R` | Extracts posterior distributions from fit objects; computes ψ̄ₜ (mean occupancy by year), λ (trend ratio), and parameter tables. Outputs `res_summary.rds` for both analyses.
`01a_SummaryFormat.R` | Combines full and sensitivity `res_summary.rds` into a single `all_res.rds` mega-list for use by figure scripts.
`occ_model_royle_onlypsi.jags` | JAGS model file. Dynamic (multi-season) single-species occupancy model following the Royle-Kéry framework with colonization/extinction dynamics.
`archive/` | Previous JAGS diagnostic and figure scripts from the 2016–2022 analysis.
