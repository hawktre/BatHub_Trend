# Code/Archive

Archived scripts from earlier iterations of this analysis. These are **not part of the active pipeline** and are kept for reference only.

## Contents

| Folder/File | Description |
|-------------|-------------|
| `stan/` | Earlier dynamic occupancy model implemented in Stan (replaced by JAGS) |
| `spOccupancy/` | Earlier spatial occupancy model implemented with the `spOccupancy` package (replaced by JAGS) |
| `spOccupancy/laci_redo_2019/` | Single-species re-analysis for LACI using 2019 data only |
| `TrendModel.R` | Standalone trend estimation script from the original analysis approach |
| `MultiOcc_AutoLogit_Bats_VaguePriors.txt` | Early JAGS model file draft using vague priors |

## Why Archived
The Stan and spOccupancy implementations were replaced with the JAGS Royle-Kéry dynamic occupancy model (`Code/jags/`) because it offered more flexibility for multi-species occupancy tracking without spatial random effects and better fit the NABat reporting framework.
