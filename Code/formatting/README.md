# Code/formatting

Data extraction, cleaning, and formatting scripts. Run these in order before any modeling.

## Execution Order

| Step | File | Input | Output | Description |
|------|------|-------|--------|-------------|
| 1 | `00a_DetectionDataFormat.R` | `DataRaw/tables/tbl*.csv`, `tlu*.csv` | `DataProcessed/detections/detections_to2024.rds` | Joins database tables; standardizes clutter percent and water body indicator; drops records missing those covariates |
| 2 | `00b_BatGrid_Covariates.R` | `DataRaw/covariates/NABat_grid_covariates/`, `DataRaw/batgrid/`, `DataRaw/covariates/LandFire/` | `DataProcessed/occurrence/batgrid_covars.shp` | Collects grid-level landscape covariates (forest, elevation, cliff/canyon) for all PNW NABat 10 km cells |
| 3 | `download_daymet.R` | `DataRaw/covariates/daymet/daymet_batch.csv` | `DataRaw/covariates/daymet/all_daymet.rds` | Downloads DAYMET daily weather (tmin, prcp, dayl, vp) for each site-night using the `daymetr` package. Run this **after** step 4 writes the batch file, or re-run 4 once the batch file exists. |
| 4 | `00b_DetectionDataFormat.R` | `detections_to2024.rds`, `DataRaw/tables/calls_*.csv`, `DataRaw/covariates/daymet/all_daymet.rds` | `DataProcessed/detections/detections_formatted_2016-2024.rds` | Joins raw acoustic detections to deployment metadata; standardizes species codes; pivots to wide (presence/absence per species); joins DAYMET weather data |
| 5 | `01a_detections_modprep.R` | `detections_formatted_2016-2024.rds` | `DataProcessed/detections/nw_nights_to2024.rds` | Keeps first-night detection per site-year; drops rows missing covariates; recodes clutter as ordered factor |
| 6 | `01b_occurrence_modprep.R` | `batgrid_covars.shp`, `nw_nights_to2024.rds` | `DataProcessed/occurrence/nw_grid_shp_to2024.rds` | Builds sampling history by year for each 10 km grid cell and merges with landscape covariates |

## Other Files

| File | Description |
|------|-------------|
| `get_daymet.py` | Alternative Python script to pull DAYMET data via `pydaymet`. Use if `daymetr` is unavailable. Reads `daymet_batch.csv`, outputs `all_daymet.csv`. |

## Notes
- All paths use the `here` package relative to `BatHub_Trend.Rproj`.
- Steps 3 and 4 are co-dependent: step 4 writes `daymet_batch.csv`, step 3 reads it. If starting fresh, run through step 4 up to the `write.csv()` call, then run step 3, then complete step 4.
- DAYMET coverage is CONUS only; all sites are within OR/WA/ID so this is not an issue.
