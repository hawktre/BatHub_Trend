This folder contains all the code.

Details about the subfolders in this folder:

Folder | Description
---|---------------------------------------------------------------------
 `formatting/` | Data extraction, cleaning, and formatting scripts. Run these first. See `formatting/README.md` for execution order.
 `jags/` | JAGS occupancy model fitting and posterior summary scripts. See `jags/ReadMe.md`.
 `Figures_Tables/` | Figure and table generation scripts. Run after the full JAGS pipeline.
 `Archive/` | Archived code from earlier iterations of this analysis (spOccupancy and Stan).

## Supporting Files

File | Description
---|---------------------------------------------------------------------
 `download_daymet.R` | Downloads DAYMET daily climate data via `daymetr::download_daymet_batch()`. Run after `formatting/00b_DetectionDataFormat.R` writes the batch CSV.
