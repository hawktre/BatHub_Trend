# BatHub_Trend Pipeline Cleanup Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Document every folder with a README.md, fix all broken file paths in R/Python scripts, and verify all data sources have corresponding retrieval code.

**Architecture:** The pipeline is a sequential R-based workflow: raw DB tables → formatted detections → DAYMET integration → grid covariates → JAGS occupancy modeling → posterior summaries → figures. All paths use the `here` package relative to the `.Rproj` root.

**Tech Stack:** R (tidyverse, here, sf, terra, jagsUI), Python (pydaymet), JAGS

---

## Broken Paths & Code Issues Found

### 1. `Code/formatting/00b_DetectionDataFormat.R` (untracked, active version)
- **Line 116–118**: References undefined variable `test$replicate` — leftover debug code
- **Line 137**: `left_join(daymet_clean, by = c("LocationName" = "site", "Night" = "date"))` — wrong column names; after `clean_names()` in 00a, columns are `location_name` and `night`

### 2. `Code/jags/00_jagsMod.R`
- **Line 72**: `drop_na(DEM_max, p_forest, precip, cliff_canyon)` — column is `cliff_cover`, not `cliff_canyon`
- **Line 99**: `readRDS(here("DataProcessed/occurrence/nw_grid_shp.rds"))` — file is `nw_grid_shp_to2024.rds`
- **Lines 161–163**: Hardcoded year selection stops at `'2022'`; 2023 and 2024 are missing

### 3. `Code/formatting/00b_BatGrid_Covariates.R`
- **Lines 110–114**: Leftover debug code referencing `bat.dat` (undefined) after the final `write_sf()` call

### 4. `Code/formatting/get_daymet.py`
- **Line 6**: `pd.read_csv(,` — missing filename argument (should be `daymet_batch.csv` path)

### 5. `Code/formatting/Untitled.R`
- Stale draft of `00b_DetectionDataFormat.R` that reads from `all_detections_join.csv` (no pipeline step produces this file). Move to `Code/Archive/`.

## Data Source Orphan Check
- `DataRaw/covariates/DEM/` — folder exists but `DEM_max` comes from `NABat_grid_covariates` shapefile; verify DEM folder contents aren't referenced anywhere

## Missing README Files
Code/formatting/, Code/Archive/, Code/Figures_Tables/,
DataRaw/tables/, DataRaw/covariates/, DataRaw/batgrid/,
DataProcessed/detections/, DataProcessed/occurrence/, DataProcessed/results/

---

### Task 1: Write Plan (this document)
**Done.**

---

### Task 2: Add READMEs to Code subdirectories

**Files:**
- Create: `Code/formatting/README.md`
- Create: `Code/Archive/README.md`
- Create: `Code/Figures_Tables/README.md`
- Update: `Code/ReadMe.md`
- Update: `Code/jags/ReadMe.md`

---

### Task 3: Add READMEs to DataRaw subdirectories

**Files:**
- Create: `DataRaw/tables/README.md`
- Create: `DataRaw/covariates/README.md`
- Create: `DataRaw/batgrid/README.md`

---

### Task 4: Add READMEs to DataProcessed subdirectories

**Files:**
- Create: `DataProcessed/detections/README.md`
- Create: `DataProcessed/occurrence/README.md`
- Create: `DataProcessed/results/README.md`

---

### Task 5: Fix `Code/formatting/00b_DetectionDataFormat.R`

Remove lines 116–118 (debug code).
Fix line 137 join column names: `location_name`/`night`.

---

### Task 6: Fix `Code/jags/00_jagsMod.R`

- Line 72: `cliff_canyon` → `cliff_cover`
- Line 99: `nw_grid_shp.rds` → `nw_grid_shp_to2024.rds`
- Lines 161–163: extend year columns through `'2024'`

---

### Task 7: Clean `Code/formatting/00b_BatGrid_Covariates.R`

Remove lines 110–114 (dead debug code after write_sf).

---

### Task 8: Fix `Code/formatting/get_daymet.py`

Add correct filename to `pd.read_csv()` call.

---

### Task 9: Archive `Code/formatting/Untitled.R`

Move to `Code/Archive/Untitled_00b_draft.R` with note at top.
