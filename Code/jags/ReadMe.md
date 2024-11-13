# JAGS code
This folder contains all the code relevant to model-fitting and summarizing results  
  
Details about the files in this folder:
  
File | Description
---|---------------------------------------------------------------------
00_jagsPrep.R | Code required to prepare and run the data analysis for each species (WARNING: Takes ~2 hours to run on my machine).
01_jagsMod_sensitivity.R | Mirrored code from `00_jagsPrep.R`, but only runs the analysis for Oregon and Washington code. 
02_jagsMod_diagnostics.R | Results summary and model diagnostics
