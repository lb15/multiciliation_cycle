# multiciliation_cycle

Scripts for manuscript titled "An alternative cell cycle coordinates
multiciliated cell differentiation".

Raw and processed data is deposited on GEO at GSE228110.

R objects referenced in the integrated analysis scripts are posted on Zenodo at DOI: 10.5281/zenodo.10892464.

Individual datasets were analyzed with scripts at github.com/lb15/autoSeurat. Folders for each individual dataset contain parameters.csv file for autoSeurat and scripts for doublet identification. 

Individual seurat objects were integrated with seurat_integration.R.

Integrated dataset analyses are in mcc_timecourse_analysis.R, E2f7_analysis.R, and Ribo_v1_analysis.R scripts.

CUT&RUN datasets analyzed with scripts at github.com/lb15/autoCutandRun. Figures produced with script CUTRUN.R.
