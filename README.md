# multiciliation_cycle

Scripts for manuscript titled "A cell cycle variant coordinates multiciliated cell differentiation".

Raw and processed data is deposited on GEO at GSE228110.

Individual datasets were analyzed with scripts at github.com/lb15/autoSeurat. Folders for each individual dataset contain parameters.csv file for autoSeurat and scripts for doublet identification. 

Individual seurat objects were integrated with seurat_integration.R.

Integrated dataset analyses are in mcc_timecourse_analysis.R and E2f7_analysis.R scripts.
