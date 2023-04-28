## check doublets
library(Seurat)
library(dplyr)

### remove clusters with > 15% doublets

####################### WT A ###################################

wt_a=readRDS("WT_A/v3/basic_analysis/WT_A_v3_basic.rds")
wt_a_dubs=read.csv("WT_A/v3/WT_A_scrublet_08perc_thresh0.30_output_table.csv")
identical(rownames(wt_a@meta.data),wt_a_dubs$cell_barcodes)
wt_a$predicted_doublet<- wt_a_dubs$predicted_doublet
wt_a$doublet_score<-wt_a_dubs$doublet_score

DimPlot(wt_a, group.by="predicted_doublet")
DimPlot(wt_a,label=T)

FeaturePlot(wt_a, "doublet_score",min.cutoff = 0.25)

prop.table(table(wt_a$seurat_clusters,wt_a$predicted_doublet),1)
DimPlot(wt_a,label=T)

clus14_vs_clus4=FindMarkers(wt_a,assay="RNA",ident.1 = 14,ident.2 = 4)
View(clus14_vs_clus4)

## remove doublets and cluster 12 and 16

wt_a_dubs$cluster = wt_a$seurat_clusters
wt_a_dubs$predicted_doublet[wt_a_dubs$cluster == 12] <-"True"
wt_a_dubs$predicted_doublet[wt_a_dubs$cluster == 16] <-"True"

prop.table(table(wt_a_dubs$cluster,wt_a_dubs$predicted_doublet),1)

write.csv(wt_a_dubs, file="WT_A/v3/WT_A_v3_predicted_dubs_v2_clus12_16.csv")

######################### WT B ###################################
wt_b=readRDS("WT_B/v3/basic_analysis/WT_B_v3_basic.rds")
wt_b_dubs=read.csv("WT_B/v3/WT_B_scrublet_08perc_thresh0.30_output_table.csv")
identical(rownames(wt_b@meta.data),wt_b_dubs$cell_barcodes)
wt_b$predicted_doublet = wt_b_dubs$predicted_doublet
wt_b$doublet_score=wt_b_dubs$doublet_score

DimPlot(wt_b, group.by="predicted_doublet")
FeaturePlot(wt_b, feature="doublet_score")
DimPlot(wt_b,label=T)

prop.table(table(wt_b$seurat_clusters,wt_b$predicted_doublet),1)

## remove doublets and cluster 9

wt_b_dubs$cluster = wt_b$seurat_clusters
wt_b_dubs$predicted_doublet[wt_b_dubs$cluster == 9] <-"True"

prop.table(table(wt_b_dubs$cluster,wt_b_dubs$predicted_doublet),1)

write.csv(wt_b_dubs, file="WT_B/v3/WT_B_v3_predicted_dubs_v2_clus9.csv")

######################### Hom A ###############################

hom_a=readRDS("Hom_A/v3/basic_analysis/Hom_A_v3_basic.rds")
hom_a_dubs=read.csv("Hom_A/v3/Hom_A_scrublet_08perc_thresh0.28_output_table.csv")

identical(rownames(hom_a@meta.data),hom_a_dubs$cell_barcodes)
hom_a$predicted_doublet = hom_a_dubs$predicted_doublet
hom_a$doublet_score = hom_a_dubs$doublet_score

DimPlot(hom_a, group.by="predicted_doublet")
FeaturePlot(hom_a, features="doublet_score")

DimPlot(hom_a,label=T)

prop.table(table(hom_a$seurat_clusters,hom_a$predicted_doublet),1)

## remove doublets and cluster 9,17

hom_a_dubs$cluster = hom_a$seurat_clusters
hom_a_dubs$predicted_doublet[hom_a_dubs$cluster == 9] <-"True"
hom_a_dubs$predicted_doublet[hom_a_dubs$cluster == 17] <-"True"

prop.table(table(hom_a_dubs$cluster,hom_a_dubs$predicted_doublet),1)

write.csv(hom_a_dubs, file="Hom_A/v3/Hom_A_v3_predicted_dubs_v2_clus9_17.csv")

######################### Hom B #################################

hom_b=readRDS("Hom_B/v3/basic_analysis/Hom_B_v3_basic.rds")
hom_b_dubs=read.csv("Hom_B/v3/Hom_B_scrublet_08perc_thresh0.30_output_table.csv")
identical(rownames(hom_b@meta.data),hom_b_dubs$cell_barcodes)
hom_b$predicted_doublet = hom_b_dubs$predicted_doublet
hom_b$doublet_score = hom_b_dubs$doublet_score

DimPlot(hom_b, group.by="predicted_doublet")
DimPlot(hom_b, label=T)

prop.table(table(hom_b$seurat_clusters,hom_b$predicted_doublet),1)

## remove doublets and cluster 13, 15

hom_b_dubs$cluster = hom_b$seurat_clusters
hom_b_dubs$predicted_doublet[hom_b_dubs$cluster == 15] <-"True"
hom_b_dubs$predicted_doublet[hom_b_dubs$cluster == 13] <-"True"

prop.table(table(hom_b_dubs$cluster,hom_b_dubs$predicted_doublet),1)

write.csv(hom_b_dubs, file="Hom_B/v3/Hom_B_v3_predicted_dubs_v2_clus15_13.csv")
