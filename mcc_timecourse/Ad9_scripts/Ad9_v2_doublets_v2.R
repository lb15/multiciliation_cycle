### Ad9 v2
library(Seurat)
library(dplyr)

setwd("Ad9")

seur=readRDS("v2/basic_analysis/Ad9_v2_basic.rds")

dubs = read.csv("v2/Ad9_v2_scrublet_12perc_thresh0.25_output_table.csv")

identical(rownames(seur@meta.data),dubs$cell_barcodes)

seur$doublet_score <- dubs$doublet_score
seur$predicted_doublet <- dubs$predicted_doublet

FeaturePlot(seur, "doublet_score")

pdf("v2/Ad9_v2_predicted_doublet.pdf",height=8,width=11,useDingbats = F)
DimPlot(seur, group.by="predicted_doublet", pt.size=2)
dev.off()


DimPlot(seur,label=T)

perc_doublets=prop.table(table(seur$seurat_clusters,seur$predicted_doublet),1)

## remove whole clusters if majority are doublets
## cluster 10, 13 in resolution 0.8.

clus10=rownames(filter(seur@meta.data, seurat_clusters==10))

## make predicted_doublets say TRUE for all of cluster 10

dubs$cluster = seur$seurat_clusters
dubs$predicted_doublet[dubs$cluster == 10] <-"True"
dubs$predicted_doublet[dubs$cluster == 13] <-"True"

seur$predicted_doublet <- dubs$predicted_doublet

DimPlot(seur, group.by="predicted_doublet")

write.csv(dubs, file="v2/Ad9_v2_predicted_dubs_v2_clus10_clus13.csv")
