### Ad9 v2
library(Seurat)
library(dplyr)

setwd("Ad36")

seur=readRDS("v2/basic_analysis/Ad36_v2_basic.rds")

dubs = read.csv("v2/Ad36_v2_scrublet_12perc_thresh0.28_output_table.csv")

identical(rownames(seur@meta.data),dubs$cell_barcodes)

seur$doublet_score <- dubs$doublet_score
seur$predicted_doublet <- dubs$predicted_doublet

FeaturePlot(seur, "doublet_score")

pdf("v2/Ad36_v2_predicted_doublet.pdf",height=8,width=11,useDingbats = F)
DimPlot(seur, group.by="predicted_doublet", pt.size=1)
dev.off()

DimPlot(seur,label=T)

perc_doublets=prop.table(table(seur$seurat_clusters,seur$predicted_doublet),1)

## remove whole clusters if majority are doublets
## cluster 8 in resolution 0.8 has 45% doublets. Looks similar to previous datasets with doublet clusters and has a spread that indicates doublets - going to conservatively remove whole cluster even though < 50% doublets.

clus8=rownames(filter(seur@meta.data, seurat_clusters==8))

## make predicted_doublets say TRUE for all of cluster 10

dubs$cluster = seur$seurat_clusters
dubs$predicted_doublet[dubs$cluster == 8] <-"True"

seur$predicted_doublet <- dubs$predicted_doublet

pdf("v2/Ad36_v2_dubs_v2_clus8.pdf",heigh=8,width=11,useDingbats = F)
DimPlot(seur, group.by="predicted_doublet")
dev.off()

write.csv(dubs, file="v2/Ad36_v2_predicted_dubs_v2_clus8.csv")
