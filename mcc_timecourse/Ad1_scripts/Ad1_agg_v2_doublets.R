library(Seurat)

setwd("Ad1_agg")

seur= readRDS("/Volumes/Reiterlab_3/mcc_timecourse/Ad1_agg/v2/basic_analysis/Ad1_agg_v2_basic.rds")

dubs = read.csv("v2/Ad1_agg_v2_scrublet_8perc_thresh0.32_output_table.csv")

identical(rownames(seur@meta.data),dubs$cell_barcodes)

seur$doublet_score <- dubs$doublet_score
seur$predicted_doublet <- dubs$predicted_doublet

FeaturePlot(seur, "doublet_score")

pdf("v2/Ad1_agg_v2_predicted_doublet.pdf",height=8,width=11,useDingbats = F)
DimPlot(seur, group.by="predicted_doublet", pt.size=1,order=T)
dev.off()

DimPlot(seur,label=T)


perc_doublets=prop.table(table(seur$seurat_clusters,seur$predicted_doublet),1)
