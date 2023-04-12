seur=readRDS("v4/basic_analysis/Ad1_agg_soupx_v4_basic.rds")
DimPlot(seur)


ad1_cells=sample(rownames(seur@meta.data),.5 * length(rownames(seur@meta.data)))
ad1_subset = subset(seur, cells=ad1_cells)

DimPlot(ad1_subset)

saveRDS(ad1_subset, "v4/basic_analysis/Ad1_agg_soupx_v4_basic_subset.rds")

seur=readRDS("v3/mcc_timecourse_v3.rds")

DefaultAssay(seur) <- "RNA"

FeaturePlot(seur, c("Myb","Mycl","Mcidas"))

DimPlot(seur, label=T)

sub = subset(seur, idents = c(0:4,6:11))
DimPlot(sub)
FEaturePlot(sub, c("Mycl","Myb"))
