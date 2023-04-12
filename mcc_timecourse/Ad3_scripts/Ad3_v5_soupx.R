## analysis for Ad3 with updated scrublet

library(Seurat)
library(dplyr)
library(DoubletFinder)
library(ggplot2)

## SoupX

out_folder ="/Volumes/LB_ReiterLab_2/10X_080519/Ad3/analysis/outs/"

sc=load10X(out_folder)
sc = autoEstCont(sc)
out = adjustCounts(sc)

dir.create("Ad3_soupx")
writeMM(out, file="Ad3_soupx/matrix.mtx")

## read in Soupx corrected data

data_dir ="Ad3_soupx/"
basename="Ad3_soupx"
folder="v5"

data = Read10X(data.dir=data_dir)

seur=CreateSeuratObject(counts=data, project="Ad3",min.cells=3,min.features=200)
seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^mt-")

seur <- subset(seur,nFeature_RNA > 3000 & percent.mt < 15)
ncount_mean = mean(seur$nCount_RNA)
std_count = sd(seur$nCount_RNA)
ncount_hi = ncount_mean + 2*std_count

seur <- subset(seur,nCount_RNA < ncount_hi)

png(paste0(folder,"/",basename,"_",folder,"_violinQC_filtered.png"),height = 800,width=1100)
VlnPlot(seur,features = c("nFeature_RNA","nCount_RNA","percent.mt"),group.by="orig.ident")
dev.off()

seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seur)
seur <- ScaleData(seur, features = all.genes)

seur <- RunPCA(seur, features = VariableFeatures(object = seur))

png(paste0(folder,"/",basename,"_",folder,"_elbow.png"),height = 800,width=1100)
ElbowPlot(seur)
dev.off()

seur <- FindNeighbors(seur, dims = 1:15)
seur <- FindClusters(seur, resolution = 0.3)

seur <- RunUMAP(seur, dims = 1:15)
DimPlot(seur, reduction = "umap")

## remove doublets
dubs = read.csv("v5/Ad3_v5_scrublet_12perc_thresh0.38_output_table_rerun.csv")

table(dubs$predicted_doublet)

cells_dubs = dubs$cell_barcodes[dubs$predicted_doublet == "True"]

seur$scrublet <- "Singlet"
seur$scrublet[match(cells_dubs,rownames(seur@meta.data))] <- "Doublet"
DimPlot(seur, group.by="scrublet")

Idents(seur) <- "scrublet"
seur_sub = subset(seur, idents = "Singlet")

seur_sub <- NormalizeData(seur_sub, normalization.method = "LogNormalize", scale.factor = 10000)
seur_sub <- FindVariableFeatures(seur_sub, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seur_sub)
seur_sub <- ScaleData(seur_sub, features = all.genes)

seur_sub <- RunPCA(seur_sub, features = VariableFeatures(object = seur_sub))

ElbowPlot(seur_sub)

seur_sub <- FindNeighbors(seur_sub, dims = 1:15)
seur_sub <- FindClusters(seur_sub, resolution = 0.3)

seur_sub <- RunUMAP(seur_sub, dims = 1:15)
DimPlot(seur_sub, reduction = "umap")

FeaturePlot(seur_sub, "Mcidas")

saveRDS(seur_sub, file="v5/Ad3_soupx_v5.rds")
