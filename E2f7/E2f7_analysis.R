#### Analysis script for integrated E2f7 WT and Hom dataset.

library(Seurat)
library(ggplot2)
library(viridis)

###### ####   color palette ##################

mcc_color <- colorRampPalette(c("lightskyblue2", "royalblue4"))

colors=c("#009245","#F7A12D","darkgoldenrod","indianred",mcc_color(5)[5],"#996699","plum","#998675",mcc_color(5)[1],mcc_color(5)[3])

mcc_clean_col=c(colors[1],mcc_color(4)[c(4,2,3,1)])

################### OBJECTS ##############################
setwd("E2f7_em3_merge/v2/")

## created from seurat.integration.R script and updated in script below
seur=readRDS("E2f7_em3_merge_v2.rds")

## generated in this script below
mccs_clean=readRDS("mccs_clean/E2f7_em3_mccs_clean.rds")


DimPlot(seur, label=T)

DefaultAssay(seur) <- "integrated"

seur=FindNeighbors(seur, dims=1:15)
seur=RunUMAP(seur, dims=1:15)
seur=FindClusters(seur, dims=1:15)

ElbowPlot(seur, ndims=40)

DimPlot(seur, label=T)


#### label with mcc_timecourse labels

mcc_timecourse=readRDS("/Volumes/Reiterlab_3/mcc_timecourse/v3/mcc_timecourse_v3_13dim.rds")
DimPlot(mcc_timecourse, label=T,group.by="integrated_snn_res.0.3")

anchors <- FindTransferAnchors(reference = mcc_timecourse, query = seur,
                                        dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = mcc_timecourse$integrated_snn_res.0.3,dims = 1:20)

write.csv(predictions, file="E2f7_em3_merge/v2/E2f7_em3_merge_predictions_mcc_timecourse.csv")


###################### TRICYCLE ######################
library(tricycle)

sc=as.SingleCellExperiment(seur,assay = "RNA")

sc=project_cycle_space(sc, species = "mouse",gname.type = "SYMBOL")
sc=estimate_cycle_position(sc)

data=as.data.frame(reducedDim(sc,type = "tricycleEmbedding"))
lab= colData(sc)$tricyclePosition

data$tricyclePosition <- lab

ggplot(data, aes(x=PC1, y=PC2,color=tricyclePosition))+
        geom_point()+
        scale_color_gradientn(colors=c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"))

top2a.idx <- which(rownames(sc) == 'Top2a')
fit.l <- fit_periodic_loess(sc$tricyclePosition,
                            assay(sc, 'logcounts')[top2a.idx,],
                            plot = TRUE,
                            x_lab = "Cell cycle position \u03b8", y_lab = "log2(Top2a)",
                            fig.title = paste0("Expression of Top2a along \u03b8 (n=",
                                               ncol(sc), ")"))
names(fit.l)
fit.l$fig + theme_bw(base_size = 14)

identical(colnames(seur),colnames(sc))
seur$tricyclePosition <- sc$tricyclePosition

FeaturePlot(seur, "tricyclePosition")+scale_color_gradientn(colors=c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"))

write.csv(seur@meta.data[,"tricyclePosition",drop=F],file="E2f7_em3_merge/v2/E2f7_em3_tricyclePosition.csv")

#tricyclepos = read.csv("E2f7_em3_merge/v2/E2f7_em3_tricyclePosition.csv",row.names=1)

identical(rownames(tricyclepos),rownames(seur@meta.data))
seur$tricyclePosition=tricyclepos$tricyclePosition

saveRDS(seur, file="E2f7_em3_merge/v2/E2f7_em3_merge_v2.rds")

seur$orig.ident <- factor(seur$orig.ident,levels=c("WT_A","WT_B","Hom_A","Hom_B"))

pdf("E2f7_em3_merge/v2/E2f7_em3_merge_v2_clusters_pred_split.pdf",height=2,width=11,useDingbats = F)
DimPlot(seur, group.by="timecourse_pred",cols=colors,pt.size=0.5,split.by="orig.ident")&theme_void() &scale_x_reverse()&theme(strip.text = element_blank(),title = element_blank(),legend.position = "none")
dev.off()

DimPlot(seur, group.by="timecourse_pred",cols=colors,pt.size=0.5,split.by="orig.ident")

###################### SUBSET MULTICILIATED CELLS ################

predictions=read.csv("E2f7_em3_merge/v2/E2f7_em3_merge_predictions_mcc_timecourse.csv",row.names=1)

seur$timecourse_pred = predictions$predicted.id

pdf("E2f7_em3_merge/v2/E2f7_em3_v2_predictions.pdf",height=8,width=11,useDingbats = F)
DimPlot(seur, group.by="timecourse_pred",label=T,pt.size=1)+theme_void()
dev.off()


Idents(seur) <- "timecourse_pred"
mccs = subset(seur , idents=c(4,8,9,0))

DimPlot(mccs, label=T)

DefaultAssay(mccs) <- "integrated"
mccs = ScaleData(mccs)
mccs=RunPCA(mccs, npcs = 40)
ElbowPlot(mccs, ndims = 40)
mccs=RunUMAP(mccs,dims = 1:13)
mccs=FindNeighbors(mccs, dims = 1:13)
mccs=FindClusters(mccs,resolution =0.4)
DimPlot(mccs, label=T)

dir.create("E2f7_em3_merge/v2/mcc_subset")
saveRDS(mccs,file="E2f7_em3_merge/v2/E2f7_em3_mcc_subset.rds")
mccs=readRDS("E2f7_em3_merge/v2/E2f7_em3_mcc_subset.rds")

## add wild-type IDs

wt_mccs = readRDS("/Volumes/Reiterlab_3/mcc_timecourse/v3/mccs_clean_subset/mccs_clean.rds")
DimPlot(wt_mccs, label=T,group.by="integrated_snn_res.0.2")
DefaultAssay(wt_mccs) <- "RNA"
wt_mccs = FindVariableFeatures(wt_mccs)
DefaultAssay(mccs) <- "RNA"
anchors <- FindTransferAnchors(reference = wt_mccs, query = mccs,
                               dims = 1:20, reference.reduction = "pca")

predictions <- TransferData(anchorset = anchors, refdata = wt_mccs$integrated_snn_res.0.2,dims = 1:20)

identical(rownames(mccs@meta.data),rownames(predictions))

mccs$mcc_predictions=predictions$predicted.id
mccs$pred_0 = predictions$prediction.score.0

pdf("E2f7_em3_merge/v2/mcc_subset/E2f7_em3_mccs_res0.4_clusters.pdf",height=8,width=11,useDingbats = F)
DimPlot(mccs, label=T,group.by="integrated_snn_res.0.4",pt.size=1)+theme_void()
dev.off()
pdf("E2f7_em3_merge/v2/mcc_subset/E2f7_em3_mccs_predictions.pdf",height=8,width=11,useDingbats = F)
DimPlot(mccs, group.by="mcc_predictions",label=T,pt.size=1)+theme_void()
dev.off()

FeaturePlot(mccs, features=c("pred_1"))
mccs$clus1_pred=predictions$prediction.score.1
FeaturePlot(mccs, features="clus1_pred")

Idents(mccs) <- "integrated_snn_res.0.4"
marks=FindAllMarkers(mccs, assay="RNA",only.pos=T)

write.csv(marks, file="E2f7_em3_merge/v2/mcc_subset/E2f7_em3_mcc_subset_res0.4_markers.csv")
weird_genes = c("Areg",
                "Nupr1",
                "Crip1",
                "Clca3b","Agr2",
                "Cyp2a5",
                "Cadm2",
                "Ces1f","Nrxn1",
                "Gnat3",
                "St18",
                "Ascl1")
png("E2f7_em3_merge/v2/mcc_subset/E2f7_em3_mcc_weirdgenes.png",height=800,width=1100)
FeaturePlot(mccs, weird_genes)&theme_void()&theme(legend.position="none")
dev.off()
pdf("E2f7_em3_merge/v2/mcc_subset/E2f7_em3_mcc_maturemcc_genes.pdf")
FeaturePlot(mccs, c("AU040972","Tmem212","Clca3b","Sec14l3"))&theme_void()&theme(legend.position="none")
dev.off()

png("E2f7_em3_merge/v2/mcc_subset/E2f7_em3_mcc_moreweird_genes.png",height=400,width=800)
FeaturePlot(mccs, c("S100g","Itln1"))&theme_void()&theme(legend.position="none")
dev.off()


#### remove cluster 6 and 8, which I believe could contain doublets similar to the groups removed in the wild-type timecourse dataset.
## and 10 which looks like neuroendocrine
Idents(mccs) <- "integrated_snn_res.0.4"
DimPlot(mccs, label=T)
mccs_sub = subset(mccs, idents=c(0:5,7:8))
DimPlot(mccs_sub, label=T)

DefaultAssay(mccs_sub) <- "integrated"
mccs_sub = ScaleData(mccs_sub)
mccs_sub=RunPCA(mccs_sub, npcs = 40)
ElbowPlot(mccs_sub, ndims = 40)
mccs_sub=FindNeighbors(mccs_sub, dims = 1:11)
mccs_sub=RunUMAP(mccs_sub,dims = 1:11)
mccs_sub=FindClusters(mccs_sub,resolution =0.3)

pdf("E2f7_em3_merge/v2/mcc_subset/E2f7_em3_mccs_clean_res0.3_clusters.pdf",height=8,width=11,useDingbats = F)
DimPlot(mccs_sub, label=T,pt.size=1)+theme_void()
dev.off()

sub_marks=FindAllMarkers(mccs_sub, assay="RNA",idents="integrated_snn_res.0.3",only.pos = T)

write.csv(sub_marks, file="E2f7_em3_merge/v2/mcc_subset/E2f7_em3_mccs_clean_res0.3_markers.csv")

saveRDS(mccs_sub, file="E2f7_em3_merge/v2/mcc_subset/E2f7_em3_mccs.rds")
mccs_sub=readRDS("E2f7_em3_merge/v2/mcc_subset/E2f7_em3_mccs.rds")

##### round two of clean-up. remove cluster 5 (appears to be remains of side population) and cluster 6 (which still seems like remaining transitory/dying/doublet of mature MCC and other genes)

DefaultAssay(mccs_sub) <- "integrated"
DimPlot(mccs_sub, label=T,group.by="integrated_snn_res.0.3")

mccs_clean=subset(mccs_sub, idents=c(0:4))
DimPlot(mccs_clean)
mccs_clean = ScaleData(mccs_clean)
mccs_clean=RunPCA(mccs_clean, npcs = 40)
ElbowPlot(mccs_clean, ndims = 40)
mccs_clean=FindNeighbors(mccs_clean, dims = 1:11)
mccs_clean=RunUMAP(mccs_clean,dims = 1:11)
mccs_clean=FindClusters(mccs_clean,resolution =0.3)

dir.create("E2f7_em3_merge/v2/mccs_clean")
setwd("E2f7_em3_merge/v2")
pdf("mccs_clean/E2f7_em3_mcc_clean_res0.3_clusters.pdf",height=8,width=11,useDingbats = F)
DimPlot(mccs_clean,label=T, pt.size=1) + theme_void()
dev.off()

DefaultAssay(mccs_clean) <- "RNA"
DefaultAssay(wt_mccs) <- "RNA"
wt_mccs=FindVariableFeatures(wt_mccs)
mccs_clean=FindVariableFeatures(mccs_clean)
anchors_clean <- FindTransferAnchors(reference = wt_mccs, query = mccs_clean,
                               dims = 1:15, reference.reduction = "pca")
predictions_clean <- TransferData(anchorset = anchors_clean, refdata = wt_mccs$integrated_snn_res.0.2,dims = 1:15)

mccs_clean$mccs_pred = predictions_clean$predicted.id
DimPlot(mccs_clean, group.by="mccs_pred")
DimPlot(mccs_clean, label=T,split.by="orig.ident",group.by="mccs_pred")
pdf("mccs_clean/E2f7_em3_mccs_clean_predictions.pdf",height=8,width=11,useDingbats = F)
DimPlot(mccs_clean, group.by="mccs_pred",label=T,cols=c(colors[1],mcc_color(5)[c(5,3,1)]),pt.size=1)
dev.off()

saveRDS(mccs_clean, file="mccs_clean/E2f7_em3_mccs_clean.rds")

mccs_clean_col=c(colors[1],mcc_color(5)[c(5,3,1)])


######################## PSEUDOTIME ########################

library(monocle3)
library(SeuratWrappers)

mccs_clean=readRDS("mccs_clean/E2f7_em3_mccs_clean.rds")

cds=as.cell_data_set(mccs_clean)
cds=cluster_cells(cds)
cds=learn_graph(cds,use_partition = F)

pdf("mccs_clean/E2f7_em3_mccs_clean_pseudotime.pdf",height=8,width=11,useDingbats = F)
plot_cells(cds,cell_size = 1,trajectory_graph_segment_size = 2)
dev.off()

cds=order_cells(cds)

cds_sub=choose_graph_segments(cds,clear_cds = F)


plot_cells(cds_sub,cell_size = 1)


cds_sub=order_cells(cds_sub)


plot_cells(cds_sub, color_cells_by = "pseudotime")
ps_df=as.data.frame(pseudotime(cds_sub))
colnames(ps_df)="pseudotime"
ps_df$pseudo_bin = cut_interval(ps_df$pseudotime, n=20,labels=F)

write.csv(ps_df, file="mccs_clean/E2f7_em3_mcc_clean_pseudotime_values.csv")
saveRDS(cds, file="mccs_clean/E2f7_em3_mcc_clean_full_monocle_object.rds")
saveRDS(cds_sub, file="mccs_clean/E2f7_em3_mcc_clean_sub_monocle_object.rds")

mccs_clean$pseudotime <- NA
mccs_clean$pseudotime[match(rownames(ps_df),rownames(mccs_clean@meta.data))] <- ps_df$pseudotime

mccs_clean$pseudo_bin <- NA
mccs_clean$pseudo_bin[match(rownames(ps_df),rownames(mccs_clean@meta.data))] <- ps_df$pseudo_bin

DimPlot(mccs_clean, group.by="pseudo_bin")
pdf("mccs_clean/E2f7_em3_mccs_clean_pseudotime.pdf",height=8,width=11,useDingbats = F)
FeaturePlot(mccs_clean, features="pseudotime",pt.size=1)+scale_color_viridis()
dev.off()


##################### SEURAT CELL CYCLING SCORING ##############
library(UCell)
library(ggpubr)
library(dplyr)
library(scales)

cc_s=read.csv("~/Documents/Reiter_Seq/s_phase_human_mgi_genes.csv")
cc_g2m=read.csv("~/Documents/Reiter_Seq/g2m_phase_human_mgi_genes.csv")

sgenes=cc_s$MGI.symbol
g2mgenes=cc_g2m$MGI.symbol

mccs_clean = AddModuleScore_UCell(mccs_clean,features = list(sgenes,g2mgenes),assay = "RNA") 
mccs_clean$orig.ident=factor(mccs_clean$orig.ident,levels=c("WT_A","WT_B","Hom_A","Hom_B"))

VlnPlot(mccs_clean, features=c("signature_1_UCell","signature_2_UCell"),group.by="orig.ident")

s=ggplot(mccs_clean@meta.data %>% filter(integrated_snn_res.0.3 %in% c(1:4)),aes(x=orig.ident,y=signature_1_UCell))+
        geom_violin()+ylab("S Score")+
        theme_classic()+
        scale_y_continuous(labels=label_number(accuracy=.01))+
        theme(axis.text=element_text(size=12,color="black"),axis.title.x = element_blank())

g2m=ggplot(mccs_clean@meta.data %>% filter(integrated_snn_res.0.3 %in% c(1:4)),aes(x=orig.ident,y=signature_2_UCell))+
        geom_violin()+ylab("G2M Score")+
        theme_classic()+
        scale_y_continuous(labels=label_number(accuracy=.01))+
        theme(axis.text=element_text(size = 12,color="black"),axis.title.x = element_blank())


ggsave(ggarrange(s,g2m,nrow=2), file = "mccs_clean/E2f7_em3_mccs_clean_S_G2M_score_violin.pdf",width=6,height=6)

mccs_clean$group<-factor(mccs_clean$group,levels=c("WT","Hom"))

pdf("mccs_clean/E2f7_em3_mccs_clean_S_G2M_featureplots.pdf",height=6,width=8,useDingbats = F)
FeaturePlot(mccs_clean, features=c("signature_1_UCell","signature_2_UCell"),split.by="group",max.cutoff = "q99",pt.size=1) & theme_void() & scale_y_reverse() &scale_color_viridis() & theme(legend.position="none",strip.background = element_blank(),strip.text.x  = element_text(size=0)) & ggtitle("")
dev.off()


pdf("mccs_clean/E2f7_em3_mccs_clean_S_G2M_featureplots_legends.pdf",height=6,width=8,useDingbats = F)
FeaturePlot(mccs_clean, features=c("signature_1_UCell","signature_2_UCell"),split.by="group",max.cutoff = "q99",pt.size=1,raster = T) & theme_void() & scale_y_reverse() &scale_color_viridis() & theme(strip.background = element_blank(),strip.text.x  = element_text(size=0)) & ggtitle("")
dev.off()



########################  DE ANALYSIS ####################### 
##will use pseudobulk methods
library(DESeq2)
library(Matrix.utils)
library(trqwe)
library(purrr)
library(apeglm)
library(tibble)
library(Seurat)
library(dplyr)

setwd("E2f7_em3_merge/v2")
mccs_clean=readRDS("mccs_clean/E2f7_em3_mccs_clean.rds")

mccs_clean$orig.ident<-factor(mccs_clean$orig.ident)
mccs_clean$group <- "WT"
mccs_clean$group[grep("Hom*",mccs_clean$orig.ident)]<-"Hom"
table(mccs_clean$group)
sce=as.SingleCellExperiment(mccs_clean,assay = "RNA")
head(colData(sce))

### sample_id = orig.ident or 4 samples
### cluster_id = use integrated_snn_res.0.3

kids <- purrr::set_names(levels(sce$integrated_snn_res.0.3))
# Total number of clusters
nk <- length(kids)
# Named vector of sample names
sids <- purrr::set_names(levels(sce$orig.ident))
# Total number of samples 
ns <- length(sids)
ns


n_cells <- as.numeric(table(sce$orig.ident))
m <- match(sids, sce$orig.ident)
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
        select(c("orig.ident","group","n_cells"))
ei

# Aggregate the counts per sample_id and cluster_id
# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("integrated_snn_res.0.3", "orig.ident")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 
class(pb)
dim(pb)
pb[1:6,1:6]

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), `[`, 1)
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
        lapply(function(u) 
                set_colnames(t(u), 
                             stringr::str_remove(rownames(u), "^[0-9]_")))
str(pb)

#Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
        pb[[x]] %>%
                colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
        unlist()
# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
        rep(names(pb)[x], 
            each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
        unlist()
# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)
colnames(gg_df) <- c("cluster_id","orig.ident")
gg_df <- left_join(gg_df, ei[, c("orig.ident", "group")]) 


metadata <- gg_df %>%
        dplyr::select(cluster_id, orig.ident, group) 

metadata$cluster_id <- factor(metadata$cluster_id)
clusters <- levels(metadata$cluster_id)
clusters

## write function for DESeq on each cluster

deseq_clusters = function(pb, cluster_metadata, clus){
        cluster_metadata <- metadata[which(metadata$cluster_id == clus), ]
        head(cluster_metadata)

        # Assign the rownames of the metadata to be the sample IDs
        rownames(cluster_metadata) <- cluster_metadata$orig.ident
        head(cluster_metadata)

# Subset the counts to only clus
        counts <- pb[[match(clus,clusters)]]

        cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

        # Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
        all(rownames(cluster_metadata) == colnames(cluster_counts))   

        dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group)
        dds$group=relevel(dds$group, ref = "WT")
        
        # Transform counts for data visualization
        rld <- rlog(dds, blind=TRUE)

        # Plot PCA
        png(paste0("mccs_clean/E2f7_em3_mcc_clean_DESeq_PCA_clus",clus,".png"),height=400,width=400)
        print(DESeq2::plotPCA(rld, intgroup = "group"))
        dev.off()
        
        dds <- DESeq(dds)
        
        png(paste0("mccs_clean/E2f7_em3_mcc_clean_DESeq_DispEsts_clus",clus,".png"),height=400,width=400)
        plotDispEsts(dds)
        dev.off()

        # Output results of Wald test for contrast for WT vs. Hom
        cluster_metadata$group <- factor(cluster_metadata$group)
        levels(cluster_metadata$group)[2]
        levels(cluster_metadata$group)[1]

        contrast <- c("group", levels(cluster_metadata$group)[1], levels(cluster_metadata$group)[2])

        #dds$group=relevel(dds$group, ref = "WT")
        # resultsNames(dds)

        res2 <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

        res <- lfcShrink(dds, 
                 coef=2,
                 res=res2,lfcThreshold = 0.58)
        
        res$padj <- res2$padj
        
        res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
        
        write.csv(res_tbl,paste0("mccs_clean/E2f7_em3_mcc_clean_DESeq_genes_Hom_vs_WT_clus",clus,".csv"))
        # Set thresholds
        padj_cutoff <- 0.05
        
        # Subset the significant results
        sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
                dplyr::arrange(padj)
        
        write.csv(sig_res,paste0("mccs_clean/E2f7_em3_mcc_clean_DESeq_siggenes_Hom_vs_WT_clus",clus,".csv"))



}

for(x in clusters){
        print(x)
        deseq_clusters(pb, cluster_metadata, x)
}


## what p-value threshold to choose?
VlnPlot(mccs_clean, features=c("Mcm3","Hells"),group.by="integrated_snn_res.0.3",split.by="orig.ident")

library(EnhancedVolcano)
clusters=c("0","1","2","3","4")
for(clus in clusters){
        clus1 =read.csv(paste0("mccs_clean/E2f7_em3_mcc_clean_DESeq_genes_Hom_vs_WT_clus",clus,".csv"),row.names=1)
        clus1_filt=clus1 %>% filter(!is.na(log2FoldChange))
        
        png(paste0("mccs_clean/E2f7_em3_mcc_clean_clus",clus,"_enhancedvolcano.png"),height=400,width=600)
        print(EnhancedVolcano(clus1_filt, 
                        lab=clus1_filt$gene,
                        x="log2FoldChange",drawConnectors = T,
                        y="svalue",labSize = 6,pointSize = 4,
                        pCutoff = .005,FCcutoff = log2(1.5)))
        dev.off()

}

## make master list of DE genes
### read in unfiltered list, combine

clus0=read.csv("mccs_clean/E2f7_em3_mcc_clean_DESeq_genes_Hom_vs_WT_clus0.csv",row.names = 2)
clus1=read.csv("mccs_clean/E2f7_em3_mcc_clean_DESeq_genes_Hom_vs_WT_clus1.csv",row.names=2)
clus2=read.csv("mccs_clean/E2f7_em3_mcc_clean_DESeq_genes_Hom_vs_WT_clus2.csv",row.names=2)
clus3=read.csv("mccs_clean/E2f7_em3_mcc_clean_DESeq_genes_Hom_vs_WT_clus3.csv",row.names=2)
clus4=read.csv("mccs_clean/E2f7_em3_mcc_clean_DESeq_genes_Hom_vs_WT_clus4.csv",row.names=2)

colnames(clus0)=paste0("clus0_",colnames(clus0))
colnames(clus1)=paste0("clus1_",colnames(clus1))
colnames(clus2)=paste0("clus2_",colnames(clus2))
colnames(clus3)=paste0("clus3_",colnames(clus3))
colnames(clus4)=paste0("clus4_",colnames(clus4))

full_dat=cbind(clus0%>% select(-"clus0_X"),clus1%>% select(-"clus1_X"))
full_dat=cbind(full_dat,clus2%>% select(-"clus2_X"))
full_dat=cbind(full_dat, clus3%>% select(-"clus3_X"))
full_dat=cbind(full_dat,clus4%>% select(-"clus4_X"))


write.csv(full_dat, file="mccs_clean/E2f7_em3_mccs_clean_DESeq_genes_Hom_Vs_WT_allclusters.csv")
full_dat=read.csv("mccs_clean/E2f7_em3_mccs_clean_DESeq_genes_Hom_Vs_WT_allclusters.csv",row.names=1)

#### rename the clusters to match the Figure. Should have done this before I ran all these analyses, but alas. I jsut renamed in the final excel document.
##old clus id = new clus id
## cluster 0 = cluster a
## cluster 4=cluster b
## cluster 2=cluster c
## cluster 3=cluster d
## cluster 1 = cluster e


clus0_sig = clus0 %>% filter(abs(log2FoldChange) > log2(2) & padj < 0.05)
clus1_sig = clus1 %>% filter(abs(log2FoldChange) > log2(2) & padj < 0.05)
clus2_sig = clus2 %>% filter(abs(log2FoldChange) > log2(2) & padj < 0.05)
clus3_sig = clus3 %>% filter(abs(log2FoldChange) > log2(2) & padj < 0.05)
clus4_sig = clus4 %>% filter(abs(log2FoldChange) > log2(2) & padj < 0.05)

sig_genes =unique(c(rownames(clus0_sig),rownames(clus1_sig),rownames(clus2_sig),rownames(clus3_sig),rownames(clus4_sig)))

full_dat_sig=full_dat[match(sig_genes, rownames(full_dat)),]

### also supplementary Table #3
cilia=read.csv("~/Documents/Reiter_Seq/CuratedCiliaAndBasalBodyGenes_18.csv",row.names=1)

cilia_sig=cilia[cilia$Gene.name..Mus.musculus. %in% rownames(full_dat_sig),]

sig_dat=full_dat_sig[match(cilia_sig$Gene.name..Mus.musculus.,rownames(full_dat_sig)),]
###################### FIGURE MAKING ########################


### UMAP
setwd("E2f7_em3_merge/v2")
mccs_clean=readRDS(file="mccs_clean/E2f7_em3_mccs_clean.rds")

mccs_clean$group <- "WT"
mccs_clean$group[grep("Hom*",mccs_clean$orig.ident)]<-"Hom"

pdf("mccs_clean/E2f7_em3_mccs_clean_clusters.pdf",height=8,width=11,useDingbats = F)
DimPlot(mccs_clean, label=F,pt.size = 1,group.by = "integrated_snn_res.0.3",cols = c(colors[1],mcc_color(6)[c(6,2,4,1)])) + theme_void()+scale_y_reverse()+ggtitle("")+theme(legend.position = "none",strip.text  = element_blank())
dev.off()

pdf("mccs_clean/E2f7_em3_mccs_clean_clusters_split.pdf",height=4,width=11,useDingbats = F)
DimPlot(mccs_clean, label=F,pt.size = 1,group.by = "integrated_snn_res.0.3",cols = c(colors[1],mcc_color(6)[c(6,2,4,1)]),split.by="geno") & theme_void()&scale_y_reverse()&ggtitle("")&theme(legend.position = "none",strip.text  = element_blank())
dev.off()


#### Pseudotime

m=FeaturePlot(mccs_clean,"pseudotime",pt.size=1)+scale_y_reverse()+scale_color_viridis(direction=-1,option = "C")+theme_void()+theme(title=element_blank())

pdf("mccs_clean/E2f7_em3_mccs_clean_pseudotime.pdf",height=8,width=11,useDingbats = F)
FeaturePlot(mccs_clean,"pseudotime",pt.size=1)+scale_y_reverse()+scale_color_viridis(direction=-1,option = "C")+theme_void()+theme(legend.position = "none",title=element_blank())
dev.off()

leg=get_legend(m)

pdf("mccs_clean/E2f7_em3_mccs_clean_pseudotime_legend.pdf",height=2,width=2,useDingbats = F)
grid.newpage()
grid.draw(leg)
dev.off()


#### Enhanced volcano
library(EnhancedVolcano)
clus2=read.csv("mccs_clean/E2f7_em3_mcc_clean_DESeq_genes_Hom_vs_WT_clus2.csv",row.names=2)

clus_fc = clus2 %>% filter(abs(log2FoldChange) > log2(1.5) & padj < 1e-05)

clus2$FC1.5_padj0.00001 = "NO"
clus2$FC1.5_padj0.00001[match(rownames(clus_fc),rownames(clus2))] = "YES"

write.csv(clus2,"mccs_clean/E2f7_em3_mcc_clean_DESeq_genes_Hom_vs_WT_clus2_rawdata.csv")

gene_labels=c("Gins2","Hells","Mcm4","Mcm7","Stmn1","Pcna","Mcm3","Cdt1","E2f1","Mcm5")

pdf("mccs_clean/E2f7_em3_mccs_clean_DESeq_volcano.pdf",height=6,width=9,useDingbats = F)
EnhancedVolcano(clus2,
                lab=rownames(clus2),FCcutoff = log2(1.5),pCutoff = 1e-05,selectLab = gene_labels,
                x = "log2FoldChange",y = "padj",
                subtitle = "",titleLabSize = 0,subtitleLabSize = 1,captionLabSize = 0,
                pointSize = 4,colAlpha = 1,
                col = c("grey30","grey30","grey30","#1b75ba"),
                drawConnectors = T,
                gridlines.major = F,gridlines.minor = F)
dev.off()


###################### SUBSET CYCLING BASAL CELLS ################

predictions=read.csv("E2f7_em3_merge/v2/E2f7_em3_merge_predictions_mcc_timecourse.csv",row.names=1)

seur$timecourse_pred = predictions$predicted.id

DimPlot(seur, group.by="timecourse_pred")
Idents(seur)<-"timecourse_pred"
cycling=subset(seur, idents=c(1,5,6))
DimPlot(cycling)

DefaultAssay(cycling) <- "integrated"

cycling= ScaleData(cycling)
cycling=RunPCA(cycling, npcs = 40)
ElbowPlot(cycling, ndims = 40)
cycling=RunUMAP(cycling,dims = 1:13)
cycling=FindNeighbors(cycling, dims = 1:13)
cycling=FindClusters(cycling,resolution =0.4)
DimPlot(cycling, label=T)

DefaultAssay(cycling) <- "RNA"
FeaturePlot(cycling, "Krt5",min.cutoff = 0.5)

k5=cycling["Krt5",]
k5_names=names(which(k5[["RNA"]]@data[1,] >0.5))

basal_cycling =subset(cycling, cells=k5_names)
DimPlot(basal_cycling)
FeaturePlot(basal_cycling,"Krt5")

DefaultAssay(basal_cycling) <- "integrated"

basal_cycling <- ScaleData(basal_cycling)
basal_cycling <- RunPCA(basal_cycling, npcs = 40)
ElbowPlot(basal_cycling, ndims = 40)

basal_cycling <- RunUMAP(basal_cycling, dims = 1:13)
basal_cycling <- FindNeighbors(basal_cycling, dims=1:13)
basal_cycling <- FindClusters(basal_cycling, resolution = 0.3)

DimPlot(basal_cycling)

tricycle=read.csv("E2f7_em3_tricyclePosition.csv",row.names=1)

tri_basal=tricycle[match(rownames(basal_cycling@meta.data),rownames(tricycle)),,drop=F]
identical(rownames(tri_basal),rownames(basal_cycling@meta.data))

basal_cycling$tricyclePosition = tri_basal$tricyclePosition
FeaturePlot(basal_cycling,"tricyclePosition")+scale_color_gradientn(colors=c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"))

FeaturePlot(basal_cycling,"Krt5")

basal_cycling$tricyclePhase <- "G1/G0"
basal_cycling$tricyclePhase[basal_cycling$tricyclePosition > pi/2 & basal_cycling$tricyclePosition < pi] <- "S"
basal_cycling$tricyclePhase[basal_cycling$tricyclePosition >=pi & basal_cycling$tricyclePosition < pi*7/4] <- "G2M"

seur_basal=seur@meta.data[match(rownames(basal_cycling@meta.data),rownames(seur@meta.data)),]
identical(rownames(seur_basal),rownames(basal_cycling@meta.data))

basal_cycling$timecourse_pred <- seur_basal$timecourse_pred

############ BASAL CELL CLUSTERS, TRICYCLE ########
pdf("cycling_basal/E2f7_em3_cycling_basal_clusters.pdf",height=8,width=11,useDingbats = F)
DimPlot(basal_cycling, pt.size=1.5,group.by="timecourse_pred",cols = colors[c(2,6,7)]) + theme_void()+theme(legend.position = "none",title=element_blank())+scale_x_reverse()
dev.off()

basal_cycling$geno <-factor(basal_cycling$geno,levels=c("WT","Hom"))

pdf("cycling_basal/E2f7_em3_cycling_basal_clusters_geno_split.pdf",height=4,width=11,useDingbats = F)
DimPlot(basal_cycling, pt.size=1.5,group.by="timecourse_pred",split.by="geno",cols = colors[c(2,6,7)]) + theme_void()+theme(legend.position = "none",title=element_blank(),strip.text = element_blank())+scale_x_reverse()
dev.off()

pdf("cycling_basal/E2f7_em3_cycling_basal_tricycle.pdf",height=4,width=11,useDingbats = F)
FeaturePlot(basal_cycling,"tricyclePosition",split.by="geno",pt.size=1.5)&scale_color_gradientn(colors=c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9")) & theme_void() & theme(title=element_blank(),legend.position = "none")
dev.off()

dir.create("E2f7_em3_merge/v2/cycling_basal")
saveRDS(basal_cycling,file="cycling_basal/E2f7_em3_cycling_basal_subset.rds")

##################### COMPARE MITOTIC AND MCC ####################
tricycle=read.csv("E2f7_em3_tricyclePosition.csv",row.names=1)

basal_cycling=readRDS("cycling_basal/E2f7_em3_cycling_basal_subset.rds")

tri_mcc=tricycle[match(rownames(mccs_clean@meta.data),rownames(tricycle)),,drop=F]
identical(rownames(tri_mcc),rownames(mccs_clean@meta.data))

mccs_clean$tricyclePosition = tri_mcc$tricyclePosition

FeaturePlot(mccs_clean,"tricyclePosition",split.by="orig.ident")&scale_color_gradientn(colors=c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"))

mccs_clean$tricyclePhase <- "G1/G0"
mccs_clean$tricyclePhase[mccs_clean$tricyclePosition > pi/2 & mccs_clean$tricyclePosition < pi] <- "S"
mccs_clean$tricyclePhase[mccs_clean$tricyclePosition >=pi & mccs_clean$tricyclePosition < pi*7/4] <- "G2M"

library(UCell)
cc_s=read.csv("~/Documents/Reiter_Seq/s_phase_human_mgi_genes.csv")
cc_g2m=read.csv("~/Documents/Reiter_Seq/g2m_phase_human_mgi_genes.csv")

sgenes=cc_s$MGI.symbol
g2mgenes=cc_g2m$MGI.symbol
DefaultAssay(basal_cycling)<-"RNA"
basal_cycling=AddModuleScore_UCell(basal_cycling, features = list(sgenes,g2mgenes),assay="RNA")
mccs_clean=AddModuleScore_UCell(mccs_clean, features = list(sgenes,g2mgenes),assay="RNA")

basal_cycling$group <- "Basal"
basal_cycling$geno <-"WT"
basal_cycling$geno[grep("Hom*",basal_cycling$orig.ident)]<-"Hom"

mccs_clean$group <- "MCC"
mccs_clean$geno<-"WT"
mccs_clean$geno[grep("Hom*",mccs_clean$orig.ident)]<-"Hom"

basal_meta=basal_cycling@meta.data[,c("signature_1_UCell","signature_2_UCell","orig.ident","group","geno","tricyclePosition","timecourse_pred")]
colnames(basal_meta)<-c("signature_1_UCell","signature_2_UCell","orig.ident","group","geno","tricyclePosition","clusters")
basal_meta$clusters <- paste0("basal_",basal_meta$clusters)
basal_meta$clusters <- factor(basal_meta$clusters,levels=c("basal_1","basal_5","basal_6"))

mccs_meta=mccs_clean@meta.data[,c("signature_1_UCell","signature_2_UCell","orig.ident","group","geno","tricyclePosition","integrated_snn_res.0.3")]

colnames(mccs_meta)=c("signature_1_UCell","signature_2_UCell","orig.ident","group","geno","tricyclePosition","clusters")
mccs_meta$clusters <- paste0("MCC_",mccs_meta$clusters)
mccs_meta$clusters<-factor(mccs_meta$clusters,levels = c("MCC_0","MCC_4","MCC_2","MCC_3","MCC_1"))

total_dat=rbind(basal_meta,mccs_meta)
total_dat$geno<-factor(total_dat$geno,levels=c("WT","Hom"))

pdf("S_Score_basal_MCC_clusters.pdf",height=4,width=8,useDingbats = F)
ggplot(total_dat, aes(x=clusters,y=signature_1_UCell,color=geno,fill=clusters,pattern=geno))+
        geom_boxplot_pattern(color="black",
                             pattern_fill = "black",
                             pattern_angle = 45,
                             pattern_density = 0.1,
                             pattern_spacing = 0.025,
                             pattern_key_scale_factor = 0.6)+
        scale_pattern_manual(values = c(WT = "none", Hom = "stripe")) +
        theme_classic()+
        scale_color_manual(values=c("black","black"))+
        scale_fill_manual(values=c(colors[c(2,6,7,1)],mcc_color(6)[c(1,2,4,6)]))+
        xlab("Clusters")+
        ylab("S Score")+
        theme(axis.text = element_text(size=14,color="black"),axis.title = element_text(size=16))
dev.off()

pdf("G2M_Score_basal_MCC_clusters.pdf",height=4,width=8,useDingbats = F)
ggplot(total_dat, aes(x=clusters,y=signature_2_UCell,color=geno,fill=clusters,pattern=geno))+
        geom_boxplot_pattern(color="black",
                             pattern_fill = "black",
                             pattern_angle = 45,
                             pattern_density = 0.1,
                             pattern_spacing = 0.025,
                             pattern_key_scale_factor = 0.6)+
        scale_pattern_manual(values = c(WT = "none", Hom = "stripe")) +
        theme_classic()+
        scale_color_manual(values=c("black","black"))+
        scale_fill_manual(values=c(colors[c(2,6,7,1)],mcc_color(6)[c(1,2,4,6)]))+
        xlab("Clusters")+
        ylab("G2M Score")+
        theme(axis.text = element_text(size=14,color="black"),axis.title = element_text(size=16))
dev.off()

### print out rawdata

colnames(total_dat) <- c("S Score","G2/M Score","Replicate","Celltype","Genotype","tricyclePosition","Cluster Name")
total_dat$Cluster <- "Basal stem"
total_dat$Cluster[which(total_dat$`Cluster Name` == "basal_5")] <- "Proliferating basal stem (S)"
total_dat$Cluster[which(total_dat$`Cluster Name` == "basal_6")] <- "Proliferating basal stem (G2/M)"
total_dat$Cluster[which(total_dat$`Cluster Name` == "MCC_0")] <- "A"
total_dat$Cluster[which(total_dat$`Cluster Name` == "MCC_4")] <- "B"
total_dat$Cluster[which(total_dat$`Cluster Name` == "MCC_2")] <- "C"
total_dat$Cluster[which(total_dat$`Cluster Name` == "MCC_3")] <- "D"
total_dat$Cluster[which(total_dat$`Cluster Name` == "MCC_1")] <- "E"

write.csv(total_dat[,c(1:5,8)],file="G2M_S_Score_basal_MCC_rawdata.csv")

mccs_clean$geno<-factor(mccs_clean$geno,levels=c("WT","Hom"))

mc1=ggplot(mccs_clean@meta.data %>% filter(!is.na(mccs_clean$pseudotime)), aes(x=pseudotime,y=signature_1_UCell,linetype=geno))+
               geom_smooth(color="black")+
        theme_classic()+
        xlab("Pseudotime")+
        ylab("S Score")+
        theme(axis.text=element_text(size=13,color="black"),axis.title = element_text(size=14),legend.position = "none")

mc2=ggplot(mccs_clean@meta.data %>% filter(!is.na(mccs_clean$pseudotime)), aes(x=pseudotime,y=signature_2_UCell,linetype=geno))+
        geom_smooth(color="black")+
        theme_classic()+
        xlab("Pseudotime")+
        ylab("G2M Score")+
        theme(axis.text=element_text(size=13,color="black"),axis.title = element_text(size=14),legend.position = "none")


dots=ggplot(mccs_clean@meta.data %>% filter(!is.na(mccs_clean$pseudotime)), aes(x=pseudotime,y=1.00,color=integrated_snn_res.0.3))+
        geom_jitter(width=0.005)+
        scale_color_manual(values=c(mcc_clean_col))+
        theme_classic()+
        scale_y_continuous(
                labels = scales::number_format(accuracy = 0.01))+
        xlab("Pseudotime")+
        theme(axis.text=element_text(size=13,color="black"),axis.title = element_text(size=14),legend.position = "none")

pdf("mccs_clean/mccs_clean_pseudo_vs_S_or_G2M_score.pdf",height=4,width=4.5,useDingbats = F)
ggarrange(mc1,mc2,ncol=1)
dev.off()

pdf("mccs_clean/mccs_clean_pseudo_vs_S_or_G2M_score_dots.pdf",height=1,width=4.5,useDingbats = F)
dots
dev.off()

########### Cluster S and G2M Scores #############
basal_cycling$timecourse_pred <- factor(basal_cycling$timecourse_pred,levels=c(1,5,6))
ggplot(basal_cycling@meta.data, aes(x=timecourse_pred,y=signature_1_UCell,color=geno))+
        geom_boxplot()


means_S_G2M=mccs_clean@meta.data %>% group_by(orig.ident,integrated_snn_res.0.3) %>%
        summarise(mean_sig1=mean(signature_1_UCell),
                  mean_sig2= mean(signature_2_UCell),
                  .groups = 'drop') 

write.csv(means_S_G2M, file="mccs_clean/E2f7_em3_means_S_G2M_multiciliated_clusters.csv")

means_basal_S_G2M=basal_cycling@meta.data %>% group_by(orig.ident,timecourse_pred) %>% 
        summarise(mean_sig1=mean(signature_1_UCell),
                  mean_sig2= mean(signature_2_UCell),
                  .groups = 'drop') 

write.csv(means_basal_S_G2M, file="cycling_basal/E2f7_em3_means_S_G2M_basal_clusters.csv")


############# PROPORTION TESTING ############
library(speckle)
mccs_clean=readRDS("mccs_clean/E2f7_em3_mccs_clean.rds")
basal_cycling=readRDS("cycling_basal/E2f7_em3_cycling_basal_subset.rds")
DimPlot(mccs_clean, label=T,group.by="integrated_snn_res.0.3")

mccs_clean$geno<-"WT"
mccs_clean$geno[grep("Hom*",mccs_clean$orig.ident)]<-"Hom"
basal_cycling$geno <-"WT"
basal_cycling$geno[grep("Hom*",basal_cycling$orig.ident)]<-"Hom"

mccs_clean$clusters <- paste0("MCC_",mccs_clean$integrated_snn_res.0.3)
basal_cycling$clusters <- paste0("Basal_",basal_cycling$timecourse_pred)
DimPlot(mccs_clean, group.by="clusters",split.by="geno")
DimPlot(basal_cycling,group.by="clusters" ,split.by="geno")

dat=mccs_clean@meta.data[,c("clusters","orig.ident","geno")]
dat$geno<-factor(dat$geno,levels=c("WT","Hom"))

prop_test=propeller(clusters =dat$clusters, sample = dat$orig.ident, 
                    group =dat$geno)
prop_test_melt =reshape2::melt(prop_test)

prop_test_melt$BaselineProp.clusters<-factor(prop_test_melt$BaselineProp.clusters,levels = rev(c("MCC_0","MCC_4","MCC_2","MCC_3","MCC_1")))

##rename the cluster IDs to print out the raw data ##
rawdata=prop_test_melt[c(grep("PropMean*",prop_test_melt$variable),grep("FDR",prop_test_melt$variable)),]
rawdata$clustername <- "A"
rawdata$clustername[grep("MCC_4",rawdata$BaselineProp.clusters)] = "B"
rawdata$clustername[grep("MCC_2",rawdata$BaselineProp.clusters)] = "C"
rawdata$clustername[grep("MCC_3",rawdata$BaselineProp.clusters)] = "D"
rawdata$clustername[grep("MCC_1",rawdata$BaselineProp.clusters)] = "E"

final_raw=cbind(rawdata$clustername,rawdata[,c(2:3)])
colnames(final_raw )<-c("Cluster", "Measurement","Value")

write.csv(final_raw, file="mccs_clean/E2f7_em3_mccs_clean_proportion_rawdata.csv")

pdf("mccs_clean/E2f7_em3_mccs_clean_proportions.pdf",height=6,width=6,useDingbats = F)
ggplot(prop_test_melt[grep("PropMean*",prop_test_melt$variable),], aes(x=variable,y=value,fill=BaselineProp.clusters))+
        geom_col()+
        theme_classic()+
        scale_fill_manual(values=rev(c(colors[1],mcc_color(6)[c(1,2,4,6)])))+
        xlab("")+
        ylab("Cell type proportion")+
        theme(axis.text = element_text(size=14,color="black"),axis.title=element_text(size=16))
dev.off()

dat2=basal_cycling@meta.data[,c("clusters","orig.ident","geno")]
dat2$geno<-factor(dat2$geno,levels=c("WT","Hom"))

prop_test2=propeller(clusters =dat2$clusters, sample = dat2$orig.ident, 
                    group =dat2$geno)
prop_test_melt2 =reshape2::melt(prop_test2)

##rename the cluster IDs to print out the raw data ##
rawdata2=prop_test_melt2[c(grep("PropMean*",prop_test_melt2$variable),grep("FDR",prop_test_melt2$variable)),]
rawdata2$clustername <- "Basal Stem"
rawdata2$clustername[grep("Basal_5",rawdata2$BaselineProp.clusters)] = "Proliferating Basal Stem (S)"
rawdata2$clustername[grep("Basal_6",rawdata2$BaselineProp.clusters)] = "Proliferating Basal Stem (G2/M)"


final_raw2=cbind(rawdata2$clustername,rawdata2[,c(2:3)])
colnames(final_raw2 )<-c("Cluster", "Measurement","Value")

write.csv(final_raw2, file="cycling_basal/E2f7_em3_cycling_basal_proportions_rawdata.csv")

prop_test_melt2$BaselineProp.clusters<-factor(prop_test_melt2$BaselineProp.clusters,levels = c("Basal_6","Basal_5","Basal_1"))

pdf("cycling_basal//E2f7_em3_cycling_basal_proportions.pdf",height=6,width=6,useDingbats = F)
ggplot(prop_test_melt2[grep("PropMean*",prop_test_melt2$variable),], aes(x=variable,y=value,fill=BaselineProp.clusters))+
        geom_col()+
        theme_classic()+
        scale_fill_manual(values=rev(colors[c(2,6,7)]))+
        xlab("")+
        ylab("Cell type proportion")+
        theme(axis.text = element_text(size=12,color="black"),axis.title=element_text(size=14))
dev.off()

########### Plot focused datasets back onto original UMAP########

seur$group <- NA
seur$group[match(rownames(mccs_clean@meta.data),rownames(seur@meta.data))]<-"MCC"
seur$group[match(rownames(basal_cycling@meta.data),rownames(seur@meta.data))]<-"Basal"

pdf("E2f7_em3_v2_highlight_basal_MCCs.pdf",height=8,width=11,useDingbats = F)
DimPlot(seur, group.by="group",pt.size=1,cols=c(colors[2],mcc_color(6)[3]),na.value = "grey30")+theme_void()+theme(title=element_blank(),legend.position = "none")+scale_x_reverse()
dev.off()

m=DimPlot(seur, group.by="group",pt.size=1,cols=c(colors[2],mcc_color(6)[3]),na.value = "grey30")+theme_void()+theme(title=element_blank())+scale_x_reverse()

leg=get_legend(m)

pdf("E2f7_em3_v2_highlight_basal_MCCs_legend.pdf",height=2,width=4,useDingbats = F)
grid.newpage()
grid.draw(leg)
dev.off()

################# GENE EXPRESSION OVER PSEUDOTIME ###########################
DefaultAssay(mccs_clean) <- "RNA"
mccs_clean$geno <- "WT"
mccs_clean$geno[grep("Hom*",mccs_clean$orig.ident)]<-"Hom"

##print as 1 plot
gene_dat=FetchData(mccs_clean, c("Mcm5","Gins2"))

## use this list for loop
gene_dat=FetchData(mccs_clean, c("Espl1","Pttg1","Ccnb1","Rad21","Cdk1","Ccne1","Ccna1", "Ccna2","Nacpd2","Mycl","Deup1","Plk1","Plk4","Ift88","Ift20","Dnah5","Dnah9","Ccp110","Cep164", "Cep83","Stil","Foxj1","Foxn4", "Myb","Mcidas","Gmnc","Rfx2","Rfx3","Cdkn1a","Ccno","E2f7","E2f8","Kif23",
"Racgap1","Ect2", "Prc1"))

##more looping
cilia=read.csv("~/Documents/Reiter_Seq/CuratedCiliaAndBasalBodyGenes_18.csv")
gene_list=cilia$Gene.name..Mus.musculus.
gene_dat=FetchData(mccs_clean, gene_list)

gene_dat_2=cbind(gene_dat, mccs_clean@meta.data[,c("geno","integrated_snn_res.0.3","pseudotime")])

gene_dat_filt=gene_dat_2 %>% filter(!is.na(pseudotime))
melted=reshape2::melt(gene_dat_filt,id.vars=c("pseudotime","geno","integrated_snn_res.0.3"))

melted$geno<-factor(melted$geno,levels=c("WT","Hom"))

pdf("mccs_clean/mccs_clean_pseudotime_vs_gene_lineplots.pdf",height=8,width=6,useDingbats = F)
ggplot(melted, aes(x=pseudotime,y=value,linetype=geno,color=geno))+
        geom_smooth()+
        scale_color_manual(values=c("grey47",mcc_color(5)[3]))+
        facet_wrap(~variable,scales = "free",ncol = 1)+
        theme_classic()+
        xlab("Pseudotime")+
        ylab("Expression")+
        theme(axis.text=element_text(size=10,color="black"),axis.title=element_text(size=12),strip.background = element_blank())
        #geom_jitter(aes(x=pseudotime, y=1, color=integrated_snn_res.0.3),height=0.05)
dev.off()


### loop to print many of these plots but as single genes
dir.create("mccs_clean/line_plots")
for(x in unique(melted$variable)){
        pdf(paste0("mccs_clean/line_plots/mccs_clean_pseudotime_vs_",x,"_lineplots.pdf"),height=4,width=6,useDingbats = F)
        print(ggplot(melted %>% filter(variable == x), aes(x=pseudotime,y=value,linetype=geno,color=geno,fill=geno))+
                geom_smooth()+
                facet_wrap(~variable,scales = "free")+
                theme_classic()+
                xlab("Pseudotime")+
                ylab("Expression")+
                scale_color_manual(values=c("grey30","#1b75ba"))+
                scale_fill_manual(values=c("grey30","#1b75ba"))+
                theme(axis.text=element_text(size=10,color="black"),axis.title=element_text(size=12),strip.background = element_blank()))
        dev.off()
}

gene_dat_3=cbind(gene_dat, mccs_clean@meta.data[,c("orig.ident","integrated_snn_res.0.3","pseudotime")])

gene_dat_filt_3=gene_dat_3 %>% filter(!is.na(pseudotime))
melted=reshape2::melt(gene_dat_filt_3,id.vars=c("pseudotime","orig.ident","integrated_snn_res.0.3"))


### loop to print many of these plots but as single genes

for(x in unique(melted$variable)){
        pdf(paste0("mccs_clean/line_plots/mccs_clean_pseudotime_vs_",x,"_lineplots_replicates.pdf"),height=4,width=6,useDingbats = F)
        print(ggplot(melted %>% filter(variable == x), aes(x=pseudotime,y=value,color=orig.ident,fill=orig.ident))+
                      geom_smooth()+
                      facet_wrap(~variable,scales = "free")+
                      theme_classic()+
                      xlab("Pseudotime")+
                      ylab("Expression")+
                      scale_color_manual(values=c("#1b75ba","#1b75ba","grey30","grey30"))+
                      scale_fill_manual(values=c("#1b75ba","#1b75ba","grey30","grey30"))+
                      theme(axis.text=element_text(size=10,color="black"),axis.title=element_text(size=12),strip.background = element_blank()))
        dev.off()
}

############## Ciliary TF, Basal body ,and ciliogenesis scores #########
library(UCell)
library(dplyr)
library(ggpubr)
library(scales)

### also info in Supplementary Table 4, printed below
cilia=read.csv("~/Documents/Reiter_Seq/CuratedCiliaAndBasalBodyGenes_18.csv")

cilia$Merged_Function2[cilia$Molecular.Function.or.Complex == "Ciliary TF"]="Ciliary TF"

ciliogen=cilia %>% filter(Merged_Function2 %in% c("Ciliogenesis","Cilia Function/Motility"))
basalbody=cilia %>% filter(Merged_Function2 %in% c("Basal body maturation","Centriole amplification"))

tfs=cilia %>% filter(Molecular.Function.or.Complex %in% c("Ciliary TF"))

seur=AddModuleScore_UCell(seur, features = list(ciliogen=ciliogen$Gene.name..Mus.musculus.,basalbody=basalbody$Gene.name..Mus.musculus.,tfs=tfs$Gene.name..Mus.musculus.),assay="RNA")

seur_mccs = seur@meta.data[match(rownames(mccs_clean@meta.data),rownames(seur@meta.data)),]
identical(rownames(seur_mccs),rownames(mccs_clean@meta.data))

mccs_clean$ciliogen_score <- seur_mccs$ciliogen
mccs_clean$basalbody_score <- seur_mccs$basalbody
mccs_clean$tfs_score <- seur_mccs$tfs_UCell

mccs_clean$geno <- "WT"
mccs_clean$geno[grep("Hom",mccs_clean$orig.ident)]<-"Hom"
mccs_clean$geno <- factor(mccs_clean$geno,levels=c("WT","Hom"))

FeaturePlot(mccs_clean, "ciliogen_score",split.by="geno")

p1=ggplot(mccs_clean@meta.data%>% filter(!is.na(mccs_clean$pseudotime)), aes(x=pseudotime,y=ciliogen_score,linetype=geno))+
        geom_smooth(color="black")+
        theme_classic()+
        xlab("Pseudotime")+
        ylab("Ciliogenesis score")+
        scale_y_continuous(labels = label_number(accuracy = 0.1))+
        theme(axis.text=element_text(size=13,color="black"),axis.title = element_text(size=14),legend.position = "none")

p2=ggplot(mccs_clean@meta.data%>% filter(!is.na(mccs_clean$pseudotime)), aes(x=pseudotime,y=basalbody_score,linetype=geno))+
        geom_smooth(color="black")+
        theme_classic()+
        xlab("Pseudotime")+
        ylab("Basal body score")+
        scale_y_continuous(labels = label_number(accuracy = 0.1))+
        theme(axis.text=element_text(size=13,color="black"),axis.title = element_text(size=14),legend.position = "none")

p3=ggplot(mccs_clean@meta.data%>% filter(!is.na(mccs_clean$pseudotime)), aes(x=pseudotime,y=tfs_score,linetype=geno))+
        geom_smooth(color="black")+
        theme_classic()+
        xlab("Pseudotime")+
        ylab("Ciliary TF score")+
        scale_y_continuous(labels = label_number(accuracy = 0.1))+
        theme(axis.text=element_text(size=13,color="black"),axis.title = element_text(size=14),legend.position = "none")

pdf("mccs_clean/mccs_clean_pseudo_vs_ciliogenesis_or_basalbody_or_tfs_score.pdf",height=5,width=4.5,useDingbats = F)
ggarrange(p1,p2,p3,ncol=1)
dev.off()

mccs_pseudo=mccs_clean@meta.data%>% filter(!is.na(mccs_clean$pseudotime))

v1=VlnPlot(subset(mccs_clean, cells=rownames(mccs_pseudo)), features="ciliogen_score",group.by="geno",pt.size = 0,cols=c("grey60","skyblue1"))+ggtitle("")+ylab("Ciliogenesis Score")+theme(axis.title.x = element_blank(),legend.position = "none",axis.text.x = element_text(angle=0,hjust = 0.5))

v2=VlnPlot(subset(mccs_clean, cells=rownames(mccs_pseudo)), features="basalbody_score",group.by="geno",pt.size = 0,cols=c("grey60","skyblue1"),)+ggtitle("")+ylab("Basal body Score")+theme(axis.title.x = element_blank(),legend.position = "none",axis.text.x = element_text(angle=0,hjust = 0.5))

v3=VlnPlot(subset(mccs_clean, cells=rownames(mccs_pseudo)), features="tfs_score",group.by="geno",pt.size = 0,cols=c("grey60","skyblue1"),)+ggtitle("")+ylab("Ciliary TF score")+theme(axis.title.x = element_blank(),legend.position = "none",axis.text.x = element_text(angle=0,hjust = 0.5))

pdf("mccs_clean/mccs_clean_ciliogenesis_basalbody_tf_score_violin.pdf",height=8,width=5,useDingbats = F)
ggarrange(v1,v2,v3,ncol = 1)
dev.off()

means_scores=mccs_clean@meta.data %>%  group_by(orig.ident) %>%
        summarise(mean_ciliogen=mean(ciliogen_score),
                  mean_basalbody_score= mean(basalbody_score),
                  .groups = 'drop') 

write.csv(means_scores,file="mccs_clean/mccs_clean_ciliogen_basalbody_means.csv")

means_scores=mccs_pseudo %>% group_by(orig.ident) %>%
        summarise(mean_ciliogen=mean(ciliogen_score),
                  mean_basalbody_score= mean(basalbody_score),
                  mean_tfs_score=mean(tfs_score),
                  .groups = 'drop') 

write.csv(means_scores,file="mccs_clean/mccs_clean_ciliogen_basalbody_tfs_means_only_pseudotime.csv")

##write out supp table
supp=cilia %>% select(Gene.Name..Homo.sapiens.,Gene.name..Mus.musculus.,Synonyms,Merged_Function2)
colnames(supp) <- c("Homo_sapiens","Mus_Musculus","Synonyms","Function")
write.csv(supp, file="~/Box Sync/E2f7_paper/Supplementary_Table4.csv")

##############################   HEATMAP for DNA Synthesis and Cytoskeletal Gene Expression   ###################################### 
library(dplyr)
dna_syn=read.csv("~/Box Sync/E2f7_paper/DNA Synthesis Genes 3.csv")
deseq=read.csv("~/Box Sync/E2f7_paper/E2f7_DESeq_padj0.05_FC1.5_new_oldCUTRUN_WESTENDORP_genelist.csv")

dna_syn_in = dna_syn %>% filter(dna_syn$Gene %in% deseq$Column1)

deseq$DNA_Synthesis <- NA
deseq$DNA_Synthesis[match(dna_syn_in$Gene,deseq$Column1)]<- dna_syn_in$Function.L1
deseq_dna = deseq %>% filter(!is.na(DNA_Synthesis))
genes_to_test=c(deseq_dna$Column1)

mccs_clean$group <- "WT"
mccs_clean$group[grep("Hom*",mccs_clean$orig.ident)]<-"Hom"

mccs_clean$group_pseudo_bin <- paste0(mccs_clean$group,mccs_clean$pseudo_bin)
pseudo_cells=mccs_clean@meta.data[!is.na(mccs_clean$pseudotime),]

mat=AverageExpression(subset(mccs_clean,cells = rownames(pseudo_cells)), group.by = "group_pseudo_bin",features = genes_to_test,assay="RNA")
library(UCell)

cc_s=read.csv("~/Documents/Reiter_Seq/s_phase_human_mgi_genes.csv")
cc_g2m=read.csv("~/Documents/Reiter_Seq/g2m_phase_human_mgi_genes.csv")

sgenes=cc_s$MGI.symbol
g2mgenes=cc_g2m$MGI.symbol

mccs_clean = AddModuleScore_UCell(mccs_clean,features = list(sgenes,g2mgenes),assay = "RNA") 

library(pheatmap)

col_order = c(paste("WT",c(1:20),sep=""),paste("Hom",c(1:20),sep=""))
mat_order=mat$RNA[,match(col_order,colnames(mat$RNA))]


##normalize min max
library(scales)
cc_scores=mccs_clean@meta.data %>% group_by(group_pseudo_bin) %>%
        summarise_at(vars(c(signature_1_UCell,signature_2_UCell)),              
                     list(name = mean)) 
cc_scores=cc_scores[match(col_order,cc_scores$group_pseudo_bin),]
cc_scores$s_minmax =rescale(cc_scores$signature_1_UCell_name)
cc_scores$g2m_minmax =rescale(cc_scores$signature_2_UCell_name)
cc_scores=as.data.frame(cc_scores)
rownames(cc_scores) <- cc_scores$group_pseudo_bin
cc_scores$Pseudo_bin <- c(1:20,1:20)
cc_scores$Dataset <- "E2f7 WT"
cc_scores$Dataset[21:40]<-"E2f7 Hom"

dna_syn_in=dna_syn[dna_syn$Gene %in% rownames(mat_order),]
row_anno = as.data.frame(rownames(mat_order))
rownames(row_anno)<-row_anno$`rownames(mat_order)`
row_anno$Function <- NA
row_anno$Function[match(dna_syn_in$Gene,rownames(row_anno))] <- dna_syn_in$Function.L1

col_anno=cc_scores[,c("Dataset","Pseudo_bin","s_minmax","g2m_minmax")]
colnames(col_anno) <- c("Dataset","Pseudotime Bin", "S score", "G2M score")

vir_col=viridis_pal()

colors_anno=rev(vir_col(20))
names(colors_anno) <- unique(col_anno$`Pseudotime Bin`)
colors_anno = list(`Pseudotime Bin`=colors_anno)

color_fnc=colorRampPalette(c("#CCFFFF","#330066"))
color_fnc2=colorRampPalette(c("#CCFFFF","#006600"))
col_scores=color_fnc(20)

colors_anno[["S score"]]<-col_scores
colors_anno[["G2M score"]]<-color_fnc2(20)

anno_row_color=c("palegreen"  ,"lightskyblue3",   "lavenderblush2",  "royalblue1" ,     "olivedrab4"  ,   "deeppink2"  ,     "darkolivegreen1" ,"seagreen3"   , "blue3"  ,   "goldenrod4" ,"mediumpurple3"  , "tomato" ,  "lightsteelblue1" ,"turquoise1"  ,  "deeppink" , "paleturquoise3" , "brown1" )[1:13]

names(anno_row_color) <- unique(dna_syn_in$Function.L1)

colors_anno[["Function"]]<-anno_row_color

data_colors =c("#666666","#6666CC")
names(data_colors) <-c("E2f7 WT","E2f7 Hom")

colors_anno[["Dataset"]] <- data_colors

pdf("~/Box Sync/E2f7_paper/E2f7_WT_vs_Hom_DE_DNAsyngenes_heatmap.pdf",height=8,width=11,useDingbats = F)
p=pheatmap(mat_order,cluster_cols = F,scale="row",cutree_rows = 5,gaps_col = 20,annotation_row = row_anno[,"Function",drop=F],annotation_col = col_anno,annotation_colors = colors_anno,show_colnames = F)
dev.off()

### group by functional category

row_anno_order = row_anno[order(row_anno$Function),]

mat_order_fnc = mat_order[match(rownames(row_anno_order),rownames(mat_order)),]

anno_row_color=c("palegreen"  ,"lightskyblue3",   "lavenderblush2",  "royalblue1" ,     "olivedrab4"  ,   "deeppink2"  ,     "darkolivegreen1" ,"seagreen3"   , "blue3"  ,   "goldenrod4" ,"mediumpurple3"  , "tomato" ,  "lightsteelblue1" ,"turquoise1"  ,  "deeppink" , "paleturquoise3" , "brown1" )[1:7]

names(anno_row_color) <- unique(row_anno_order$Function)

colors_anno[["Function"]]<-anno_row_color

pdf("~/Box Sync/E2f7_paper/E2f7_WT_vs_Hom_DE_DNAsyngenes_heatmap_orderfnccat.pdf",height=8,width=11,useDingbats = F)
p=pheatmap(mat_order_fnc,cluster_cols = F,cluster_rows = F,scale="row",cutree_rows = 5,gaps_col = 20,gaps_row=c(),annotation_row = row_anno_order[,"Function",drop=F],annotation_col = col_anno,annotation_colors = colors_anno,show_colnames = F)
dev.off()

################### Cytoskeletal Heatmap ###################
### cytoskeleton

cyto=read.csv("~/Downloads/gene sets for cytoskeleton.gmt",sep="\t",header = F)
cyto_list=cyto[,3:length(cyto[1,])]

cyto_vec=NULL
for(i in 1:length(cyto_list$V3)){
        cyto_i=unlist((cyto_list[i,,drop=T]))
        cyto_vec=c(cyto_vec, cyto_i)
}

cyto_un=cyto_vec=unique(cyto_vec)
cyto_final=cyto_un[!cyto_un == ""]

write.csv(cyto_final, file="cytoskeleton_GO.csv")
## use gprofiler to add human IDs

cyto_mouse = read.csv("cytoskeleton_GO.csv")
mat=AverageExpression(mccs_clean, group.by = "group_pseudo_bin",features =cyto_mouse$ortholog_name,assay="RNA")

dat_melt  =reshape2::melt(mat)

threshed_genes=NULL
thresh=0.5
for(x in unique(dat_melt$Var1)){
        dat_x=dat_melt %>% filter(Var1 == x)
        bs_dat=dat_x[grep("WT*",dat_x$Var2),]
        bs_sum= sum( bs_dat$value > thresh)
        mc_dat=dat_x[grep("Hom*",dat_x$Var2),]
        mc_sum= sum( mc_dat$value > thresh)
        if(bs_sum | mc_sum > 1){
                threshed_genes <- c( threshed_genes,x)
        }
}
mat=AverageExpression(mccs_clean, group.by = "group_pseudo_bin",features =threshed_genes,assay="RNA")

col_order = c(paste("WT",c(1:20),sep=""),paste("Hom",c(1:20),sep=""))
mat_order=mat$RNA[,match(col_order,colnames(mat$RNA))]

intersect(deseq$Column1,rownames(mat_order))

p=pheatmap(mat_order,cluster_cols = F,scale="row",cutree_rows = 4,gaps_col = 20)

## deseq only in cluster b and c

deseq_bc=deseq %>% filter(abs(clus_b_log2FoldChange)> log2(1.5) | abs(clus_c_log2FoldChange) > log2(1.5) )

deseq_bc_in=deseq_bc[deseq_bc$Column1 %in% rownames(mat_order),]
deseq_in=deseq[deseq$Column1 %in% rownames(mat_order),]

row_anno = as.data.frame(rownames(mat_order))
rownames(row_anno)<-row_anno$`rownames(mat_order)`
row_anno$Deseq <- ""
row_anno$Deseq[match(deseq_in$Column1,rownames(row_anno))] <- "DE"
row_anno$Deseq[match(deseq_bc_in$Column1,rownames(row_anno))] <- "DE cluster b or c"

de_colors=c("white","red","blue")
names(de_colors)<- c("","DE","DE cluster b or c")
colors_anno[["Deseq"]]<-de_colors

pdf("~/Box Sync/E2f7_paper/E2f7_WT_vs_Hom_cytoskeleton_heatmap.pdf",height=30,width=11,useDingbats = F)
pheatmap(mat_order,cluster_cols=F,scale="row",cutree_rows = 9,gaps_col=20,annotation_row = row_anno[,"Deseq",drop=F],annotation_col=col_anno,annotation_colors = colors_anno,show_colnames = F)
dev.off()

###################### Make combined heatmap for DNA Synthesis and Cytoskeletal Genes #######
## Fig. S13a
## use plot from DNA Synthesis heatmap and add on DE cytoskeletal genes from above
row_anno_order ## this is order of genes in prior plot

cyto = c("Diaph3","Kif23","Stmn1","Tubg1","Incenp","Nrg1")
cyto_df=as.data.frame(cyto)
cyto_df$Function <- "Cyto"

colnames(cyto_df) <- colnames(row_anno_order)

all_dat = rbind(row_anno_order, cyto_df)
colnames(all_dat) <-c("Gene","Function")

### resuse heatmap code from DNA synthesis seciton but add cyto category
mat2=AverageExpression(subset(mccs_clean,cells = rownames(pseudo_cells)), group.by = "group_pseudo_bin",features = all_dat$Gene,assay="RNA")
col_order = c(paste("WT",c(1:20),sep=""),paste("Hom",c(1:20),sep=""))
mat_order2=mat2$RNA[,match(col_order,colnames(mat2$RNA))]

gene_order = c("E2f1","E2f3","Ccne1","Ccne2","Orc1","Cdt1","Mcm2","Mcm3","Mcm4","Mcm5","Mcm6","Mcm7","Mcm10","Cdc45","Gins2","Gins4","Pola1","Pola2","Prim2","Rpa1","Rpa2","Pold1","Pold2","Pold3","Pole4","Pcna","Rfc1","Rfc2","Rfc3","Rfc4","Lig1","Rnaseh2a","Rnaseh2b","Slbp","Rbbp4","Nap1l1","Chaf1a","Chaf1b","Diaph3","Kif23","Stmn1","Tubg1","Incenp","Nrg1")
gaps= c(4,12,16,21,30,33,38)

new_row_anno = as.data.frame(rownames(mat_order2))
colnames(new_row_anno) <- "Gene"
new_row_anno$Function <- "Cyto"
new_row_anno$Function[match(rownames(row_anno),new_row_anno$Gene)] <- row_anno$Function

new_row_anno_order=new_row_anno[match(gene_order, new_row_anno$Gene),]
rownames(new_row_anno_order) <- new_row_anno_order$Gene

mat_order_fnc2 = mat_order2[match(gene_order,rownames(mat_order2)),]

anno_row_color2=c("palegreen"  ,"lightskyblue3",   "lavenderblush2",  "royalblue1" ,     "olivedrab4"  ,   "deeppink2"  ,     "darkolivegreen1" ,"seagreen3"   , "blue3"  ,   "goldenrod4" ,"mediumpurple3"  , "tomato" ,  "lightsteelblue1" ,"turquoise1"  ,  "deeppink" , "paleturquoise3" , "brown1" )[1:8]

names(anno_row_color2) <- unique(new_row_anno_order$Function)

colors_anno[["Function"]]<-anno_row_color2

## set scale to prevent outliers from dominating scale

pdf("~/Box Sync/E2f7_paper/E2f7_WT_vs_Hom_DE_DNAsyngenes_cyto_heatmap_ordered.pdf",height=8,width=7,useDingbats = F)
p=pheatmap(mat_order_fnc2,cluster_cols = F,cluster_rows = F,scale="row",cutree_rows = 5,gaps_col = 20,gaps_row=gaps,annotation_row = new_row_anno_order[,"Function",drop=F],annotation_col = col_anno,annotation_colors = colors_anno,show_colnames = F,legend_breaks = c(-3.5,-1.5,0,1.5,3.5))
dev.off()

####### print count and metadata for GEO ########

seur=readRDS("E2f7_em3_merge/v2/E2f7_em3_merge_v2.rds")
DimPlot(seur, label=T, group.by="timecourse_pred")

mccs_clean=readRDS("E2f7_em3_merge/v2/mccs_clean/E2f7_em3_mccs_clean.rds")
basal_cycling=readRDS("E2f7_em3_merge/v2/cycling_basal/E2f7_em3_cycling_basal_subset.rds")

### need MCC or Basal cell subset, pseudotime values for MCCs, S and G2M scores
library(dplyr)
seur$subset_group <- NA
seur$subset_group[match(rownames(mccs_clean@meta.data),rownames(seur@meta.data))]<-"Multiciliated"
seur$subset_group[match(rownames(basal_cycling@meta.data),rownames(seur@meta.data))]<-"Basal"

DimPlot(seur, group.by="subset_group")
meta=seur@meta.data
meta$cluster_names <- meta$timecourse_pred

DimPlot(seur, group.by="timecourse_pred",label=T)+scale_y_reverse()
meta=meta %>% mutate(cluster_names=recode(cluster_names,
                  `1`="Basal stem",
                  `2`="Basal to intermediate",
                  `3`="Secretory",
                  `4`="Multiciliated 3",
                  `5`="S",
                  `6`="G2M",
                  `7`="Tnfrsf12a+",
                  `8`="Multiciliated 1",
                  `9`="Multiciliated 2",
                  `0`="Intermediate"
                  ))

seur@meta.data=meta
DimPlot(seur, group.by="cluster_names",label=T)

Idents(seur) <- "cluster_names"
DefaultAssay(seur) <- "RNA"

marks=FindAllMarkers(seur, assay = "RNA",only.pos = T)
write.csv(marks, file="E2f7_merge_v2_markers.csv")
write.csv(marks, file="~/Box Sync/E2f7_paper/Supplementary_Table3.csv")

## add Multiciliated pseidotime values and clusters
seur$multiciliated_pseudotime <- NA
seur$multiciliated_pseudotime[match(rownames(mccs_clean@meta.data),rownames(seur@meta.data))]<-mccs_clean$pseudotime

FeaturePlot(seur, "multiciliated_pseudotime")

seur$multiciliated_clusters <- NA
seur$multiciliated_clusters[match(rownames(mccs_clean@meta.data),rownames(seur@meta.data))]<- as.character(mccs_clean$integrated_snn_res.0.3)
seur$multiciliated_cluster_names <- seur$multiciliated_clusters
seur@meta.data=seur@meta.data %>% mutate(multiciliated_cluster_names=recode(multiciliated_cluster_names,
                                                             `0`="Intermediate",
                                                             `1`="Multiciliated 4",
                                                             `2`="Multciliated 2",
                                                             `3`="Multciliated 3",
                                                             `4`="Multiciliated 1"))
## add S and G2M scores calculated by UCell
seur$basal_S_score <- NA
seur$basal_S_score[match(rownames(basal_cycling@meta.data),rownames(seur@meta.data))]<-basal_cycling$signature_1_UCell
seur$basal_G2M_score <- NA
seur$basal_G2M_score[match(rownames(basal_cycling@meta.data),rownames(seur@meta.data))]<-basal_cycling$signature_2_UCell

FeaturePlot(seur, "basal_S_score")

seur$multiciliated_S_score <- NA
seur$multiciliated_S_score[match(rownames(mccs_clean@meta.data),rownames(seur@meta.data))]<-mccs_clean$signature_1_UCell
seur$multiciliated_G2M_score <- NA
seur$multiciliated_G2M_score[match(rownames(mccs_clean@meta.data),rownames(seur@meta.data))]<-mccs_clean$signature_2_UCell

FeaturePlot(seur, "multiciliated_S_score",split.by="orig.ident")

#### Add UMAP for original, subsetted data

full_umap =Embeddings(seur, "umap")
mcc_umap=Embeddings(mccs_clean,"umap")
basal_umap=Embeddings(basal_cycling,"umap")

seur$UMAP_1<-full_umap[,c("UMAP_1")]
seur$UMAP_2<-full_umap[,c("UMAP_2")]

seur$multiciliated_UMAP_1 <-NA
seur$multiciliated_UMAP_1[match(rownames(mcc_umap),rownames(seur@meta.data))] <- mcc_umap[,c("UMAP_1")]
seur$multiciliated_UMAP_2 <-NA
seur$multiciliated_UMAP_2[match(rownames(mcc_umap),rownames(seur@meta.data))] <-mcc_umap[,c("UMAP_2")]

seur$basal_UMAP_1<-NA
seur$basal_UMAP_1[match(rownames(basal_umap),rownames(seur@meta.data))] <- basal_umap[,c("UMAP_1")]
seur$basal_UMAP_2<-NA
seur$basal_UMAP_2[match(rownames(basal_umap),rownames(seur@meta.data))] <- basal_umap[,c("UMAP_2")]

columns_to_keep=c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","UMAP_1","UMAP_2","timecourse_pred","cluster_names","tricyclePosition","subset_group","multiciliated_clusters","multiciliated_cluster_names","multiciliated_UMAP_1","multiciliated_UMAP_2","multiciliated_pseudotime","multiciliated_S_score","multiciliated_G2M_score","basal_UMAP_1","basal_UMAP_2","basal_S_score","basal_G2M_score")

for_printing=seur@meta.data[,columns_to_keep]

write.csv(for_printing, file="E2f7_GEO_files/E2f7_merge_metadata.csv")
library(Matrix)

writeMM(seur[["RNA"]]@counts, file="E2f7_GEO_files/E2f7_merge_counts.mtx")

## cellxgene

new_meta=read.csv("E2f7_GEO_files/temp/E2f7_merge_metadata.csv",row.names = 1)
identical(rownames(new_meta),rownames(seur@meta.data))

new_meta$assay_ontology_term_id <- "EFO:0009922"
new_meta$cell_type_ontology_term_id <- "CL:0000307"
new_meta$development_stage_ontology_term_id <-"MmusDv:0000110"
new_meta$disease_ontology_term_id <- "PATO:0000461"
new_meta$donor_id <- "na"
new_meta$is_primary_data <-TRUE
new_meta$organism_ontology_term_id <-"NCBITaxon:10090"
new_meta$self_reported_ethnicity_ontology_term_id <- "na"
new_meta$sex_ontology_term_id <- "unknown"
new_meta$suspension_type <- "cell"
new_meta$tissue_type <- "cell culture"
new_meta$tissue_ontology_term_id <- "UBERON:0001901"

sub_temp=seur
sub_temp$UMAP_1 = sub_temp$UMAP_1*-1

library(sceasy)
library(reticulate)
use_condaenv('sceasy')
dir.create("cellxgene")
sub_temp@meta.data = new_meta
sceasy::convertFormat(sub_temp, assay='RNA', from="seurat", to="anndata", main_layer='data', transfer_layers='counts', drop_single_values=FALSE, outFile='cellxgene/E2f7_AD07_full.h5ad')

## print out raw cellranger counts for cell x gene
library(Matrix)
meta_final=new_meta

wtA_names = rownames(meta_final)[meta_final$orig.ident == "WT_A"]
wtB_names=rownames(meta_final)[meta_final$orig.ident == "WT_B"]
homA_names=rownames(meta_final)[meta_final$orig.ident == "Hom_A"]
homB_names=rownames(meta_final)[meta_final$orig.ident == "Hom_B"]


wtA=ReadMtx("/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/WT_A_matrix.mtx.gz",features="/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/WT_A_features.tsv.gz",cells="/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/WT_A_barcodes.tsv.gz")
colnames(wtA)=paste0(colnames(wtA),"_1")
wtA_sub = wtA[,match(wtA_names,colnames(wtA))]


wtB=ReadMtx("/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/WT_B_matrix.mtx.gz",features="/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/WT_B_features.tsv.gz",cells="/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/WT_B_barcodes.tsv.gz")
colnames(wtB)=paste0(colnames(wtB),"_2")
wtB_sub = wtB[,match(wtB_names,colnames(wtB))]


homA=ReadMtx("/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/Hom_A_matrix.mtx.gz",features="/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/Hom_A_features.tsv.gz",cells="/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/Hom_A_barcodes.tsv.gz")
colnames(homA)=paste0(colnames(homA),"_3")
homA_sub = homA[,match(homA_names,colnames(homA))]

homB=ReadMtx("/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/Hom_B_matrix.mtx.gz",features="/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/Hom_B_features.tsv.gz",cells="/Volumes/Reiterlab_3/E2f7_em3/E2f7_em3/E2f7_GEO_files/Hom_B_barcodes.tsv.gz")
colnames(homB)=paste0(colnames(homB),"_4")
homB_sub = homB[,match(homB_names,colnames(homB))]

full=cbind(wtA_sub,wtB_sub,homA_sub,homB_sub)

write(rownames(full), file="cellxgene/E2f7_AD07_full_combined_genes.tsv")
write(colnames(full), file="cellxgene/E2f7_AD07_full_combined_barcodes.tsv")
identical(colnames(full),rownames(meta_final))

writeMM(full, "cellxgene/E2f7_AD07_full_combined_rawCR_counts.mtx")

### focused datasets

mcc_umap=Embeddings(mccs_clean,"umap")
mccs_clean$UMAP_1<-mcc_umap[,c("UMAP_1")]
mccs_clean$UMAP_2<-mcc_umap[,c("UMAP_2")]*-1
mccs_clean@meta.data=mccs_clean@meta.data %>% mutate(multiciliated_cluster_names=recode(integrated_snn_res.0.3,
                                                                            `0`="Intermediate",
                                                                            `1`="Multiciliated 4",
                                                                            `2`="Multciliated 2",
                                                                            `3`="Multciliated 3",
                                                                            `4`="Multiciliated 1"))

mccs_clean$geno <- "WT"
mccs_clean$geno[grep("Hom*",mccs_clean$orig.ident)] <- "Hom"

columns_to_keep=c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","geno","UMAP_1","UMAP_2","multiciliated_cluster_names","tricyclePosition","pseudotime","signature_1_UCell","signature_2_UCell")

new_meta=mccs_clean@meta.data[,columns_to_keep]
colnames(new_meta)=c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Genotype","UMAP_1","UMAP_2","multiciliated_cluster_names","tricyclePosition","pseudotime","S Score","G2M Score")
## cellxgene

identical(rownames(new_meta),rownames(mccs_clean@meta.data))

new_meta$assay_ontology_term_id <- "EFO:0009922"
new_meta$cell_type_ontology_term_id <- "CL:0000307"
new_meta$development_stage_ontology_term_id <-"MmusDv:0000110"
new_meta$disease_ontology_term_id <- "PATO:0000461"
new_meta$donor_id <- "na"
new_meta$is_primary_data <-FALSE
new_meta$organism_ontology_term_id <-"NCBITaxon:10090"
new_meta$self_reported_ethnicity_ontology_term_id <- "na"
new_meta$sex_ontology_term_id <- "unknown"
new_meta$suspension_type <- "cell"
new_meta$tissue_type <- "cell culture"
new_meta$tissue_ontology_term_id <- "UBERON:0001901"

sub_temp=mccs_clean

library(sceasy)
library(reticulate)
use_condaenv('sceasy')

sub_temp@meta.data = new_meta
sceasy::convertFormat(sub_temp, assay='RNA', from="seurat", to="anndata", main_layer='data', transfer_layers='counts', drop_single_values=FALSE, outFile='cellxgene/E2f7_AD07_MCCfocused.h5ad')


basal_umap=Embeddings(basal_cycling,"umap")
basal_cycling$UMAP_1<-basal_umap[,c("UMAP_1")]*-1
basal_cycling$UMAP_2<-basal_umap[,c("UMAP_2")]

basal_cycling@meta.data=basal_cycling@meta.data %>% mutate(basal_cluster_names=recode(timecourse_pred,
                                                                                        `1`="Basal stem",
                                                                                        `5`="Proliferating basal stem (S)",
                                                                                        `6`="Proliferating basal stem (G2/M)"))

columns_to_keep=c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","geno","UMAP_1","UMAP_2","basal_cluster_names","tricyclePosition","signature_1_UCell","signature_2_UCell")

new_meta=basal_cycling@meta.data[,columns_to_keep]
colnames(new_meta)=c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Genotype","UMAP_1","UMAP_2","basal_cluster_names","tricyclePosition","S Score","G2M Score")

identical(rownames(new_meta),rownames(basal_cycling@meta.data))

new_meta$assay_ontology_term_id <- "EFO:0009922"
new_meta$cell_type_ontology_term_id <- "CL:0000307"
new_meta$development_stage_ontology_term_id <-"MmusDv:0000110"
new_meta$disease_ontology_term_id <- "PATO:0000461"
new_meta$donor_id <- "na"
new_meta$is_primary_data <-FALSE
new_meta$organism_ontology_term_id <-"NCBITaxon:10090"
new_meta$self_reported_ethnicity_ontology_term_id <- "na"
new_meta$sex_ontology_term_id <- "unknown"
new_meta$suspension_type <- "cell"
new_meta$tissue_type <- "cell culture"
new_meta$tissue_ontology_term_id <- "UBERON:0001901"

sub_temp=basal_cycling

library(sceasy)
library(reticulate)
use_condaenv('sceasy')

sub_temp@meta.data = new_meta
sceasy::convertFormat(sub_temp, assay='RNA', from="seurat", to="anndata", main_layer='data', transfer_layers='counts', drop_single_values=FALSE, outFile='cellxgene/E2f7_AD07_Basalfocused.h5ad')

