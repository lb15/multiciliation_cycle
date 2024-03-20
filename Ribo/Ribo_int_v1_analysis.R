library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

setwd("/Volumes/ReiterLab4/Ribo_AD03/Ribo_int/")

seur=readRDS("/Volumes/ReiterLab4/Ribo_AD03/Ribo_int/v1/Ribo_int_v1.rds")


#### COLORS #####
mcc_color <- colorRampPalette(c("lightskyblue2", "royalblue4"))


colors=c("#009245","#F7A12D","indianred",mcc_color(5)[3],"darkgoldenrod","#998675",mcc_color(5)[1],"#996699","plum",mcc_color(5)[5],"grey17")

basal_to_mcc_col=c(colors[1:2],mcc_color(4)[c(4,2,3,1)])
pseudo_mcc=c(colors[1:2],"indianred",mcc_color(5)[c(3,1,5)])

tricolors=c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9")

#######

png("v1/Ribo_int_v1_umap.png",height=800,width=1100)
DimPlot(seur,label=T,group.by="integrated_snn_res.0.3",label.size = 6,pt.size = 1.5,cols = colors)
dev.off()

png("v1/Ribo_int_v1_umap_split.png",height=1100,width=1100)
DimPlot(seur,split.by = "orig.ident",label=T,pt.size=1,label.size=2,ncol = 3,cols=colors)
dev.off()

DimPlot(seur, group.by="integrated_snn_res.0.3",label=T)
tab=melt(as.data.frame(prop.table(table(seur$orig.ident,seur$integrated_snn_res.0.3),1)*100))
tab$treatment <- "DMSO"
tab$treatment[grep("Ribo*",tab$Var1)]<-"Ribo"
tab$replicate <- "A"
tab$replicate[grep("*B",tab$Var1)] <- "B"
tab$replicate[grep("*C",tab$Var1)]<-"C"


write.csv(tab, file="v1/Ribo_int_v1_res0.3_cluster_perc.csv")

DefaultAssay(seur) <- "RNA"
markers=FindAllMarkers(seur, assay="RNA",group.by="integrated_snn_res.0.3")
write.csv(markers, file="v1/Ribo_int_v1_res0.3_markers.csv")

###### LABEL TRANSFER #####
##ran on wynton, import predictions

pred=read.csv("v1/Ribo_int_v1_mcc_timecourse_v3_13dim_predictions.csv")

identical(pred$X,rownames(seur@meta.data))

seur$LT_clusters <- pred$predicted.id
seur$treatment<-"DMSO"
seur$treatment[grep("Ribo*",seur$orig.ident)]<-"Ribo"

png("v1/Ribo_int_v1_mcc_timecourse_predictions_umap.png",height=800,width=1100)
DimPlot(seur, group.by="LT_clusters",pt.size=1.5,label.size = 6,label=T)
dev.off()

DimPlot(seur, split.by="treatment")
tab=melt(as.data.frame(prop.table(table(seur$orig.ident,seur$LT_clusters),1)*100))
tab$treatment <- "DMSO"
tab$treatment[grep("Ribo*",tab$Var1)]<-"Ribo"
tab$replicate <- "A"
tab$replicate[grep("*B",tab$Var1)] <- "B"
tab$replicate[grep("*C",tab$Var1)]<-"C"

pdf("v1/Ribo_int_v1_cluster_percentage.pdf",height=8,width=11,useDingbats = F)
ggplot(tab, aes(x=treatment, y=value))+
        geom_line(aes(group=replicate))+
        geom_point(aes(color=Var1))+
        facet_wrap(~Var2,scales = "free")+
        scale_color_manual(values=c(rep("grey39",3),rep("magenta",3)))+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust=1),axis.title.y = element_text(size=12))+
        xlab("")+
        ylab("Percentage of Cells")
dev.off()

DefaultAssay(seur) <- "RNA"

FeaturePlot(seur, "Gmnc",split.by="orig.ident")


####### TRICYCLE ########

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


plot_ccposition_den(sc$tricyclePosition,
                    sc$orig.ident, 'orig.ident',
                    bw = 10, fig.title = "Kernel density of \u03b8") +
        theme_bw(base_size = 14)

#### place values back into Seurat object
identical(colnames(seur),colnames(sc))
seur$tricyclePosition <- sc$tricyclePosition

png("v1/Ribo_int_v1_tricycle_split.png",height=800,width=1100)
FeaturePlot(seur, "tricyclePosition",pt.size=1.5,split.by = "orig.ident")+
        patchwork::plot_layout(ncol = 3, nrow = 2) &        
        scale_color_gradientn(colors=c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9")) &
        theme_void()&
        theme(legend.position = "none")
dev.off()

png("v1/Ribo_int_v1_tricycle.png",height=800,width=1100)
FeaturePlot(seur, "tricyclePosition",pt.size=1.5,raster=T) &        
        scale_color_gradientn(colors=c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9")) &
        theme_void()&
        theme(legend.position = "none")
dev.off()


seur$tricycle_bin = cut_interval(seur$tricyclePosition, n=20,labels=F)
DimPlot(seur, group.by="tricycle_bin",split.by = "orig.ident")

#### SEURAT CELL CYCLE #####

library(UCell)

cc_s=read.csv("~/Documents/Reiter_Seq/s_phase_human_mgi_genes.csv")
cc_g2m=read.csv("~/Documents/Reiter_Seq/g2m_phase_human_mgi_genes.csv")

sgenes=cc_s$MGI.symbol
g2mgenes=cc_g2m$MGI.symbol

seur = AddModuleScore_UCell(seur,features = list(sgenes,g2mgenes),assay = "RNA") 
seur$orig.ident=factor(seur$orig.ident,levels=c("DMSOA","DMSOB","DMSOC","RiboA","RiboB","RiboC"))

saveRDS(seur,"/Volumes/ReiterLab4/Ribo_AD03/Ribo_int/v1/Ribo_int_v1.rds")

VlnPlot(seur, features=c("signature_1_UCell","signature_2_UCell"),group.by="orig.ident",pt.size = 0)

png("v1/Ribo_int_v1_seurat_cc_scores.png",height=400,width=800)
FeaturePlot(seur, c("signature_1_UCell","signature_2_UCell"),order=T,pt.size=1.5)&
        theme_void()&
        scale_color_gradient(low = "grey17",high="cyan3")
dev.off()

png("v1/Ribo_int_v1_seurat_cc_scores_max0.2.png",height=400,width=800)
FeaturePlot(seur, c("signature_1_UCell","signature_2_UCell"),order=T,pt.size=1.5,max.cutoff = 0.2,)&
        theme_void()&
        scale_color_gradient(low = "grey17",high="cyan3")
dev.off()


g2m=ggplot(seur@meta.data,aes(x=pseudotime,y=signature_2_UCell))+
        ylab("G2M Score")+
        geom_point(aes(color=integrated_snn_res.0.3))+
        scale_color_manual(values=c(colors))+
        geom_smooth(color="black")+
        theme_classic()+
        theme(axis.text=element_text(size=12,color="black"),axis.title.x = element_blank())+
        facet_wrap(~treatment)

s=ggplot(seur@meta.data,aes(x=pseudotime,y=signature_1_UCell))+
        ylab("S Score")+
        geom_point(aes(color=integrated_snn_res.0.3))+
        scale_color_manual(values=c(colors))+
        geom_smooth(color="black")+
        theme_classic()+
        theme(axis.text=element_text(size=12,color="black"),axis.title.x = element_blank())+
        facet_wrap(~treatment)

png("v1/Ribo_int_v1_sphase_score_vs_pseudotime.png", height=400,width=600)
s
dev.off()

png("v1/Ribo_int_v1_G2Mphase_score_vs_pseudotime.png", height=400,width=600)
g2m
dev.off()

smcc=ggplot(seur@meta.data,aes(x=pseudotime_mcc,y=signature_1_UCell))+
        ylab("S Score")+
        xlab("Pseudotime")+
        geom_point(aes(color=integrated_snn_res.0.3))+
        scale_color_manual(values=c(colors))+
        geom_smooth(color="black")+
        theme_classic()+
        theme(axis.text=element_text(size=12,color="black"),axis.title.x = element_blank())+
        facet_wrap(~treatment)
        
g2mmcc=ggplot(seur@meta.data,aes(x=pseudotime_mcc,y=signature_2_UCell))+
        ylab("G2M Score")+
        xlab("Pseudotime")+
        geom_point(aes(color=integrated_snn_res.0.3))+
        scale_color_manual(values=c(colors))+
        geom_smooth(color="black")+
        theme_classic()+
        theme(axis.text=element_text(size=12,color="black"),axis.title.x = element_blank())+
        facet_wrap(~treatment)        

g2m_auc=ggplot(seur@meta.data,aes(x=pseudotime_mcc,y=signature_2_UCell))+
        ylab("G2M Score")+
        xlab("Pseudotime")+
        geom_smooth(color="black")+
        theme_classic()+
        theme(axis.text=element_text(size=12,color="black"),axis.title.x = element_blank())+
        facet_wrap(~orig.ident) 

ggp_data <- ggplot_build(g2m_auc)$data[[1]]

s_auc=ggplot(seur@meta.data,aes(x=pseudotime_mcc,y=signature_1_UCell))+
        ylab("S Score")+
        xlab("Pseudotime")+
        geom_smooth(color="black")+
        theme_classic()+
        theme(axis.text=element_text(size=12,color="black"),axis.title.x = element_blank())+
        facet_wrap(~orig.ident) 

ggp_data_s <- ggplot_build(s_auc)$data[[1]]

g2m_aucs = NULL

for(x in 1:6){
        dat= ggp_data[ggp_data$PANEL == x,]
        area = MESS::auc(dat$x, dat$y)
        g2m_aucs = c(g2m_aucs, area)
}

s_aucs=NULL
for(x in 1:6){
        dat= ggp_data_s[ggp_data_s$PANEL == x,]
        area = MESS::auc(dat$x, dat$y)
        s_aucs = c(s_aucs, area)
}

write.csv(s_aucs, file="v1/Ribo_int_v1_s_score_vs_pseudotime_mcc_AUC.csv")
write.csv(g2m_aucs, file="v1/Ribo_int_v1_g2m_score_vs_pseudotime_mcc_AUC.csv")

ggsave(ggarrange(s,g2m,nrow=2), file = "mccs_clean/E2f7_em3_mccs_clean_S_G2M_score_violin.pdf",width=6,height=6)

png("v1/Ribo_int_v1_sphase_score_vs_pseudotime_mcc.png", height=400,width=600)
smcc
dev.off()

png("v1/Ribo_int_v1_G2Mphase_score_vs_pseudotime_mcc.png", height=400,width=600)
g2mmcc
dev.off()

sdat=seur@meta.data %>% group_by(integrated_snn_res.0.3,orig.ident) %>% summarise(s_score=mean(signature_1_UCell))
sdat$treatment = "DMSO"
sdat$treatment[grep("Ribo*",sdat$orig.ident)]<-"Ribo"
ggplot(sdat, aes(x=integrated_snn_res.0.3, y=s_score, color=treatment))+
        geom_boxplot()

g2mdat=seur@meta.data %>% group_by(integrated_snn_res.0.3,orig.ident) %>% summarise(g2m_score=mean(signature_2_UCell))
g2mdat$treatment = "DMSO"
g2mdat$treatment[grep("Ribo*",g2mdat$orig.ident)]<-"Ribo"
ggplot(g2mdat, aes(x=integrated_snn_res.0.3, y=g2m_score, color=treatment))+
        geom_boxplot()

write.csv(sdat, file="v1/Ribo_int_v1_sphase_meanscores_res0.3.csv")
write.csv(g2mdat, file="v1/Ribo_int_v1_g2mphase_meanscores_res0.3.csv")

########## Scores over pseudotime #######
## make 10 equal range sized bins across pseudotime
seur$pseudotime_mcc_bin = cut_interval(seur$pseudotime_mcc, n=10,labels=F)

DimPlot(seur, group.by = "pseudotime_mcc_bin")

sdat_pseudo=seur@meta.data %>% 
        filter(!is.na(pseudotime_mcc_bin)) %>% 
        group_by(pseudotime_mcc_bin,orig.ident) %>% 
        summarise(s_score=mean(signature_1_UCell))
sdat_pseudo$treatment = "DMSO"
sdat_pseudo$treatment[grep("Ribo*",sdat_pseudo$orig.ident)]<-"Ribo"
ggplot(sdat_pseudo, aes(x=as.factor(pseudotime_mcc_bin), y=s_score, color=treatment))+
        geom_boxplot()

g2mdat_pseudo=seur@meta.data %>% 
        filter(!is.na(pseudotime_mcc_bin)) %>% 
        group_by(pseudotime_mcc_bin,orig.ident) %>% 
        summarise(g2m_score=mean(signature_2_UCell))
g2mdat_pseudo$treatment = "DMSO"
g2mdat_pseudo$treatment[grep("Ribo*",g2mdat_pseudo$orig.ident)]<-"Ribo"
ggplot(g2mdat_pseudo, aes(x=as.factor(pseudotime_mcc_bin), y=g2m_score, color=treatment))+
        geom_boxplot()

write.csv(sdat_pseudo, file="v1/Ribo_int_v1_sphase_meanscores_pseudotime_mcc_bin.csv")
write.csv(g2mdat_pseudo, file="v1/Ribo_int_v1_g2mphase_meanscores_pseudotime_mcc_bin.csv")

####### cell cycle module scoring ######
library(stringr)
library(UCell)
cc_mods=read.table("v1/Cyclebase_Genes.txt",header=T)
cc_mods$Gene_MGI = str_to_title(cc_mods$Gene)

g1=cc_mods$Gene_MGI[cc_mods$Stage == "G1"]
s=cc_mods$Gene_MGI[cc_mods$Stage == "S"]
g2=cc_mods$Gene_MGI[cc_mods$Stage == "G2"]
m=cc_mods$Gene_MGI[cc_mods$Stage == "M"]

seur = AddModuleScore_UCell(seur,features = list(g1score=g1,sscore=s,g2score=g2,mscore=m),assay = "RNA")


png("v1/Ribo_int_g1score.png",height=800,width=1100)
FeaturePlot(seur,"g1score_UCell",pt.size=1.5,max.cutoff = "q99",order=T) + scale_color_viridis()
dev.off()
png("v1/Ribo_int_sscore.png",height=800,width=1100)
FeaturePlot(seur,"sscore_UCell",max.cutoff = "q99",pt.size = 1.5,order=T) + scale_color_viridis()
dev.off()

png("v1/Ribo_int_g2mscore.png",height=800,width=1100)
FeaturePlot(seur,"g2score_UCell",max.cutoff = "q99",pt.size = 1.5,order=T) + scale_color_viridis()
dev.off()
png("v1/Ribo_int_mscore.png",height=800,width=1100)
FeaturePlot(seur,"mscore_UCell",max.cutoff = "q99",pt.size = 1.5,order=T) + scale_color_viridis()
dev.off()

####### gene plot over psuedotime ######

meta=melt(seur@meta.data[,c("pseudotime_mcc","g1score_UCell","sscore_UCell","g2score_UCell","mscore_UCell","treatment")],id.vars = c("pseudotime_mcc","treatment"))

pdf("v1/Ribo_int_v1_g1_s_g2_m_scores_pseudotime.pdf",height = 11,width=8,useDingbats = F)
ggplot(meta, aes(x=pseudotime_mcc,y=value,color=treatment))+
        geom_smooth()+
        facet_wrap(~variable,ncol = 1,scales="free")+
        theme_classic()+
        ylab("Score")+
        xlab("Pseudotime")
dev.off()

DefaultAssay(seur) <- "RNA"
gene_dat=FetchData(seur,unique(c(g1,s)))
#gene_dat=FetchData(seur, c("Myb","Foxj1"))

gene_dat=cbind(gene_dat, seur@meta.data[,c("treatment","integrated_snn_res.0.3","pseudotime_mcc")])

gene_dat_filt=gene_dat %>% filter(!is.na(pseudotime_mcc))
melted=reshape2::melt(gene_dat_filt,id.vars=c("pseudotime_mcc","treatment","integrated_snn_res.0.3"))

dir.create("v1/gene_v_pseudotime")
#for(x in c("Myb","Foxj1")){ 
for(x in unique(c(g1,s))){
       pdf(paste0("v1/gene_v_pseudotime/Ribo_int_v1_",x,"_v_pseudotime.pdf"),height=4,width=5.5,useDingbats = F)
        print(ggplot(melted %>% filter(variable==x), aes(x=pseudotime_mcc,y=value,group=treatment,color=treatment))+
                geom_smooth()+
                theme_classic()+
                ggtitle(x)+
                ylab("Expression level")+
                xlab("Pseudotime"))
        dev.off()
}

##plots cluster id bars

pdf("v1/gene_v_pseudotime/example_cluster_id.pdf",height=4,width=5.5,useDingbats = F)
ggplot(melted %>% filter(variable=="Myb"), aes(x=pseudotime_mcc,y=value,group=treatment,color=treatment))+
        geom_smooth()+
        theme_classic()+
        ggtitle(x)+
        ylab("Expression level")+
        xlab("Pseudotime")+
        geom_point(aes(x=pseudotime_mcc,y=3,color=integrated_snn_res.0.3),position = "jitter")+
        scale_color_manual(values=c(pseudo_mcc,"grey67","black"))
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

setwd("/Volumes/ReiterLab4/Ribo_AD03/Ribo_int/")

seur=readRDS("/Volumes/ReiterLab4/Ribo_AD03/Ribo_int/v1/Ribo_int_v1.rds")

DimPlot(seur, label=T,group.by="integrated_snn_res.0.3")

seur$orig.ident<-factor(seur$orig.ident)
seur$group <- "DMSO"
seur$group[grep("Ribo*",seur$orig.ident)]<-"Ribo"
table(seur$group)

### need to round SoupX counts for integer input for DESeq2
mat=round(GetAssayData(seur, slot='counts',assay = "RNA"))
new_seur=SetAssayData(seur,slot = 'counts',assay="RNA",new.data = mat)

sce=as.SingleCellExperiment(new_seur,assay = "RNA")
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
                             stringr::str_remove(rownames(u), "^[0-9]+_")))
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

deseq_clusters = function(pb, clus){
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
        dds$group=relevel(dds$group, ref = "DMSO")
        
        # Transform counts for data visualization
        rld <- rlog(dds, blind=TRUE)
        
        # Plot PCA
        png(paste0("v1/Ribo_int_v1_DESeq_PCA_clus",clus,".png"),height=400,width=400)
        print(DESeq2::plotPCA(rld, intgroup = "group"))
        dev.off()
        
        dds <- DESeq(dds)
        
        png(paste0("v1/Ribo_int_v1__DESeq_DispEsts_clus",clus,".png"),height=400,width=400)
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
        
        print(resultsNames(dds))
        res <- lfcShrink(dds, 
                         coef=2,
                         lfcThreshold = 0.58)
        
        res$padj <- res2$padj
        
        res_tbl <- res %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()
        
        write.csv(res_tbl,paste0("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus",clus,".csv"))
        # Set thresholds
        padj_cutoff <- 0.05
        
        # Subset the significant results
        sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
                dplyr::arrange(padj)
        
        write.csv(sig_res,paste0("v1/Ribo_int_v1_DESeq_siggenes_Ribo_vs_DMSO_clus",clus,".csv"))
        
        
        
}

for(x in clusters){
        print(x)
        deseq_clusters(pb, x)
}



library(EnhancedVolcano)
clusters=c("0","1","2","3","4","5","6","7","8","9","10")
for(clus in clusters){
        clus1 =read.csv(paste0("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus",clus,".csv"),row.names=1)
        clus1_filt=clus1 %>% filter(!is.na(log2FoldChange))
        
        png(paste0("v1/Ribo_int_v1_clean_clus",clus,"_enhancedvolcano.png"),height=400,width=600)
        print(EnhancedVolcano(clus1_filt, 
                              lab=clus1_filt$gene,
                              x="log2FoldChange",drawConnectors = T,
                              y="padj",labSize = 6,pointSize = 4,
                              pCutoff = .0005,FCcutoff = log2(1.5)))
        dev.off()
        
}



## make master list of DE genes
### read in unfiltered list, combine

clus0=read.csv("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus0.csv",row.names = 2)
clus1=read.csv("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus1.csv",row.names=2)
clus2=read.csv("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus2.csv",row.names=2)
clus3=read.csv("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus3.csv",row.names=2)
clus4=read.csv("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus4.csv",row.names=2)
clus5=read.csv("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus5.csv",row.names=2)
clus6=read.csv("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus6.csv",row.names=2)
clus7=read.csv("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus7.csv",row.names=2)
clus8=read.csv("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus8.csv",row.names=2)
clus9=read.csv("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus9.csv",row.names=2)
clus10=read.csv("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus10.csv",row.names=2)



find_sig_genes=function(clus){
        clus_sig = clus %>% filter(abs(log2FoldChange) > log2(1.5) & padj < .0005) 
        clus_sig$dir <- "DOWN"
        clus_sig$dir[clus_sig$log2FoldChange > 0] <- "UP"
        print(paste0("finished cluster "))
        return(clus_sig)
}

to_test= list(clus0,clus1,clus2,clus3,clus4,clus5,clus6,clus8)
sig_genes=NULL
x=1
for(i in to_test){
        genes=find_sig_genes(i)
        sig_genes[[x]]<-genes
        x=x+1
}

df=NULL
for(x in 1:8){
        genez=table(sig_genes[[x]]$dir)
        df=rbind(df,genez)
}

rownames(df) <- c("clus0","clus1","clus2","clus3","clus4","clus5","clus6","clus8")

write.csv(df, file="v1/Ribo_int_siggenes_updowntotals_padj0.0005_logFC2.csv")

library(UpSetR)

listInput <- list(clus0 = rownames(sig_genes[[1]])[sig_genes[[1]]$dir == "DOWN"],
                  clus1 = rownames(sig_genes[[2]])[sig_genes[[2]]$dir == "DOWN"],
                  clus2 = rownames(sig_genes[[3]])[sig_genes[[3]]$dir == "DOWN"],
                  clus3 = rownames(sig_genes[[4]])[sig_genes[[4]]$dir == "DOWN"],
                  #clus4 = rownames(sig_genes[[5]])[sig_genes[[5]]$dir == "UP"],
                  #clus5 = rownames(sig_genes[[6]])[sig_genes[[6]]$dir == "UP"],
                  clus6 = rownames(sig_genes[[7]])[sig_genes[[7]]$dir == "DOWN"]
                  #clus8 = rownames(sig_genes[[8]])[sig_genes[[8]]$dir == "UP"]
                  )

pdf("v1/Ribo_int_v1_upset_UP_logFC2_padj0005.pdf",height=4,width=6,useDingbats = F)
upset(fromList(listInput),order.by="freq")
dev.off()

dir.create("v1/pseudotime_mcc_vs_DEgenes")

gene_list = unique(c(rownames(sig_genes[[1]]),rownames(sig_genes[[2]]),rownames(sig_genes[[4]]),rownames(sig_genes[[7]])))

seur_dat=FetchData(seur, gene_list)
for(gene in gene_list){
       
        tmp=seur_dat[,gene,drop=F]
        colnames(tmp)<- "GOI"
        tmp$pseudotime_mcc = seur$pseudotime_mcc
        tmp$treatment<-seur$treatment
        tmp$clusters <- (seur$integrated_snn_res.0.3)
        
        
        colors2=c("#009245","#F7A12D","grey17","indianred",mcc_color(5)[3],"darkgoldenrod","#998675",mcc_color(5)[1],"#996699","plum",mcc_color(5)[5])
        
        pdf(paste0("v1/pseudotime_mcc_vs_DEgenes/Ribo_int_v1_mccs_subset_pseudotime_vs_DEgenes_",gene,".pdf"),height=4,width=6,useDingbats = F)
        print(ggplot(tmp, aes(x=pseudotime_mcc,y=GOI))+
                      #geom_point(aes(color=clusters),alpha=0.2)+
                      geom_smooth(aes(color=treatment))+
                      scale_color_manual(values=c("#585858","#A00000"))+
                      theme_classic()+
                      xlab("Pseudotime")+
                      ylab(gene)
        )
        dev.off()
}


## convert to cluster names for supp table
colnames(clus0)=paste0("Intermediate_",colnames(clus0))
colnames(clus1)=paste0("Basal_stem_",colnames(clus1))
colnames(clus2)=paste0("Secretory_",colnames(clus2))
colnames(clus3)=paste0("Multiciliated_2_",colnames(clus3))
colnames(clus4)=paste0("Tnfrsf12a_",colnames(clus4))
colnames(clus5)=paste0("Nupr1_",colnames(clus5))
colnames(clus6)=paste0("Multiciliated_1_",colnames(clus6))
colnames(clus7)=paste0("Proliferating_",colnames(clus7))
colnames(clus8)=paste0("Neuroendocrine_",colnames(clus8))
colnames(clus9)=paste0("Multiciliated_2_",colnames(clus9))
colnames(clus10)=paste0("Krt14_basal_stem_",colnames(clus10))

full_dat=cbind(clus0%>% select(-"Intermediate_X"),clus1%>% select(-"Basal_stem_X"))
full_dat=cbind(full_dat,clus2%>% select(-"Secretory_X"))
full_dat=cbind(full_dat, clus3%>% select(-"Multiciliated_2_X"))
full_dat=cbind(full_dat,clus4%>% select(-"Tnfrsf12a_X"))
full_dat=cbind(full_dat,clus5%>% select(-"Nupr1_X"))
full_dat=cbind(full_dat,clus6%>% select(-"Multiciliated_1_X"))
full_dat=cbind(full_dat,clus7%>% select(-"Proliferating_X"))
full_dat=cbind(full_dat,clus8%>% select(-"Neuroendocrine_X"))
full_dat=cbind(full_dat,clus9%>% select(-"Multiciliated_2_X"))
full_dat=cbind(full_dat,clus10%>% select(-"Krt14_basal_stem_X"))


write.csv(full_dat, file="v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_allclusters_named.csv")



#### Cluster Abundance - https://bioconductor.org/books/3.13/OSCA.multisample/differential-abundance.html

abundances <- table(seur$integrated_snn_res.0.3,seur$orig.ident)
abundances <- unclass(abundances)
abundances 

library(edgeR)
sample_dat=as.data.frame(unclass(table(seur$orig.ident)))
colnames(sample_dat) <-"cell_num"
sample_dat$group <- "A"
sample_dat$group[grep("*B",rownames(sample_dat))]<-"B"
sample_dat$group[grep("*C",rownames(sample_dat))]<-"C"

sample_dat$treatment <- "DMSO"
sample_dat$treatment[grep("Ribo*",rownames(sample_dat))]<-"Ribo"

y.ab <- DGEList(abundances, samples=sample_dat)
keep <- filterByExpr(y.ab, group=y.ab$samples$treatment)
y.ab <- y.ab[keep,]
summary(keep)

design <- model.matrix(~factor(group) + factor(treatment), y.ab$samples)
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)

write.csv(res$table,file="v1/Ribo_int_1_DA_testing.csv")

### include assumption that most labels do not change to combat composition effects

y.ab2 <- calcNormFactors(y.ab)
y.ab2$samples$norm.factors

y.ab2 <- estimateDisp(y.ab2, design, trend="none")
fit.ab2 <- glmQLFit(y.ab2, design, robust=TRUE, abundance.trend=FALSE)
res2 <- glmQLFTest(fit.ab2, coef=ncol(design))
topTags(res2, n=10)

write.csv(topTags(res2,n=11), file="v1/Ribo_int_1_DA_calcNormFac_testing.csv")


dat=unclass(topTags(res2,n=11))$table
dat$cluster <- rownames(dat)

dat_gg=melt(dat,id.vars = c("FDR","logFC"),measure.vars = "cluster")
dat_gg$value <- factor((dat_gg$value),levels=c(0:11))

pdf("v1/Ribo_int_v1_DA_plot.pdf",height=3,width=4,useDingbats = F)
ggplot(dat_gg, aes(x=logFC,y=FDR, color=value))+
        geom_vline(xintercept = log2(1.5),linetype=3)+
        geom_vline(xintercept = -log2(1.5),linetype=3)+
        geom_hline(yintercept = 0.05,linetype=3)+
        geom_point(size=4)+
        scale_color_manual(values=colors)+
        xlim(c(-7,7))+
        theme_classic()+
        xlab("Log2 Fold Change")+
        ylab("FDR")+
        theme(axis.text=element_text(size=12,color="black"))
dev.off()

##### plot log fold changes of proportion ratios for visualization of above results
tab <- table(seur$orig.ident, seur$integrated_snn_res.0.3)
props <- tab/rowSums(tab)
props.pseudo <- t((tab + 0.5)/rowSums(tab + 0.5))

props.pseudo=data.frame(unclass(props.pseudo))

props.pseudo$Ribo_vs_DMSO_A <- props.pseudo$RiboA/props.pseudo$DMSOA
props.pseudo$Ribo_vs_DMSO_B <- props.pseudo$RiboB/props.pseudo$DMSOB
props.pseudo$Ribo_vs_DMSO_C <- props.pseudo$RiboC/props.pseudo$DMSOC

props.pseudo$clusters=rownames(props.pseudo)

props.pseudo=props.pseudo %>% mutate(cluster_names=recode(clusters,
                                                          `1`="Basal stem",
                                                          `2`="Secretory",
                                                          `3`="Multiciliated 2",
                                                          `4`="Tnfrsf12a+",
                                                          `5`="Nupr1+",
                                                          `6`="Multiciliated 1",
                                                          `7`="Proliferating",
                                                          `8`="Neuroendocrine",
                                                          `9`="Multiciliated 3",
                                                          `0`="Intermediate",
                                                          `10`="Krt14+ basal stem"
))
melt_prop=reshape2::melt(props.pseudo[,c(10,grep("Ribo_vs_DMSO*",colnames(props.pseudo)))])

melt_prop$log2FC = log2(melt_prop$value)

melt_prop=melt_prop %>% mutate(cluster_names=recode(clusters,
                                                    `1`="Basal stem",
                                                    `2`="Secretory",
                                                    `3`="Multiciliated 2",
                                                    `4`="Tnfrsf12a+",
                                                    `5`="Nupr1+",
                                                    `6`="Multiciliated 1",
                                                    `7`="Proliferating",
                                                    `8`="Neuroendocrine",
                                                    `9`="Multiciliated 3",
                                                    `0`="Intermediate",
                                                    `10`="Krt14+ basal stem"
))


melt_prop_agg= melt_prop %>% 
        group_by(clusters) %>%
        summarise(mean_log2FC=mean(log2FC),sd=sd(log2FC),flcse=sd/sqrt(n()))

melt_prop_agg=melt_prop_agg %>% mutate(cluster_names=recode(clusters,
                                                            `1`="Basal stem",
                                                            `2`="Secretory",
                                                            `3`="Multiciliated 2",
                                                            `4`="Tnfrsf12a+",
                                                            `5`="Nupr1+",
                                                            `6`="Multiciliated 1",
                                                            `7`="Proliferating",
                                                            `8`="Neuroendocrine",
                                                            `9`="Multiciliated 3",
                                                            `0`="Intermediate",
                                                            `10`="Krt14+ basal stem"
))

write.csv(props.pseudo, file="v1/Ribo_int_v1_props_pseudo.csv")
write.csv(melt_prop_agg, file="v1/Ribo_int_v1_props_se_plot.csv")
pdf("v1/Ribo_int_v1_Ribo_DMSO_proportions_speckle_sebars.pdf",height=6,width=6,useDingbats = F)

pdf("v1/Ribo_int_v1_celltypeproporation_reps_mean_se.pdf",height=8,width=8,useDingbats = F)        
ggplot() +
        geom_bar(data = melt_prop_agg, aes(x = clusters, y = mean_log2FC, fill = clusters), stat = "identity") +
        geom_point(data = melt_prop, aes(x = clusters, y = log2FC, group = clusters), position = position_dodge2(width = 0.9), size = 4, color = "black") +
        geom_errorbar(data = melt_prop_agg, aes(x = clusters, ymin = mean_log2FC - flcse, ymax = mean_log2FC + flcse), linewidth = 1, width = 0.7, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0) +
        scale_fill_manual(values=colors[c(1:2,11,3:10)]) +
        xlab("") +
        theme_classic() +
        theme(axis.text = element_text(color = "black", size = 14))
dev.off()


sub_clus_melt = melt_prop %>% filter(cluster_names %in% c("Basal stem", "Intermediate","Multiciliated 1","Multiciliated 2","Multiciliated 3","Secretory"))
sub_clus_agg = melt_prop_agg %>% filter(cluster_names %in% c("Basal stem", "Intermediate","Multiciliated 1","Multiciliated 2","Multiciliated 3","Secretory"))

pdf("v1/Ribo_int_v1_celltypeproporation_reps_mean_se.pdf",height=8,width=8,useDingbats = F)        
ggplot() +
        geom_bar(data = sub_clus_agg, aes(x = cluster_names, y = mean_log2FC, fill = clusters), stat = "identity") +
        geom_point(data = sub_clus_melt, aes(x = cluster_names, y = log2FC, group = clusters), position = position_dodge2(width = 0.9), size = 4, color = "black") +
        geom_errorbar(data = sub_clus_agg, aes(x = cluster_names, ymin = mean_log2FC - flcse, ymax = mean_log2FC + flcse), linewidth = 1, width = 0.7, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0) +
        scale_fill_manual(values=colors[c(1:2,3,4,7,10)]) +
        xlab("") +
        theme_classic() +
        theme(axis.text = element_text(color = "black", size = 14))
dev.off()

ggplot(melt_prop_agg, aes(x=clusters,y=mean_log2FC,fill=clusters))+
        geom_bar(stat="identity")+
        geom_errorbar(aes(x=clusters, ymin=mean_log2FC-flcse, ymax=mean_log2FC+flcse)) +theme_classic()+
        scale_fill_manual(values=colors[c(1:2,11,3:10)])+
        xlab("Cluster")+
        ylab("Log2 of Ribo/DMSO cluster proportion")
dev.off()  

#### SPECKLE ####
library(speckle)
speck = propeller(cluster=seur$integrated_snn_res.0.3,group=seur$treatment,sample=seur$orig.ident)

write.csv(speck, file="v1/Ribo_int_v1_speckle_DA.cssv")

###### PSEUDOTIME ######
library(monocle3)
library(SeuratWrappers)

### use a subset of data for just the main cells that are connected.

Idents(seur)<- "integrated_snn_res.0.3"
seur_clean=subset(seur,idents=c(0:3,6,9))
DimPlot(seur_clean)

DefaultAssay(seur_clean) <- "RNA"
cds=as.cell_data_set(seur_clean)
cds=cluster_cells(cds)
cds=learn_graph(cds,close_loop = F,use_partition = F,learn_graph_control =list( ncenter=300))

plot_cells(cds,cell_size = 1,trajectory_graph_segment_size = 2,color_cells_by = "ident")


cds=order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime")
ps_df=as.data.frame(pseudotime(cds))

cds_mcc=choose_graph_segments(cds,clear_cds = F)
cds_sec=choose_graph_segments(cds,clear_cds = F)

saveRDS(cds, file="v1/Ribo_int_v1_monocle.rds")
cds=readRDS("v1/Ribo_int_v1_monocle.rds")

colnames(ps_df)="pseudotime"
ps_df$pseudo_bin = cut_interval(ps_df$pseudotime, n=20,labels=F)
ps_df$pseudotime_mcc = NA
ps_df$pseudotime_mcc[match(names(pseudotime(cds_mcc)),rownames(ps_df))] <- pseudotime(cds_mcc)

ps_df$pseudotime_sec = NA
ps_df$pseudotime_sec[match(names(pseudotime(cds_sec)),rownames(ps_df))] <- pseudotime(cds_sec)

seur$pseudotime <-NA
seur$pseudotime_mcc <-NA
seur$pseudotime_sec <- NA

seur$pseudotime[match(rownames(ps_df),rownames(seur@meta.data))]<-ps_df$`pseudotime(cds)`
seur$pseudotime_mcc[match(rownames(ps_df),rownames(seur@meta.data))]<-ps_df$pseudotime_mcc
seur$pseudotime_sec[match(rownames(ps_df),rownames(seur@meta.data))]<-ps_df$pseudotime_sec


png("v1/Ribo_int_v1_pseudotime.png",height=400,width=500)
FeaturePlot(seur, "pseudotime") + scale_color_viridis(option="B")
dev.off()

png("v1/Ribo_int_v1_pseudotime_mcc.png",height=400,width=500)
FeaturePlot(seur, "pseudotime_mcc") + scale_color_viridis(option="B")
dev.off()

png("v1/Ribo_int_v1_pseudotime_sec.png",height=400,width=500)
FeaturePlot(seur, "pseudotime_sec") + scale_color_viridis(option="B")
dev.off()

write.csv(ps_df, file="v1/Ribo_int_v1_clean_pseudotime_values.csv")
saveRDS(cds, file="v1/Ribo_int_v1_clean_full_monocle_object.rds")
cds=readRDS(file="v1/Ribo_int_v1_clean_full_monocle_object.rds")
saveRDS(seur, file="v1/Ribo_int_v1.rds")

png("v1/Ribo_int_v1_tricycle_vs_pseudotime_mcc.png",height=400,width=600)
ggplot(seur@meta.data, aes(x=pseudotime_mcc,y=tricyclePosition,color=integrated_snn_res.0.3))+
        geom_point()+
        facet_wrap(~treatment)+
        scale_color_manual(values=colors)+
        theme_classic()+
        xlab("Pseudotime")+
        ylab("Tricycle Position")
dev.off()

png("v1/Ribo_int_v1_tricycle_vs_pseudotime_sec.png",height=400,width=600)
ggplot(seur@meta.data, aes(x=pseudotime_sec,y=tricyclePosition,color=integrated_snn_res.0.3))+
        geom_point()+
        facet_wrap(~treatment)+
        scale_color_manual(values=colors)+
        theme_classic()+
        xlab("Pseudotime")+
        ylab("Tricycle Position")
dev.off()

############ SUBSET Basal, int, and MCC clusters ############
Idents(seur) <- "integrated_snn_res.0.3"
DimPlot(seur, label=T)
mccs = subset(seur, idents=c(1,0,3,6,9))
DimPlot(mccs, label=T)

DefaultAssay(mccs) <- "integrated"
mccs = ScaleData(mccs)
mccs=RunPCA(mccs, npcs = 40)
ElbowPlot(mccs, ndims = 40)
mccs=RunUMAP(mccs,dims = 1:13)
mccs=FindNeighbors(mccs, dims = 1:13)
mccs=FindClusters(mccs,resolution =0.2)
DimPlot(mccs, label=T,split.by = "treatment")

new_names <- c("Int","Basal","MCC 2","MCC 1","MCC 3")
names(new_names) <- levels(mccs)

mccs <- RenameIdents(object = mccs, new_names)
DimPlot(mccs, label = TRUE)
mccs[["res0.2_names"]] <- Idents(mccs)

dir.create("v1/mccs_subset")
saveRDS(mccs, file="v1/mccs_subset/Ribo_int_v1_mccs_subset.rds")
mccs=readRDS("v1/mccs_subset/Ribo_int_v1_mccs_subset.rds")

png("v1/mccs_subset/Ribo_int_v1_mcc_subset_clusters.png",width=800,height=800)
DimPlot(mccs, label=T,label.size = 12,pt.size = 1.5,group.by="res0.2_names",cols = c(basal_to_mcc_col[c(1,2,5,4,3)]))
dev.off()

png("v1/mccs_subset/Ribo_int_v1_mcc_subset_tricycle.png",width=800,height=800)
FeaturePlot(mccs, "tricyclePosition",label=T,label.size = 12,pt.size = 1.5,cols = tricolors)
dev.off()

####### DA #########

abundances <- table(mccs$res0.2_names,mccs$orig.ident)
abundances <- unclass(abundances)
abundances 

library(edgeR)
sample_dat=as.data.frame(unclass(table(mccs$orig.ident)))
colnames(sample_dat) <-"cell_num"
sample_dat$group <- "A"
sample_dat$group[grep("*B",rownames(sample_dat))]<-"B"
sample_dat$group[grep("*C",rownames(sample_dat))]<-"C"

sample_dat$treatment <- "DMSO"
sample_dat$treatment[grep("Ribo*",rownames(sample_dat))]<-"Ribo"

y.ab <- DGEList(abundances, samples=sample_dat)
keep <- filterByExpr(y.ab, group=y.ab$samples$treatment)
y.ab <- y.ab[keep,]
summary(keep)

design <- model.matrix(~factor(group) + factor(treatment), y.ab$samples)
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)

write.csv(res$table,file="v1/mccs_subset/Ribo_int_v1_mcc_subset_DA_testing.csv")

### include assumption that most labels do not change to combat composition effects

y.ab2 <- calcNormFactors(y.ab)
y.ab2$samples$norm.factors

y.ab2 <- estimateDisp(y.ab2, design, trend="none")
fit.ab2 <- glmQLFit(y.ab2, design, robust=TRUE, abundance.trend=FALSE)
res2 <- glmQLFTest(fit.ab2, coef=ncol(design))
topTags(res2, n=10)

write.csv(res2$table, file="v1/mccs_subset/Ribo_int_1_mcc_subset_DA_calcNormFac_testing.csv")

dat=unclass(topTags(res2,n=11))$table
dat$cluster <- rownames(dat)

dat_gg=melt(dat,id.vars = c("FDR","logFC"),measure.vars = "cluster")
dat_gg$value <- factor((dat_gg$value))

pdf("v1/mccs_subset/Ribo_int_v1_mcc_subset_DA_plot.pdf",height=3,width=4,useDingbats = F)
ggplot(dat_gg, aes(x=logFC,y=FDR, color=value))+
        geom_vline(xintercept = log2(2),linetype=3)+
        geom_vline(xintercept = -log2(2),linetype=3)+
        geom_hline(yintercept = 0.05,linetype=3)+
        geom_point(size=3)+
        scale_color_manual(values=c(basal_to_mcc_col[c(2,1,4,5,3)]))+
        xlim(c(-7,7))+
        theme_classic()+
        xlab("Log2 Fold Change")+
        ylab("FDR")+
        theme(axis.text=element_text(size=12,color="black"))
dev.off()


#### tricycle abundances ####


mccs$tricyclePhase <- "G1/G0"
mccs$tricyclePhase[mccs$tricyclePosition > pi/2 & mccs$tricyclePosition < pi] <- "S"
mccs$tricyclePhase[mccs$tricyclePosition >=pi & mccs$tricyclePosition < pi*7/4] <- "G2M"


DimPlot(mccs, group.by = "tricyclePhase")

abundances <- table(mccs$tricyclePhase,mccs$orig.ident)
abundances <- unclass(abundances)
abundances 

library(edgeR)
sample_dat=as.data.frame(unclass(table(mccs$orig.ident)))
colnames(sample_dat) <-"cell_num"
sample_dat$group <- "A"
sample_dat$group[grep("*B",rownames(sample_dat))]<-"B"
sample_dat$group[grep("*C",rownames(sample_dat))]<-"C"

sample_dat$treatment <- "DMSO"
sample_dat$treatment[grep("Ribo*",rownames(sample_dat))]<-"Ribo"

y.ab <- DGEList(abundances, samples=sample_dat)
keep <- filterByExpr(y.ab, group=y.ab$samples$treatment)
y.ab <- y.ab[keep,]
summary(keep)

design <- model.matrix(~factor(group) + factor(treatment), y.ab$samples)
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)

write.csv(res$table,file="v1/mccs_subset/Ribo_int_v1_mcc_subset_DA_tricyclePhase_testing.csv")

### include assumption that most labels do not change to combat composition effects

y.ab2 <- calcNormFactors(y.ab)
y.ab2$samples$norm.factors

y.ab2 <- estimateDisp(y.ab2, design, trend="none")
fit.ab2 <- glmQLFit(y.ab2, design, robust=TRUE, abundance.trend=FALSE)
res2 <- glmQLFTest(fit.ab2, coef=ncol(design))
topTags(res2, n=10)

write.csv(res2$table, file="v1/mccs_subset/Ribo_int_1_mcc_subset_DA_calcNormFac_tricyclePhase_testing.csv")

dat=unclass(topTags(res2,n=11))$table
dat$cluster <- rownames(dat)

dat_gg=melt(dat,id.vars = c("FDR","logFC"),measure.vars = "cluster")
dat_gg$value <- factor((dat_gg$value))

binned_col=colorRampPalette(tricolors)(20)

pdf("v1/mccs_subset/Ribo_int_v1_mcc_subset_DA_tricycle_plot.pdf",height=3,width=4,useDingbats = F)
ggplot(dat_gg, aes(x=logFC,y=FDR, color=value))+
        geom_vline(xintercept = log2(2),linetype=3)+
        geom_vline(xintercept = -log2(2),linetype=3)+
        geom_hline(yintercept = 0.05,linetype=3)+
        geom_point(size=3)+
        scale_color_manual(values=c(binned_col[1],binned_col[14],binned_col[8]))+
        xlim(c(-7,7))+
        theme_classic()+
        xlab("Log2 Fold Change")+
        ylab("FDR")+
        theme(axis.text=element_text(size=12,color="black"))
dev.off()

#### SUBSET only MCCs and Int


sub_mccs = subset(mccs, idents=c("Int","MCC 1","MCC 2","MCC 3"))
DimPlot(sub_mccs, label=T)

DefaultAssay(sub_mccs) <- "integrated"
sub_mccs = ScaleData(sub_mccs)
sub_mccs=RunPCA(sub_mccs, npcs = 40)
ElbowPlot(sub_mccs, ndims = 40)
sub_mccs=RunUMAP(sub_mccs,dims = 1:9)
sub_mccs=FindNeighbors(sub_mccs, dims = 1:9)
sub_mccs=FindClusters(sub_mccs,resolution =0.1)
DimPlot(sub_mccs, label=T)


new_names <- c("Int","MCC 2","MCC 1","MCC 3")
names(new_names) <- levels(sub_mccs)

sub_mccs <- RenameIdents(object = sub_mccs, new_names)
DimPlot(sub_mccs, label = TRUE)
sub_mccs[["res0.1_names"]] <- Idents(sub_mccs)

dir.create("v1/mccs_int_subset")
saveRDS(sub_mccs, file="v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset.rds")
sub_mccs=readRDS("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset.rds")

png("v1/mccs_int_subset/Ribo_int_v1_mcc_int_subset_clusters.png",width=800,height=800)
DimPlot(sub_mccs, label=T,label.size = 12,pt.size = 1.5,group.by="res0.1_names",cols = c(basal_to_mcc_col[c(1,5,4,3)]))
dev.off()

png("v1/mccs_subset/Ribo_int_v1_mcc_subset_tricycle.png",width=800,height=800)
FeaturePlot(mccs, "tricyclePosition",label=T,label.size = 12,pt.size = 1.5,cols = tricolors)
dev.off()

abundances <- table(sub_mccs$tricyclePhase,sub_mccs$orig.ident)
abundances <- unclass(abundances)
abundances 


library(edgeR)
sample_dat=as.data.frame(unclass(table(sub_mccs$orig.ident)))
colnames(sample_dat) <-"cell_num"
sample_dat$group <- "A"
sample_dat$group[grep("*B",rownames(sample_dat))]<-"B"
sample_dat$group[grep("*C",rownames(sample_dat))]<-"C"

sample_dat$treatment <- "DMSO"
sample_dat$treatment[grep("Ribo*",rownames(sample_dat))]<-"Ribo"

y.ab <- DGEList(abundances, samples=sample_dat)
keep <- filterByExpr(y.ab, group=y.ab$samples$treatment)
y.ab <- y.ab[keep,]
summary(keep)

design <- model.matrix(~factor(group) + factor(treatment), y.ab$samples)
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)

write.csv(res$table,file="v1/mccs_int_subset/Ribo_int_v1_mcc_int_subset_DA__Tricycle_testing.csv")


dat=unclass(topTags(res,n=11))$table
dat$cluster <- rownames(dat)

dat_gg=melt(dat,id.vars = c("FDR","logFC"),measure.vars = "cluster")
dat_gg$value <- factor((dat_gg$value))

pdf("v1/mccs_int_subset/Ribo_int_v1_mcc_int_subset_DA_tricycle_plot.pdf",height=3,width=4,useDingbats = F)
ggplot(dat_gg, aes(x=logFC,y=FDR, color=value))+
        geom_vline(xintercept = log2(1.5),linetype=3)+
        geom_vline(xintercept = -log2(1.5),linetype=3)+
        geom_hline(yintercept = 0.05,linetype=3)+
        geom_point(size=4)+
        scale_color_manual(values=c(binned_col[1],binned_col[14],binned_col[8]))+
        xlim(c(-7,7))+
        theme_classic()+
        xlab("Log2 Fold Change")+
        ylab("FDR")+
        theme(axis.text=element_text(size=12,color="black"))
dev.off()

png("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_tricyclePhase.png",height=800,width=800)
DimPlot(sub_mccs, pt.size=2.5,group.by="tricyclePhase",cols = c(binned_col[1],binned_col[14],binned_col[8]))
dev.off()

#### plot replicate log fold changes

tab <- table(sub_mccs$orig.ident, sub_mccs$tricyclePhase)
props <- tab/rowSums(tab)
props.pseudo <- t((tab + 0.5)/rowSums(tab + 0.5))

props.pseudo=data.frame(unclass(props.pseudo))

props.pseudo$Ribo_vs_DMSO_A <- props.pseudo$RiboA/props.pseudo$DMSOA
props.pseudo$Ribo_vs_DMSO_B <- props.pseudo$RiboB/props.pseudo$DMSOB
props.pseudo$Ribo_vs_DMSO_C <- props.pseudo$RiboC/props.pseudo$DMSOC

props.pseudo$clusters=rownames(props.pseudo)

props.pseudo$clusters = factor(props.pseudo$clusters, levels=c("G1/G0","S","G2M"))
write.csv(props.pseudo, file="v1/mccs_int_subset/Ribo_int_v1_mccs_sub_tricyclePhase_proportions_manual_reps.csv")

melt_prop=reshape2::melt(props.pseudo[,c(10,grep("Ribo_vs_DMSO*",colnames(props.pseudo)))])

melt_prop$log2FC = log2(melt_prop$value)
melt_prop$rep <- "A"
melt_prop$rep[grep("Ribo_vs_DMSO_B",melt_prop$variable)] <- "B"
melt_prop$rep[grep("Ribo_vs_DMSO_C",melt_prop$variable)] <- "C"

melt_prop_agg= melt_prop %>% 
        group_by(clusters) %>%
        summarise(mean_log2FC=mean(log2FC),sd=sd(log2FC),flcse=sd/sqrt(n()))


pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_sub_tricycleproporation_reps_mean_se.pdf",height=8,width=6,useDingbats = F)        
ggplot() +
        geom_bar(data = melt_prop_agg, aes(x = clusters, y = mean_log2FC, fill = clusters), stat = "identity") +
        geom_point(data = melt_prop, aes(x = clusters, y = log2FC, group = clusters), position = position_dodge2(width = 0.9), size = 4, color = "black") +
        geom_errorbar(data = melt_prop_agg, aes(x = clusters, ymin = mean_log2FC - flcse, ymax = mean_log2FC + flcse), linewidth = 1, width = 0.7, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0) +
        scale_fill_manual(values = alpha(c("darkred", "darkblue", "darkgreen"), alpha = 0.8)) +
        xlab("") +
        theme_classic() +
        theme(axis.text = element_text(color = "black", size = 14))
dev.off()

write.csv(melt_prop_agg, file="v1/mccs_int_subset/Ribo_int_v1_mccs_sub_tricyclePhase_proportions_manual_sebars.csv")

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_sub_tricyclePhase_proportions_speckle_sebars.pdf",height=6,width=6,useDingbats = F)
ggplot(melt_prop_agg, aes(x=clusters,y=mean_log2FC,fill=clusters))+
        geom_bar(stat="identity")+
        geom_errorbar(aes(x=clusters, ymin=mean_log2FC-flcse, ymax=mean_log2FC+flcse)) +theme_classic()+
        scale_fill_manual(values=c(binned_col[1],binned_col[8],binned_col[14]))+
        xlab("Cluster")+
        ylab("Log2 of Ribo/DMSO cluster proportion")
dev.off()  



library(speckle)
speck = propeller(cluster=sub_mccs$tricyclePhase,group=sub_mccs$treatment,sample=sub_mccs$orig.ident)

write.csv(speck, "v1/mccs_int_subset/Ribo_int_v1_mccs_sub_tricyclePhase_proportions_speckle.csv")


#### plot S and G2M scores in TricyclePhase

meta=sub_mccs@meta.data

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_int_S_tricyclePhase_violin.pdf",height=6,width=6,useDingbats = F)
ggplot(meta %>% filter(tricyclePhase == "S"), aes(x=treatment,y=signature_1_UCell,fill=treatment))+
        geom_violin()+
        scale_fill_manual(values=c("grey30","indianred3"))+
        xlab("")+
        ylab("S Score")+
        theme_classic()
dev.off()

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_int_G2M_tricyclePhase_violin.pdf",height=6,width=6,useDingbats = F)
ggplot(meta %>% filter(tricyclePhase == "G2M"), aes(x=treatment,y=signature_2_UCell,fill=treatment))+
        geom_violin() +
        scale_fill_manual(values=c("grey30","indianred3"))+
        xlab("")+
        ylab("G2/M Score")+
        theme_classic()
dev.off()

meta_score_sum = meta %>% group_by(orig.ident,tricyclePhase) %>%
        summarise(mean_signature_1_UCell =mean(signature_1_UCell))

write.csv(meta_score_sum, file="v1/mccs_int_subset/Ribo_int_v1_mccs_int_S_score_tricyclePhase.csv")

meta_score_sum_2 = meta %>% group_by(orig.ident,tricyclePhase) %>%
        summarise(mean_signature_1_UCell =mean(signature_2_UCell))

write.csv(meta_score_sum_2, file="v1/mccs_int_subset/Ribo_int_v1_mccs_int_G2M_score_tricyclePhase.csv")

full_umap =Embeddings(sub_mccs, "umap")

meta$UMAP_1<-full_umap[,c("UMAP_1")]
meta$UMAP_2<-full_umap[,c("UMAP_2")]

ggplot(meta %>% filter(tricyclePhase =="S"), aes(x=UMAP_1,y=UMAP_2,color=signature_1_UCell))+
        geom_point()+
        facet_wrap(~treatment)+
        theme_classic()+
        scale_color_viridis()

##### MCC INT PSEUDOTIME

library(monocle3)
library(SeuratWrappers)

DimPlot(sub_mccs,label=T)

DefaultAssay(sub_mccs) <- "RNA"
cds=as.cell_data_set(sub_mccs)
cds=cluster_cells(cds)
cds=learn_graph(cds,close_loop = F,use_partition = F,learn_graph_control =list( ncenter=200))
#l)

plot_cells(cds,cell_size = 1,trajectory_graph_segment_size = 2)


cds=order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime")

saveRDS(cds, file="v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_monocle.rds")
cds=readRDS("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_monocle.rds")
ps_df=as.data.frame(pseudotime(cds))


sub_mccs$pseudotime <-NA
sub_mccs$pseudotime[match(rownames(ps_df),rownames(sub_mccs@meta.data))]<-ps_df$`pseudotime(cds)`
DefaultAssay(sub_mccs) <- "RNA"
dat = FetchData(sub_mccs, c("Mycl","Gmnc","Mcidas","E2f7","Foxj1","Myb","Cav1"))

dat$treatment <- sub_mccs$treatment
dat$res0.1_names <- sub_mccs$res0.1_names
dat$pseudotime <- sub_mccs$pseudotime
dat$rep <- sub_mccs$orig.ident

ggplot(dat, aes(x=pseudotime,y=Mycl))+
        geom_point(aes(color=res0.1_names),alpha=0.2)+
        geom_smooth(aes(color=treatment))


#################### PLOTS #######################
## full dataset
pdf("v1/Ribo_int_v1_umap.pdf",height=8,width=11,useDingbats = F)
DimPlot(seur, pt.size=1,group.by="integrated_snn_res.0.3",cols=colors)+
        theme_void()+
        ggtitle("")+
        NoLegend()
dev.off()

pdf("v1/Ribo_int_v1_umap_clusterlabels.pdf",height=8,width=11,useDingbats = F)
DimPlot(seur, pt.size=1,group.by="integrated_snn_res.0.3",cols=colors,label=T,label.size = 12)+
        theme_void()+
        ggtitle("")+
        NoLegend()
dev.off()

pdf("v1/Ribo_int_v1_umap_split.pdf",height=4,width=11,useDingbats = F)
DimPlot(seur, pt.size=0.5,group.by="integrated_snn_res.0.3",cols=colors,split.by = "treatment")+
        theme_void()+
        ggtitle("")+
        NoLegend()
dev.off()

pdf("v1/Ribo_int_v1_trajectory.pdf",height=8,width=11,useDingbats = F)
plot_cells(cds,cell_size = 1,trajectory_graph_segment_size = 2,color_cells_by = "ident",label_cell_groups = F,label_roots = F,label_leaves = F,label_branch_points = F) + scale_color_manual(values = colors[c(1:4,7,10)])+ theme_void() + theme(legend.position = "none")
dev.off()

pdf("v1/Ribo_int_v1_pseudotime_mcc.pdf",height=8,width=10,useDingbats = F)
FeaturePlot(seur, "pseudotime_mcc",pt.size=1) + scale_color_viridis(option = "plasma",direction=-1) +theme_void()+ ggtitle("")
dev.off()

pdf("v1/Ribo_int_v1_tricycle_split.pdf",height=4,width=11,useDingbats = F)
FeaturePlot(seur, "tricyclePosition",pt.size=0.5,cols=tricolors,split.by = "treatment")&
        theme_void()&
        NoLegend()
dev.off()

seur$MCC_dataset = FALSE
seur$MCC_dataset[match(rownames(sub_mccs@meta.data),rownames(seur@meta.data))]<-TRUE

pdf("v1/Ribo_int_v1_MCCs_onumap.pdf",height=8,width=11,useDingbats = F)
DimPlot(seur, group.by = "MCC_dataset",cols = c("grey30",mcc_color(6)[3]),pt.size=1)+
        theme_void()+
        ggtitle("")
dev.off()

write.csv(cbind(c(0:10),c("Intermediate","Basal stem","Secretory","Multiciliated 2","Tnfrsf12a+","Nupr1+","Multiciliated 1","Proliferating","Neuroendocrine","Multiciliated 3","Krt14+ basal stem")), "v1/Ribo_int_res0.3_to_clusternames.csv")

## mcc_int dataset

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_umap.pdf",height=8,width=11,useDingbats = F)
DimPlot(sub_mccs, pt.size=1, group.by="res0.1_names",cols=c(basal_to_mcc_col[c(1,5,4,3)])) +
        theme_void()+
        ggtitle("")+
        NoLegend()
dev.off()

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_umap_split.pdf",height=4,width=11,useDingbats = F)
DimPlot(sub_mccs, pt.size=1, group.by="res0.1_names",cols=c(basal_to_mcc_col[c(1,5,4,3)]),split.by = "treatment") +
        theme_void()+
        ggtitle("")+
        NoLegend()
dev.off()

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_umap_trajectory.pdf",height=4,width=11,useDingbats = F)
plot_cells(cds, color_cells_by = "ident",trajectory_graph_segment_size =1.5,label_groups_by_cluster = F,label_cell_groups = F,label_roots = F,label_leaves = F) + scale_color_manual(values=c(basal_to_mcc_col[c(1,5,4,3)])) + theme_void()
dev.off()

binned_col=colorRampPalette(tricolors)(20)

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_triyclePhase_split.pdf",height=4,width=11,useDingbats = F)
DimPlot(sub_mccs, pt.size=1, group.by="tricyclePhase",cols=c(binned_col[1],binned_col[14],binned_col[8]),split.by="treatment") +
        theme_void()+
        ggtitle("")+
        NoLegend()
dev.off()

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_triyclePhase_split_withlegend.pdf",height=4,width=11,useDingbats = F)
DimPlot(sub_mccs, pt.size=1, group.by="tricyclePhase",cols=c(binned_col[1],binned_col[14],binned_col[8]),split.by="treatment") +
        theme_void()+
        ggtitle("")
dev.off()

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_pseudotime.pdf",height=8,width=11,useDingbats = F)
FeaturePlot(sub_mccs, "pseudotime",pt.size=1) +
        scale_color_viridis(direction=-1,option = "C")+
        theme_void()+
        ggtitle("")+
        NoLegend()
dev.off()

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_pseudotime_withlegend.pdf",height=8,width=11,useDingbats = F)
FeaturePlot(sub_mccs, "pseudotime",pt.size=1) +
        scale_color_viridis(direction=-1,option = "C")+
        theme_void()+
        ggtitle("")
dev.off()

library(EnhancedVolcano)
clusters=c("0","1","2","3","4","5","6","8")

for(clus in clusters){
        clus1 =read.csv(paste0("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus",clus,".csv"),row.names=1)
        clus1_filt=clus1 %>% filter(!is.na(log2FoldChange))
        
        pdf(paste0("v1/Ribo_int_v1_clus",clus,"_enhancedvolcano.pdf"),height=6,width=6)
        print(EnhancedVolcano(clus1_filt, 
                              lab=clus1_filt$gene,
                              x="log2FoldChange",drawConnectors = F,
                              y="padj",labSize = 0,pointSize = 4,
                              pCutoff = .0005,FCcutoff = log2(1.5),
                              col=c("grey30","grey30","grey30","indianred3"),
                              title=paste0("Cluster ",clus),
                              subtitle=""))
        dev.off()
        
}

clus0 =read.csv(paste0("v1/Ribo_int_v1_DESeq_genes_Ribo_vs_DMSO_clus0.csv"),row.names=1)
clus0_filt=clus0 %>% filter(!is.na(log2FoldChange))
clus0_sig=clus0 %>% filter(log2FoldChange < log2(1.5) & padj <0.0005)
pdf(paste0("v1/Ribo_int_v1_clus0_labels_enhancedvolcano.pdf"),height=6,width=6)
print(EnhancedVolcano(clus0_filt, 
                      lab=clus0_filt$gene,
                      selectLab = c('Gmnc','Mycl','Mcidas','Ccno','Foxj1','Myb',"Gins2"),
                      x="log2FoldChange",drawConnectors = T,
                      y="padj",pointSize = 4,
                      pCutoff = .0005,FCcutoff = log2(1.5),
                      col=c("grey30","grey30","grey30","indianred3"),
                      title=paste0("Cluster 0"),
                      subtitle=""))
dev.off()

##### mccs only genes vs. pseudotime

sub_mccs=readRDS("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset.rds")
cds=readRDS("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_monocle.rds")
ps_df=as.data.frame(pseudotime(cds))

sub_mccs$pseudotime_submccs <-NA
sub_mccs$pseudotime_submccs[match(rownames(ps_df),rownames(sub_mccs@meta.data))]<-ps_df$`pseudotime(cds)`

meta=melt(sub_mccs@meta.data[,c("pseudotime_submccs","signature_1_UCell","signature_2_UCell","res0.1_names","treatment","orig.ident")],id.vars = c("pseudotime_submccs","treatment","res0.1_names","orig.ident"))

meta$score = "S"
meta$score[meta$variable == "signature_2_UCell"] <- "G2M"
meta$score=factor(meta$score,levels=c("S","G2M"))

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_int_s_g2m_scores_pseudotime_clusterids.pdf",height = 11,width=8,useDingbats = F)
ggplot(meta, aes(x=pseudotime_submccs,y=value,color=treatment))+
        geom_smooth()+
        facet_wrap(~score,ncol = 1,scales="free")+
        theme_classic()+
        ylab("Score")+
        xlab("Pseudotime")+
        geom_jitter(aes(x=pseudotime_submccs,y=0.1,color=res0.1_names),height=0.005)+
        scale_color_manual(values=c("black",basal_to_mcc_col[c(1,4,5,3)],"grey47"))
dev.off()

pdf("v1/mccs_int_subset/Ribo_int_v1_mccs_int_s_g2m_scores_pseudotime.pdf",height = 11,width=8,useDingbats = F)
ggplot(meta, aes(x=pseudotime_submccs,y=value,color=treatment))+
        geom_smooth()+
        facet_wrap(~score,ncol = 1,scales="free")+
        theme_classic()+
        ylab("Score")+
        xlab("Pseudotime")+
        scale_color_manual(values=c("grey30","indianred3"))
dev.off()


######## AUC scores ########

ggp_s=ggplot(meta %>% filter(score == "S"), aes(x=pseudotime_submccs,y=value,color=treatment))+
        geom_smooth()+
        theme_classic()+
        ylab("Score")+
        xlab("Pseudotime")+
        scale_color_manual(values=c("grey30","indianred3"))+
        facet_wrap(~orig.ident) 

ggp_data_s <- ggplot_build(ggp_s)$data[[1]]


s_aucs=NULL
for(x in 1:6){
        dat= ggp_data_s[ggp_data_s$PANEL == x,]
        area = MESS::auc(dat$x, dat$y)
        s_aucs = c(s_aucs, area)
}

ggp_g2m=ggplot(meta %>% filter(score == "G2M"), aes(x=pseudotime_submccs,y=value,color=treatment))+
        geom_smooth()+
        theme_classic()+
        ylab("Score")+
        xlab("Pseudotime")+
        scale_color_manual(values=c("grey30","indianred3"))+
        facet_wrap(~orig.ident) 

ggp_data_g2m <- ggplot_build(ggp_g2m)$data[[1]]


g2m_aucs=NULL
for(x in 1:6){
        dat= ggp_data_g2m[ggp_data_g2m$PANEL == x,]
        area = MESS::auc(dat$x, dat$y)
        g2m_aucs = c(g2m_aucs, area)
}

write.csv(s_aucs, file="v1/mccs_int_subset/Ribo_int_v1_mccs_int_s_score_vs_pseudotime_mcc_AUC.csv")
write.csv(g2m_aucs, file="v1/mccs_int_subset/Ribo_int_v1_mccs_int_g2m_score_vs_pseudotime_mcc_AUC.csv")

####### gene v pseudotime mccs int subset ######

DefaultAssay(sub_mccs) <- "RNA"
gene_dat=FetchData(sub_mccs, c("Gmnc","Mycl","Ccne1"))
#gene_dat=FetchData(seur, c("Myb","Foxj1"))

gene_dat=cbind(gene_dat, sub_mccs@meta.data[,c("treatment","res0.1_names","pseudotime_submccs")])

melted=reshape2::melt(gene_dat,id.vars=c("pseudotime_submccs","treatment","res0.1_names"))



pdf(paste0("v1/mccs_int_subset/Ribo_int_v1_mccs_int_gmnc_mycl_ccne1_v_pseudotime.pdf"),height=4,width=8,useDingbats = F)
ggplot(melted, aes(x=pseudotime_submccs,y=value,group=treatment,color=treatment))+
              geom_smooth()+
              theme_classic()+
              ylab("Expression level")+
              xlab("Pseudotime")+
        facet_wrap(~variable,scales = "free")+
        theme_classic()+
        scale_color_manual(values=c("grey30","indianred3"))
dev.off()

pdf(paste0("v1/mccs_int_subset/Ribo_int_v1_mccs_int_gmnc_mycl_ccne1_v_pseudotime_clusterids.pdf"),height=4,width=8,useDingbats = F)
ggplot(melted, aes(x=pseudotime_submccs,y=value,group=treatment,color=treatment))+
        geom_smooth()+
        theme_classic()+
        ylab("Expression level")+
        xlab("Pseudotime")+
        facet_wrap(~variable,scales = "free")+
        theme_classic()+
        geom_jitter(aes(x=pseudotime_submccs, y=1,color=res0.1_names),height=0.05)+
        scale_color_manual(values=c("grey30",basal_to_mcc_col[c(1,4,5,3)],"indianred3"))
dev.off()


### bin over pseudotime

#sub_mccs$pseudo_submccs_bin <- cut_interval(sub_mccs$pseudotime_submccs,n=10,labels=F)
sub_mccs$pseudo_submccs_bin <- cut_number(sub_mccs$pseudotime_submccs,n=20,labels=F)
sub_mccs$pseudo_bin_treat <- paste(sub_mccs$treatment,sub_mccs$pseudo_submccs_bin,sep="_")
Idents(sub_mccs) <- "pseudo_bin_treat"
avg_exp=AverageExpression(sub_mccs, assays = "RNA",features=c("Gmnc","Mycl","Ccne1"))

dat=melt(avg_exp$RNA)
dat$treatment <- "DMSO"
dat$treatment[grep("Ribo*",dat$Var2)]<-"Ribo"
dat$bin <- str_replace(dat$Var2,"DMSO_","")
dat$bin <- str_replace(dat$bin,"Ribo_","")
dat$bin = as.numeric(dat$bin)

ggplot(dat, aes(x=bin,y=value,color=treatment,group=treatment))+
        geom_line()+
        facet_wrap(~Var1, scales="free")
####### gene v pseudotime full dataset ######


DefaultAssay(seur) <- "RNA"
gene_dat=FetchData(seur,unique(c(g1,s)))
#gene_dat=FetchData(seur, c("Myb","Foxj1"))

gene_dat=cbind(gene_dat, seur@meta.data[,c("treatment","integrated_snn_res.0.3","pseudotime_mcc")])

gene_dat_filt=gene_dat %>% filter(!is.na(pseudotime_mcc))
melted=reshape2::melt(gene_dat_filt,id.vars=c("pseudotime_mcc","treatment","integrated_snn_res.0.3"))

dir.create("v1/gene_v_pseudotime")
#for(x in c("Myb","Foxj1")){ 
for(x in unique(c(g1,s))){
        pdf(paste0("v1/gene_v_pseudotime/Ribo_int_v1_",x,"_v_pseudotime.pdf"),height=4,width=5.5,useDingbats = F)
        print(ggplot(melted %>% filter(variable==x), aes(x=pseudotime_mcc,y=value,group=treatment,color=treatment))+
                      geom_smooth()+
                      theme_classic()+
                      ggtitle(x)+
                      ylab("Expression level")+
                      xlab("Pseudotime"))
        dev.off()
}

### plot average expression of replicates from DESEQ hits

seur$cluster_rep = paste(seur$integrated_snn_res.0.3,seur$orig.ident,sep="_")

dat=AverageExpression(seur, assays="RNA",features=c("Mycl","Gmnc","Mcidas","Myb","E2f7"),group.by="cluster_rep")

library(pheatmap)
PurpleAndYellow()

pdf("v1/Ribo_int_v1_clus0_heatmap.pdf",height=6,width=6,useDingbats = F)
pheatmap(dat$RNA[,grep("^0_*",colnames(dat$RNA))],scale="row",color = PurpleAndYellow(50))
dev.off()

pdf("v1/Ribo_int_v1_clus0_heatmap.pdf",height=6,width=6,useDingbats = F)
pheatmap(dat$RNA[,grep("^6_*",colnames(dat$RNA))],scale="row",color = PurpleAndYellow(50))
dev.off()
                 
DefaultAssay(seur) <- "RNA"
gene_list=c("Mycl","Ccne1","Gmnc","Mcidas","Myb","E2f7")
for(gene in gene_list){
        pdf(paste0("v1/Ribo_int_v1_",gene,"_split.pdf"),height=4,width=11,useDingbats = F)
        print(FeaturePlot(seur,pt.size=0.5, gene,split.by="treatment",order=T,cols = c("grey17","cyan")) & theme_void())
        dev.off()
}


### GEO METADATA ####
DimPlot(seur, group.by="integrated_snn_res.0.3",label=T)
DimPlot(seur, group.by="LT_clusters",label=T)
meta=seur@meta.data
meta$cluster_names=meta$integrated_snn_res.0.3
meta=meta %>% mutate(cluster_names=recode(cluster_names,
                                          `1`="Basal stem",
                                          `2`="Secretory",
                                          `3`="Multiciliated 2",
                                          `4`="Tnfrsf12a+",
                                          `5`="Nupr1+",
                                          `6`="Multiciliated 1",
                                          `7`="Proliferating",
                                          `8`="Neuroendocrine",
                                          `9`="Multiciliated 3",
                                          `0`="Intermediate",
                                          `10`="Krt14+ basal stem"
))

seur@meta.data=meta
DimPlot(seur, group.by = "cluster_names")

full_umap =Embeddings(seur, "umap")

seur$UMAP_1<-full_umap[,c("UMAP_1")]
seur$UMAP_2<-full_umap[,c("UMAP_2")]

columns_to_keep=c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","UMAP_1","UMAP_2","LT_clusters","integrated_snn_res.0.3","cluster_names","tricyclePosition","signature_1_UCell","signature_2_UCell","pseudotime","pseudotime_mcc","pseudotime_sec")

for_printing=seur@meta.data[,columns_to_keep]
colnames(for_printing)<-c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","UMAP_1","UMAP_2","labeltransfer_clusters","integrated_snn_res.0.3","cluster_names","tricyclePosition","S_score_UCell","G2M_score_UCell","pseudotime","pseudotime_mcc","pseudotime_sec")
write.csv(for_printing, file="v1/GEO/Ribo_int_metadata.csv")
library(Matrix)

writeMM(seur[["RNA"]]@counts, file="v1/GEO/Ribo_int_counts.mtx")


#### MARKERS ####

DefaultAssay(seur) <- "RNA"
Idents(seur) <- "cluster_names"
marks=FindAllMarkers(seur, assay="RNA",only.pos = T)

marks_ord=marks %>% group_by(cluster) %>% arrange(desc(avg_log2FC),.by_group = T)

write.csv(marks_ord, file="v1/Ribo_int_v1_res0.3_cluster_names_markers.csv")

############# CELL x GENE ##################

## add ontology terms
new_meta$assay_ontology_term_id <- "EFO:0009922"
new_meta$cell_type_ontology_term_id <- "CL:0000307"
new_meta$development_stage_ontology_term_id <-"MmusDv:0000110"
new_meta$disease_ontology_term_id <- "PATO:0000461"
new_meta$donor_id <- "na"
new_meta$is_primary_data <- TRUE
new_meta$organism_ontology_term_id <-"NCBITaxon:10090"
new_meta$self_reported_ethnicity_ontology_term_id <- "na"
new_meta$sex_ontology_term_id <- "unknown"
new_meta$suspension_type <- "cell"
new_meta$tissue_type <- "cell culture"
new_meta$tissue_ontology_term_id <- "UBERON:0001901"

seur_temp=seur

library(sceasy)
library(reticulate)
use_condaenv('sceasy')
dir.create("v1/cellxgene")
seur_temp@meta.data = new_meta
sceasy::convertFormat(seur_temp, assay='RNA', from="seurat", to="anndata", main_layer='data', transfer_layers='counts', drop_single_values=FALSE, outFile='v1/cellxgene/Ribo_AD03_full.h5ad')

## print out raw cellranger counts for cell x gene
library(Matrix)
meta_final=read.csv("/Volumes/ReiterLab4/Ribo_AD03/Ribo_int/v1/GEO/Ribo_int_metadata.csv",row.names=1)

dmsoA_names = rownames(meta_final)[meta_final$orig.ident == "DMSOA"]
dmsoB_names=rownames(meta_final)[meta_final$orig.ident == "DMSOB"]
dmsoC_names=rownames(meta_final)[meta_final$orig.ident == "DMSOC"]
riboA_names=rownames(meta_final)[meta_final$orig.ident == "RiboA"]
riboB_names=rownames(meta_final)[meta_final$orig.ident == "RiboB"]
riboC_names=rownames(meta_final)[meta_final$orig.ident == "RiboC"]

dmsoA=ReadMtx("/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/DMSOA_matrix.mtx.gz",features="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/DMSOA_features.tsv.gz",cells="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/DMSOA_barcodes.tsv.gz")
colnames(dmsoA)=paste0(colnames(dmsoA),"_1")
dmsoA_sub = dmsoA[,match(dmsoA_names,colnames(dmsoA ))]

dmsoB=ReadMtx("/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/DMSOB_matrix.mtx.gz",features="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/DMSOB_features.tsv.gz",cells="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/DMSOB_barcodes.tsv.gz")
colnames(dmsoB)=paste0(colnames(dmsoB),"_2")
dmsoB_sub = dmsoB[,match(dmsoB_names,colnames(dmsoB))]

dmsoC=ReadMtx("/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/DMSOC_matrix.mtx.gz",features="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/DMSOC_features.tsv.gz",cells="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/DMSOC_barcodes.tsv.gz")
colnames(dmsoC)=paste0(colnames(dmsoC),"_3")
dmsoC_sub = dmsoC[,match(dmsoC_names,colnames(dmsoC))]

riboA=ReadMtx("/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/RiboA_matrix.mtx.gz",features="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/RiboA_features.tsv.gz",cells="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/RiboA_barcodes.tsv.gz")
colnames(riboA)=paste0(colnames(riboA),"_4")
riboA_sub = riboA[,match(riboA_names,colnames(riboA))]

riboB=ReadMtx("/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/RiboB_matrix.mtx.gz",features="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/RiboB_features.tsv.gz",cells="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/RiboB_barcodes.tsv.gz")
colnames(riboB)=paste0(colnames(riboB),"_5")
riboB_sub = riboB[,match(riboB_names,colnames(riboB))]

riboC=ReadMtx("/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/RiboC_matrix.mtx.gz",features="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/RiboC_features.tsv.gz",cells="/Volumes/Reiterlab4/Ribo_AD03/Ribo_int/v1/GEO/RiboC_barcodes.tsv.gz")
colnames(riboC)=paste0(colnames(riboC),"_6")
riboC_sub = riboC[,match(riboC_names,colnames(riboC))]

full=cbind(dmsoA_sub,dmsoB_sub,dmsoC_sub,riboA_sub,riboB_sub,riboC_sub)

write(rownames(full), file="v1/cellxgene/Ribo_AD03_full_combined_genes.tsv")
write(colnames(full), file="v1/cellxgene/Ribo_AD03_full_combined_barcodes.tsv")
identical(colnames(full),rownames(meta_final))

writeMM(full, "v1/cellxgene/Ribo_AD03_full_combined_rawCR_counts.mtx")

### focused dataset
sub_mccs=readRDS("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset.rds")
DimPlot(sub_mccs, group.by = "res0.2_names")

full_umap =Embeddings(sub_mccs, "umap")
full_umap[,c("UMAP_2")] <-full_umap[,c("UMAP_2")]*-1

ggplot(as.data.frame(full_umap), aes(x=UMAP_1,y=UMAP_2))+
        geom_point()

sub_mccs$UMAP_1<-full_umap[,c("UMAP_1")]
sub_mccs$UMAP_2<-full_umap[,c("UMAP_2")]

cds=readRDS("v1/mccs_int_subset/Ribo_int_v1_mccs_int_subset_monocle.rds")
ps_df=as.data.frame(pseudotime(cds))

sub_mccs$pseudotime_submccs <-NA
sub_mccs$pseudotime_submccs[match(rownames(ps_df),rownames(sub_mccs@meta.data))]<-ps_df$`pseudotime(cds)`

columns_to_keep=c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","UMAP_1","UMAP_2","res0.1_names","treatment","tricyclePosition","tricyclePhase","signature_1_UCell","signature_2_UCell","pseudotime_submccs")

for_printing=sub_mccs@meta.data[,columns_to_keep]
colnames(for_printing)<-c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","UMAP_1","UMAP_2","cluster_names","treatment","tricyclePosition","tricyclePhase" ,"S_score_UCell","G2M_score_UCell","pseudotime")

new_meta=for_printing
new_meta$assay_ontology_term_id <- "EFO:0009922"
new_meta$cell_type_ontology_term_id <- "CL:0000307"
new_meta$development_stage_ontology_term_id <-"MmusDv:0000110"
new_meta$disease_ontology_term_id <- "PATO:0000461"
new_meta$donor_id <- "na"
new_meta$is_primary_data <- FALSE
new_meta$organism_ontology_term_id <-"NCBITaxon:10090"
new_meta$self_reported_ethnicity_ontology_term_id <- "na"
new_meta$sex_ontology_term_id <- "unknown"
new_meta$suspension_type <- "cell"
new_meta$tissue_type <- "cell culture"
new_meta$tissue_ontology_term_id <- "UBERON:0001901"

sub_temp=sub_mccs

library(sceasy)
library(reticulate)
use_condaenv('sceasy')

sub_temp@meta.data = new_meta
sceasy::convertFormat(sub_temp, assay='RNA', from="seurat", to="anndata", main_layer='data', transfer_layers='counts', drop_single_values=FALSE, outFile='v1/cellxgene/Ribo_AD03_focused.h5ad')

