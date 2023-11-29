
###################### NEW CUT&RUN CR57_E2f7 #################################

### peak list is also Supplementary Table S8

setwd("/Volumes/Reiterlab_3/CR57_E2f7/")
e2f7_peaks=read.csv("E2f7_A_vs_nlsGFP_A_dupmark_peaks_vs_E2f7_B_vs_nlsGFP_B_dupmark_peaks/common_unique_E2f7_A_vs_nlsGFP_A_dupmark_peaks_E2f7_B_vs_nlsGFP_B_dupmark_peaks_chipseeker_peakannotations.csv",row.names=1)

e2f7_common= e2f7_peaks[e2f7_peaks$V11 == "common",]

## any gene near e2f7 peaks

e2f7_genes=unique(e2f7_common$SYMBOL)

e2f7_flank=e2f7_common$gene_flank_symbol

library(stringr)

flanks=unlist(strsplit(e2f7_flank,split =" ; "))
flanks_clean=unique(flanks[!is.na(flanks)])

total_genes=unique(c(e2f7_genes,flanks_clean))


west_peaks=read.csv("~/Box Sync/E2f7_paper/E2f7_peaks_Westendorp2012.csv")
inter_west=intersect(toupper(total_genes),west_peaks$Gene.symbol)

###### try with only promoter peaks (within 5 kb of TSS)

promoter_e2f7=e2f7_common[abs(e2f7_common$distanceToTSS) < 5000,]

prom_genes=unique(promoter_e2f7$SYMBOL)
#prom_flank=unlist(strsplit(promoter_e2f7$gene_flank_symbol,split =" ; "))
#prom_flank_clean=unique(prom_flank[!is.na(prom_flank)])

inter_west_prom=intersect(toupper(prom_genes),west_peaks$Gene.symbol)


#### eluer diagrams

library(eulerr)

fit1 <- euler(c("CUTRUN" = 95, "Westendorp" = 746,
                "CUTRUN&Westendorp" = 27))

error_plot(fit1)
plot(fit1,quantities = T,main = "All intersecting genes")


fit2 <- euler(c("CUTRUN" = length(prom_genes), "Westendorp" = 746,
                "CUTRUN&Westendorp" = length(inter_west_prom)))
plot(fit2,quantities = T,main = "CUT&RUN promoter intersecting genes")

#### WESTENDORP TOP GENES ###

west_top=read.csv("~/Box Sync/E2f7_paper/Westendorp_E2f7_toptarget.csv",header=F)

intersect(west_top$V1,toupper(total_genes))
intersect(west_top$V1,toupper(prom_genes))

fit3= euler(c("CUTRUN" = length(prom_genes), "Westendorp" = length(west_top$V1),
              "CUTRUN&Westendorp" = length(intersect(west_top$V1,toupper(prom_genes)))))

plot(fit3,quantities = T)

### fold change in b/c clusters
#### DE genes between E2f7 WT and MUT is now updated Supplementary Table7
deseq=read.csv("~/Box Sync/E2f7_paper/Supplementary_Table5.csv")
deseq_a = deseq[,c(1,grep("clus_a_*",colnames(deseq)))]
deseq_a_filt = deseq_a %>% filter(clus_a_padj < .05 & abs(clus_a_log2FoldChange) > log2(1.5))

deseq_b = deseq[,c(1,grep("clus_b_*",colnames(deseq)))]
deseq_b_filt = deseq_b %>% filter(clus_b_padj < .05 & abs(clus_b_log2FoldChange) > log2(1.5))

deseq_c = deseq[,c(1,grep("clus_c_*",colnames(deseq)))]
deseq_c_filt = deseq_c %>% filter(clus_c_padj < .05 & abs(clus_c_log2FoldChange) > log2(1.5))

deseq_d = deseq[,c(1,grep("clus_d_*",colnames(deseq)))]
deseq_d_filt = deseq_d %>% filter(clus_d_padj < .05 & abs(clus_d_log2FoldChange) > log2(1.5))

deseq$ClusterA_DE <- "N"
deseq$ClusterA_DE[match(deseq_a_filt$Column1,deseq$Column1)] <- "Y"               
deseq$ClusterB_DE <- "N"
deseq$ClusterB_DE[match(deseq_b_filt$Column1,deseq$Column1)] <- "Y"    

deseq$ClusterC_DE <- "N"
deseq$ClusterC_DE[match(deseq_c_filt$Column1,deseq$Column1)] <- "Y"    

deseq$ClusterD_DE <- "N"
deseq$ClusterD_DE[match(deseq_d_filt$Column1,deseq$Column1)] <- "Y"  

equal_Y <- function(x){
        x == "Y"
}

deseq_thresh=deseq %>% filter(if_any(contains("_DE"), equal_Y))

write.csv(deseq_thresh,"~/Box Sync/E2f7_paper/E2f7_DE_0.05padj_1.5FC_cluster_thresh.csv")

deseq_thresh$Westendorp <- NA
deseq_thresh$Westendorp[na.omit(match(west_peaks$Gene.symbol, toupper(deseq_thresh$Column1)))]<-"Westendorp peak"
deseq_thresh$Westendorp[na.omit(match(west_top$V1,toupper(deseq_thresh$Column1)))]<-"Westendorp TOP"

deseq_thresh$CUTRUN <- NA
deseq_thresh$CUTRUN[na.omit(match(total_genes,deseq_thresh$Column1))]<-"CUT&RUN peak"
deseq_thresh$CUTRUN[na.omit(match(prom_genes,deseq_thresh$Column1))]<-"CUT&RUN promoter peak"

table(deseq_thresh$Westendorp,deseq_thresh$CUTRUN)

## add old CUT&RUN
old_CR=read.csv("/Volumes/Reiterlab_3/CR45/CR45/CR45_E2f7GFP_A_vs_CR45_nlsGFP_A_dupmark_q05_peaks_vs_CR45_E2f7GFP_B_vs_CR45_nlsGFP_B_dupmark_q05_peaks/common_unique_CR45_E2f7GFP_A_vs_CR45_nlsGFP_A_dupmark_q05_peaks_CR45_E2f7GFP_B_vs_CR45_nlsGFP_B_dupmark_q05_peaks_chipseeker_peakannotations.csv")

old_common= old_CR[old_CR$V11 == "common",]
old_genes=unique(old_common$SYMBOL)

old_flank=old_common$gene_flank_symbol

library(stringr)

old_flanks=unlist(strsplit(old_flank,split =" ; "))
old_flanks_clean=unique(old_flanks[!is.na(old_flanks)])

old_total_genes=unique(c(old_genes,old_flanks_clean))

old_promoter_e2f7=old_common[abs(old_common$distanceToTSS) < 5000,]

old_prom_genes=unique(old_promoter_e2f7$SYMBOL)

deseq_thresh$OLD_CUTRUN <- NA
deseq_thresh$OLD_CUTRUN[na.omit(match(old_total_genes,deseq_thresh$Column1))]<-"CUT&RUN peak"
deseq_thresh$OLD_CUTRUN[na.omit(match(old_prom_genes,deseq_thresh$Column1))]<-"CUT&RUN promoter peak"

final_deseq=deseq_thresh %>% select(-contains("baseMean")) %>% select(-contains("lfcSE")) %>% select(-contains("svalue"))

write.csv(final_deseq,"~/Box Sync/E2f7_paper/E2f7_DESeq_padj0.05_FC1.5_new_oldCUTRUN_WESTENDORP_genelist.csv")
final_deseq=read.csv("~/Box Sync/E2f7_paper/E2f7_DESeq_padj0.05_FC1.5_new_oldCUTRUN_WESTENDORP_genelist.csv",row.names=1)

### add Westendorp microarray data

west_de=read.csv("~/Box Sync/E2f7_paper/westendorp_microarray_gProfiler_hsapiens_mmusculus_10-25-2023_1-50-52 PM.csv")
west_de_clean=west_de[!(west_de$ortholog_name == "N/A"),]


### also make list without DE thresholding (to see all peaks)
deseq$Westendorp <- NA
deseq$Westendorp[na.omit(match(west_peaks$Gene.symbol, toupper(deseq$Column1)))]<-"Westendorp peak"
deseq$Westendorp[na.omit(match(west_top$V1,toupper(deseq$Column1)))]<-"Westendorp TOP"

deseq$CUTRUN <- NA
deseq$CUTRUN[na.omit(match(total_genes,deseq$Column1))]<-"CUT&RUN peak"
deseq$CUTRUN[na.omit(match(prom_genes,deseq$Column1))]<-"CUT&RUN promoter peak"

deseq$OLD_CUTRUN <- NA
deseq$OLD_CUTRUN[na.omit(match(old_total_genes,deseq$Column1))]<-"CUT&RUN peak"
deseq$OLD_CUTRUN[na.omit(match(old_prom_genes,deseq$Column1))]<-"CUT&RUN promoter peak"

final_nothresh=deseq %>% select(-contains("baseMean")) %>% select(-contains("lfcSE")) %>% select(-contains("svalue"))

write.csv(final_nothresh, file="~/Box Sync/E2f7_paper/E2f7_DESeq_allgenes_new_oldCUTRUN_WESTENDORP_genelist.csv")

#### CUT&RUN figures #####

##peak distribution

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(reshape2)
library(ggplot2)
txdb=TxDb.Mmusculus.UCSC.mm10.knownGene

makeGRangesFromDataFrame(e2f7_common)

peakAnnoList=annotatePeak(makeGRangesFromDataFrame(e2f7_common), tssRegion = c(-5000,5000),TxDb=txdb)

annopeaks=peakAnnoList@anno

annopeaks$peak_distr <- annopeaks$annotation
annopeaks$peak_distr[grep("Promoter*",annopeaks$annotation)]<- "Promoter"
annopeaks$peak_distr[grep("Intron*",annopeaks$annotation)]<- "Intron"
annopeaks$peak_distr[grep("Exon*",annopeaks$annotation)]<- "Exon"

anno_tab = melt(prop.table(table(annopeaks$peak_distr))*100)

pdf("~/Box Sync/E2f7_paper/new_CUTRUN_E2f7_peak_distribution.pdf",height=6,width=4,useDingbats = F)
ggplot(anno_tab,aes(x=1,y=value,fill=Var1))+geom_col() + 
        theme_classic() + 
        ylab("% of Peaks")+ 
        scale_fill_manual(values=rev(c("skyblue","plum","indianred","darkgreen","darkgoldenrod")))+
        theme(axis.title.x  = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_text(size=14),axis.text.y=element_text(size=12))
dev.off()

library(eulerr)

## compare promoter peaks of our CUT&RUN with westendorp peaks by promoters

table(deseq$CUTRUN)
table(deseq$Westendorp)
fit1 <- euler(c("CUTRUN" = length(prom_genes), "Westendorp" = length(west_peaks$Gene.symbol),
                "CUTRUN&Westendorp" = length(intersect(toupper(prom_genes),west_peaks$Gene.symbol))))

error_plot(fit1)
plot(fit1,quantities = T,main = "Genes within 5 kb of promoter")

fit2 <- euler(c("CUTRUN" = length(prom_genes), "Westendorp" = length(west_top$V1),
                "CUTRUN&Westendorp" = length(intersect(toupper(prom_genes),west_top$V1))))

plot(fit2,quantities = T,main = "Genes in Westendorp TOP peak list")

x=list(CUTRUN=toupper(prom_genes),
       Westendorp=west_peaks$Gene.symbol)

write.csv(prom_genes, file="CUTRUN_Westendorp_venn.csv")

y=list(CUTRUN=toupper(prom_genes),
       Westendorp=west_top$V1)

write.csv(west_top$V1, file="Westendorp_top.csv")

z=list(CUTRUN=prom_genes,
       DESEQ=final_deseq$Column1)


deseq_bc=final_deseq %>% filter(ClusterB_DE == "Y" | ClusterC_DE == "Y")
a=list(CUTRUN=prom_genes,
       DESEQ=deseq_bc$Column1)

write.csv(deseq_bc$Column1, file="DE_E2f7_clusterb_c_venn.csv")

b=list(DESEQ=toupper(deseq_bc$Column1), 
       WestendorpTOP=west_top$V1)

c=list(DESEQ=deseq_bc$Column1,
       WestendorpMA=west_de_clean$ortholog_name)

library(VennDiagram) 
venn.diagram(x, fill = c("skyblue", "gray67"), 
             alpha = c(0.5, 0.5), lwd =0, "~/Box Sync/E2f7_paper/CUTRUN_Westendorp_promoter.tiff")

venn.diagram(y, fill = c("skyblue", "gray67"), 
             alpha = c(0.5, 0.5), lwd =0, "~/Box Sync/E2f7_paper/CUTRUN_Westendorp_top.tiff")

venn.diagram(z,fill = c("skyblue", "gray67"), 
             alpha = c(0.5, 0.5), lwd =0, "~/Box Sync/E2f7_paper/CUTRUN_Deseq.tiff")

venn.diagram(a,fill = c("skyblue", "gray67"), 
             alpha = c(0.5, 0.5), lwd =0, "~/Box Sync/E2f7_paper/CUTRUN_Deseq_clusBC.tiff")

venn.diagram(b,fill = c("skyblue", "gray67"), 
             alpha = c(0.5, 0.5), lwd =0, "~/Box Sync/E2f7_paper/Deseq_clusBC_Westendorp_top.tiff")

venn.diagram(c,fill = c("skyblue", "gray67"), 
             alpha = c(0.5, 0.5), lwd =0, "~/Box Sync/E2f7_paper/Deseq_clusBC_Westendorp_microarray.tiff")


##### tracks
setwd("/Volumes/Reiterlab_3/CR57_E2f7/")
library(trackplot)
library(data.table)
library(ChIPseeker)

## picked out regions I want to display
regions=read.csv("E2f7_CUTRUNpromoter_genelist_bc.csv",header=F)

e2f7_a="/Volumes/Reiterlab_3/CR57_E2f7/E2f7_A/dup.marked.120bp/E2f7_A_henikoff_dupmark_120bp_RPGC.bw"
e2f7_b="/Volumes/Reiterlab_3/CR57_E2f7/E2f7_B/dup.marked.120bp/E2f7_B_henikoff_dupmark_120bp_RPGC.bw"
nls_a="/Volumes/Reiterlab_3/CR57_E2f7/nlsGFP_A/dup.marked.120bp/nlsGFP_A_henikoff_dupmark_120bp_RPGC.bw"
nls_b="/Volumes/Reiterlab_3/CR57_E2f7/nlsGFP_B/dup.marked.120bp/nlsGFP_B_henikoff_dupmark_120bp_RPGC.bw"

bigWigs = c(e2f7_a,e2f7_b,nls_a,nls_b)

for(i in 1:length(regions$V1)){
        track_data = track_extract(bigWigs = bigWigs, loci = regions$V2[i])
        #track_sum=track_summarize(track_data, c("E2f7GFP","E2f7GFP","nlsGFP","nlsGFP"))
        
        pdf(paste0("CR57_",regions$V1[i],"replicates.pdf"),height=6,width=11,useDingbats = F)
        track_plot(summary_list = track_data, draw_gene_track = T, gene_model="~/Documents/Reiter_Seq/mm10.refGene.gtf.gz", isGTF = T,groupScaleByCondition = T ,gene_fsize = 2,build = "mm10")
        dev.off()
}