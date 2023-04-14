
library(Seurat)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
print(args)
obj_list=c()

for(i in 1:(length(args)-2)){
	print(args[i])
	obj = readRDS(args[i])
	obj_list = c(obj_list,obj)
}
print(obj_list)
folder=args[length(args)-1]
version = args[length(args)]


setwd("/wynton/group/reiter/lauren/")
dir.create(folder)
dir.create(paste0(folder,"/",version,"/"))
setwd(paste0(folder,"/",version,"/"))

write.table(args, file=paste0(folder,"_list_of_args.txt"))


for (i in 1:length(obj_list)) {
    obj_list[[i]] <- NormalizeData(obj_list[[i]], verbose = FALSE)
    obj_list[[i]] <- FindVariableFeatures(obj_list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}

anchors = FindIntegrationAnchors(obj_list, dims=1:30)
integrated = IntegrateData(anchors, dims=1:30)

DefaultAssay(integrated)<- "integrated"

integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs=40)
integrated <- RunUMAP(integrated, reduction="pca",dims=1:20)

integrated = FindNeighbors(integrated, dims=1:20)
integrated = FindClusters(integrated, resolution = 0.3)

pdf(paste0(folder,"_",version,"_clusters.pdf"),height=8.5,width=11)
DimPlot(integrated, pt.size=2,label=T)
dev.off()

pdf(paste0(folder,"_",version,"_batch.pdf"),height=8.5,width=11)
DimPlot(integrated, pt.size=2,group.by="orig.ident",label=T)
dev.off()

saveRDS(integrated, file=paste0(folder,"_",version,".rds"))
