library(Seurat)
#normal_merged data
AT2_1000<-subset_AT2[,sample(ncol(subset_AT2),1000)] #ncol(subset_AT2)=16609
AT1_1000<-subset_AT1[,sample(ncol(subset_AT1),1000)] #ncol(subset_AT1)=1029
Basal_cell_1000<-subset_Basal_cell[,sample(ncol(subset_Basal_cell),1000)] #ncol(subset_Basal_cell)= 5350
normal_merged<-merge(x=AT2_1000,y=list(AT1_1000,subset_Club_cell,Basal_cell_1000,subset_Ciliated_cell)) 
normal_merged%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
normal_merged@meta.data[,"protocol"]<-"Normal"

#Seurat v3 ; data integration & alignment
test<-merge(x=normal_merged,y=LUAD_merged)
normal_LUAD.list<-SplitObject(object = test,split.by = "protocol")
for (i in 1:length(x=normal_LUAD.list)) {
  normal_LUAD.list[[i]]<-NormalizeData(object=normal_LUAD.list[[i]],verbose = FALSE)
  normal_LUAD.list[[i]]<-FindVariableFeatures(object = normal_LUAD.list[[i]],selection.method = "vst", 
      nfeatures=2000, verbose = FALSE)
}

reference.list<-normal_LUAD.list[c("Normal","LUAD")]
data.anchors<-FindIntegrationAnchors(object.list = normal_LUAD.list,dims=1:30,k.filter = 100) #문제 k.filter를 default 200에서 100으로 낮춤
data.integrated<-IntegrateData(anchorset = data.anchors, dims= 1:30)

#Visualization
library(ggplot2)
library(cowplot)
DefaultAssay(object = data.integrated)<-"integrated"
data.integrated<-ScaleData(object = data.integrated, verbose = FALSE)
data.integrated<-RunPCA(object = data.integrated, npcs = 30, verbose= FALSE)
data.integrated<-RunUMAP(object = data.integrated, reduction = "pca", dims = 1:30) #pip install 필수 

TSNE1_normal_merged<-DimPlot(object=data.integrated, reduction = "pca", group.by = "protocol")
TSNE2_normal_merged<-DimPlot(object=data.integrated, reduction = "pca", group.by = "protocol")
UMAP1_normal_merged_LUAD<-DimPlot(object=data.integrated, reduction = "umap", group.by = "protocol")
UMAP2_normal_merged_LUAD<-DimPlot(object=data.integrated, reduction = "umap", group.by = "Clustername_cell",label = TURE,repel = TRUE)+NoLegend()
plot_grid(p1,p2)








