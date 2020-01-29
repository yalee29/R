#data upload (H358/LC_PT_45/LC_PT_45_Re/LUAD_total)
#total file은 bulk seq까지 다 포함됨.
library(Seurat)
library(dplyr) 
H358_total<-read.csv("./LUAD/H358_bulk_sc.csv",header=TRUE)
LUAD_total<-read.csv("./LUAD/GSE69405_PROCESSED_GENE_TPM_ALL.csv",header=TRUE)
LC_PT_45_total<-read.csv("./LUAD/LC-PT-45_bulk_sc.csv",header=TRUE)
LC_PT_45_Re_total<-read.csv("./LUAD/LC-PT-45_Re_bulk_sc.csv")

#SC data 모으기(bulk seq data 제외); H358_sc & LC_PT_sc
#gene_name column은 character로 지정 
H358_sc<-H358_total%>%select(-gene_id,-H358_Pooled)
LC_PT_45_sc<-LC_PT_45_total%>%select(-gene_id,-.PT.45_Pooled)
LC_PT_45_Re_sc<-LC_PT_45_Re_total%>%select(-gene_id,-LC.PT.45.Re_Pooled)
H358_sc$gene_name<-as.character(H358_sc$gene_name)
LC_PT_45_sc$gene_name<-as.character(LC_PT_45_Re_sc$gene_name)
LC_PT_45_Re_sc$gene_name<-as.character(LC_PT_45_Re_sc$gene_name)

#중복되는 gene name은 더해서 하나로 합치기 
#my_collapse_df(df,by) 함수 생성
my_collapse_df <- function(df, by) {
  
  by <- unique(by)
  
  df_return <- data.frame(gene_names=by)
  df_return[,colnames(df)[-1]] <- NA
  
  for(g in by) {
    print(g)
    df_temp  <- df[df$gene_name==g, ]
    df_return[df_return$gene_name==g, 2:ncol(df_return)] <- colSums(df_temp[,-1])
    
  }
  
  return(df_return)
}

H358_sc_merged <- my_collapse_df(H358_sc, H358_sc$gene_name)
LC_PT_45_sc_merged<-my_collapse_df(LC_PT_45_sc,LC_PT_45_sc$gene_name)
LC_PT_45_Re_sc_merged<-my_collapse_df(LC_PT_45_Re_sc,LC_PT_45_Re_sc$gene_name)

#gene_name열을 rowname으로 지정
H358_sc_merged$gene_names<-as.character(H358_sc_merged$gene_names)
rownames(H358_sc_merged)<-H358_sc_merged[,1]
H358_sc_merged<-select(H358_sc_merged,-gene_names)
LC_PT_45_sc_merged$gene_names<-as.character(LC_PT_45_sc_merged$gene_names)
rownames(LC_PT_45_sc_merged)<-LC_PT_45_sc_merged[,1]
LC_PT_45_sc_merged<-select(LC_PT_45_sc_merged,-gene_names)
LC_PT_45_Re_sc_merged$gene_names<-as.character(LC_PT_45_Re_sc_merged$gene_names)
rownames(LC_PT_45_Re_sc_merged)<-LC_PT_45_Re_sc_merged[,1]
LC_PT_45_Re_sc_merged<-select(LC_PT_45_Re_sc_merged,-gene_names)

rownames(H358_sc_merged)<-make.unique(rownames(H358_sc_merged))
rownames(LC_PT_45_sc_merged)<-make.unique(rownames(LC_PT_45_sc_merged))
rownames(LC_PT_45_Re_sc_merged)<-make.unique(rownames(LC_PT_45_Re_sc_merged))

#Set up Seurat Object for cancer data
H358_sc_merged<-CreateSeuratObject(count=H358_sc_merged)
LC_PT_45_sc_merged<-CreateSeuratObject(count=LC_PT_45_sc_merged)
LC_PT_45_Re_sc_merged<-CreateSeuratObject(count=LC_PT_45_Re_sc_merged)
LC_PT_sc_merged<-merge(x=LC_PT_45_sc_merged,y=LC_PT_45_Re_sc_merged)
LUAD_merged<-merge(x=H358_sc_merged,y=list(LC_PT_45_sc_merged,LC_PT_45_Re_sc_merged))

LC_PT_sc_merged%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
LUAD_merged%>%NormalizeData%>%ScaleData%>%FindVariableFeatures

#Set up normal cell data according to cell type(-immune cell,endothelial ell,fibroblast)/ Immune cell data 
data.frame(donors@assays$RNA@data) #donors data 확인 (data.frame 형태로!)
subset_donors <- SubsetData(donors, cells=names(Idents(donors))[as.character(Idents(donors)) %in% c("Basal cell", "AT2", "AT1","Ciliated cell", "Club cell")])
subset_AT2<-subset(x=donors,idents = "AT2")
subset_AT1<-subset(x=donors,idents="AT1")
subset_Club_cell<-subset(x=donors,idents = "Club cell")
subset_Basal_cell<-subset(x=donors,idents="Basal cell")
subset_Ciliated_cell<-subset(x=donors,idents="Ciliated cell")
subset_Immune_cell<-subset(x=donors,idents=c("monocyte","Dendritic cell","B cell","Alveolar macrophage"))
subset_donors%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
subset_AT2%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
subset_AT1%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
subset_Club_cell%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
subset_Basal_cell%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
subset_Ciliated_cell%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
subset_Immune_cell%>%NormalizeData%>%ScaleData%>%FindVariableFeatures

AT2_300<-subset_AT2[,sample(ncol(subset_AT2),300)]
AT1_300<-subset_AT1[,sample(ncol(subset_AT1),300)]
Club_cell_300<-subset_Club_cell[,sample(ncol(subset_Club_cell),300)]
Basal_cell_300<-subset_Basal_cell[,sample(ncol(subset_Basal_cell),300)]
Ciliated_cell_300<-subset_Ciliated_cell[,sample(ncol(subset_Ciliated_cell),300)]
Immune_cell_300<-subset_Immune_cell[,sample(ncol(subset_Immune_cell),300)]
AT2_300%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
AT1_300%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
Club_cell_300%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
Ciliated_cell_300%>%NormalizeData%>%ScaleData%>%FindVariableFeatures
Immune_cell_300%>%NormalizeData%>%ScaleData%>%FindVariableFeatures

#CCA on subset_donors(normal_merged/AT1/AT2/Club/Basal/Ciliated/Immune) ↔ LUAD_merged & visualization
library(cowplot)
hvg.union<-union(x=HVFInfo(object=normal_100_merged),y=HVFInfo(object=LUAD_merged))
normal_100_merged@meta.data[,"protocol"]<-"normal"
LUAD_merged@meta.data[,"protocol"]<-"LUAD"
normal_100_CCA<-RunCCA(object1 = normal_100_merged, object2 = LUAD_merged)
normal_LUAD_plot<-DimPlot(object = normal_100_CCA,reduction = "cca", return=TRUE)
normal_LUAD_vln<-VlnPlot(object = normal_100_CCA,features = "CC_1",group.by = "protocol")
normal_LUAD_heat<-DimHeatmap(object = normal_100_CCA,reduction = "cca",cells = 120,dims = 1:9,balanced = TRUE)

hvg.union<-union(x=HVFInfo(object=AT1_300),y=HVFInfo(object=LUAD_merged))
AT1_300@meta.data[,"protocol"]<-"AT1"
LUAD_merged@meta.data[,"protocol"]<-"LUAD"
AT1_300_CCA<-RunCCA(object1 = AT1_300,object2 = LUAD_merged)
AT1_LUAD_plot<-DimPlot(object = AT1_300_CCA,reduction = "cca", return=TRUE)
AT1_LUAD_vln<-VlnPlot(object = AT1_300_CCA,features = "CC_1",group.by = "protocol")
AT1_LUAD_heat<-DimHeatmap(object = AT1_300_CCA,reduction = "cca",cells = 120,dims = 1:9,balanced = TRUE)

#CCA on subset_immune(AM/BC/Mono/DC) ↔ LUAD_merged & visualization
subset_Alveolar_macrophage<-subset(x=donors,idents="Alveolar macrophage")
subset_B_cell<-subset(x=donors,idents = "B cell")
subset_Monocyte<-subset(x=donors, idents = "monocyte")
subset_Dendritic_cell<-subset(x=donors, idents = "Dendritic cell")
subset_Alveolar_macrophage%>%NormalizeData%>%FindVariableFeatures%>%ScaleData
subset_B_cell%>%NormalizeData%>%FindVariableFeatures%>%ScaleData
subset_Monocyte%>%NormalizeData%>%FindVariableFeatures%>%ScaleData
subset_Dendritic_cell%>%NormalizeData%>%FindVariableFeatures%>%ScaleData

hvg.union<-union(x=HVFInfo(object=LC_PT_sc_merged),y=HVFInfo(object=AM300))
LC_PT_sc_merged@meta.data[,"protocol"]<-"LUAD"
AM300@meta.data[,"protocol"]<-"Alveolar macrophage"
LC_AM300_CCA<-RunCCA(object1 = LC_PT_sc_merged, object2 = subset_Alveolar_macrophage)
LC_AM300_plot<-DimPlot(object = LC_AM300_CCA,reduction = "cca",return=TRUE)
LC_AM300_vln<-VlnPlot(object = LC_AM300_CCA,features = "CC_1",group.by = "protocol")
LC_AM300_heat<-DimHeatmap(object = LC_AM300_CCA,reduction = "cca",cells = 200,dims = 1:9,balanced = TRUE)

# LUAD_merged clustering -> 결과 LC_PT_sc와 H358은 확실하게 clustering
LUAD_merged<-FindNeighbors(object=LUAD_merged, dims = 1:50)
LUAD_merged<-FindClusters(object=LUAD_merged,dims=1:50,resolution=0.7,save.SNN=T,force.recalc=T)
LUAD_merged<-RunTSNE(object = LUAD_merged,dims=1:10,do.fast=T,check_duplicates=F)
DimPlot(object=LUAD_merged,reduction = "tsne")
LC_PT_sc_merged<-RunPCA(object=LC_PT_sc_merged, npcs = 50)
ElbowPlot(LC_PT_sc_merged, ndims = 50)
LC_PT_sc_merged<-FindNeighbors(object=LC_PT_sc_merged,dims = 1:50)
LC_PT_sc_merged<-FindClusters(object=LC_PT_sc_merged,dims=1:50,resolution=0.7,save.SNN=T,force.recalc=T)
LC_PT_sc_merged<-RunTSNE(object=LC_PT_sc_merged,dims = 1:50,do.fast=T)
DimPlot(object=LC_PT_sc_merged,reduction="tsne")


