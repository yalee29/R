# download & merge.h5 file (donor1~donor8)
donor1<-Read10X_h5("./Donor_single_cell/GSM3489182_Donor_01_filtered_gene_bc_matrices_h5.h5")
donor2<-Read10X_h5("./Donor_single_cell/GSM3489185_Donor_02_filtered_gene_bc_matrices_h5.h5")
donor3<-Read10X_h5("./Donor_single_cell/GSM3489187_Donor_03_filtered_gene_bc_matrices_h5.h5")
donor4<-Read10X_h5("./Donor_single_cell/GSM3489189_Donor_04_filtered_gene_bc_matrices_h5.h5")
donor5<-Read10X_h5("./Donor_single_cell/GSM3489191_Donor_05_filtered_gene_bc_matrices_h5.h5")
donor6<-Read10X_h5("./Donor_single_cell/GSM3489193_Donor_06_filtered_gene_bc_matrices_h5.h5")
donor7<-Read10X_h5("./Donor_single_cell/GSM3489195_Donor_07_filtered_gene_bc_matrices_h5.h5")
donor8<-Read10X_h5("./Donor_single_cell/GSM3489197_Donor_08_filtered_gene_bc_matrices_h5.h5")
rownames(donor1)<-make.unique(rownames(donor1))
rownames(donor2)<-make.unique(rownames(donor2))
rownames(donor3)<-make.unique(rownames(donor3))
rownames(donor4)<-make.unique(rownames(donor4))
rownames(donor5)<-make.unique(rownames(donor5))
rownames(donor6)<-make.unique(rownames(donor6))
rownames(donor7)<-make.unique(rownames(donor7))
rownames(donor8)<-make.unique(rownames(donor8))
colnames(donor1)=paste("D1",colnames(donor1))
colnames(donor2)=paste("D2",colnames(donor2))
colnames(donor3)=paste("D3",colnames(donor3))
colnames(donor4)=paste("D4",colnames(donor4))
colnames(donor5)=paste("D5",colnames(donor5))
colnames(donor6)=paste("D6",colnames(donor6))
colnames(donor7)=paste("D7",colnames(donor7))
colnames(donor8)=paste("D8",colnames(donor8))
donor1<-CreateSeuratObject(count=donor1)
donor2<-CreateSeuratObject(count=donor2)
donor3<-CreateSeuratObject(count=donor3)
donor4<-CreateSeuratObject(count=donor4)
donor5<-CreateSeuratObject(count=donor5)
donor6<-CreateSeuratObject(count=donor6)
donor7<-CreateSeuratObject(count=donor7)
donor8<-CreateSeuratObject(count=donor8)
donors<-merge(x=donor1,y=list(donor2,donor3,donor4,donor5,donor6,donor7,donor8))
human_ENSEMBL_genename<-read_tsv("./Donor_single_cell/human_ENSEMBL_genename.tsv")

# sample(donors) clustering 1 (nPCA↑, default → RunTSNE dims 1:10)
donors<-NormalizeData(object=donors)
donors<-FindVariableFeatures(object=donors,x.low.cutoff=0.0125,x.high.cutoff=3,y.cutoff=0.5)
donors<-ScaleData(object=donors,do.par=T,num.cores=2)
donors<-RunPCA(object=donors, npcs = 50)
ElbowPlot(donors, ndims = 50)
donors<-FindNeighbors(object=donors, dims = 1:50)
donors<-FindClusters(object=donors,dims=1:50,resolution=0.7,save.SNN=T,force.recalc=T)
donors<-RunTSNE(object = donors,h)
DimPlot(object=donors,reduction = "tsne")

# 박성열선생님따라하기 → tSNE plot with cluster number
library(tidyverse)
library(dplyr)
library(tibble)
tsne_dt<-donors@reductions$tsne@cell.embeddings
tsne_dt<-tsne_dt%>%as.data.frame%>%rownames_to_column('cell_id')
dim(tsne_dt)
ident_dt<-Idents(donors)%>%as.data.frame%>%rownames_to_column('cell_id')
colnames(ident_dt)<-c('cell_id','cluster')
tsne_dt<-left_join(tsne_dt,ident_dt)
pos_dt<-tsne_dt%>%group_by(cluster)%>%summarise(med1=median(tSNE_1),med2=median(tSNE_2))  
ggplot(tsne_dt,aes(x=tSNE_1,y=tSNE_2))+
  geom_point(aes(color=cluster),alpha=0.5)+
  geom_text(data=pos_dt,aes(x=med1,y=med2,label=cluster))

#Find Marker gene_1 #그룹간 비교(박성열선생님) 
#9번 cluster와 나머지 cluster 비교
#찾은 marker는 marker_dt에 저장
group1 <- c('9')
group2 <- NULL  
marker_dt <- FindMarkers(donors, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_name') %>% as.tibble()
cur_markers <- marker_dt$gene_name
g1_markers <- marker_dt %>% filter(avg_logFC >= 0) %>% .$gene_name
g2_markers <- marker_dt %>% filter(avg_logFC < 0) %>% .$gene_name

cell1 <- names(Idents(donors)[Idents(donors) %in% group1])
cell2 <- names(Idents(donors)[Idents(donors) %in% group2])
tsne_dt <- donors@reductions$tsne@cell.embeddings %>% as.data.frame() %>% rownames_to_column('cell_id')
tsne_dt <- tsne_dt %>% mutate(group = ifelse(cell_id %in% cell1, paste(group1,collapse=','), ifelse(cell_id %in% cell2, paste(group2, collapse=','),'other')))
ggplot(tsne_dt, aes(tSNE_1, tSNE_2))+
  geom_point(aes(color=group), alpha=0.5)

FeaturePlot(object = donors, reduction = 'tsne',features = g2_markers[1:12], ncol=4)
FeaturePlot(object = cdonors, reduction = 'tsne',features = g2_markers[1:12], ncol=4)

#Finde Marker Gene_2 
cluster.markers<-FindAllMarkers(object=donors,only.pos = TRUE,min.pct = 0.5,logfc.threshold = 0.5,max.cells.per.ident = 200)
cluster_marker_group<-cluster.markers%>%group_by(cluster)

#FeaturePlot 
cluster0.markers<-subset(x=cluster.markers,cluster==0)
top10_cluster1.markers<-head(x=cluster1.markers,n=10)
FeaturePlot(object=donors,features=c("S100A14","SLPI","RPS2","SFTPC","B2M","FTL","HLA-DQA1","C1QA","PGC","CAV1","SRGN","CSTB","C1QA","TIMP3","CAPS","TOP2A","JCHAIN","GPX3","RPS27A","KIAA0101","RPS2","AGER"))

#Re-clustering (cluster 합치기)
current.cluster.ids<-c(0:29)
new.cluster.ids<-c("Basal cell","AT2","AT2","Alveolar macrophage","AT2","Alveolar macrophage","Alveolar macrophage","AT2","AT2","Monocyte","Alveolar macrophage","AT2","Dendritic cell","AT2","Alveolar macrophage","Alveolar macrophage","AT1","Alveolar macrophage","Endothelial cell","Basal cell","Ciliated cell","Club cell","Alveolar macrophage","B cell","Club cell","Basal cell","Fibroblast","Alveolar macrophage","Endothelial cell","Alveolar macrophage")
donors@active.ident<-plyr::mapvalues(x=donors@active.ident,from = current.cluster.ids,to=new.cluster.ids)
TSNEPlot(object=donors,labels=TRUE,pt.size=0.5)

#Stash identities for later
donors<-StashIdent(object=donors,save.name="??")  