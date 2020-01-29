#Mouse embryo data loading & processing (extract foregut/ gene name>ENSEMBL)
#E8.25_foregut & E9.5_11.5 %E14.5 % P6~P10(1/2/3)
E8.25<-read_tsv("./Mouse_embryo/E8.25_E-MTAB-6153.tsv")
E9.5_11.5<-read_csv("./Mouse_embryo/GSE87038_E9.5_E10.5_E11.5_counts_matrix.csv")
E14.5<-read.csv("./Mouse_embryo/GSM2906422_FetalLung_dge.txt", sep = ' ', stringsAsFactors = FALSE, na.strings = "")
P6_10_1<-read.csv("./Mouse_embryo/GSM2906429_Lung1_dge.txt", sep = ' ', stringsAsFactors = FALSE, na.strings = "")
P6_10_2<-read.csv("./Mouse_embryo/GSM2906430_Lung2_dge.txt", sep = ' ', stringsAsFactors = FALSE, na.strings = "")
P6_10_3<-read.csv("./Mouse_embryo/GSM2906431_Lung3_dge.txt", sep = ' ', stringsAsFactors = FALSE, na.strings = "")
E8.25_foregut<-E8.25%>%select(gene_id, starts_with('foregut')) #ncol(E8.25_foregut)=185
E8.25_foregut_ensembl<-E8.25%>%select(gene_id, starts_with('foregut')) #ncol(E8.25_foregut)=185 
P1_1<-read_tsv("./projects/LUAD_scRNA/mouse_embryo/yan/GSM3464101_Batch1_out_gene_exon_tagged.dge.txt")
P1_2<-read_tsv("./projects/LUAD_scRNA/mouse_embryo/yan/GSM3464102_Batch2_out_gene_exon_tagged.dge.txt")

gene_name_ensembl_mouse<-read_tsv("./Mouse_embryo/gene_name_ensembl_mouse.tsv")
gene_name_ensembl_mouse<- gene_name_ensembl_mouse[!duplicated(gene_name_ensembl_mouse[1]),] #remove duplicated gene name

E9.5_11.5_ensembl<-merge(E9.5_11.5,gene_name_ensembl_mouse,by.x = "Gene",by.y = "Gene name")
E9.5_11.5_ensembl<-E9.5_11.5_ensembl[,-1] # remove mouse gene name # Gene stable ID in [,225]
E14.5<-rownames_to_column(E14.5,"gene_name")
E14.5_ensembl<-E14.5_ensembl[,-1]
P6_10_1<-P6_10_1%>%rownames_to_column("gene_name")
P6_10_1_ensembl<-merge(P6_10_1,gene_name_ensembl_mouse,by.x = "gene_name",by.y = "Gene name")
P6_10_1_ensembl<-P6_10_1_ensembl[,-1]
P6_10_2<-P6_10_2%>%rownames_to_column("gene_name")
P6_10_2_ensembl<-merge(P6_10_2,gene_name_ensembl_mouse,by.x = "gene_name",by.y = "Gene name")
P6_10_2_ensembl<-P6_10_2_ensembl[,-1]
P6_10_3<-P6_10_3%>%rownames_to_column("gene_name")
P6_10_3_ensembl<-merge(P6_10_3,gene_name_ensembl_mouse,by.x = "gene_name",by.y = "Gene name")
P6_10_3_ensembl<-P6_10_3_ensembl[,-1]
P1_1_ensembl<-merge(P1_1, gene_name_ensembl_mouse, by.x = "GENE", by.y = "Gene name")
P1_1_ensembl<-P1_1_ensmebl[,-1]
names(P1_1_ensembl)[ncol(P1_1_ensembl)] <- c("Gene stable ID")
P1_2_ensembl<-merge(P1_2, gene_name_ensembl_mouse, by.x = "GENE", by.y = "Gene name")
P1_2_ensembl<-P1_2_ensembl[,-1]
names(P1_2_ensembl)[ncol(P1_2_ensembl)] <- c("Gene stable ID")

ensembl_mouse_human<-read_tsv("./Mouse_embryo/ensembl_mouse_human.tsv") #나중에 human과 비교해야하니깐!

# E14.5 cell 수(8396)↓ by KNN pooling & random sampling 
# 1. knn pooling 
knn <- function(A,k) { # column-wisely similar columns were summed
  i = nrow(A)
  j = ncol(A)
  B = matrix(0, nrow = i, ncol = j)
  D = cor(A)
  for(m in 1:j){
    r = rank(-D[,m],ties.method = "first") <= k + 1
    B[,m] = rowSums(A[,r,drop = F])
  }
  colnames(B) <- colnames(A)
  rownames(B) <- rownames(A)
  B
}

E14.5_1k<-E14.5_ensembl 
rownames(E14.5_1k) = E14.5_ensembl$`Gene stable ID` #마지막 column이 gene id인 상태 → rowname으로 지정, column 삭제
E14.5_1k<-E14.5_1k[,-ncol(E14.5_1k)]

E14.5_1k<-knn(E14.5_1k,10) #k 값은 임의로 10정도로 줌. data가 많지 않으면 5

# 2. random sampling (8396 → 4000)
set.seed(42)
sample<-sample(ncol(E14.5_1k),4000)
E14.5_1k_test<-E14.5_1k[,c(sample)]
E14.5_1k<-E14.5_1k_test
rm(E14.5_1k_test)

#P6_10_1,2,3 
P6to10.1<-P6_10_1_ensembl
P6to10.2<-P6_10_2_ensembl
P6to10.3<-P6_10_3_ensembl
rownames(P6to10.1)=P6_10_1_ensembl$'Gene stable ID'
rownames(P6to10.2)=P6_10_2_ensembl$'Gene stable ID'
rownames(P6to10.3)=P6_10_3_ensembl$'Gene stable ID'
P6to10.1<-P6to10.1[,-ncol(P6to10.1)]
P6to10.2<-P6to10.2[,-ncol(P6to10.2)]
P6to10.3<-P6to10.3[,-ncol(P6to10.3)]

P6to10.1<-knn(P6to10.1,5)
set.seed(42)
sample_1<-sample(ncol(P6to10.1),500)
P6to10.1_test<-P6to10.1[,c(sample_1)]
P6to10.1<-P6to10.1_test
rm(P6to10.1_test)
P6to10.2<-knn(P6to10.2,5)
set.seed(42)
sample_2<-sample(ncol(P6to10.2),500)
P6to10.2_test<-P6to10.2[,c(sample_2)]
P6to10.2<-P6to10.2_test
rm(P6to10.2_test)
P6to10.3<-knn(P6to10.3,5)
set.seed(42)
sample_3<-sample(ncol(P6to10.3),500)
P6to10.3_test<-P6to10.3[,c(sample_3)]
P6to10.3<-P6to10.3_test
rm(P6to10.3_test)

# Embryo datas 공통된 genes로 합치기(cbind) & batch vector 부여
union_gene_vector = intersect(E8.25_foregut_ensembl$gene_id, E9.5_11.5_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,E14.5_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,P6_10_2_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,P6_10_1_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,P6_10_3_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,P1_1_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,P1_2_ensembl$`Gene stable ID`)

A = E8.25_foregut_ensembl[match(union_gene_vector, E8.25_foregut_ensembl$gene_id),-1]
colnames(A) = paste0("E8.25_", colnames(A))
B = E9.5_11.5_ensembl[match(union_gene_vector, E9.5_11.5_ensembl$`Gene stable ID`),-ncol(E9.5_11.5_ensembl)]
colnames(B) = paste0("E9.5to11.5_", colnames(B))
C = E14.5_1k[match(union_gene_vector, rownames(E14.5_1k)), ]
colnames(C) = paste0("E14.5_", colnames(C))
D = P6_10_1_ensembl[match(union_gene_vector,P6_10_1_ensembl$`Gene stable ID`),-ncol(P6_10_1_ensembl)]
colnames(D) = paste0("P6to10.1_", colnames(D))
G = P6_10_2_ensembl[match(union_gene_vector,P6_10_2_ensembl$`Gene stable ID`),-ncol(P6_10_2_ensembl)]
colnames(G) = paste0("P6to10.2_", colnames(G))
H = P6_10_3_ensembl[match(union_gene_vector,P6_10_3_ensembl$`Gene stable ID`),-ncol(P6_10_3_ensembl)]
colnames(H) = paste0("P6to10.3_", colnames(H))
I = P1_1_ensembl[match(union_gene_vector,P1_1_ensembl$`Gene stable ID`),-ncol(P1_1_ensembl)]
colnames(I) = paste0("P1.1_", colnames(I))
J = P1_2_ensembl[match(union_gene_vector,P1_2_ensembl$`Gene stable ID`),-ncol(P1_2_ensembl)]
colnames(J) = paste0("P1.2_", colnames(J))

E = cbind(A,B,C,D,G,H,I,J)
rownames(E) = union_gene_vector
table(duplicated(t(E)))
E = E[,!duplicated(t(E))]
table(duplicated(t(E)))
batch_vector = stringr::str_replace(colnames(E), "_.*","")

# Seurat으로 processing & batch
E_seurat = CreateSeuratObject(E)
E_seurat$batch = batch_vector
E_seurat <-NormalizeData(E_seurat)
E_seurat <-FindVariableFeatures(E_seurat)
E_seurat <-ScaleData(E_seurat)
E_seurat <- RunPCA(E_seurat, npcs = 50)
E_seurat <- JackStraw(E_seurat,dims = 50)
E_seurat <- ScoreJackStraw(E_seurat, dims = 1:50)
JackStrawPlot(E_seurat, dims = 1:50)
ElbowPlot(E_seurat, ndims = 50)
E_seurat <- FindNeighbors(E_seurat, dims = 1:30)
E_seurat <- FindClusters(E_seurat, resolution = 0.7)
E_seurat <- RunTSNE(E_seurat, dims = 1:30)
E_seurat <- RunUMAP(E_seurat, dims = 1:30)
DimPlot(E_seurat, reduction = "tsne")
DimPlot(E_seurat, reduction = "tsne", group.by = "batch")
DimPlot(E_seurat, reduction = "umap")
DimPlot(E_seurat, reduction = "umap", group.by = "batch")

#bbknn 실행 
#'Cannot allocate memory'라는 error 뜨면 node memory 부족 or environmnet 정리하기 
source("~kjyi/Projects/cav/bbknn.R")
E_seurat <- E_seurat %>% RunBBKNN(dims.use = 1:20, trim = 15 ,batch.key = "batch", python.path = "/home/users/yunah1029/anaconda3/bin/python")
E_seurat_bb = E_seurat
umap.emb = read_tsv("umap.tmp.tsv", col_names = F) %>% as.matrix()
plot(umap.emb)
colnames(umap.emb) = colnames(E_seurat_bb@reductions$umap@cell.embeddings)
rownames(umap.emb) = rownames(E_seurat@reductions$umap@cell.embeddings)
E_seurat_bb@reductions$umap@cell.embeddings = as.matrix(umap.emb)
E_seurat_bb@reductions$umap@cell.embeddings %>% head
E_seurat@reductions$umap@cell.embeddings %>% head

DimPlot(E_seurat_bb, reduction = "umap", group.by = "batch")
DimPlot(E_seurat_bb, reduction = "umap")

#Clustering
# 박성열선생님따라하기 → UMAP with cluster number
library(tidyverse)
library(dplyr)
library(tibble)
umap_dt<-E_seurat_bb@reductions$umap@cell.embeddings
umap_dt<-umap_dt%>%as.data.frame%>%rownames_to_column('cell_id')
dim(umap_dt)
ident_dt<-Idents(E_seurat_bb)%>%as.data.frame%>%rownames_to_column('cell_id')
colnames(ident_dt)<-c('cell_id','cluster')
umap_dt<-left_join(umap_dt,ident_dt)
umap_dt<-umap_dt[is.na(umap_dt$cluster)==FALSE,] #bbknn이후에 Idents()%>%as.data.frame()%>%nrow 하면 원래의 cell 수보다 감소한 상태. 따라서 NA가 생김. 
pos_dt<-umap_dt%>%group_by(cluster)%>%summarise(med1=median(UMAP_1),med2=median(UMAP_2))  
ggplot(umap_dt,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=cluster),alpha=0.5)+
  geom_text(data=pos_dt,aes(x=med1,y=med2,label=cluster))+
  theme_classic()

#Find Marker gene_1 #그룹간 비교(박성열선생님) 
#group1 cluster와 나머지 cluster 비교
group1 <- c('10')
group2 <- NULL  
marker_dt <- FindMarkers(E_seurat_bb, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
cur_markers <- marker_dt$gene_id
g1_markers <- marker_dt %>% filter(avg_logFC >= 0) %>% .$gene_id
g2_markers <- marker_dt %>% filter(avg_logFC < 0) %>% .$gene_id

cell1 <- names(Idents(E_seurat_bb)[Idents(E_seurat_bb) %in% group1])
cell2 <- names(Idents(E_seurat_bb)[Idents(E_seurat_bb) %in% group2])
umap_dt <- E_seurat_bb@reductions$umap@cell.embeddings %>% as.data.frame() %>% rownames_to_column('cell_id')
umap_dt <- umap_dt %>% mutate(group = ifelse(cell_id %in% cell1, paste(group1,collapse=','), ifelse(cell_id %in% cell2, paste(group2, collapse=','),'other')))

ggplot(umap_dt, aes(UMAP_1, UMAP_2))+
  geom_point(aes(color=group), alpha=0.5)

#Re-clustering (cluster 합치기)
current.cluster.ids<-c(0:19)
new.cluster.ids<-c("AT2","Immune","Pericyte","Immune","Pericyte","Pericyte","Foregut","Immune","Immune","Mesenchymal","Endothelial","Immune","Epithelial","Clara","Fibroblast","Immune","Ciliated","AT1","Immune","Immune")
E_seurat_bb@active.ident<-plyr::mapvalues(x=E_seurat_bb@active.ident,from = current.cluster.ids,to=new.cluster.ids)
UMAPPlot(object=E_seurat_bb,labels=TRUE,pt.size=0.5)

ggplot(umap_dt,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=cluster),alpha=0.5)+
  geom_text(data=pos_dt,aes(x=med1,y=med2,label=cluster))+
  theme_classic()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)



