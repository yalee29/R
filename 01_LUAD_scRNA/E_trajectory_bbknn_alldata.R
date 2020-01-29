#Mouse embryo data loading & processing (extract foregut/ gene name>ENSEMBL)
#E8.25_foregut & E9.5_11.5 & E14.5 & P6~P10(1/2/3) & collected(E14.5/E16.5/E18.5/PND107)
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
collected<-read_tsv("./projects/LUAD_scRNA/mouse_embryo/treulein/GSE52583/collected.txt")
LungMap_E18.5 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_E18.5.xlsx")
LungMap_P1 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_P01.xlsx")
LungMap_P3 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_P03.xlsx")
LungMap_P7 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_P07.xlsx")
LungMap_P10 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_P10.xlsx")
LungMap_P14 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_P14.xlsx")
LungMap_E18.5 <- LungMap_E18.5[-1,] #LungMap data는 1행에 필요없는 내용 있음. 
LungMap_P1 <- LungMap_P1[-1,]
LungMap_P3 <- LungMap_P3[-1,]
LungMap_P7 <- LungMap_P7[-1,]
LungMap_P10 <- LungMap_P10[-1,]
LungMap_P14 <- LungMap_P14[-1,]
LungMap_E18.5 <- rename(LungMap_E18.5, gene_name = ...1) 
LungMap_P1 <- rename(LungMap_P1, gene_name = ...1) 
LungMap_P3 <- rename(LungMap_P3, gene_name = ...1) 
LungMap_P7 <- rename(LungMap_P7, gene_name = ...1) 
LungMap_P10 <- rename(LungMap_P10, gene_name = ...1)
LungMap_P14 <- rename(LungMap_P14, gene_name = "Type") 

#mouse gene_name ↔ gene ensembl table (reference)
gene_name_ensembl_mouse<-read_tsv("./Mouse_embryo/gene_name_ensembl_mouse.tsv")
gene_name_ensembl_mouse<- gene_name_ensembl_mouse[!duplicated(gene_name_ensembl_mouse[1]),] #remove duplicated gene name

#gene name > gene_ensembl
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
collected_ensembl<-merge(collected, gene_name_ensembl_mouse, by.x = "gene_symbol", by.y = "Gene name")
collected_ensembl<-collected_ensembl[,-1]
names(collected_ensembl)[ncol(collected_ensembl)] <- c("Gene stable ID")
LungMap_E18.5<-merge(LungMap_E18.5, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_E18.5_ensembl<-LungMap_E18.5[,-1]
LungMap_P1<-merge(LungMap_P1, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_P1_ensembl<-LungMap_P1[,-1]
LungMap_P3<-merge(LungMap_P3, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_P3_ensembl<-LungMap_P3[,-1]
LungMap_P7<-merge(LungMap_P7, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_P7_ensembl<-LungMap_P7[,-1]
LungMap_P10<-merge(LungMap_P10, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_P10_ensembl<-LungMap_P10[,-1]
LungMap_P14<-merge(LungMap_P14, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_P14_ensembl<-LungMap_P14[,-1]


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
union_gene_vector = intersect(union_gene_vector,collected_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_E18.5_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_P1_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_P3_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_P7_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_P10_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_P14_ensembl$`Gene stable ID`)

A = E8.25_foregut_ensembl[match(union_gene_vector, E8.25_foregut_ensembl$gene_id),-1]
colnames(A) = paste0("E8.25_", colnames(A))
B = E9.5_11.5_ensembl[match(union_gene_vector, E9.5_11.5_ensembl$`Gene stable ID`),-ncol(E9.5_11.5_ensembl)]
colnames(B) = paste0("E9.5to11.5_", colnames(B))
C = E14.5_ensembl[match(union_gene_vector, E14.5_ensembl$'Gene stable ID'), -ncol(E14.5_ensembl)]
colnames(C) = paste0("E14.5.1_", colnames(C))
D = P6_10_1_ensembl[match(union_gene_vector,P6_10_1_ensembl$`Gene stable ID`),-ncol(P6_10_1_ensembl)]
colnames(D) = paste0("P6to10.1_", colnames(D))
G = P6_10_2_ensembl[match(union_gene_vector,P6_10_2_ensembl$`Gene stable ID`),-ncol(P6_10_2_ensembl)]
colnames(G) = paste0("P6to10.2_", colnames(G))
H = P6_10_3_ensembl[match(union_gene_vector,P6_10_3_ensembl$`Gene stable ID`),-ncol(P6_10_3_ensembl)]
colnames(H) = paste0("P6to10.3_", colnames(H))
I = P1_1_ensembl[match(union_gene_vector,P1_1_ensembl$`Gene stable ID`),-ncol(P1_1_ensembl)]
colnames(I) = paste0("PND1.1_", colnames(I))
J = P1_2_ensembl[match(union_gene_vector,P1_2_ensembl$`Gene stable ID`),-ncol(P1_2_ensembl)]
colnames(J) = paste0("PND1.2_", colnames(J))
L = LungMap_E18.5_ensembl[match(union_gene_vector,LungMap_E18.5_ensembl$`Gene stable ID`),-ncol(LungMap_E18.5_ensembl)]
colnames(L) = paste0("LungMap.E18.5_", colnames(L))
M = LungMap_P1_ensembl[match(union_gene_vector,LungMap_P1_ensembl$`Gene stable ID`),-ncol(LungMap_P1_ensembl)]
colnames(M) = paste0("LungMap.PND1_", colnames(M))
N = LungMap_P3_ensembl[match(union_gene_vector,LungMap_P3_ensembl$`Gene stable ID`),-ncol(LungMap_P3_ensembl)]
colnames(N) = paste0("LungMap.PND3_", colnames(N))
O = LungMap_P7_ensembl[match(union_gene_vector,LungMap_P7_ensembl$`Gene stable ID`),-ncol(LungMap_P7_ensembl)]
colnames(O) = paste0("LungMap.PND7_", colnames(O))
P = LungMap_P10_ensembl[match(union_gene_vector,LungMap_P10_ensembl$`Gene stable ID`),-ncol(LungMap_P10_ensembl)]
colnames(P) = paste0("LungMap.PND10_", colnames(P))
Q = LungMap_P14_ensembl[match(union_gene_vector,LungMap_P14_ensembl$`Gene stable ID`),-ncol(LungMap_P14_ensembl)]
colnames(Q) = paste0("LungMap.PND14_", colnames(Q))

# downsampling → 500 cell
C_500<-knn(C,10)
D_500<-knn(D,10)
G_500<-knn(G,10)
H_500<-knn(H,10)
I_500<-knn(I,10)
J_500<-knn(J,10)

sample<-sample(ncol(C_500),500)
C_500<-C_500[,c(sample)]
sample<-sample(ncol(D_500),500)
D_500<-D_500[,c(sample)]
sample<-sample(ncol(G_500),500)
G_500<-G_500[,c(sample)]
sample<-sample(ncol(H_500),500)
H_500<-H_500[,c(sample)]
sample<-sample(ncol(I_500),500)
I_500<-I_500[,c(sample)]
sample<-sample(ncol(J_500),500)
J_500<-J_500[,c(sample)]

#collected_ensembl 자료는 E8.5~P15까지 data가 다 섞여있기 때문에 추가적인 과정 필요. 
K = collected_ensembl[match(union_gene_vector,collected_ensembl$`Gene stable ID`),-ncol(collected_ensembl)]
K1<-K[1:80]
K2<-K[81:125]
K3<-K[126:171]
K4<-K[172:198]
colnames(K1) = paste0("E18.5_", colnames(K1))
colnames(K2) = paste0("E14.5.2_", colnames(K2))
colnames(K3) = paste0("P15_", colnames(K3))
colnames(K4) = paste0("E16.5_", colnames(K4))
K = cbind(K1, K2, K3, K4)

E = cbind(A,B,C_500,D_500,G_500,H_500,I_500,J_500,K,L,M,N,O,P,Q)
rownames(E) = union_gene_vector
table(duplicated(t(E)))
E = E[,!duplicated(t(E))]
table(duplicated(t(E)))
batch_vector = stringr::str_replace(colnames(E), "_.*","")

# Seurat으로 processing & batch
E_seurat1 = CreateSeuratObject(E)
E_seurat1$batch = batch_vector
E_seurat1 <-NormalizeData(E_seurat1)
E_seurat1 <-FindVariableFeatures(E_seurat1)
E_seurat1 <-ScaleData(E_seurat1)
E_seurat1 <- RunPCA(E_seurat1, npcs = 50)
E_seurat1 <- JackStraw(E_seurat1,dims = 50)
E_seurat1 <- ScoreJackStraw(E_seurat1, dims = 1:50)
JackStrawPlot(E_seurat1, dims = 1:50)
ElbowPlot(E_seurat, ndims = 50)
E_seurat1 <- FindNeighbors(E_seurat1, dims = 1:50)
E_seurat1 <- FindClusters(E_seurat1, resolution = 0.7)
E_seurat1 <- RunTSNE(E_seurat1, dims = 1:50)
E_seurat1 <- RunUMAP(E_seurat1, dims = 1:50)
DimPlot(E_seurat1, reduction = "tsne")
DimPlot(E_seurat1, reduction = "tsne", group.by = "batch")
DimPlot(E_seurat1, reduction = "umap")
DimPlot(E_seurat1, reduction = "umap", group.by = "batch")

#bbknn1 실행 
#bbknn1 함수등록
RunBBKNN <- function(object,   
                     batch.key = "batch", # slot of object that store batch information
                     dims.use = 1:50,
                     neighbors_within_batch = 6,
                     trim = NULL,
                     reduction.use = "pca",
                     clustering=FALSE,
                     python.path = "/home/users/yunah1029/anaconda3/bin/python") {
  
  write_tsv(as.data.frame(object@reductions[[reduction.use]]@cell.embeddings[,dims.use]),"pca.tmp.tsv")
  write_tsv(data.frame(object[[batch.key]]), "batch.tmp.tsv")
  trim = ifelse(is.null(trim), "", paste0(", trim =", trim))
  write_lines(paste0("
import numpy as np
import pandas as pd
import scanpy.api as sc
import anndata
import bbknn
import os
from scipy import sparse
X=pd.read_csv('pca.tmp.tsv',sep='\t')
obs=pd.read_csv('batch.tmp.tsv',sep='\t')
adata=anndata.AnnData(X=X.values, obs=obs)
sc.tl.pca(adata)
adata.obsm.X_pca = X.values[:,0:(adata.obsm.X_pca[1].size)]
bbknn.bbknn(adata, batch_key = '",batch.key,"', neighbors_within_batch = ",neighbors_within_batch, trim,")
sc.tl.umap(adata)
np.savetxt('umap.tmp.tsv',adata.obsm.X_umap,delimiter='\t')
",ifelse(!clustering,"","sc.tl.louvain(adata, resolution=2)
adata.obs.to_csv('louvain.tmp.csv')
")),"bbknn.tmp.py")
  pythonexitstatus <- system(paste(python.path, "bbknn.tmp.py"))
  if(pythonexitstatus != 0) {stop()}
  umap <- read_tsv("umap.tmp.tsv",col_names=c("BBKNN_1","BBKNN_2")) %>%
    as.data.frame %>% as.matrix
  
  rownames(umap) <- rownames(object@reductions[[reduction.use]]@cell.embeddings)
  object[["bbknn"]] <- CreateDimReducObject(embeddings = umap, assay = DefaultAssay(object), key = "BBKNN_")
  if(clustering){
    louvain <- read_csv("louvain.tmp.csv")$louvain
    object$louvain <- read_csv("louvain.tmp.csv")$louvain
  }
  system("rm -f louvain.tmp.csv umap.tmp.tsv pca.tmp.tsv batch.tmp.tsv bbknn.tmp.py")
  object
}

#RunBBKNN + 여러조건으로 한꺼번에 돌려서 pdf에 저장 → 원하는 plot, 조건 찾기.
E_seurat1_bb <- E_seurat1 %>% RunBBKNN(dims.use = 1:50,
                                  neighbors_within_batch = 14,
                                  trim = 40,
                                  batch.key = "batch",
                                  python.path = "/home/users/yunah1029/anaconda3/bin/python")

DimPlot(E_seurat1_bb, reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch")+ggtitle("bbknn1; k=14, trim=40")

pdf("tmp.pdf", 6,5)
for(k in c(3, 6, 10, 15, 20)){
  for(trim in (k*length(unique(E_seurat$batch))*c(0, .9,.8,.7,.6,.5,.4)) %>% round){
    if(trim == 0) { trim = NULL}
    print(
    E_seurat %>% RunBBKNN(dims.use = 1:30,
                          neighbors_within_batch = k,
                          trim = trim,
                          batch.key = "batch",
                          python.path = "/home/users/yunah1029/anaconda3/bin/python") %>%
      DimPlot(reduction = "bbknn") %>% 
      LabelClusters(id="ident") +
      ggtitle(paste0("k=",k,",trim=",trim,",every cells"))
    )
  }
}
E_seurat_epi <- subset(E_seurat, idents = c(14,21,23,12,22,10,24,5,1))
for(k in c(3, 6, 10, 15, 20)){
  for(trim in (k*length(unique(E_seurat$batch))*c(0, .9,.8,.7,.6,.5,.4,.3)) %>% round){
    if(trim == 0) { trim = NULL}
    print(
    E_seurat_epi %>% RunBBKNN(dims.use = 1:30,
                          neighbors_within_batch = k,
                          trim = trim,
                          batch.key = "batch",
                          python.path = "/home/users/yunah1029/anaconda3/bin/python") %>%
      DimPlot(reduction = "bbknn") %>%
      LabelClusters(id="ident") +
      ggtitle(paste0("k=",k,",trim=",trim,",endoderm,epithelial cells"))
    )
  }
}
dev.off()

#원하는 plot, 조건 ; E_seurat_epi에서 k=3, trim=18
E_seurat_epi<- E_seurat_epi %>% RunBBKNN(dims.use = 1:30,
                                         neighbors_within_batch = 3,
                                         trim = 18,
                                         batch.key = "batch",
                                         python.path = "/home/users/yunah1029/anaconda3/bin/python") 

test_ver2 <- E_seurat_epi[,E_seurat_epi$batch %in% c("P6to10.2","PND1.1","E18.5","P15","E14.5.2","E16.5","E8.25")]
test_ver2 <- test_ver2 %>% RunBBKNN(dims.use = 1:30,
                                       neighbors_within_batch = 3,
                                       trim = 18,
                                       batch.key = "batch",
                                       python.path = "/home/users/yunah1029/anaconda3/bin/python")
DimPlot(test_ver2, reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; k=3, trim=18, epi")


#Find Marker genes
group1 <- c('12')
group2 <- NULL  
marker_ver2 <- FindMarkers(test_ver2, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
marker_ver2_12 <- marker_dt

#UMAP with batch
plot <- DimPlot(E_seurat_epi, reduction = "bbknn", group.by = "batch") 
LabelClusters(plot = plot, id = "batch") + ggtitle(paste0("k= 3 /","trim= 18 /","endoderm & epithelial cells"))
#UMAP with cluster number
plot_1 <- DimPlot(E_seurat_epi, reduction = "bbknn", group.by = "ident") 
LabelClusters(plot = plot_1, id = "ident")

#Find Marker gene_1 #그룹간 비교(박성열선생님) 
#group1 cluster와 나머지 cluster 비교
group1 <- c('1')
group2 <- NULL  
marker_dt <- FindMarkers(E_seurat, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
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




#bbknn2 test with epithelial cells

P6to10.1_down <- E_seurat_epi[,E_seurat_epi$batch %in% "P6to10.1"]
table(test)

unique(test$batch)

test <- test %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()
test_bbknn <- test %>%
  FindBatchBalancedNeighbors(batch.key = "batch",
                             neighbors_within_batch = 3,
                             prune.SNN = 2/15) %>%
  FindClusters(graph.name = "bbknn_snn") %>%
  RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_")
DimPlot(test_bbknn, reduction = "bbknn", group.by = "batch", cells = WhichCells(test_bbknn, idents = c(0,1,2,3))) %>% LabelClusters(id = "batch") 

#cf with bbknn1
test_standard <- test %>%
  FindNeighbors(prune.SNN = 2/15) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)

DimPlot(test_standard, reduction = "umap", group.by = "batch")


table(test_standard$batch)




