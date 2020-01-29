# Seurat으로 processing & batch
E_seurat1 = CreateSeuratObject(E2)
E_seurat1$batch = batch_vector
E_seurat1 <-NormalizeData(E_seurat1)
E_seurat1 <-FindVariableFeatures(E_seurat1)
E_seurat1 <-ScaleData(E_seurat1)
E_seurat1 <- RunPCA(E_seurat1, npcs = 50)
E_seurat1 <- JackStraw(E_seurat1,dims = 50)
E_seurat1 <- ScoreJackStraw(E_seurat1, dims = 1:50)
JackStrawPlot(E_seurat1, dims = 1:50)
ElbowPlot(E_seurat1, ndims = 50)
E_seurat1 <- FindNeighbors(E_seurat1, dims = 1:40)
E_seurat1 <- FindClusters(E_seurat1, resolution = 0.7)
E_seurat1 <- RunTSNE(E_seurat1, dims = 1:40)
E_seurat1 <- RunUMAP(E_seurat1, dims = 1:40)
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

#bbknn1_total cell
E_seurat1_bb<- E_seurat1 %>% RunBBKNN(dims.use = 1:30,
                                      neighbors_within_batch =  4,
                                      trim = 30,
                                      batch.key = "batch",
                                      python.path = "/home/users/yunah1029/anaconda3/bin/python") 
DimPlot(E_seurat1_bb, reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; k=4, trim=30, all")
DimPlot(E_seurat1_bb, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; k=4, trim=30, all")

group1 <- c('21','3','16','12')
group2 <- NULL  
marker_dt <- FindMarkers(E_seurat1_bb, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()

new.cluster.ids <- c("Epithelial", "Immune", "Immune", "Mesenchymal", "Immune", "Immune", "Mesenchymal", "Immune", "Immune","Mesenchymal", "Immune", "Epithelial", "Mesenchymal","Mesenchymal", "Immune","Epithelial","Mesenchymal","Mesenchymal","Epithelial","Fetal lung","Immune","Mesenchymal","Epithelial","Immune","Mesenchymal","Fetal lung","Immune","Immune") 
names(new.cluster.ids) <- levels(E_seurat1_bb)
E_seurat1_bb <- RenameIdents(E_seurat1_bb, new.cluster.ids)
DimPlot(E_seurat1_bb, reduction = "bbknn", label = TRUE, pt.size = 0.5) 

#bbknn1_Epithelial cells
E_seurat1_epi<-subset(E_seurat1_bb, idents = c("Epithelial", "Fetal lung"))

E_seurat1_epi <-NormalizeData(E_seurat1_epi) #생략.?
E_seurat1_epi <-FindVariableFeatures(E_seurat1_epi)
E_seurat1_epi <-ScaleData(E_seurat1_epi)
E_seurat1_epi <- RunPCA(E_seurat1_epi, npcs = 50)
E_seurat1_epi <- JackStraw(E_seurat1_epi,dims = 50)
E_seurat1_epi <- ScoreJackStraw(E_seurat1_epi, dims = 1:50)
JackStrawPlot(E_seurat1_epi, dims = 1:50)
ElbowPlot(E_seurat1_epi, ndims = 50)

E_seurat1_epi <- FindNeighbors(E_seurat1_epi, dims = 1:15)
E_seurat1_epi <- FindClusters(E_seurat1_epi, resolution = 0.7)
E_seurat1_epi <- RunTSNE(E_seurat1_epi, dims = 1:15)
E_seurat1_epi <- RunUMAP(E_seurat1_epi, dims = 1:15)

E_seurat1_epi<- E_seurat1_epi %>% RunBBKNN(dims.use = 1:15,
                                      neighbors_within_batch = 5,
                                      trim = 20,
                                      batch.key = "batch",
                                      python.path = "/home/users/yunah1029/anaconda3/bin/python") 
DimPlot(E_seurat1_epi, reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; k=5, trim=20, epi")
DimPlot(E_seurat1_epi, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; k=3, trim=30, epi")

pdf("bbknn1_test.pdf", 6,5)
for(k in c(3:10)){
  for(trim in c(20,25,30,35,40,45)){
    print(
      E_seurat1_epi %>% RunBBKNN(dims.use = 1:15,
                                 neighbors_within_batch = k,
                                 trim = trim,
                                 batch.key = "batch",
                                 python.path = "/home/users/yunah1029/anaconda3/bin/python") %>%  
        DimPlot(reduction = "bbknn", group.by = "batch") %>% 
        LabelClusters(id="batch") + ggtitle(paste0("bbknn1; neighbors=",k,",trim=",trim,",epi"))
    )
  }
}
dev.off()





#cluster 찾기
group1 <- c('11')
group2 <- NULL  
marker_dt <- FindMarkers(E_seurat1_epi, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
marker_11<-marker_dt

FeaturePlot(E_seurat1_epi, "ENSMUSG00000020932", reduction = "bbknn")

#re-clustering lung_epi
new.cluster.ids <- c("AT2","AT2","AT2","Clara","AT2","Fetal lung","AT1","AT1","Ciliated","AT2","Fetal lung","AT1","AT2") 
names(new.cluster.ids) <- levels(E_seurat1_epi)
E_seurat1_epi <- RenameIdents(E_seurat1_epi, new.cluster.ids)
DimPlot(E_seurat1_epi, reduction = "bbknn", label = TRUE, pt.size = 0.5) + ggtitle("bbknn1; k=4, trim=30, epi")

#bbknn2; epi 
E_seurat2_epi<-subset(E_seurat1_bb, idents = c("Epithelial", "Fetal lung"))
E_seurat2_epi<- E_seurat2_epi%>% NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()
  
E_seurat2_epi_bb<-E_seurat2_epi%>%FindBatchBalancedNeighbors(batch.key = "batch",
                                                     neighbors_within_batch = 3,
                                                     prune.SNN = 1/15)   
E_seurat2_epi_bb<-E_seurat2_epi_bb %>% FindClusters(graph.name = "bbknn_snn", resolution = 1.2) 
E_seurat2_epi_bb<-E_seurat2_epi_bb %>% RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_") 

DimPlot(E_seurat2_epi_bb, reduction = "bbknn", group.by = "batch")%>%LabelClusters(id = "batch") + ggtitle("bbknn2; k=3, prune= 1/15, resolution=1.2, all")

DimPlot(E_seurat2_epi_bb, reduction = "bbknn", group.by = "ident")%>%LabelClusters(id = "ident") + ggtitle("bbknn2; k=5, prune= 1/15, resolution=1.2, all")

#bbknn2_up; epi
E_seurat2_up_epi<-subset(E_seurat1_bb, idents = c("Epithelial", "Fetal lung"))
E_seurat2_up_epi<- E_seurat2_up_epi%>% NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

E_seurat2_up_epi_bb<- E_seurat2_up_epi %>% FindBatchBalancedNeighbors(batch.key = "batch",
                                                           neighbors_within_own = 4,
                                                           neighbors_within_batch = 3,
                                                           prune = 1/15) 
E_seurat2_up_epi_bb <- E_seurat2_up_epi_bb %>%  FindClusters(graph.name = "bbknn_snn", resolution = 1) 
E_seurat2_up_epi_bb <- E_seurat2_up_epi_bb %>%  RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_")

DimPlot(E_seurat2_up_epi_bb, reduction = "bbknn", group.by = "batch")%>%LabelClusters(id = "batch") + ggtitle("bbknn2_up;k=4(3), prune= 1/15, resolution = 1,epi(b1)")
DimPlot(E_seurat2_up_epi_bb, reduction = "bbknn", group.by = "ident")%>%LabelClusters(id = "ident") + ggtitle("bbknn2_up;k=4(3), prune= 1/15, resolution = 1.2,epi(b1)")












test_ver2 <- E_seurat1_bb[,E_seurat1_bb$batch %in% c("P6to10.2","PND1.1","E18.5","P15","E14.5.2","E16.5","E8.25")]
test_ver2 <- test_ver2 %>% RunBBKNN(dims.use = 1:30,
                                    neighbors_within_batch = 3,
                                    trim = 18,
                                    batch.key = "batch",
                                    python.path = "/home/users/yunah1029/anaconda3/bin/python")
DimPlot(test_ver2, reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; k=3, trim=18, epi")

