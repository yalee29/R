E8.25_foregut_ensembl[1:3,1:3]
E9.5_11.5_ensembl[1:3,(ncol(E9.5_11.5_ensembl)-3):ncol(E9.5_11.5_ensembl)]


union_gene_vector = intersect(E8.25_foregut_ensembl$gene_id, E9.5_11.5_ensembl$`Gene stable ID`)


A = E8.25_foregut_ensembl[match(union_gene_vector, E8.25_foregut_ensembl$gene_id),-1]
colnames(A) = paste0("E8.25_", colnames(A))
B = E9.5_11.5_ensembl[match(union_gene_vector, E9.5_11.5_ensembl$`Gene stable ID`),-ncol(E9.5_11.5_ensembl)]
colnames(B) = paste0("E9.5to11.5_", colnames(B))
C = cbind(A,B)
rownames(C) = union_gene_vector
table(duplicated(t(C)))
table(!duplicated(t(C)))
C = C[,!duplicated(t(C))]
batch_vector = stringr::str_replace(colnames(C), "_.*","")

C_seurat = CreateSeuratObject(C)

C_seurat$batch = batch_vector

C_seurat <- C_seurat %>% 
  NormalizeData() %>% 
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  FindClusters() %>%
  RunTSNE() %>%
  RunUMAP(dims=1:30)

source("~kjyi/Projects/cav/bbknn.R")
C_seurat <- C_seurat %>% RunBBKNN(dims.use = 1:10, batch.key = "batch", python.path = "/home/users/yunah1029/anaconda3/bin/python")
umap.emb = read_tsv("umap.tmp.tsv", col_names = F) %>% as.data.frame
plot(umap.emb)


DimPlot(C_seurat, reduction = "tsne")
DimPlot(C_seurat, reduction = "umap", group.by = "batch")
C_seurat_bb = C_seurat 
colnames(umap.emb) = colnames(C_seurat_bb@reductions$umap@cell.embeddings)
rownames(umap.emb) = rownames(C_seurat_bb@reductions$umap@cell.embeddings)
C_seurat_bb@reductions$umap@cell.embeddings = as.matrix(umap.emb)
C_seurat_bb@reductions$umap@cell.embeddings %>% head
C_seurat@reductions$umap@cell.embeddings %>% head

DimPlot(C_seurat, reduction = "umap", group.by = "batch")
DimPlot(C_seurat_bb, reduction = "umap", group.by = "batch")




