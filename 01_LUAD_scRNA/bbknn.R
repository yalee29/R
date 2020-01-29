#bbknn 실행 
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
