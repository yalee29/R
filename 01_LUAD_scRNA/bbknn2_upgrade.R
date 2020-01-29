# bbknn2_upgrade; 같은 batch, 다른 batch 각각 조절 
# bbknn less ugly version

library(Seurat)
library(tidyverse)

# k neighbors in the same batch of a cell,  k' from each other batches

.batchrank <- function(v,b,i=1, d=0) {
  r = v
  for (k in unique(b)) {
    if (k == b[i]) {
      r[b == k] = rank(v[b == k])
    } else {
      r[b == k] = rank(v[b == k]) + d
    }
  }
  r
}

.batchrankmat <- function(E, b, d = 0) {
  D = 1-.sparse.cor(E)
  R = D
  for (i in 1:ncol(R)) {
    R[i,] = .batchrank(D[,i],b,i,d)
  }
  R
}

.sparse.cor <- function(x){
  n <- nrow(x)
  cMeans <- Matrix::colMeans(x)
  cSums <- Matrix::colSums(x)
  covmat <- Matrix::tcrossprod(cMeans, (-2*cSums+n*cMeans))
  crossp <- as.matrix(Matrix::crossprod(x))
  covmat <- covmat+crossp
  sdvec <- sqrt(Matrix::diag(covmat)) # standard deviations of columns
  covmat/crossprod(t(sdvec)) # correlation matrix
}

FindBatchBalancedNeighbors <- function(object,
                                       batch.key = "batch",
                                       neighbors_within_own = 10,
                                       neighbors_within_batch = 6,
                                       ...) {
  E <- GetAssayData(object)[VariableFeatures(object),]
  b <- object@meta.data[[batch.key]]
  k <- neighbors_within_batch * (length(unique(b)) - 1) + neighbors_within_own
  d <- neighbors_within_own - neighbors_within_batch
  R <- .batchrankmat(E, b, d)
  O <- Seurat::FindNeighbors(R,
                             distance.matrix = T,
                             k.param = k,
                             force.recalc = T, ...)
  object@graphs$bbknn_nn = O$nn
  object@graphs$bbknn_snn = O$snn
  return(object)
}

cat('object_standard <- object %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  FindClusters() %>%
  RunUMAP()

object_bbknn <- object %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindBatchBalancedNeighbors(batch.key = "batch",
                             neighbors_within_own = 10,
                             neighbors_within_batch = 6) %>%
  FindClusters(graph.name = "bbknn_snn") %>%
  RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_")
')

