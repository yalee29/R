#bbknn2; R에서 bbknn
library(Seurat)
library(tidyverse)

.batchrank <- function(v,b) {
  r = v
  for (k in unique(b)) {
    r[b == k] = rank(v[b == k])
  }
  r
}

.batchrankmat <- function(E, b) {
  D = 1-.sparse.cor(E) # 1-cor(as.matrix(E))
  R = D
  for (i in 1:ncol(R)) {
    R[i,] = .batchrank(D[,i],b)
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
                                       neighbors_within_batch = 6,
                                       ...) {
  E <- GetAssayData(object)[VariableFeatures(object),]
  b <- object@meta.data[[batch.key]]
  k <- neighbors_within_batch * length(unique(b))
  R <- .batchrankmat(E, b)
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
                             neighbors_within_batch = 6) %>%
  FindClusters(graph.name = "bbknn_snn") %>%
  RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_")
')
