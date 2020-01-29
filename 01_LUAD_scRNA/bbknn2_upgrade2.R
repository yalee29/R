# bbknn2 upgrade ver2; UMAP 문제 해결
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
  require(progress)
  D = 1-.sparse.cor(E)
  R = D
  pb <- progress_bar$new(total = ncol(R))
  for (i in 1:ncol(R)) {
    R[i,] = .batchrank(D[,i],b,i,d)
    pb$tick()
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
                                       check.dup = T,
                                       seed = 42,
                                       rm.dup = T,
                                       ...) {
  set.seed(seed)
  E <- GetAssayData(object)[VariableFeatures(object),]
  b <- object@meta.data[[batch.key]]
  k <- neighbors_within_batch * (length(unique(b)) - 1) + neighbors_within_own
  message("k is ",k)
  d <- neighbors_within_own - neighbors_within_batch
  message("Calculating Pearson correlation coefficients for all pairs of cells ...")
  R <- .batchrankmat(E, b, d)
  O <- Seurat::FindNeighbors(R,
                             distance.matrix = T,
                             k.param = k,
                             force.recalc = T, ...)
  if(check.dup) {
    message("Checking duplicates in batch balanced SNN matrix ...")
    dup <- duplicated(as.matrix(O$snn))
    if (any(dup)) {
      message(sum(dup), " duplicated cells in batch balanced SNN matrix is detected.")
      if(rm.dup) {
        message("Duplicated cells will be removed")
        object <- object[,-dup]
        object@graphs$bbknn_nn = O$nn[-dup,-dup]
        object@graphs$bbknn_snn = Seurat::as.Graph(O$snn[-dup,-dup])
        return(object)
      } else {
        message("Random small noise will be added to these cells")
        
        dup_cols <- O$snn[,which(dup)]
        non_zero_dup_vals <- dup_cols[dup_cols>0&dup_cols<1]
        noise <- rnorm(length(non_zero_dup_vals), mean = 0, sd = 0.0001)
        dup_cols_with_noise <- dup_cols
        dup_cols_with_noise[dup_cols>0&dup_cols<1] <- non_zero_dup_vals + noise
        dup_cols_with_noise[dup_cols_with_noise<0] <- 0
        dup_cols_with_noise[dup_cols_with_noise>1] <- 1
        O$snn[,which(dup)] <- dup_cols_with_noise
        O$snn[which(dup),] <- Matrix::t(dup_cols_with_noise)
        
        object@graphs$bbknn_nn = O$nn
        object@graphs$bbknn_snn = Seurat::as.Graph(O$snn)
        return(object)        
      }
    } else {
      message("No duplicates!")
      object@graphs$bbknn_nn = O$nn
      object@graphs$bbknn_snn = Seurat::as.Graph(O$snn)
      return(object)
    }
  } else {
    object@graphs$bbknn_nn = O$nn
    object@graphs$bbknn_snn = Seurat::as.Graph(O$snn)
    return(object)
  }
}



RunUMAP.graph.mod <- function(
  object, # graph (for bbknn-snn)
  assay = NULL,
  n.components = 2L,
  metric = "correlation",
  n.epochs = 0L,
  learning.rate = 1,
  min.dist = 0.3,
  spread = 1,
  repulsion.strength = 1,
  negative.sample.rate = 5L,
  a = NULL,
  b = NULL,
  seed.use = 42L,
  metric.kwds = NULL,
  verbose = TRUE,
  reduction.key = 'BBKNN_',
  init = "random",
  ...
) {
  require(reticulate)
  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
  }
  if (!py_module_available(module = 'numpy')) {
    stop("Cannot find numpy, please install through pip (e.g. pip install numpy).")
  }
  if (!py_module_available(module = 'sklearn')) {
    stop("Cannot find sklearn, please install through pip (e.g. pip install scikit-learn).")
  }
  if (!py_module_available(module = 'scipy')) {
    stop("Cannot find scipy, please install through pip (e.g. pip install scipy).")
  }
  set.seed(seed.use)
  np <- import("numpy", delay_load = TRUE)
  sp <- import("scipy", delay_load = TRUE)
  sklearn <- import("sklearn", delay_load = TRUE)
  umap <- import("umap", delay_load = TRUE)
  diag(x = object) <- 0
  data <- object
  object <- sp$sparse$coo_matrix(arg1 = object)
  ab.params <- umap$umap_$find_ab_params(spread = spread, min_dist = min.dist)
  a <- a %||% ab.params[[1]]
  b <- b %||% ab.params[[2]]
  n.epochs <- n.epochs %||% 0L
  random.state <- sklearn$utils$check_random_state(seed = as.integer(x = seed.use))
  embeddings <- umap$umap_$simplicial_set_embedding(
    data = data,
    graph = object,
    n_components = n.components,
    initial_alpha = learning.rate,
    a = a,
    b = b,
    gamma = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    n_epochs = as.integer(x = n.epochs),
    random_state = random.state,
    init = init,
    metric = metric,
    metric_kwds = metric.kwds,
    verbose = verbose
  )
  rownames(x = embeddings) <- colnames(x = data)
  colnames(x = embeddings) <- paste0("BBKNN_", 1:n.components)
  # center the embeddings on zero
  embeddings <- scale(x = embeddings, scale = FALSE)
  umap <- CreateDimReducObject(embeddings = embeddings, key = reduction.key, assay = assay)
  return(umap)
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