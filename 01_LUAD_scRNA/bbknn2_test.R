#결과 확인 
test <- E_seurat_epi[,E_seurat_epi$batch %in% c("P6to10.2","PND1.1","E18.5","P15","E14.5.2","E16.5","E8.25")]
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

#FindMarker
group1 <- c('0')
group2 <- NULL  
marker_dt <- FindMarkers(test_bbknn, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
marker_0 <- marker_dt



