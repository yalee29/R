#Trouble shooting

#I. PCC 확인
a <- Embeddings(E_seurat_epi,reduction = "bbknn") %>% as.data.frame
plot(a)
tmp_1 <- gatepoints::fhs(a)
tmp_2 <- gatepoints::fhs(a)
tmp_3 <- gatepoints::fhs(a)
tmp_4 <- which(E_seurat_epi$batch == "E18.5")

TMP_1 <- E_seurat_epi[E_seurat_epi@assays$RNA@var.features,tmp_1] %>% GetAssayData() %>% Matrix::rowMeans()
TMP_2 <- E_seurat_epi[E_seurat_epi@assays$RNA@var.features,tmp_2] %>% GetAssayData() %>% Matrix::rowMeans()
TMP_3 <- E_seurat_epi[E_seurat_epi@assays$RNA@var.features,tmp_3] %>% GetAssayData() %>% Matrix::rowMeans()
TMP_4 <- E_seurat_epi[E_seurat_epi@assays$RNA@var.features,tmp_4] %>% GetAssayData() %>% Matrix::rowMeans()

# 1 = embryo
# 2 = adult at1
# 3 = fetus
cor(cbind(TMP_1, TMP_2, TMP_3, TMP_4))

#II. cell cycle 보정
subset_E14.5.1 <- subset(E_seurat_epi, batch = "E14.5.1")
subset_E14.5.2 <- subset(E_seurat_epi, batch = "E14.5.2")
subset_E16.5 <- subset(E_seurat_epi, batch = "E16.5")
subset_P6_AT1<- subset(E_seurat_epi, idents = "22")
subset_PND1_AT1<-subset(E_seurat_epi, idents = "10")
subset_E18.5_AT1<-subset(E_seurat_epi, idents = "24")
cell_cycle_data <- merge(x= subset_E14.5.2, y= list(subset_E14.5.1,subset_E16.5, subset_E18.5_AT1, subset_P6_AT1, subset_PND1_AT1))
cell_cycle_data<-cell_cycle_data@assays$RNA@counts %>% as.data.frame() %>% rownames_to_column("gene_id")
cell_cycle_gene<-merge(cell_cycle_data,gene_name_ensembl_mouse,by.x = "gene_id",by.y = "Gene stable ID")
cell_cycle_gene<-cell_cycle_gene[,-1]
rownames(cell_cycle_gene) <- cell_cycle_gene$`Gene name`
cell_cycle_gene<-cell_cycle_gene[,-ncol(cell_cycle_gene)]
batch_vector_AT1 = stringr::str_replace(colnames(cell_cycle_gene), "_.*","")

cell_cycle_gene = CreateSeuratObject(cell_cycle_gene)
cell_cycle_gene$batch = batch_vector_AT1

rm(cell_cycle_data)

#cell cycle gene load
data(cc.genes)
s.gene <- cc.genes$s.genes %>% tolower() %>% {substr(.,1,1) <- toupper(substr(.,1,1));.}
g2m.gene <- cc.genes$g2m.genes %>% tolower() %>% {substr(.,1,1) <- toupper(substr(.,1,1));.}

#방법1
object_1 <- cell_cycle_gene%>% NormalizeData() %>% 
  FindVariableFeatures() %>%
  CellCycleScoring(s.features = s.gene, g2m.features = g2m.gene, set.ident = TRUE) %>% 
  ScaleData(vars.to.regress = c("S.Score","G2M.Score"), features = rownames(object_1)) %>% 
  RunPCA(features = VariableFeatures(object_1)) %>%
  FindNeighbors() %>%
  FindClusters() %>%
  RunTSNE()

#방법2
object_2 <- cell_cycle_data%>% NormalizeData() %>%
  FindVariableFeatures() %>%
  CellCycleScoring(s.features = s.gene, g2m.features = g2m.gene, set.ident = "TRUE") %>%
  {.$CC.Difference <- .$S.score - .$G2M.Score;.} %>%
  ScaleData(vars.to.regress = c("CC.Diference"), features = rownames(.)) %>%
  FindNeighbors() %>%
  FindClusters() %>%
  RunTSNE()

#III. marker gene 확인 
group1 <- c('14','22')
group2 <- NULL  
marker_dt <- FindMarkers(E_seurat_epi, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
