#For extracting neuroendocrine cell and basal cells

#1. GSE1003354
##Data load(GSE1003354)
prox_epi<-read_tsv("GSE103354_Trachea_droplet_UMIcounts.tsv")
prox_epi <- prox_epi %>% column_to_rownames(var = "M1_GCTTGAGAAGTCGT_Club")
head(prox_epi[,1:10])
prox_epi <- prox_epi %>%  rownames_to_column(var = "Gene")
head(prox_epi[,1:10])
head(gene_name_ensembl_mouse)

prox_epi_ensembl <- merge(prox_epi, gene_name_ensembl_mouse, by.x = "Gene", by.y = "Gene name")
head(prox_epi_ensembl[,1:10])
dim(prox_epi)
dim(prox_epi_ensembl) #(17701genes/ 7192 cells)
head(prox_epi_ensembl[ncol(prox_epi_ensembl)])
prox_epi_ensembl <- prox_epi_ensembl[,-1]
prox_epi_ensembl <- prox_epi_ensembl %>%  column_to_rownames(var = "Gene stable ID")

##Create Seurat Object
prox_epi_seurat = CreateSeuratObject(prox_epi_ensembl)
prox_epi_seurat$batch = "P6to12"
prox_epi_seurat <- prox_epi_seurat %>% NormalizeData()
prox_epi_seurat <- prox_epi_seurat %>% FindVariableFeatures()
prox_epi_seurat <- prox_epi_seurat %>% ScaleData()
dim(prox_epi_seurat@assays$RNA@scale.data)
prox_epi_seurat <- prox_epi_seurat %>% RunPCA(npcs = 50)
ElbowPlot(prox_epi_seurat,ndims = 50 )
prox_epi_seurat <- prox_epi_seurat %>%FindNeighbors(dims = 1:30)
prox_epi_seurat <- prox_epi_seurat %>% FindClusters(dims= 1:30, resolution = 0.7)
prox_epi_seurat <- prox_epi_seurat %>% RunTSNE(dims = 1:30)
prox_epi_seurat <- prox_epi_seurat %>% RunUMAP(dims = 1:30)
DimPlot(prox_epi_seurat, reduction = "tsne", group.by = "ident", label = TRUE)
DimPlot(prox_epi_seurat, reduction = "umap", group.by = "ident", label = TRUE)

##Find Neuroendocrine & basal cell by known markers (refer to article)
###neuroendocrine
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000021194","ENSMUSG00000030669", "ENSMUSG00000021194")) #68cells/ cluster 12
##basal cells
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000028435", "ENSMUSG00000061527", "ENSMUSG00000022510")) #
###Ionocyte
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000035951", "ENSMUSG00000047861")) #17 cells
###Tuft
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000027797")) #124
###Goblet
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000030954")) #59
###total
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000021194","ENSMUSG00000028435", "ENSMUSG00000047861","ENSMUSG00000027797","ENSMUSG00000030954")) #Basal/PNEC/Ionocyte/Tuft/Goblet
new.cluster <- c("Basal","Basal","Basal","Clara","Basal","Clara","Clara","Clara","Ciliated","Clara","Basal","Tuft","PNEC","Goblet","Ionocyte")
names(new.cluster) <- levels(prox_epi_seurat)
prox_epi_seurat <- RenameIdents(prox_epi_seurat, new.cluster)
DimPlot(prox_epi_seurat, reduction = "umap", group.by = "ident", label = TRUE)

prox_epi.markers <- FindAllMarkers(prox_epi_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

##basal downsampling to 500 cells (수가 너무 많음)
basal <- prox_epi_seurat[,prox_epi_seurat@active.ident == "Basal"]
random_basal<-sample(colnames(basal), 300, replace = FALSE)
length(unique(random_basal))
basal <- basal[,colnames(basal) %in% random_basal]  
dim(basal@assays$RNA@counts)
##pnec
pnec <- prox_epi_seurat[,(prox_epi_seurat@active.ident == "PNEC")]

##merge basal + pnec
basal_pnec <- merge(basal, pnec)

#2. GSE102580
##GSE102580 data load
id_55 <- read_tsv("./projects/01_LUAD_scRNA/mouse_traj/PNEC_basal_scRNA(mouse)/GSE102580/GSM2741420_Nov2016_2_Uninjured_mouse_id_55.counts.tsv")
id_56 <- read_tsv("./projects/01_LUAD_scRNA/mouse_traj/PNEC_basal_scRNA(mouse)/GSE102580/GSM2741421_Nov2016_15_uninjured_mouse_id_56.counts.tsv")
id_57 <- read_tsv("./projects/01_LUAD_scRNA/mouse_traj/PNEC_basal_scRNA(mouse)/GSE102580/GSM2741422_Nov2016_10_uninjured_mouse_id_57.counts.tsv")
dim(id_56)
dim(id_57)
id_55 <- id_55 %>% t()
id_56 <- id_56 %>% t()
id_57 <- id_57 %>% t()
colnames(id_55) <- as.character(unlist(id_55[1,]))
colnames(id_55) = paste0("P6to8.1_", colnames(id_55))
id_55 = id_55[-1, ]
colnames(id_56) <- as.character(unlist(id_56[1,]))
colnames(id_56) = paste0("P6to8.2_",colnames(id_56))
id_56 = id_56[-1, ]
colnames(id_57) <- as.character(unlist(id_57[1,]))
colnames(id_57) = paste0("P6to8.3_",colnames(id_57))
id_57 = id_57[-1, ]
id_55 <- id_55 %>% as.data.frame()
id_56 <- id_56 %>% as.data.frame()
id_57 <- id_57 %>% as.data.frame()
head(id_55[,1:10])
head(id_56[,1:10])
head(id_57[,1:10])

##Quality control 
id_s <- CreateSeuratObject(counts = id_57) #여기에 input만 바꿔서 넣어줄것
id_s[["mt"]] <- PercentageFeatureSet(id_s, pattern = "^mt-")
VlnPlot(id_s, features = c("nFeature_RNA", "nCount_RNA", "mt"), ncol = 3)
plot1 <- FeatureScatter(id_s, feature1 = "nCount_RNA", feature2 = "mt")
plot2 <- FeatureScatter(id_s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
id_55_s <- subset(id_s, subset = nFeature_RNA < 3000 & mt < 20) #1376
id_56_s <- subset(id_s, subset = nFeature_RNA < 3000 & mt < 20) #2745
id_57_s <- subset(id_s, subset = nFeature_RNA < 3000 & mt < 20) #2196
dim(id_55_s)
dim(id_56_s)
dim(id_57_s)

##merge gse102580 datas to gse2
gse2 <-merge(id_55_s, id_56_s)
gse2 <- merge(gse2, id_57_s)
dim(gse2)

##gse2 gene name -> mouse ensembl
gse2 <- gse2@assays$RNA@counts %>% as.data.frame()
dim(gse2)
gse2 <- gse2 %>% rownames_to_column(var = "Gene")
head(gse2[,1:10])
gse2 <- merge(gse2, gene_name_ensembl_mouse, by.x = "Gene", by.y = "Gene name")
dim(gse2)
head(gse2[,1:10])
head(gse2[ncol(gse2)])
gse2 <- gse2[,-1]
gse2 <- gse2 %>% column_to_rownames(var = "Gene stable ID")
head(gse2[,1:10])
batch_vector = stringr::str_replace(colnames(gse2), "_.*","")

##merge gse1(gse103354) data + gse10280 data
prox_pnec_basal <- gse2 %>% CreateSeuratObject()
prox_pnec_basal$batch = batch_vector
prox_pnec_basal <- merge(prox_pnec_basal, basal_pnec)
dim(prox_pnec_basal)

rm(id_55)
rm(id_56)
rm(id_57)
rm(id_55_s)
rm(id_56_s)
rm(id_57_s)
rm(gse2)

##bbknn on prox_pnec_basal(batch; P6to12/ P6to8.1/ P6to8.2/ P6to8.3) -> PNEC, basal cell 최대한 추출 -> E traj와 다시 bbknn으로 합쳐줄것. 
prox_pnec_basal <- prox_pnec_basal %>% NormalizeData()
prox_pnec_basal <- prox_pnec_basal %>% FindVariableFeatures()
prox_pnec_basal <- prox_pnec_basal %>% ScaleData()
prox_pnec_basal <- prox_pnec_basal %>% RunPCA(npcs = 50)
ElbowPlot(prox_pnec_basal, ndims = 50)
prox_pnec_basal <- prox_pnec_basal %>% FindNeighbors(dims = 1:30)
prox_pnec_basal <- prox_pnec_basal %>% FindClusters(resolution = 0.7)
prox_pnec_basal <- prox_pnec_basal %>% RunUMAP(dims = 1:30)
DimPlot(prox_pnec_basal, reduction = "umap", group.by = "batch", label = TRUE)
DimPlot(prox_pnec_basal, reduction = "umap", group.by = "ident", label = TRUE)
FeaturePlot(prox_pnec_basal, reduction = "umap", c("ENSMUSG00000021194","ENSMUSG00000030669", "ENSMUSG00000021194"), max.cutoff = 5) 
FeaturePlot(prox_pnec_basal, reduction = "umap", c("ENSMUSG00000028435", "ENSMUSG00000061527", "ENSMUSG00000022510"), max.cutoff = 5) 

##basal, pnec subsetting
basal <- subset(prox_pnec_basal, idents = c("5", "3"))
pnec <- subset(prox_pnec_basal, idents = "10")
dim(basal) #1301
dim(pnec) #112

##basal down sampling
random_basal<-sample(colnames(basal), 300, replace = FALSE)
length(unique(random_basal))
basal <- basal[,colnames(basal) %in% random_basal]  
dim(basal@assays$RNA@counts)

##pnec + basal -> pnec_basal 
pnec_basal <- merge(pnec, basal)
pnec_basal <-NormalizeData(pnec_basal)
pnec_basal <-FindVariableFeatures(pnec_basal)
pnec_basal <-ScaleData(pnec_basal)
pnec_basal <- RunPCA(pnec_basal, npcs = 50)
ElbowPlot(pnec_basal, ndims = 50)
pnec_basal <- FindNeighbors(pnec_basal, dims = 1:20)
pnec_basal <- FindClusters(pnec_basal, resolution = 0.7) #resolution 변경한 것임. 
pnec_basal <- RunUMAP(pnec_basal, dims = 1:20)
DimPlot(pnec_basal, reduction = "umap", group.by = "batch")
FeaturePlot(pnec_basal, reduction = "umap", c("ENSMUSG00000021194","ENSMUSG00000030669", "ENSMUSG00000021194")) #68cells/ cluster 12
FeaturePlot(pnec_basal, reduction = "umap", c("ENSMUSG00000028435", "ENSMUSG00000061527", "ENSMUSG00000022510")) #

##Doublet check for pnec_basal
library(DoubletFinder)
### pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.pnec_basal <- paramSweep_v3(pnec_basal, PCs = 1:20, sct = FALSE)
sweep.stats_pnec_basal <- summarizeSweep(sweep.res.pnec_basal, GT = FALSE)
bcmvn_pnec_basal <- find.pK(sweep.stats_pnec_basal)

### Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotation <- pnec_basal@active.ident
homotypic.prop <- modelHomotypic(annotation)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*length(colnames(pnec_basal)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

### Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
pnec_basal <- doubletFinder_v3(pnec_basal, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
pnec_basal <- doubletFinder_v3(pnec_basal, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

DimPlot(pnec_basal, reduction = "umap", group.by = "batch")
DimPlot(pnec_basal, reduction = "umap", group.by = "DF.classifications_0.25_0.09_31")

# 3.E_traj + pnec_basal
test_merge <- merge(E_seurat1_ad_epi_bb, pnec_basal)
test_merge <- test_merge %>%  NormalizeData()
test_merge <- test_merge %>% FindVariableFeatures()
test_merge <- test_merge %>% ScaleData(features = rownames(test_merge@assays$RNA@counts))
test_merge <- test_merge %>% RunPCA(npcs = 50)
ElbowPlot(test_merge, ndims = 50)
test_merge <- FindNeighbors(test_merge, dims = 1:30)
test_merge <- FindClusters(test_merge, resolution = 1)
test_merge <- RunUMAP(test_merge, dims = 1:30)

pdf("test_merge2.pdf", 7,5)
for(k in c(2:7)){
  for(trim in c(25,30,35,40,45)){
    print(
      test_merge %>% RunBBKNN(dims.use = 1:30,
                              neighbors_within_batch = k,
                              trim = trim,
                              batch.key = "batch2",
                              python.path = "/home/users/yunah1029/anaconda3/bin/python") %>%  
        DimPlot(reduction = "bbknn", group.by = "batch") %>% 
        LabelClusters(id="batch") + ggtitle(paste0("bbknn1; neighbors=",k,",trim=",trim,",add_PNEC&Basal"))
    )
  }
}
dev.off()

test_merge<- test_merge %>% RunBBKNN(dims.use = 1:30,
                                     neighbors_within_batch =  5,
                                     trim=33,
                                     batch.key = "batch",
                                     python.path = "/home/users/yunah1029/anaconda3/bin/python") 
DimPlot(test_merge, reduction = "bbknn", group.by = "batch") + ggtitle("bbknn1; dims=30, k=5, trim=33, epi_PNEC&Basal")
DimPlot(test_merge, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=30, k=5, trim=33, epi_PNEC&Basal")

FeaturePlot(test_merge[,test_merge$batch %in% c("E8.25","E9.5to11.5","E14.5.2","E16.5","E18.5","LungMap.E18.5")], c("ENSMUSG00000056370"), reduction = "bbknn")
FeaturePlot(test_merge[,test_merge$batch %in% c("E8.25","E9.5to11.5","E14.5.2","E16.5","E18.5","LungMap.E18.5")], reduction = "bbknn", c("ENSMUSG00000037846","ENSMUSG00000000216")) #
FeaturePlot(test_merge[,test_merge$batch %in% c("E8.25","E9.5to11.5","E14.5.2","E16.5","E18.5")], c("ENSMUSG00000042784","ENSMUSG00000028583"), reduction = "bbknn")
FeaturePlot(test_merge, reduction = "umap", c("ENSMUSG00000021194","ENSMUSG00000030669", "ENSMUSG00000021194"), reduction = "bbknn") #68cells/ cluster 12

## Ident 부여
plot <- DimPlot(test_merge[,test_merge$batch %in% c("E8.25","E9.5to11.5","E14.5.2","E16.5","E18.5","LungMap.E18.5")], reduction = "bbknn", group.by = "ident") + ggtitle("bbknn1; dims=30, k=5, trim=33, epi_PNEC&Basal")
plot <- DimPlot(test_merge, reduction = "bbknn", group.by = "ident") + ggtitle("bbknn1; dims=30, k=5, trim=33, epi_PNEC&Basal")
plot <-FeaturePlot(test_merge[,test_merge$batch %in% c("E8.25","E9.5to11.5","E14.5.2","E16.5","E18.5","LungMap.E18.5")], reduction = "bbknn", c("ENSMUSG00000037846")) #
plot
select.fetal <- CellSelector(plot = plot)
Idents(test_merge, cells = select.fetal) <- "Fetal_AT1"

test_merge <- RenameIdents(test_merge, '0' = 'AT2', '1' = 'AT2', '2' = 'AT2', '3' = 'AT1', '4' = 'AT1', '5' = 'Clara', '6' = 'Ciliated', '7' = 'AT2', '8' = 'Basal', '10' = 'AT2', '11' = 'AT1', '12' = 'Basal', '13' = 'PNEC', '14' = 'AT2','15' = 'AT2', '17' = 'AT1')
DimPlot(test_merge, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=30, k=5, trim=33, epi_PNEC&Basal")
test_merge <- RenameIdents(test_merge, '9' = 'Fetal lung', '16' = 'Fetal lung', '18' = 'Fetal lung') 
DimPlot(test_merge, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=30, k=5, trim=33, epi_PNEC&Basal")
markers.all <- FindAllMarkers(test_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Gene set 설정; 최종 = 792 genes
ad_top100 <- markers.all %>% group_by(cluster) %>% top_n(n = 120, wt = avg_logFC) #840
DoHeatmap(test_merge, features = ad_top100$gene) 
ad_top100_rm_nonepi<-setdiff(ad_top100$gene, nonepi_marker_dt$gene_name) #658

##epi marker + lung subtype markers = final_genes 
ad_fianl_genes <- union(ad_top100_rm_nonepi, cur_markers$gene_name ) #859
ad_final_genes<-setdiff(ad_fianl_genes, luad_genes) #luad_gene들은 제거/ 그대로 유지.

##Final_gene to human ensembl; 813
ad_final_genes_df <- ad_final_genes%>%as.data.frame()
colnames(ad_final_genes_df) <- "epi_mouse"
ad_final<-merge(ad_final_genes_df, ensembl_mouse_human, by.x = "epi_mouse", by.y = "mouse_ensembl")
table(duplicated(ad_final$human_ensembl))
ad_final<- ad_final[!duplicated(ad_final[2]),]
head(ad_final)
names(ad_final)[1]<-"mouse_ensembl"

##+final gene 수정; 792
## cell cycle related gene 21개 + luad related gene 1개 제거
ad_final_test <- ad_final %>% filter(human_ensembl %in% setdiff(ad_final$human_ensembl, c("ENSG00000118971","ENSG00000039068","ENSG00000198668","ENSG00000135446","ENSG00000105976","ENSG00000171848","ENSG00000073350","ENSG00000100714","ENSG00000150753","ENSG00000111057","ENSG00000148834","ENSG00000144381","ENSG00000155368","ENSG00000173207","ENSG00000010278","ENSG00000075131","ENSG00000129757","ENSG00000105968","ENSG00000188486","ENSG00000160752","ENSG00000110092","ENSG00000147883")))
ad_final<- ad_final_test

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#single cell data processing -> with 600 genes (expression data)
mouse_epi_traj<-test_merge@assays$RNA@data %>%as.data.frame() # cell; 6315 / gene; 28695
mouse_epi_traj<-mouse_epi_traj%>%rownames_to_column()
names(mouse_epi_traj)[1] <- c("ensembl")
mouse_epi_600 <- merge(mouse_epi_traj, ad_final, by.x = "ensembl", by.y = "mouse_ensembl") #792genes/ 6315cells
mouse_epi_600 <- mouse_epi_600[,-1]

#mouse traj knn pooling (sparcity 보정)
mouse_epi_600 <- mouse_epi_600%>%column_to_rownames(var = "human_ensembl")
mouse_epi_600_knn <- knn(mouse_epi_600, 5)

##Plot
mouse_epi_traj<-test_merge
group1 <- c('Fetal lung')
group2 <- c('AT1')
group3 <- c('AT2')
group4 <- c('Ciliated')
group5 <- c('Clara')
group6 <- c('Fetal_AT1')
group7 <- c('Fetal_AT2')
group8 <- c('Fetal_Clara')
group9 <- c('Fetal_Ciliated')
group10 <- c('Basal')
group11 <- c('PNEC')
cell1 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group1])
cell2 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group2])
cell3 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group3])
cell4 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group4])
cell5 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group5])
cell6 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group6])
cell7 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group7])
cell8 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group8])
cell9 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group9])
cell10 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group10])
cell11 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group11])
bbknn_dt <- mouse_epi_traj@reductions$bbknn@cell.embeddings %>% as.data.frame() %>% rownames_to_column('cell_id')
bbknn_dt <- bbknn_dt %>% mutate(group = ifelse(cell_id %in% cell1, paste(group1,collapse=','), 
                                               ifelse(cell_id %in% cell2, paste(group2, collapse=','),
                                                      ifelse(cell_id %in% cell3, paste(group3, collapse=','),
                                                             ifelse(cell_id %in% cell4, paste(group4, collapse=','),
                                                                    ifelse(cell_id %in% cell5, paste(group5, collapse=','),
                                                                           ifelse(cell_id %in% cell6, paste(group6, collapse=','),
                                                                                  ifelse(cell_id %in% cell7, paste(group7, collapse=','),
                                                                                         ifelse(cell_id %in% cell8, paste(group8, collapse=','),
                                                                                                ifelse(cell_id %in% cell9, paste(group9, collapse=','),
                                                                                                       ifelse(cell_id %in% cell10, paste(group10, collapse=','),
                                                                                                              "PNEC")))))))))))
ggplot(bbknn_dt, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), size = 0.8, alpha = 0.7)+
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#87d4f5","#7dab78", "#f57c16",
                                 "#F0E442", "#bf046a","#a8850f","#0e5e05","#9e5416","#cf99f7")) +
  theme_classic()

#Human normal 

###-### bbknn(1)_final_troubleshooting 과 연계 -> plot
mouse_epi_600_knn <- mouse_epi_600_knn %>% as.data.frame() %>% rownames_to_column(var = "human_ensembl")
hu_AT1 <- mouse_epi_600_knn$human_ensembl %>% as.data.frame()
hu_AT1$id  <- 1:nrow(hu_AT1)
colnames(hu_AT1) <- c("gene","id")
hu_AT1 <- merge(hu_AT1, donor_AT1, by.x = "gene", by.y = "gene", all.x = TRUE, all.y = FALSE) #여기 바꾸면서 조절
hu_AT1[is.na(hu_AT1)] <- 0
sum(is.na(hu_AT1))
hu_AT1 <- hu_AT1[order(hu_AT1$id),]
rownames(hu_AT1) <-NULL
hu_AT1 <- hu_AT1 %>% as.data.frame() %>% column_to_rownames(var = "gene")
hu_AT1 <- hu_AT1[,-1]
head(rownames(hu_AT1))
head(colnames(hu_AT1))
dim(hu_AT1)
mouse_epi_600_knn <- mouse_epi_600_knn%>%column_to_rownames(var = "human_ensembl")
###

##human lung cluster마다 mouse traj에 붙여보기 (#표시한곳 조건 바꿔가면서 cluster별 그릴 수 있음)  #전과정은 bbknn(1)_final_epi.R 에 정리돼있음. 여기는 knn pooling한 data와 human noraml lung을 비교하기 위함
m_h_AT1 <- cbind(mouse_epi_600_knn, hu_AT1)
cor_m_h_AT1 <- cor(m_h_AT1, method = "pearson") #method = "spearman"/ pearson(default)
cor_m_h_AT1 <- cor_m_h_AT1 %>% as.data.frame()
cor_m_h_AT1_sort <- cor_m_h_AT1[1:6315,6316:nrow(cor_m_h_AT1)]
cor_m_h_AT1_max <- apply(cor_m_h_AT1_sort,2, max) #correlation cutoff 설정 
cor_m_h_AT1_sort <- cor_m_h_AT1_sort[,which(cor_m_h_AT1_max > 0.45)]
cor_m_h_AT1_sort_rank <- apply(cor_m_h_AT1_sort,2,function(x){order(-x)})[1:7,]%>%t()
cor_m_h_AT1_sort_rank <- cor_m_h_AT1_sort_rank %>% as.data.frame()

##
a<-Embeddings(test_merge,reduction = "bbknn") %>% as.data.frame()
cor_m_h_AT1_sort_rank <- cor_m_h_AT1_sort_rank[,1:3] %>% t()
matched <- data.frame()

for(n in 1:ncol(cor_m_h_AT1_sort_rank)){
  tmp_matched <- a[cor_m_h_AT1_sort_rank[,n],]
  tmp_matched$group <- "h_AT1" #여기도 cluster 바뀔때마다 조절
  matched<- rbind(matched, tmp_matched)
}

matched <- matched %>% rownames_to_column(var = "cell_id")

bbknn_dt_plot <- rbind(bbknn_dt, matched)
ggplot(bbknn_dt_plot, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), size = 1, alpha = 0.7)+
  #scale_shape_manual(values = c( 16, 16, 16,16, 16, 2, 17))+
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#87d4f5","#7dab78", "#f57c16",
                                 "#F0E442", "#bf046a","#a8850f","#0e5e05","#9e5416","#3e0ffa","#cf99f7"))+
  theme_classic()




#LUAD_bulk_RNAseq with 600 genes (expression data_FPKM)
dir <- ("./projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUAD_551_FPKM/")
file_list<-list.files("/home/users/yunah1029/projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUAD_551_FPKM",, recursive=TRUE, pattern="*.FPKM.txt.gz$")

tcga <- mouse_epi_600$human_ensembl %>% as.data.frame()
tcga$id  <- 1:nrow(tcga)
colnames(tcga) <- "gene"

for (file in file_list) {
  temp <- read.csv(paste0(dir, file), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = FALSE)
  colnames(temp) <- c("gene", strsplit(file, "/")[[1]][1])
  temp$gene <- sapply(strsplit(temp$gene, "[.]"), function(x) x[1])
  tcga <- merge(tcga, temp, by.x = "gene", by.y = "gene", all.x = TRUE, all.y = FALSE)
}
tcga[is.na(tcga)] <- 0
sum(is.na(tcga))

#tcga sample 별 mouse_epi_600과 cor -> 가장 가까운 cell 찾기. 
tcga<-tcga[order(tcga$"NA"), ] #mouse_epi_600의 gene 순서와 동일하게 정렬 -> 나중에 cbind
rownames(tcga) <- NULL
tcga<- tcga %>% as.data.frame() %>% column_to_rownames(var = "gene")
tcga<-tcga[,-1]
head(rownames(tcga))
head(colnames(tcga))

mouse_epi_600 <- mouse_epi_600%>%column_to_rownames(var = "human_ensembl")


