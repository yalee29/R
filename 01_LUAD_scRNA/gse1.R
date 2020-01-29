#For extracting neuroendocrine cell and basal cells
#Data load(GSE1003354)
prox_epi<-read_tsv("GSE103354_Trachea_droplet_UMIcounts.tsv")
prox_epi %<>% column_to_rownames(var = "M1_GCTTGAGAAGTCGT_Club")
head(prox_epi[,1:10])
prox_epi %<>% rownames_to_column(var = "Gene")
head(prox_epi[,1:10])
head(gene_name_ensembl_mouse)

prox_epi_ensembl <- merge(prox_epi, gene_name_ensembl_mouse, by.x = "Gene", by.y = "Gene name")
head(prox_epi_ensembl[,1:10])
dim(prox_epi)
dim(prox_epi_ensembl) #(17701genes/ 7192 cells)
head(prox_epi_ensembl[ncol(prox_epi_ensembl)])
prox_epi_ensembl <- prox_epi_ensembl[,-1]
prox_epi_ensembl %<>% column_to_rownames(var = "Gene stable ID")

#Create Seurat Object
prox_epi_seurat = CreateSeuratObject(prox_epi_ensembl)
prox_epi_seurat$batch = "P6to12"
prox_epi_seurat %<>% NormalizeData()
prox_epi_seurat %<>% FindVariableFeatures()
prox_epi_seurat %<>% ScaleData(features = rownames(prox_epi_seurat@assays$RNA@counts))
dim(prox_epi_seurat@assays$RNA@scale.data)
prox_epi_seurat %<>% RunPCA(npcs = 50)
ElbowPlot(prox_epi_seurat,ndims = 50 )
prox_epi_seurat %<>% FindNeighbors(dims = 1:30)
prox_epi_seurat %<>% FindClusters(dims= 1:30, resolution = 0.7)
prox_epi_seurat %<>% RunTSNE(dims = 1:30)
prox_epi_seurat %<>% RunUMAP(dims = 1:30)
DimPlot(prox_epi_seurat, reduction = "tsne", group.by = "ident", label = TRUE)
DimPlot(prox_epi_seurat, reduction = "umap", group.by = "ident", label = TRUE)

#Find Neuroendocrine & basal cell by known markers (refer to article)
##neuroendocrine
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000021194","ENSMUSG00000030669", "ENSMUSG00000021194")) #68cells/ cluster 12
##basal cells
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000028435", "ENSMUSG00000061527", "ENSMUSG00000022510")) #
##Ionocyte
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000035951", "ENSMUSG00000047861")) #17 cells
##Tuft
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000027797")) #124
##Goblet
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000030954")) #59
#total
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000021194","ENSMUSG00000028435", "ENSMUSG00000047861","ENSMUSG00000027797","ENSMUSG00000030954")) #Basal/PNEC/Ionocyte/Tuft/Goblet

new.cluster <- c("Basal","Basal","Basal","Clara","Basal","Clara","Clara","Clara","Ciliated","Clara","Basal","Tuft","PNEC","Goblet","Ionocyte")
names(new.cluster) <- levels(prox_epi_seurat)
prox_epi_seurat <- RenameIdents(prox_epi_seurat, new.cluster)
DimPlot(prox_epi_seurat, reduction = "umap", group.by = "ident", label = TRUE)

prox_epi.markers <- FindAllMarkers(prox_epi_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# test_merge: pnec + E_seurat1_ad_epi_bb
##basal downsampling to 500 cells
basal <- prox_epi_seurat[,prox_epi_seurat@active.ident == "Basal"]
random_basal<-sample(colnames(basal), 300, replace = FALSE)
length(unique(random_basal))
basal <- basal[,colnames(basal) %in% random_basal]  
dim(basal@assays$RNA@counts)
##pnec
pnec <- prox_epi_seurat[,(prox_epi_seurat@active.ident == "PNEC")]

##merge basal + pnec
basal_pnec <- merge(basal, pnec)
