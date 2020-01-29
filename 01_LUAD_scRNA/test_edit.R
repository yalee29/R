#Batch balancing v2: basal, PNEC 추가할 때 ciliated cell, clara cell을 함께

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
prox_epi_seurat <- prox_epi_seurat %>% RunUMAP(dims = 1:30)
DimPlot(prox_epi_seurat, reduction = "umap", group.by = "ident", label = TRUE)

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
id_55_s <- subset(id_s, subset = nFeature_RNA < 2500 & mt < 20) #1262
id_56_s <- subset(id_s, subset = nFeature_RNA < 2500 & mt < 20) #2557
id_57_s <- subset(id_s, subset = nFeature_RNA < 2500 & mt < 20) #2056
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

##gse102580 data to seurat
prox_pnec_basal <- gse2 %>% CreateSeuratObject()
prox_pnec_basal$batch = batch_vector
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
prox_pnec_basal <- prox_pnec_basal %>% FindNeighbors(dims = 1:20)
prox_pnec_basal <- prox_pnec_basal %>% FindClusters(resolution = 0.7)
prox_pnec_basal <- prox_pnec_basal %>% RunUMAP(dims = 1:20)

#3.E_seurat_ad_epi_bb + PNEC & Basal-----------------------------------------------------------------------------------------------------------------------
#GSE datas for PNEC, Basal processing
##Ident
FeaturePlot(prox_pnec_basal, reduction = "umap", c("ENSMUSG00000021194","ENSMUSG00000030669", "ENSMUSG00000020052"), max.cutoff = 5) 
FeaturePlot(prox_pnec_basal, reduction = "umap", c("ENSMUSG00000028435", "ENSMUSG00000061527", "ENSMUSG00000022510","ENSMUSG00000022676"), max.cutoff = 5) 
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000021194","ENSMUSG00000030669", "ENSMUSG00000021194","ENSMUSG00000020052"), max.cutoff = 5) 
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000028435", "ENSMUSG00000061527", "ENSMUSG00000022510","ENSMUSG00000022676"), max.cutoff = 5) 

FeaturePlot(prox_pnec_basal, reduction = "umap", c("ENSMUSG00000034227"), max.cutoff = 5) #ciliated = cluster5
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000034227"), max.cutoff = 5) #ciliated = cluster8
FeaturePlot(prox_pnec_basal, reduction = "umap", c("ENSMUSG00000024653")) #club (scgb1a1)
FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000021795","ENSMUSG00000064057","ENSMUSG00000024653"), max.cutoff = 5) #club(sftpd)

DimPlot(prox_epi_seurat, reduction = "umap")
DimPlot(prox_pnec_basal, reduction = "umap")
#plot <-FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000024653")) #
plot <- DimPlot(prox_pnec_basal, reduction = "umap")
plot
select_cluster <- CellSelector(plot = plot)
Idents(prox_pnec_basal, cells = select_cluster) <- "Clara"

nrow(prox_pnec_basal@meta.data[prox_pnec_basal@active.ident =="Basal",])
nrow(prox_epi_seurat@meta.data[prox_epi_seurat@active.ident =="PNEC",])
nrow(prox_epi_seurat@meta.data[prox_epi_seurat@active.ident =="Ciliated",])
nrow(prox_pnec_basal@meta.data[prox_pnec_basal@active.ident =="Clara",])
DimPlot(prox_epi_seurat, reduction = "umap", group.by = "ident", label = TRUE)
DimPlot(prox_pnec_basal, reduction = "umap", group.by = "ident", label = TRUE)
test1 <-E_seurat1_ad_epi_bb
test6 <- subset(prox_pnec_basal, idents = c("Basal","Clara","Ciliated","PNEC")) #PNEC = 70/ Basal = 239/ Ciliated = 215/ Clara =172
test7 <- subset(prox_epi_seurat, idents = c("Basal","Clara","Ciliated","PNEC")) #PNEC = 70/ Basal = 239/ Ciliated = 148/ Clara = 153

##Preprocessing before merging mouse traj + PNEC & Basal(gene set) -> 11528 genes
gene_set <- intersect(intersect(rownames(test1), rownames(test2)), rownames(test4)) %>% as.data.frame() #11528
colnames(gene_set) <- "gene set"
dim(gene_set)
test1 <- E_seurat1_ad_epi_bb@assays$RNA@counts[rownames(E_seurat1_ad_epi_bb) %in% gene_set$`gene set`,] %>% CreateSeuratObject()
test6 <- test6@assays$RNA@counts[rownames(test6) %in% gene_set$`gene set`,] %>% CreateSeuratObject()
test7 <- test7@assays$RNA@counts[rownames(test7) %in% gene_set$`gene set`,] %>% CreateSeuratObject()

dim(test6)
dim(test7)

#test8 = test1 +test6 + test7 (batch balancing with basal, PNEC, ciliated, clara)-----------------------------------------------------------------------------------------
batch_vector = stringr::str_replace(colnames(test1), "_.*","")
test1$batch = batch_vector
test6$batch <- "P6to8"
test7$batch <- "P6to12"
test8 <- merge(test1, merge(test6, test7))
dim(test8)
unique(test8$batch)

test8 <-NormalizeData(test8)
test8 <-FindVariableFeatures(test8)
test8 <-ScaleData(test8, features = rownames(test8@assays$RNA@counts))
test8 <- RunPCA(test8, npcs = 50)
ElbowPlot(test8, ndims = 50)
test8 <- JackStraw(test8, dims = 30)
test8 <- ScoreJackStraw(test8, dims = 1:30)
JackStrawPlot(test8, dims = 1:30)
test8 <- FindNeighbors(test8, dims = 1:20)
test8 <- FindClusters(test8, resolution = 1) #resolution 변경한 것임. 
test8 <- RunUMAP(test8, dims = 1:20)
DimPlot(test8, reduction = "umap")

pdf("test8_20_pnec.pdf", 7,5)
for(k in c(4:7)){
  for(trim in c(25,27,30,33,35,37,40,43,45,47,50)){
    print(
      test8 %>% RunBBKNN(dims.use = 1:20,
                         neighbors_within_batch = k,
                         trim = trim,
                         batch.key = "batch",
                         python.path = "/home/users/yunah1029/anaconda3/bin/python") %>%  
        DimPlot(reduction = "bbknn", group.by = "batch") %>% 
        LabelClusters(id="batch") + ggtitle(paste0("bbknn1; neighbors=",k,",trim=",trim,",add_PNEC&Basal"))
    )
  }
}
dev.off()

test9 <- test8 %>% RunBBKNN(dims.use = 1:20,
                            neighbors_within_batch = 6,
                            trim = 41,
                            batch.key = "batch",
                            python.path = "/home/users/yunah1029/anaconda3/bin/python")
DimPlot(test9, reduction = "bbknn", group.by = "batch") + ggtitle(paste0("bbknn1; neighbors=",6,",trim=",41,",add_PNEC&Basal"))

##Clustering
test8 <- RenameIdents(object = test8, "0" = "AT2", "1" = "AT2", "2" = "AT2", "3" = "AT1", "4" = "AT1", "5" = "Clara", "6" = "Basal", "7" = "Ciliated", "8" = "Basal", "10" = "AT2", "11" = "AT2", "12" = "AT2", "13" = "PNEC", "15" = "AT2" ,"17" ="AT1")
test8 <- RenameIdents(object = test8, "14" = "Fetal lung")

FeaturePlot(test8, c("ENSMUSG00000023951"), reduction = "bbknn") #AT1 "ENSMUSG00000059325","ENSMUSG00000028583","ENSMUSG00000023039","ENSMUSG00000023951"
FeaturePlot(test8, c("ENSMUSG00000041247"), reduction = "bbknn") #AT2 "ENSMUSG00000029375","ENSMUSG00000037071","ENSMUSG00000042784"
FeaturePlot(test8, c("ENSMUSG00000038791","ENSMUSG00000031722"), reduction = "bbknn", max.cutoff = 8) #Clara "ENSMUSG00000022595"
FeaturePlot(test8, c("ENSMUSG00000027676","ENSMUSG00000040703"), reduction = "bbknn") #Ciliated "ENSMUSG00000040703", "ENSMUSG00000024653" ,"ENSMUSG00000027676"
FeaturePlot(test8, reduction = "bbknn", c("ENSMUSG00000021194","ENSMUSG00000030669", "ENSMUSG00000020052"), max.cutoff = 5) #PNEC
FeaturePlot(test8, reduction = "bbknn", c("ENSMUSG00000028435","ENSMUSG00000022676"), max.cutoff = 5) #BASAL

plot <-DimPlot(test8[,test8$batch %in% c("E8.25","E9.5to11.5","E14.5.2","E16.5","E18.5","LungMap.E18.5")], reduction = "bbknn", group.by = "ident")
plot
select_cluster <- CellSelector(plot = plot)
Idents(test8, cells = select_cluster) <- "Fetal lung"
