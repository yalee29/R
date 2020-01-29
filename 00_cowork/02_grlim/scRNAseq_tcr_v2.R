library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
library(gplots)
library(cowplot)
library(ggplot2)
library(cowplot)

#Setup seurat objects 
data_dir <- './projects/00_cowork/02_krlim_scRNA/00_data/scRNAseq_microgen/02_align/10x5GEX_6/outs/filtered_feature_bc_matrix/' #dir 바꿔주면서 obj_1 ~ obj_6 생성
list.files(data_dir)
data_1 <- Read10X(data.dir = data_dir)
data_2 <- Read10X(data.dir = data_dir)
data_3 <- Read10X(data.dir = data_dir)
data_4 <- Read10X(data.dir = data_dir)
data_5 <- Read10X(data.dir = data_dir)
data_6 <- Read10X(data.dir = data_dir)
colnames(data_1) <- paste("pbmc", colnames(data_1), sep = "_")
colnames(data_2) <- paste("ovary", colnames(data_2), sep = "_")
colnames(data_3) <- paste("transverse.colon", colnames(data_3), sep = "_")
colnames(data_4) <- paste("bladder", colnames(data_4), sep = "_")
colnames(data_5) <- paste("omentum", colnames(data_5), sep = "_")
colnames(data_6) <- paste("rt.ex.ilian.lymphnode", colnames(data_6), sep = "_")
obj_1 = CreateSeuratObject(counts = data_1, min.cells = 3 ,min.features = 200)
obj_2 = CreateSeuratObject(counts = data_2, min.cells = 3 ,min.features = 200)
obj_3 = CreateSeuratObject(counts = data_3, min.cells = 3 ,min.features = 200)
obj_4 = CreateSeuratObject(counts = data_4, min.cells = 3 ,min.features = 200)
obj_5 = CreateSeuratObject(counts = data_5, min.cells = 3 ,min.features = 200)
obj_6 = CreateSeuratObject(counts = data_6, min.cells = 3 ,min.features = 200)

#Seurat object quality control (obj_1 ~ obj_6 input 바꿔가면서 ㄱ)
obj_6[["percent.mt"]] <- PercentageFeatureSet(obj_6, pattern = "^MT-")
plot1 <- FeatureScatter(obj_6, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj_6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
obj_6_qc <- subset(obj_6, subset = nFeature_RNA >200 & nFeature_RNA<7000 & percent.mt <5)
dim(obj_6)
dim(obj_6_qc)
obj_6 <- obj_6_qc
rm(obj_6_qc)

#공통된 genes로 합치기 -> obj_combined
union_gene_vector = intersect(rownames(obj_1), rownames(obj_2))
union_gene_vector = intersect(union_gene_vector,rownames(obj_3))
union_gene_vector = intersect(union_gene_vector,rownames(obj_4))
union_gene_vector = intersect(union_gene_vector,rownames(obj_5))
union_gene_vector = intersect(union_gene_vector,rownames(obj_6))
obj_1 <- obj_1[union_gene_vector,]
obj_2 <- obj_2[union_gene_vector,]
obj_3 <- obj_3[union_gene_vector,]
obj_4 <- obj_4[union_gene_vector,]
obj_5 <- obj_5[union_gene_vector,]
obj_6 <- obj_6[union_gene_vector,]
obj_1$batch <- "pbmc"
obj_2$batch <- "ovary"
obj_3$batch <- "transverse.colon"
obj_4$batch <- "bladder"
obj_5$batch <- "omentum"
obj_6$batch <- "rt.ex.ilian.lymphnode"
obj_combined <- merge(obj_1, obj_2)
obj_combined <- merge(obj_combined, obj_3)
obj_combined <- merge(obj_combined, obj_4)
obj_combined <- merge(obj_combined, obj_5)
obj_combined <- merge(obj_combined, obj_6)

#Integration by Seurat(v3) standard workflow 
obj.list <- SplitObject(obj_combined, split.by = "batch")
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE)
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], selection.method = "vst", 
                                        nfeatures = nrow(obj_combined), verbose = FALSE)
}
reference.list <- obj.list[c("pbmc", "ovary", "transverse.colon","bladder","omentum","rt.ex.ilian.lymphnode")]
obj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:40, anchor.features = nrow(obj_combined))
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:40)
DefaultAssay(obj.integrated) <- "integrated"
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE, features = rownames(obj_combined))
dim(obj.integrated@assays$integrated@scale.data)

obj.integrated <- RunPCA(obj.integrated, npcs = 100, verbose = FALSE)
ElbowPlot(obj.integrated, ndims = 100)
obj.integrated <- JackStraw(obj.integrated, dims = 50)
obj.integrated <- ScoreJackStraw(obj.integrated, dims = 1:40)
JackStrawPlot(obj.integrated, dims = 1:50)
obj.integrated <- obj.integrated %>% FindNeighbors(dims = 1:50)
obj.integrated <- obj.integrated %>% FindClusters(dims= 1:50, resolution = 1.5)
obj.integrated <- RunUMAP(obj.integrated, dims = 1:50)
DimPlot(obj.integrated, reduction = "umap", label = TRUE)
DimPlot(obj.integrated, reduction = "umap", group.by = "batch")
DimPlot(obj.integrated, reduction = "umap", split.by = "batch")

FeaturePlot(obj.integrated, c("CD3D", "CD3G"), max.cutoff = 3, min.cutoff = 0) #Tcell ,"CD3G"
FeaturePlot(obj.integrated, c("CD4","FOXP3","GATA3","TBX21"), min.cutoff = 0) #Tcell ,"CD3G"
FeaturePlot(obj.integrated, c("CD3D","CD8A","CD4","TOX"), min.cutoff = 0, max.cutoff = 5) #Tcell ,"CD3G"
FeaturePlot(obj.integrated, c("GNLY","NKG7","KLRC1"), min.cutoff = 0) #Tcell ,"CD3G"
FeaturePlot(obj.integrated, c("CD79A","MS4A1"), min.cutoff = 0, max.cutoff = 5) #Bcell
FeaturePlot(obj.integrated, c("CD11"), min.cutoff = 0) #Tcell ,"CD3G"
FeaturePlot(obj.integrated, c("AIF1","CD68","CD163"), min.cutoff = 0, max.cutoff = 5) #macrophage
FeaturePlot(obj.integrated, c("LILRA4","PTCRA","CD317","CD83"), min.cutoff = 0) #plasmacytoid dendritic
FeaturePlot(obj.integrated, c("IGHG1","IGHG3","IGLC2","IGKC"), min.cutoff = 0 ) #Bcell subtypes
FeaturePlot(obj.integrated, c("CD8A","CD8B", min.cutoff = 0))
FeaturePlot(obj.integrated, c("HESX1"), max.cutoff = 3,min.cutoff = 0)
FeaturePlot(obj.integrated, c("CTLA4","FOXP3"), max.cutoff = 5, min.cutoff = 0) #plasmacytoid dendritic
FeaturePlot(obj.integrated, c("CD11C", "CD317"), min.cutoff = 0, max.cutoff = 3)
FeaturePlot(obj.integrated, c("CD68", "CD14"), min.cutoff = 0, max.cutoff = 3) #Bcell

##1차 clustering
which(colnames(t(obj.integrated@assays$integrated@scale.data))=="CD3D")
rownames(obj.integrated@assays$integrated@scale.data)[14723]
which(colnames(t(obj.integrated@assays$integrated@scale.data))=="CD3G")
rownames(obj.integrated@assays$integrated@scale.data)[10071]
DimPlot(obj.integrated[,which(obj.integrated@assays$integrated@scale.data[14723,]>0 | obj.integrated@assays$integrated@scale.data[10071,]>0)], reduction = "umap", label = TRUE, group.by = "ident")

plot <- DimPlot(obj.integrated[,which(obj.integrated@assays$integrated@scale.data[14723,]>0 | obj.integrated@assays$integrated@scale.data[10071,]>0)], reduction = "umap", label = TRUE, group.by = "ident")
plot
select_cluster <- CellSelector(plot = plot)
Idents(obj.integrated, cells = select_cluster) <- "T cell"

which(colnames(t(obj.integrated@assays$integrated@scale.data))=="CD79A")
rownames(obj.integrated@assays$integrated@scale.data)[608]
DimPlot(obj.integrated[,which(obj.integrated@assays$integrated@scale.data[608,]>0)], reduction = "umap", label = TRUE, group.by = "ident")
plot<-DimPlot(obj.integrated[,which(obj.integrated@assays$integrated@scale.data[608,]>0)], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(obj.integrated, cells = select_cluster) <- "B cell"

DimPlot(obj.integrated[,which(!obj.integrated@active.ident %in% c("T cell", "B cell"))], reduction = "umap", label = TRUE, group.by = "ident")
obj.integrated <- RenameIdents(object = obj.integrated, "0" = "T cell", "1" = "T cell","2" = "T cell","3" = "T cell","4" = "Macrophage","5" = "B cell","6" = "T cell","7" = "T cell","8" = "B cell","9" = "B cell",
             "10" = "T cell","11" = "T cell","12" = "B cell","13" = "T cell","14" = "NK cell","15" = "T cell","16" = "Macrophage","17" = "T cell","18" = "B cell","19" = "Macrophage","20" = "T cell",
             "21" = "NK cell","22" = "NK cel","23" = "23","24" = "B cell","25" = "T cell","26" = "26","27" = "27","28" = "28","29" = "29","30" = "30","31" = "T cell")
DimPlot(obj.integrated, reduction = "umap", label = TRUE)

marker_26 <- FindMarkers(obj.integrated, ident.1 = "26", ident.2 = NULL, only.pos = TRUE)
marker_29 <- FindMarkers(obj.integrated, ident.1 = "29", ident.2 = NULL, only.pos = TRUE)


#2차 clustering
##T subset
t_integrated <- subset(obj.integrated, idents = "T cell")
DefaultAssay(t_integrated) <- "integrated"

t_integrated <- FindVariableFeatures(t_integrated)
t_integrated <- t_integrated %>% ScaleData(features = rownames(obj.integrated))
dim(t_integrated@assays$integrated@scale.data)

t_integrated <- t_integrated %>% RunPCA(npcs = 100)
ElbowPlot(t_integrated,ndims = 50)
t_integrated <- t_integrated %>%FindNeighbors(dims = 1:40)
t_integrated <- t_integrated %>% FindClusters(dims= 1:40, resolution = 1.5)
t_integrated <- t_integrated %>% RunUMAP(dims = 1:40)
DimPlot(t_integrated, reduction = "umap", label = TRUE, group.by = "ident")
#save.image("~/projects/00_cowork/02_krlim_scRNA/R.Data/scRNAseq_tcr_v2_re.RData")

DimPlot(t_integrated, reduction = "umap", split.by = "batch", label = FALSE)

FeaturePlot(t_integrated, c("CD8A"), min.cutoff = 0, max.cutoff = 2) #Tcell ,"CD3G"
FeaturePlot(t_integrated, c("CD4"), min.cutoff = 0, max.cutoff = .5)
FeaturePlot(t_integrated, c("CD4","FOXP3","GATA3","TBX21"), min.cutoff = 0, max.cutoff = 2) #Tcell ,"CD3G"
FeaturePlot(t_integrated, c("CD3D","CD8A","CD4","TOX"), min.cutoff = 0, max.cutoff = 5) #Tcell ,"CD3G"

###2.1 CD4+/ CD8+
which(colnames(t(t_integrated@assays$integrated@scale.data))=="CD4")
rownames(t_integrated@assays$integrated@scale.data)[1241]
DimPlot(t_integrated[,which(t_integrated@assays$integrated@scale.data[1241,]>0)], reduction = "umap", label = TRUE, group.by = "ident")
plot<-DimPlot(t_integrated[,which(t_integrated@assays$integrated@scale.data[1241,]>0)], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(t_integrated, cells = select_cluster) <- "CD4 T cells"
DimPlot(t_integrated, reduction = "umap", label = TRUE)

which(colnames(t(t_integrated@assays$integrated@scale.data))=="CD8A")
rownames(t_integrated@assays$integrated@scale.data)[1067]
plot<-DimPlot(t_integrated[,which(t_integrated@assays$integrated@scale.data[1067,]>0 & t_integrated@active.ident != "CD4 T cells")], reduction = "umap", label = TRUE, group.by = "ident")
plot
select_cluster <- CellSelector(plot = plot)
Idents(t_integrated, cells = select_cluster) <- "CD8 T cells"
DimPlot(t_integrated, reduction = "umap", label = TRUE)

DimPlot(t_integrated[,which(!t_integrated@active.ident %in% c("CD4 T cells", "CD8 T cells"))], reduction = "umap", label = TRUE, group.by = "ident")
t_integrated <- RenameIdents(object = t_integrated, "0" = "0", "1" = "CD4 T cells","2" = "CD8 T cells","3" = "CD4 T cells","4" = "4","5" = "CD8 T cells","6" = "6","7" = "7","8" = "CD4 T cells","9" = "9",
                               "10" = "CD4 T cells","11" = "CD4 T cells","12" = "CD8 T cells","13" = "CD8 T cells","14" = "CD8 T cells","15" = "CD8 T cells","16" = "CD8 T cells","17" = "17","18" = "CD8 T cells","19" = "CD4 T cells","20" = "20",
                               "21" = "CD4 T cells")
t_integrated <- RenameIdents(object = t_integrated, "0" = "CD4 T cells", "7" = "CD4 T cells","17" = "CD8 T cells","20" = "CD4 T cells")
DimPlot(t_integrated, reduction = "umap", label = TRUE)

plot<-DimPlot(t_integrated[,which(t_integrated@assays$integrated@scale.data[1067,]>0 & t_integrated@active.ident %in% c("4","6","9"))], reduction = "umap", label = TRUE, group.by = "ident")
plot<- FeaturePlot(t_integrated[,t_integrated@active.ident %in% c("4","6","9")], "CD8A", min.cutoff = 0, max.cutoff = 2)
plot
select_cluster <- CellSelector(plot = plot)
Idents(t_integrated, cells = select_cluster) <- "CD8 T cells"
DimPlot(t_integrated, reduction = "umap", label = TRUE)

###2.2 CD4+
cd4_integrated <- subset(t_integrated, idents = "CD4 T cells")
DefaultAssay(cd4_integrated) <- "integrated"
cd4_integrated <- FindVariableFeatures(cd4_integrated)
cd4_integrated <- cd4_integrated %>% ScaleData(features = rownames(obj.integrated))
dim(cd4_integrated@assays$integrated@scale.data)
cd4_integrated <- cd4_integrated %>% RunPCA(npcs = 100)
ElbowPlot(cd4_integrated,ndims = 50)
cd4_integrated <- cd4_integrated %>%FindNeighbors(dims = 1:40)
cd4_integrated <- cd4_integrated %>% FindClusters(dims= 1:40, resolution = 1.5)
cd4_integrated <- cd4_integrated %>% RunUMAP(dims = 1:40)
DimPlot(cd4_integrated, reduction = "umap", label = TRUE, group.by = "ident")

FeaturePlot(cd4_integrated, c("CCL20", "IFNG", "IL2RA", "LTA", "CSF2"), min.cutoff = 0, max.cutoff = 5) #Tcell ,"CD3G"
FeaturePlot(cd4_integrated, c("GNLY", "CCR7"), min.cutoff = 0, max.cutoff = 2) #Tcell ,"CD3G"
FeaturePlot(cd4_integrated, c("IL21", "LRMP", "PDCD1", "CXCR5"), min.cutoff = 0, max.cutoff = 1) #Tcell ,Tfh
FeaturePlot(cd4_integrated, c("PDCD1", "HAVCR2", "LAG3", "TOX"), min.cutoff = 0, max.cutoff = 3) #Texhausted
FeaturePlot(cd4_integrated, c("FOXP3", "IL2RA"), min.cutoff = 0, max.cutoff = 3) #Treg

which(colnames(t(cd4_integrated@assays$integrated@scale.data))=="FOXP3")
rownames(cd4_integrated@assays$integrated@scale.data)[383]
which(colnames(t(cd4_integrated@assays$integrated@scale.data))=="IL2RA")
rownames(cd4_integrated@assays$integrated@scale.data)[545]
DimPlot(cd4_integrated[,which(cd4_integrated@assays$integrated@scale.data[383,]>0|cd4_integrated@assays$integrated@scale.data[545,]>0)], reduction = "umap", label = TRUE, group.by = "ident")
plot<-DimPlot(cd4_integrated[,which(cd4_integrated@assays$integrated@scale.data[383,]>0|cd4_integrated@assays$integrated@scale.data[545,]>0)], reduction = "umap", label = TRUE, group.by = "ident")
plot
select_cluster <- CellSelector(plot = plot)
Idents(cd4_integrated, cells = select_cluster) <- "Treg cells"
DimPlot(cd4_integrated, reduction = "umap", label = TRUE)

which(colnames(t(cd4_integrated@assays$integrated@scale.data))=="PDCD1")
rownames(cd4_integrated@assays$integrated@scale.data)[1293]
which(colnames(t(cd4_integrated@assays$integrated@scale.data))=="TOX")
rownames(cd4_integrated@assays$integrated@scale.data)[2186]
DimPlot(cd4_integrated[,which((cd4_integrated@assays$integrated@scale.data[1293,]>1|cd4_integrated@assays$integrated@scale.data[2186,]>1)& cd4_integrated@active.ident != "Treg cells" )], reduction = "umap", label = TRUE, group.by = "ident")
plot<-DimPlot(cd4_integrated[,which((cd4_integrated@assays$integrated@scale.data[1293,]>1|cd4_integrated@assays$integrated@scale.data[2186,]>1)& cd4_integrated@active.ident != "Treg cells" )], reduction = "umap", label = TRUE, group.by = "ident")

plot<-DimPlot(cd4_integrated[,cd4_integrated@active.ident == "CD4 T cells"])
plot
select_cluster <- CellSelector(plot = plot)
Idents(cd4_integrated, cells = select_cluster) <- "Treg cells"
DimPlot(cd4_integrated, reduction = "umap", label = TRUE, split.by = "batch")



