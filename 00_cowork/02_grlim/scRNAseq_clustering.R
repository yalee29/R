library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
library(gplots)
library(cowplot)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

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

#각각 clustering
#obj_1: PBMC annotation ----------------------------------------------------------------------------
obj_1 <- obj_1 %>% NormalizeData()
obj_1 <- obj_1 %>% FindVariableFeatures()
obj_1 <- obj_1 %>% ScaleData()
dim(obj_1@assays$RNA@scale.data)
obj_1 <- obj_1 %>% RunPCA(npcs = 50)
ElbowPlot(obj_1,ndims = 50 )
obj_1 <- JackStraw(obj_1, dims = 50)
obj_1 <- ScoreJackStraw(obj_1, dims = 1:50)
JackStrawPlot(obj_1, dims = 1:50)
obj_1 <- obj_1 %>%FindNeighbors(dims = 1:35)
obj_1 <- obj_1 %>% FindClusters(dims= 1:35, resolution = 1)
obj_1 <- obj_1 %>% RunUMAP(dims = 1:35)
DimPlot(obj_1, label = TRUE)

##Bcell(plasma/naive/activated/naive): batch 마다 모두 똑같이 해볼것. 
FeaturePlot(object = obj_1, features = c("MZB1", "FCER2","CD27","TCL1A","MS4A1"), min.cutoff = 0, max.cutoff = 5)
#Monocyte
FeaturePlot(object = obj_1, features = c("CD14","FCGR3A", "MS4A7")) #monocyte(CD16+/CD14+) 
#M
FeaturePlot(object = obj_1, features = c("MARCO","CXCL3")) #M0
FeaturePlot(object = obj_1, features = c("CD163","MS4A6A")) #M2
FeaturePlot(object = obj_1, features = c("CCR7","BCL2A1","TNF","CD86")) #M1 (BCL2A1:M0, M1 모두발현)
#M_PTGDS
FeaturePlot(object = obj_1, features = c("PTGDS")) #M0
#NKcell
FeaturePlot(object = obj_1, features = c("NKG7","GNLY"))
#DC
FeaturePlot(object = obj_1, features = c("CD1C","FCER1A"))
#T cells
FeaturePlot(object = obj_1, features = c("CD3D","CD3G","CD4","CD8A"))
FeaturePlot(object = obj_1, features = c("CD4","CD8A"))
FeaturePlot(object = obj_1, features = c("CCR7","CTLA4","PDCD1","CXCR6"))
FeaturePlot(object = obj_1, features = c("CCR7","IL7R","S100A4"))
#Unknown
FeaturePlot(object = obj_1, features = c("KRT7"))

DimPlot(obj_3, label = TRUE)
FeaturePlot(obj_6, c("AIF1", "PTGDS","CD1C"), ncol = 3, max.cutoff = 2)
FeaturePlot(obj_6, c("MZB1","MS4A1","NKG7"), ncol = 3)
FeaturePlot(obj_6, c("CD3D", "CD8A","CD4"), ncol = 3)


obj_1 <- RenameIdents(object = obj_1, "0"="CD4+ T","1"="CD4+ T","2"="CD14+ Monocyte","3"="CD8+ T","4"="NK","5"="CD8+ T","6"="NK","7"="B","8"="CD16+ Monocyte","9"="CD8+ T","10"="CD8+ T","11"="B","12"="CD4+ T","13"="CD8+ T","14"="CD4+ T","15"="DC","16"="CD8+ T","17"="NK","18"="CD4+ T")

pbmc.markers <- FindAllMarkers(obj_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc_ann.markers <- FindAllMarkers(obj_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10_ann <- pbmc_ann.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(obj_1, features = top10_ann$gene) + NoLegend()

#obj_2: ovary annotation ----------------------------------------------------------------------
obj_2 <- obj_2 %>% NormalizeData()
obj_2 <- obj_2 %>% FindVariableFeatures()
obj_2 <- obj_2 %>% ScaleData()
dim(obj_2@assays$RNA@scale.data)
obj_2 <- obj_2 %>% RunPCA(npcs = 50)
ElbowPlot(obj_2,ndims = 50 )
obj_2 <- JackStraw(obj_2, dims = 50)
obj_2 <- ScoreJackStraw(obj_2, dims = 1:50)
JackStrawPlot(obj_2, dims = 1:50)
obj_2 <- obj_2 %>%FindNeighbors(dims = 1:50)
obj_2 <- obj_2 %>% FindClusters(dims= 1:50, resolution = 1)
obj_2 <- obj_2 %>% RunUMAP(dims = 1:50)
DimPlot(obj_2, label = TRUE)

which(colnames(t(obj_2@assays$RNA@scale.data))=="CD8A")
rownames(obj_2@assays$RNA@scale.data)[260]
plot<-DimPlot(obj_2[,which(obj_2@assays$RNA@scale.data[260,]>0 & obj_2@active.ident=="4")], reduction = "umap", label = TRUE, group.by = "ident")
plot<-DimPlot(obj_2[,which(obj_2@active.ident=="CD8+ T"| obj_2@active.ident=="4")], reduction = "umap", label = TRUE, group.by = "ident")

FeaturePlot(obj_1, c("PTGDS")) 
plot
select_cluster <- CellSelector(plot = plot)
Idents(obj_2, cells = select_cluster) <- "CD8+ T"

obj_2 <- RenameIdents(object = obj_2, "0"="M","1"="CD8+ T","2"="CD4+ T","3"="B","4"="CD4+ T","5"="M","6"="CD4+ T","7"="NK","8"="CD8+ T","9"="Plasma B","10"="DC","11"="CD4+ T","12"="Unknown","13"="M","14"="M","15"="CD8+ T","16"="CD4+ T","17"="CD8+ T","18"="Plasma B","19"="M","20"="M")
ovary.markers <- FindAllMarkers(obj_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ovary_top10 <- ovary.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
ovary_ann.markers <- FindAllMarkers(obj_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ovary_top10_ann <- ovary_ann.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(obj_2, features = ovary_top10_ann$gene) + NoLegend()

#obj_3:t_colon-----------------------------------------------------------------------
obj_3 <- obj_3 %>% NormalizeData()
obj_3 <- obj_3 %>% FindVariableFeatures()
obj_3 <- obj_3 %>% ScaleData()
dim(obj_3@assays$RNA@scale.data)
obj_3 <- obj_3 %>% RunPCA(npcs = 50)
ElbowPlot(obj_3,ndims = 50 )
obj_3 <- JackStraw(obj_3, dims = 50)
obj_3 <- ScoreJackStraw(obj_3, dims = 1:50)
JackStrawPlot(obj_3, dims = 1:50)
obj_3 <- obj_3 %>%FindNeighbors(dims = 1:45)
obj_3 <- obj_3 %>% FindClusters(dims= 1:45, resolution = 1)
obj_3 <- obj_3 %>% RunUMAP(dims = 1:45)
DimPlot(obj_3, label = TRUE)

FeaturePlot(obj_3, c("FCER1A"), max.cutoff = 1, pt.size = 0.3) 
plot<-FeaturePlot(obj_3, c("FCER1A"), max.cutoff = 1, pt.size = 0.3) 
select_cluster <- CellSelector(plot = plot)
Idents(obj_3, cells = select_cluster) <- "DC"

obj_3 <- RenameIdents(object = obj_3, "8"="CD4+ T","14"="CD4+ T","2"="CD4+ T","6"="CD4+ T","0"="CD8+ T","5"="CD8+ T","17"="CD8+ T","18"="CD8+ T","16"="NK","1"="B","4"="B","7"="B","9"="B","10"="B","13"="Plasma B","21"="Unknown")
obj_3 <- RenameIdents(object = obj_3, "19"="CD8+ T","15"="B","12"="M", "3"="CD4+ T", "11"= "CD8+ T", "20" = "M_PTGDS")

colon.markers <- FindAllMarkers(obj_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
colon_top10 <- colon.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
colon_ann.markers <- FindAllMarkers(obj_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
colon_top10_ann <- colon_ann.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(obj_3, features = colon_top10_ann$gene) + NoLegend()

#obj_4:bladder-----------------------------------------------------------------------
obj_4 <- obj_4 %>% NormalizeData()
obj_4 <- obj_4 %>% FindVariableFeatures()
obj_4 <- obj_4 %>% ScaleData()
dim(obj_4@assays$RNA@scale.data)
obj_4 <- obj_4 %>% RunPCA(npcs = 50)
ElbowPlot(obj_4,ndims = 50 )
obj_4 <- JackStraw(obj_4, dims = 50)
obj_4 <- ScoreJackStraw(obj_4, dims = 1:50)
JackStrawPlot(obj_4, dims = 1:50)
obj_4 <- obj_4 %>%FindNeighbors(dims = 1:45)
obj_4 <- obj_4 %>% FindClusters(dims= 1:45, resolution = 1)
obj_4 <- obj_4 %>% RunUMAP(dims = 1:45)
DimPlot(obj_4, label = TRUE)

FeaturePlot(obj_4, c("CD1C"), max.cutoff = 1, pt.size = 0.3) 
plot<-FeaturePlot(obj_4, c("CD1C"), max.cutoff = 1, pt.size = 0.3) 
select_cluster <- CellSelector(plot = plot)
Idents(obj_4, cells = select_cluster) <- "DC"

obj_4 <- RenameIdents(object = obj_4, "2"="CD4+ T","4"="CD4+ T","10"="CD4+ T","11"="CD4+ T","0"="CD8+ T","3"="CD8+ T","5"="CD8+ T","7"="CD8+ T","13"="CD8+ T","14"="CD8+ T","15"="CD8+ T","6"="NK","1"="B","12"="B","8"="M","16"="M","9"="Plasma B")
DimPlot(obj_4, label = TRUE)

bladder_ann.markers <- FindAllMarkers(obj_4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
bladder_top10_ann <- bladder_ann.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(obj_4, features = bladder_top10_ann$gene) + NoLegend()

#obj_5:omentum-----------------------------------------------------------------------
obj_5 <- obj_5 %>% NormalizeData()
obj_5 <- obj_5 %>% FindVariableFeatures()
obj_5 <- obj_5 %>% ScaleData()
dim(obj_5@assays$RNA@scale.data)
obj_5 <- obj_5 %>% RunPCA(npcs = 50)
ElbowPlot(obj_5,ndims = 50 )
obj_5 <- JackStraw(obj_5, dims = 50)
obj_5 <- ScoreJackStraw(obj_5, dims = 1:50)
JackStrawPlot(obj_5, dims = 1:50)
obj_5 <- obj_5 %>%FindNeighbors(dims = 1:45)
obj_5 <- obj_5 %>% FindClusters(dims= 1:45, resolution = 1)
obj_5 <- obj_5 %>% RunUMAP(dims = 1:45)
DimPlot(obj_5, label = TRUE)

DimPlot(obj_5, label = TRUE)
FeaturePlot(obj_5[,obj_5@active.ident=="13"], c("CD3D"), max.cutoff = 1, pt.size = 0.3) 
plot<-FeaturePlot(obj_5[,obj_5@active.ident=="1"], c("CD4"), max.cutoff = 1, pt.size = 0.3) 
select_cluster <- CellSelector(plot = plot)
Idents(obj_5, cells = select_cluster) <- "CD4+ T"

obj_5 <- RenameIdents(object = obj_5, "10"="B","11"="B","14"="B","5"="CD4+ T","9"="CD4+ T","6"="CD4+ T","7"="CD4+ T","2"="CD8+ T","3"="CD8+ T","4"="CD8+ T","15"="CD8+ T","12"="NK","18"="NK","17"="Plasma B","19"="M_PTGDS","22"="M","21"="M","8"="M","16"="M","20"="Unknown")
obj_5 <- RenameIdents(object = obj_5, "1"="CD8+ T", "13"="CD8+ T")
DimPlot(obj_5, label = TRUE)

omentum_ann.markers <- FindAllMarkers(obj_5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
omentum_top10_ann <- omentum_ann.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(obj_5, features = omentum_top10_ann$gene) + NoLegend()

#obj_6:rt.ex.ilian.lymphnode-----------------------------------------------------------------------
obj_6 <- obj_6 %>% NormalizeData()
obj_6 <- obj_6 %>% FindVariableFeatures()
obj_6 <- obj_6 %>% ScaleData()
dim(obj_6@assays$RNA@scale.data)
obj_6 <- obj_6 %>% RunPCA(npcs = 50)
ElbowPlot(obj_6,ndims = 50 )
obj_6 <- JackStraw(obj_6, dims = 50)
obj_6 <- ScoreJackStraw(obj_6, dims = 1:50)
JackStrawPlot(obj_6, dims = 1:50)
obj_6 <- obj_6 %>%FindNeighbors(dims = 1:45)
obj_6 <- obj_6 %>% FindClusters(dims= 1:45, resolution = 1)
obj_6 <- obj_6 %>% RunUMAP(dims = 1:45)
DimPlot(obj_6, label = TRUE)

FeaturePlot(obj_6, c(""), max.cutoff = 1, pt.size = 0.3) 
plot<-FeaturePlot(obj_6, c("KRT7"), max.cutoff = 1, pt.size = 0.3) 
select_cluster <- CellSelector(plot = plot)
Idents(obj_6, cells = select_cluster) <- "Unknown"

obj_6 <- RenameIdents(object = obj_6, "20"="PTGDS","21"="Plasma B","17"="Plasma B","13"="B","15"="B","12"="M","11"="NK","16"="CD8+ T","5"="CD8+ T","10"="CD8+ T","4"="CD4+ T","22"="CD4+ T","0"="CD4+ T","3"="CD4+ T","8"="CD4+ T","7"="CD4+ T","2"="B","9"="B","1"="B","14"="B","6"="B","18"="B","19"="B")
obj_6 <- RenameIdents(object = obj_6, "PTGDS"="M_PTGDS")
DimPlot(obj_6, label = TRUE)

omentum_ann.markers <- FindAllMarkers(obj_6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
omentum_top10_ann <- omentum_ann.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(obj_6, features = omentum_top10_ann$gene) + NoLegend()

#Integrating: Seurat v3(CCA)===================================================================
obj_combined <- merge(obj_1, merge(obj_2, merge(obj_3, merge(obj_4, merge(obj_5, obj_6)))))

obj.list <- SplitObject(obj_combined, split.by = "batch")
obj.list <- obj.list[c("pbmc", "ovary", "transverse.colon", "bladder", "omentum", "rt.ex.ilian.lymphnode")]
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE)
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)

library(ggplot2)
library(cowplot)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(obj.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:30)
DimPlot(obj.integrated, reduction = "umap", group.by = "batch", label = TRUE)
DimPlot(obj.integrated, reduction = "umap", label = TRUE)

#p2 <- DimPlot(obj.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
#              repel = TRUE) + NoLegend()
#plot_grid(p1, p2)

