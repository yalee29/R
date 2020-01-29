#read count data to SeuratObject
BN1_P <- read_csv("./projects/00_cowork/02_krlim_scRNA/00_data/gse118828/BN1_P.csv") #223 cells(26364 genes)
HG1_P <- read_csv("./projects/00_cowork/02_krlim_scRNA/00_data/gse118828/HG1_P.csv") #252
HG1_M <- read_csv("./projects/00_cowork/02_krlim_scRNA/00_data/gse118828/HG1_M.csv") #330
HG2F_P <- read_csv("./projects/00_cowork/02_krlim_scRNA/00_data/gse118828/HG2F_P.csv") #262
HG2F_M <- read_csv("./projects/00_cowork/02_krlim_scRNA/00_data/gse118828/HG2F_M.csv") #267
BN1_P <-BN1_P %>% column_to_rownames(var = "CellId") %>% t()
HG1_P <-HG1_P %>% column_to_rownames(var = "CellId") %>% t()
HG1_M <-HG1_M %>% column_to_rownames(var = "CellId") %>% t()
HG2F_P <-HG2F_P %>% column_to_rownames(var = "CellId") %>% t()
HG2F_M <-HG2F_M %>% column_to_rownames(var = "CellId") %>% t()

BN1_P <- CreateSeuratObject(BN1_P)
HG1_P <- CreateSeuratObject(HG1_P)
HG1_M <- CreateSeuratObject(HG1_M)
HG2F_P <- CreateSeuratObject(HG2F_P)
HG2F_M <- CreateSeuratObject(HG2F_M)

#
HG4_P1 <- read_csv("./projects/00_cowork/02_krlim_scRNA/00_data/gse118828/HG4_P1.csv")
HG4_P2 <- read_csv("./projects/00_cowork/02_krlim_scRNA/00_data/gse118828/HG4_P2.csv")
HG4_P3 <- read_csv("./projects/00_cowork/02_krlim_scRNA/00_data/gse118828/HG4_P3.csv")
HG4_P4 <- read_csv("./projects/00_cowork/02_krlim_scRNA/00_data/gse118828/HG4_P4.csv")
HG4_P1 <-HG4_P1 %>% column_to_rownames(var = "CellId") %>% t()
HG4_P2 <-HG4_P2 %>% column_to_rownames(var = "CellId") %>% t()
HG4_P3 <-HG4_P3 %>% column_to_rownames(var = "CellId") %>% t()
HG4_P4 <-HG4_P4 %>% column_to_rownames(var = "CellId") %>% t()

HG4_P1 <- HG4_P1 %>% CreateSeuratObject()
HG4_P2 <- HG4_P2 %>% CreateSeuratObject()
HG4_P3 <- HG4_P3%>% CreateSeuratObject()
HG4_P4 <- HG4_P4 %>% CreateSeuratObject()
HG4_P1$batch2 <- "HG4_P1"
HG4_P2$batch2 <- "HG4_P2"
HG4_P3$batch2 <- "HG4_P3"
HG4_P4$batch2 <- "HG4_P4"
HG4 <- merge(HG4_P1, c(HG4_P2,HG4_P3,HG4_P4))

#Integrating HG1_P and HG1_M
BN1_P$batch <- "BN1_P"
BN1_P<-NormalizeData(BN1_P)
BN1_P<-FindVariableFeatures(BN1_P)

HG1_P$batch <- "HG1_P"
HG1_P<-NormalizeData(HG1_P)
HG1_P <- FindVariableFeatures(HG1_P)

#Seurat-Integrating
immune.anchors <- FindIntegrationAnchors(object.list = list(BN1_P, HG1_P), dims=1:30)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims=1:30)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, features = rownames(immune.combined))
immune.combined <- RunPCA(immune.combined, npcs=100)

immune.combined <- JackStraw(immune.combined, dims = 50)
immune.combined <- ScoreJackStraw(immune.combined, dims = 1:50)
JackStrawPlot(immune.combined, dims = 1:50)
ElbowPlot(immune.combined, ndims = 20)

immune.combined <- RunUMAP(immune.combined, reduction="pca", dims=1:10)
immune.combined <- FindNeighbors(immune.combined, reduction="pca", dims=1:10)
immune.combined <- FindClusters(immune.combined, resolution=0.7)
p1<- DimPlot(immune.combined, reduction="umap", group.by="batch")
p2<- DimPlot(immune.combined, reduction="umap", group.by="seurat_clusters",label=TRUE)
plot_grid(p1, p2)
DimPlot(immune.combined, reduction = "umap", split.by = "batch")
DimPlot(immune.combined, reduction = "umap",group.by="seurat_clusters",label=TRUE )

#integrated analysis
DefaultAssay(immune.combined) <- "RNA"
FeaturePlot(immune.combined, features = c("ANGPT2","FGF18", "ITGB4","ITGB8"))
FeaturePlot(immune.combined, features = c("CD3D","CD3G","CD3E"))

markers.0 <- FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = "batch", ,verbose = FALSE)
markers.1 <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "batch", verbose = FALSE)
markers.2 <- FindConservedMarkers(immune.combined, ident.1 = 2, grouping.var = "batch", verbose = FALSE)
markers.3 <- FindConservedMarkers(immune.combined, ident.1 = 3, grouping.var = "batch", verbose = FALSE)
markers.4 <- FindConservedMarkers(immune.combined, ident.1 = 4, grouping.var = "batch", verbose = FALSE)
markers <- rbind(markers.0[1:100,],markers.1[1:100,],markers.2[1:100,],markers.3[1:100,],markers.4[1:100,])

DoHeatmap(immune.combined, features = rownames(markers))


markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(immune.combined, features = top100$gene)

FeaturePlot(immune.combined, features = c("MALAT1","YPEL5","MDK","HLA-DRA","SPARC"))



MALAT1




