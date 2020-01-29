library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
library(gplots)
library(cowplot)

#Setup control & met
data_dir <- './projects/00_cowork/jhpark_scRNA/met/outs/filtered_feature_bc_matrix'
list.files(data_dir)
data <- Read10X(data.dir = data_dir)

control = CreateSeuratObject(counts = data, min.cells = 3 ,min.features = 200)
control$batch <- "control"
control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^mt-")
plot1 <- FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
VlnPlot(control, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
control <- subset(control, subset = nFeature_RNA >200 & nFeature_RNA<7500 & percent.mt <40)
control<-NormalizeData(control)
control<-FindVariableFeatures(control)

met = CreateSeuratObject(counts = data, min.cells = 3 ,min.features = 200)
met$batch <- "metformin"
met[["percent.mt"]] <- PercentageFeatureSet(met, pattern = "^mt-")
plot1 <- FeatureScatter(met, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(met, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
VlnPlot(met, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
met <- subset(met, subset = nFeature_RNA >200 & nFeature_RNA<6000 & percent.mt <50)
met <- NormalizeData(met)
met <- FindVariableFeatures(met)

#Seurat-Integrating
immune.anchors <- FindIntegrationAnchors(object.list = list(control, met), dims=1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims=1:20)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined)
immune.combined <- RunPCA(immune.combined, npcs=100)
ElbowPlot(immune.combined, ndims = 100)
immune.combined <- RunUMAP(immune.combined, reduction="pca", dims=1:50)
immune.combined <- FindNeighbors(immune.combined, reduction="pca", dims=1:50)
immune.combined <- FindClusters(immune.combined, resolution=.5)
p1<- DimPlot(immune.combined, reduction="umap", group.by="batch")
p2<- DimPlot(immune.combined, reduction="umap", group.by="seurat_clusters",label=TRUE)
plot_grid(p1, p2)
p2

#integrated analysis
DefaultAssay(immune.combined) <- "RNA"
markers.0 <- FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = "stim", ,verbose = FALSE)
head(nk.markers)


FeaturePlot(immune.combined, features = c("Trdc"))
FeaturePlot(immune.combined, features = c("Cd3d","Cd3g","Cd3e"))

#gd_t DEG
gd_t <- subset(immune.combined, idents = c("18"))
gd_t <- ScaleData(gd_t)
Idents(gd_t) <- "batch"
avg.gd_t <- log1p(AverageExpression(gd_t, verbose=FALSE)$RNA)
avg.gd_t$gene <- rownames(avg.gd_t)
gd_t.response <- FindMarkers(gd_t, ident.1="control", ident.2="metformin")
gd_t.response <- gd_t.response %>% rownames_to_column()
gd_deg <- gd_t.response %>% filter(gd_t.response$p_val < 0.05)
head(gd_deg)
gd_deg_ordr <- gd_deg [ order(gd_deg$avg_logFC), ] 
DoHeatmap(gd_t, features = gd_deg_ordr$rowname) 

#heatmap version 2
#library(NMF) #heatmap version2
#tmp_mt <- as.matrix(gd_t@assays$RNA@data) 
#tmp_mt <- tmp_mt %>% as.data.frame %>% rownames_to_column() 
#tmp_mt_merge<- merge(tmp_mt, gd_deg_ordr, by.x = "rowname", by.y = "rowname")
#tmp_mt_merge <- tmp_mt_merge [ order(tmp_mt_merge$avg_logFC), ] 
#tmp_mt_merge<-tmp_mt_merge[,-(113:117)]
#rownames(tmp_mt_merge)<-tmp_mt_merge$rowname
#tmp_mt_merge<-tmp_mt_merge[,-1]
#head(rownames(tmp_mt_merge)) 
#gd_annotation<- gd_t@meta.data %>% rownames_to_column() %>% select(rowname, batch) %>% arrange(batch) %>% column_to_rownames('rowname')
#gd_annotation
#tmp_mt_merge <- tmp_mt_merge[,rownames(gd_annotation)]
#aheatmap(tmp_mt_merge, scale = "row", annCol = gd_annotation, Rowv = NA, Colv = NA, height = 30, cexRow = 10)



#T cells
t.combined <- immune.combined[,(immune.combined@active.ident %in% c(1,4,5,8,7,12,14))]
DimPlot(t.combined, reduction = "umap", group.by = "batch") %>% + ggtitle("T cell")
DimPlot(t.combined, reduction = "umap", group.by = "ident", label = TRUE) %>% + ggtitle("T cell")

t.combined <- ScaleData(t.combined)
t.combined <- FindVariableFeatures(t.combined)
t.combined <- RunPCA(t.combined, npcs = 50)
# t.combined <- JackStraw(t.combined, dims = 30)
# t.combined <- ScoreJackStraw(t.combined, dims = 30)
# JackStrawPlot(t.combined, dims = 1:30) 
ElbowPlot(t.combined)

t.combined <- RunUMAP(t.combined, reduction = "pca", dims = 1:20)
t.combined <- FindNeighbors(t.combined, reduction = "pca", dims = 1:20)
t.combined <- FindClusters(t.combined, resolution = 0.5)

p_t1<-DimPlot(t.combined, group.by = "batch", reduction = "umap")
p_t2<-DimPlot(t.combined, group.by = "ident", reduction = "umap", label = TRUE)
plot_grid(p_t1,p_t2)
DimPlot(t.combined, reduction = "umap", split.by = "batch")
DimPlot(t.combined, reduction = "tsne", group.by = "batch")
DimPlot(t.combined, group.by = "ident", reduction = "umap", label = TRUE)

FeaturePlot(t.combined, features = c("Trdc"))




#other_than_25 <- immune.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% filter(integrated_snn_res.2 != 25) %>% .$cell_id

trdc_dt <- immune.combined@assays$RNA@data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column('gene_name') %>% filter(gene_name == 'Trdc') %>% column_to_rownames('gene_name') %>% t() %>% as.data.frame()
ggplot(trdc_dt)+
  geom_histogram(aes(Trdc), binwidth=0.1)+
  coord_cartesian(ylim=c(0,100))

target_cells <- trdc_dt %>% rownames_to_column('cell_id') %>% filter(Trdc > 1) %>% .$cell_id
t.combined@meta.data <- t.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% mutate(trdc_info= ifelse(cell_id %in% target_cells & RNA_snn_res.0.5 == 9, 'pos9',
                                                                                                                                ifelse(cell_id %in% target_cells , 'pos','neg'))) %>% column_to_rownames('cell_id')
#trdc+ subset
trdc.combined <- t.combined[,(t.combined@meta.data$trdc_info %in% c('pos','pos9'))]
#DimPlot(trdc.combined, reduction = "umap", group.by = "batch") %>% + ggtitle("gamma delta T cell")
#DimPlot(trdc.combined, reduction = "umap", group.by = "ident", label = TRUE) %>% + ggtitle("gamma deltaT cell")

trdc.combined <- ScaleData(trdc.combined)
trdc.combined <- FindVariableFeatures(trdc.combined)
trdc.combined <- RunPCA(trdc.combined, npcs = 100)
ElbowPlot(trdc.combined, ndims = 100)
trdc.combined <- RunUMAP(trdc.combined, reduction = "pca", dims = 1:20)
trdc.combined <- FindNeighbors(trdc.combined, reduction = "pca", dims = 1:20)
trdc.combined <- FindClusters(trdc.combined, resolution = 0.8)

p_t1<-DimPlot(trdc.combined, group.by = "batch", reduction = "umap")
p_t2<-DimPlot(trdc.combined, group.by = "ident", reduction = "umap", label = TRUE)
plot_grid(p_t1,p_t2)
DimPlot(trdc.combined, reduction = "umap", split.by = "batch")
DimPlot(trdc.combined, reduction = "umap", group.by = "trdc_info")
DimPlot(trdc.combined, group.by = "ident", reduction = "umap", label = TRUE)

FeaturePlot(trdc.combined, features = c("Trdc"))
FeaturePlot(trdc.combined, features = c("Ifng", "Cd27", "Il17a","Pdcd1")) #,"Klra1"))

#gd_t DEG
gd_t_ifn <- subset(trdc.combined, idents = c("0","1"))
Idents(gd_t_ifn) <- "batch"
avg.gd_t_ifn <- log1p(AverageExpression(gd_t_ifn, verbose = FALSE)$RNA)
avg.gd_t_ifn$gene <- rownames(avg.gd_t_ifn)

gd_t_ifn.response <- FindMarkers(gd_t_ifn, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
gd_t_ifn.response <- gd_t_ifn.response %>% rownames_to_column()
gd_t_ifn_deg <- gd_t_ifn.response %>% filter(gd_t_ifn.response$p_val < 0.05)
head(gd_t_ifn.response)

p_gd_t_ifn <-ggplot(avg.gd_t_ifn, aes(control, metformin)) + geom_point() + ggtitle("gdT cell")+ coord_cartesian(xlim = c(0:6), ylim = c(0:6))
p_gd_t_ifn

gd_t_ifn_deg_ordr <- gd_t_ifn_deg [ order(gd_t_ifn_deg$avg_logFC), ] 
DoHeatmap(gd_t_ifn, features = gd_t_ifn_deg_ordr$rowname) 

#특정 유전자에 대해 그리기
tmp_dt1 <-  gd_t@assays$RNA@data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column('gene_name') %>% filter(gene_name == 'Bcl2') %>% t() %>% as.data.frame() %>% rownames_to_column("cell_id")
tmp_dt1 <- tmp_dt1[-1,]
colnames(tmp_dt1) <- c("cell_id", "Bcl2")
tmp_dt2 <- gd_t@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% select(cell_id, batch)
tmp_dt <- left_join(tmp_dt1, tmp_dt2)
tmp_dt$Bcl2 <- tmp_dt$Bcl2 %>% as.numeric()
head(tmp_dt)

library(ggsignif)
ggplot(tmp_dt, aes(x=batch, y=Bcl2))+
  geom_boxplot(fatten = NULL) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  stat_summary(fun.y=mean) +
  geom_signif(xmin=start, xmax=end, y_position = 12, annotations = '*', tip_length = c(0,0))








