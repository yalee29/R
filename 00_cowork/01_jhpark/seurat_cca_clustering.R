library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
library(gplots)
library(cowplot)

#Seurat-Integrating
#Setup control & met
data_dir <- './projects/00_cowork/jhpark_scRNA/met/outs/filtered_feature_bc_matrix'
list.files(data_dir)
data <- Read10X(data.dir = data_dir)

control = CreateSeuratObject(counts = data)
control$batch <- "control"
control<-NormalizeData(control)
control<-FindVariableFeatures(control)

met = CreateSeuratObject(counts = data)
met$batch <- "metformin"
met <- NormalizeData(met)
met <- FindVariableFeatures(met)

#Integration
immune.anchor <- FindIntegrationAnchors(object.list = list(control, met), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchor, dims = 1:20)

#Integrated analysis
#visualization and clustering
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined)
immune.combined <- RunPCA(immune.combined, npcs = 20)

 #immune.combined <- JackStraw(immune.combined, dims = 100)
 #immune.combined <- ScoreJackStraw(immune.combined, dims = 100)
 #JackStrawPlot(immune.combined, dims = 1:100) 
ElbowPlot(immune.combined, ndims = 20)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

p1<-DimPlot(immune.combined, group.by = "batch", reduction = "umap")
p2<-DimPlot(immune.combined, group.by = "ident", reduction = "umap")
plot_grid(p1,p2)
DimPlot(immune.combined, reduction = "umap", split.by = "batch")
DimPlot(immune.combined, reduction = "tsne", group.by = "batch")
DimPlot(immune.combined, group.by = "ident", reduction = "umap", label = TRUE)

#Identify conserved cell type markers
DefaultAssay(immune.combined) <- "RNA"
marker.0 <-FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = "batch") #Microglia
marker.10 <-FindConservedMarkers(immune.combined, ident.1 = 10, grouping.var = "batch") #Microglia
marker.5 <-FindConservedMarkers(immune.combined, ident.1 = 5, grouping.var = "batch") #Microglia
marker.13 <-FindConservedMarkers(immune.combined, ident.1 = 13, grouping.var = "batch") #Microglia
marker.2 <- FindConservedMarkers(immune.combined, ident.1 = 2, grouping.var = "batch")
marker.4 <- FindConservedMarkers(immune.combined, ident.1 = 4, grouping.var = "batch")#Tmem
marker.6 <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "batch")
marker.8 <- FindConservedMarkers(immune.combined, ident.1 = 8, grouping.var = "batch") #CD8+ Tc
marker.9 <- FindConservedMarkers(immune.combined, ident.1 = 9, grouping.var = "batch") #Treg
marker.11 <- FindConservedMarkers(immune.combined, ident.1 = 11, grouping.var = "batch")
marker.12 <- FindConservedMarkers(immune.combined, ident.1 = 12, grouping.var = "batch") #multi_potent,,,?????????
marker.17 <- FindConservedMarkers(immune.combined, ident.1 = 17, grouping.var = "batch")
marker.18 <- FindConservedMarkers(immune.combined, ident.1 = 18, grouping.var = "batch") #gamma delta
marker.15 <- FindConservedMarkers(immune.combined, ident.1 = 15, grouping.var = "batch") #NK
marker.3 <- FindConservedMarkers(immune.combined, ident.1 = 3, grouping.var = "batch") #NK
marker.1 <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "batch") #B cell
marker.7 <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "batch") #macrophage?
marker.14 <- FindConservedMarkers(immune.combined, ident.1 = 14, grouping.var = "batch") #Dendritic
marker.16 <- FindConservedMarkers(immune.combined, ident.1 = 16, grouping.var = "batch") #macrophage?
marker.19 <- FindConservedMarkers(immune.combined, ident.1 = 19, grouping.var = "batch") #?? Mega??
marker.20 <- FindConservedMarkers(immune.combined, ident.1 = 20, grouping.var = "batch")
marker.21 <- FindConservedMarkers(immune.combined, ident.1 = 21, grouping.var = "batch") #Basophil
marker.22 <- FindConservedMarkers(immune.combined, ident.1 = 22, grouping.var = "batch") #Fibroblast
marker.23 <- FindConservedMarkers(immune.combined, ident.1 = 23, grouping.var = "batch") #Neutrophil
marker.24 <- FindConservedMarkers(immune.combined, ident.1 = 24, grouping.var = "batch") #??

#
FeaturePlot(immune.combined, features = c("Gpr34")) #,"Tmem119","P2ry12","Cx3cr1")) #Microglia = cluster0,5,10,13
FeaturePlot(immune.combined, features = c("Cd79b")) #,"Cd79a","Cd79b","Ly6d","Iglc2","Iglc3")) #B cell = cluster 1
FeaturePlot(immune.combined, features = c("Gata2")) #,"Hdc","Gata2","Ms4a2","Mcpt8"))  # Basophil = cluster21
FeaturePlot(immune.combined, features = c("Igfbp6")) #,"Col5a2","Igfbp6","Col1a1","Col6a2")) #Fibroblast = clust22
FeaturePlot(immune.combined, features = c("Retnlg")) #,"Cxcr2","Mmp8","S100a9","S100a8")) #Neutrophil = cluster23
FeaturePlot(immune.combined, features = c("Retnlg","Cxcr2","Mmp8","Ly6g")) # = cluster24
FeaturePlot(immune.combined, features = c("Ace","Ear2","Adgre4","Treml4","Fabp4","Cd209a","Cd300e","Nxpe4")) # = 
FeaturePlot(immune.combined, features = c("Cd8b1","Cd8a","Gzmk","Cxcr6","Klrc1")) # = cluster8
FeaturePlot(immune.combined, features = c("Pclaf","Ube2c","Birc5","Cenpf","Cdca8","Ccnb2")) # = cluster11
FeaturePlot(immune.combined, features = c("Il2ra","Ctla4","Foxp3","Izumo1r","Tnfrsf4")) #Treg = cluster9
FeaturePlot(immune.combined, features = c("Klra4")) #,"Ncr1","Klra8","Klrb1f")) #NK cell = cluster3,15
FeaturePlot(immune.combined, features = c("Kmo")) #,"Flt3","Kmo","Clec9a","Cd209a")) #Dendritic = cluster14 
FeaturePlot(immune.combined, features = c("Trdc"))
FeaturePlot(immune.combined, features = c("Trbc2","Trac","Cd3g", "Cd3e","Cd3d")) # T cells

#cell 개수
nrow(immune.combined@meta.data) #Total = 11647
nrow(immune.combined@meta.data %>% filter(immune.combined@meta.data$seurat_clusters %in% c(0,5,10,13))) #Microglia = 2729
nrow(immune.combined@meta.data %>% filter(immune.combined@meta.data$seurat_clusters %in% c(2,4,6,8,9,11,17,18))) #T cell = 4649
nrow(immune.combined@meta.data %>% filter(immune.combined@meta.data$seurat_clusters %in% c(18))) #gamma delta T cell = 139 (control=75/ metformin=64)

#
#other_than_25 <- immune.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% filter(integrated_snn_res.2 != 25) %>% .$cell_id

trdc_dt <- immune.combined@assays$RNA@data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column('gene_name') %>% filter(gene_name == 'Trdc') %>% column_to_rownames('gene_name') %>% t() %>% as.data.frame()
ggplot(trdc_dt)+
  geom_histogram(aes(Trdc), binwidth=0.1)+
  coord_cartesian(ylim=c(0,100))

target_cells <- trdc_dt %>% rownames_to_column('cell_id') %>% filter(Trdc > 1) %>% .$cell_id
immune.combined@meta.data <- immune.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% mutate(trdc_info= ifelse(cell_id %in% target_cells & integrated_snn_res.0.7 == 16, 'pos25',
                                                                                                             ifelse(cell_id %in% target_cells, 'pos','neg'))) %>% column_to_rownames('cell_id')
#trdc+ subset
trdc.combined <- immune.combined[,(immune.combined@meta.data$trdc_info %in% c('pos','pos25'))]
 #DimPlot(trdc.combined, reduction = "umap", group.by = "batch") %>% + ggtitle("gamma delta T cell")
 #DimPlot(trdc.combined, reduction = "umap", group.by = "ident", label = TRUE) %>% + ggtitle("gamma deltaT cell")

trdc.combined <- ScaleData(trdc.combined)
trdc.combined <- RunPCA(trdc.combined, npcs = 100)
ElbowPlot(trdc.combined, ndims = 100)

trdc.combined <- RunUMAP(trdc.combined, reduction = "pca", dims = 1:100)
trdc.combined <- FindNeighbors(trdc.combined, reduction = "pca", dims = 1:100)
trdc.combined <- FindClusters(trdc.combined, resolution = 1.5)

p_t1<-DimPlot(trdc.combined, group.by = "batch", reduction = "umap")
p_t2<-DimPlot(trdc.combined, group.by = "ident", reduction = "umap", label = TRUE)
plot_grid(p_t1,p_t2)
DimPlot(trdc.combined, reduction = "umap", split.by = "batch")
DimPlot(trdc.combined, reduction = "tsne", group.by = "batch")
DimPlot(trdc.combined, reduction = "umap", group.by = "trdc_info")
DimPlot(trdc.combined, group.by = "ident", reduction = "umap", label = TRUE)

FeaturePlot(trdc.combined, features = c("Trdc"))
FeaturePlot(trdc.combined, features = c("Ifng", "Cd27", "Il17a","Pdcd1")) #,"Klra1"))

#gd_t DEG
gd_t <- subset(trdc.combined, idents = c("0","1","2","3","4","5","6"))
Idents(gd_t) <- "batch"
avg.gd_t <- log1p(AverageExpression(gd_t, verbose = FALSE)$RNA)
avg.gd_t$gene <- rownames(avg.gd_t)

gd_t.response <- FindMarkers(gd_t, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
gd_t.response <- gd_t.response %>% rownames_to_column()
gd_t_deg <- gd_t.response %>% filter(gd_t.response$p_val < 0.05)
head(gd_t.response)

p_gd_t <-ggplot(avg.gd_t, aes(control, metformin)) + geom_point() + ggtitle("gdT cell")+ coord_cartesian(xlim = c(0:6), ylim = c(0:6))
p_gd_t

gd_t_deg_ordr <- gd_t_deg [ order(gd_t_deg$avg_logFC), ] 
DoHeatmap(gd_t, features = gd_t_deg_ordr$rowname) 

#plot 그리기 - 특정 유전자에 대해서!
tmp_dt1 <-  gd_t@assays$RNA@data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column('gene_name') %>% filter(gene_name == 'Gzmb') %>% t() %>% as.data.frame() %>% rownames_to_column("cell_id")
tmp_dt1 <- tmp_dt1[-1,]
colnames(tmp_dt1) <- c("cell_id", "Gzmb")
tmp_dt2 <- gd_t@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% select(cell_id, batch)
tmp_dt <- left_join(tmp_dt1, tmp_dt2)
tmp_dt$Gzmb <- tmp_dt$Gzmb %>% as.numeric()
head(tmp_dt)

library(ggsignif)
ggplot(tmp_dt, aes(x=batch, y=Gzmb))+
  geom_boxplot(fatten = NULL) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  stat_summary(fun.y=mean) +
  geom_signif(xmin=start, xmax=end, y_position = 12, annotations = '*', tip_length = c(0,0))
  
  #tmp_dt1 <- trdc_dt %>% rownames_to_column('cell_id')
  #tmp_dt2 <- immune.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% select(cell_id, batch) 
  #tmp_dt <- left_join(tmp_dt1, tmp_dt2)
  #library(ggsignif)
  #ggplot(tmp_dt, aes(x=batch, y=Trdc))+
   #geom_point()+
   #geom_boxplot()+
   #geom_signif(xmin=1, xmax=2, y_position = 6, annotations = 'pvalue=0.05', tip_length = c(0,0))




#T cell subset 찾기..
t.combined <- immune.combined[,(immune.combined@active.ident %in% c(2,4,6,8,9,11,17,18))]
DimPlot(t.combined, reduction = "umap", group.by = "batch") %>% + ggtitle("T cell")
DimPlot(t.combined, reduction = "umap", group.by = "ident", label = TRUE) %>% + ggtitle("T cell")

t.combined <- ScaleData(t.combined)
t.combined <- RunPCA(t.combined, npcs = 50)
 # t.combined <- JackStraw(t.combined, dims = 30)
 # t.combined <- ScoreJackStraw(t.combined, dims = 30)
 # JackStrawPlot(t.combined, dims = 1:30) 
ElbowPlot(t.combined)

t.combined <- RunUMAP(t.combined, reduction = "pca", dims = 1:20)
t.combined <- RunTSNE(t.combined, reduction = "pca", dims = 1:20)
t.combined <- FindNeighbors(t.combined, reduction = "pca", dims = 1:20)
t.combined <- FindClusters(t.combined, resolution = 2)

p_t1<-DimPlot(t.combined, group.by = "batch", reduction = "umap")
p_t2<-DimPlot(t.combined, group.by = "ident", reduction = "umap", label = TRUE)
plot_grid(p_t1,p_t2)
DimPlot(t.combined, reduction = "umap", split.by = "batch")
DimPlot(t.combined, reduction = "tsne", group.by = "batch")
DimPlot(t.combined, group.by = "ident", reduction = "umap", label = TRUE)

t.marker.0 <-FindConservedMarkers(t.combined, ident.1 = 0, grouping.var = "batch") #Treg
t.marker.1 <-FindConservedMarkers(t.combined, ident.1 = 1, grouping.var = "batch")
t.marker.2 <-FindConservedMarkers(t.combined, ident.1 = 2, grouping.var = "batch") #??????????
t.marker.3 <-FindConservedMarkers(t.combined, ident.1 = 3, grouping.var = "batch")
t.marker.4 <-FindConservedMarkers(t.combined, ident.1 = 4, grouping.var = "batch")
t.marker.5 <-FindConservedMarkers(t.combined, ident.1 = 5, grouping.var = "batch")
t.marker.6 <-FindConservedMarkers(t.combined, ident.1 = 6, grouping.var = "batch")
t.marker.7 <-FindConservedMarkers(t.combined, ident.1 = 7, grouping.var = "batch")
t.marker.8 <-FindConservedMarkers(t.combined, ident.1 = 8, grouping.var = "batch")
t.marker.9 <-FindConservedMarkers(t.combined, ident.1 = 9, grouping.var = "batch")
t.marker.10 <-FindConservedMarkers(t.combined, ident.1 = 10, grouping.var = "batch")
t.marker.11 <-FindConservedMarkers(t.combined, ident.1 = 11, grouping.var = "batch")
t.marker.12 <-FindConservedMarkers(t.combined, ident.1 = 12, grouping.var = "batch") #gamma delta T cell (control = 57/ metformin = 51)
t.marker.13 <-FindConservedMarkers(t.combined, ident.1 = 13, grouping.var = "batch")
t.marker.14 <-FindConservedMarkers(t.combined, ident.1 = 14, grouping.var = "batch") #????????
t.marker.17 <-FindConservedMarkers(t.combined, ident.1 = 17, grouping.var = "batch") #????????


FeaturePlot(t.combined, features = c("Foxp3")) #"Izumo1r","Foxp3","Il2ra","Tnfrsf18")) # 0 = Treg
FeaturePlot(t.combined, features = c("Tcrg-C1","Trdc","Tcrg-V6","Trdv4")) #
FeaturePlot(t.combined, features = c("Cd3d", "Cd3g","Cd3e"))
FeaturePlot(t.combined, features = c("Fcrla","Pax5","Fcmr"))

#gamma delta T cell
gamma_delta <- subset(t.combined, idents = "14")
Idents(gamma_delta) <- "batch"
avg.gamma_delta <- log1p(AverageExpression(gamma_delta, verbose = FALSE)$RNA)
avg.gamma_delta$gene <- rownames(avg.gamma_delta)

#gamma delta T cell DEG
t.combined$celltype.batch <- paste(Idents(t.combined), t.combined$batch, sep = "_")
t.combined$celltype <- Idents(t.combined)
Idents(t.combined) <- "celltype.batch"
gamma_delta.response <- FindMarkers(t.combined, ident.1 = "14_control", ident.2 = "14_metformin", verbose = FALSE)
gamma_delta.response <- gamma_delta.response %>% rownames_to_column()
gd_deg <- gamma_delta.response %>% filter(gamma_delta.response$p_val < 0.05)
head(gd_deg)

 #genes.to.label = c("Hif1a") #gamma delta
 #genes.to.label = c("Igha","Penk","Gh","Ccl4","Kcnq1ot1") #T reg
p_gd <-ggplot(avg.gamma_delta, aes(control, metformin)) + geom_point() + ggtitle("gamma delta T cell")
 #p_gd <- LabelPoints(plot = p_gd, points = genes.to.label, color)
p_gd

 #genes.to.label = c("Lgals1") #gamma delta
p_gd <-ggplot(avg.gamma_delta, aes(control, metformin)) + geom_point() + ggtitle("gamma delta T cell")+ coord_cartesian(xlim = c(0:6), ylim = c(0:6))
p_gd <-LabelPoints(plot = p_gd, points = gd_deg$rowname, repel = TRUE)
p_gd <- LabelPoints(plot = p_gd, points = genes.to.label, color= "red")
p_gd

gd_deg_ordr <- gd_deg [ order(gd_deg$avg_logFC), ] 
DoHeatmap(gamma_delta, features = gd_deg_ordr$rowname) 

library(NMF)

tmp_mt <- as.matrix(gamma_delta@assays$RNA@data) 
tmp_mt <- tmp_mt %>% as.data.frame %>% rownames_to_column() 
tmp_mt_merge<- merge(tmp_mt, gd_deg_ordr, by.x = "rowname", by.y = "rowname")
tmp_mt_merge <- tmp_mt_merge [ order(tmp_mt_merge$avg_logFC), ] 
tmp_mt_merge<-tmp_mt_merge[,-(113:117)]
rownames(tmp_mt_merge)<-tmp_mt_merge$rowname
tmp_mt_merge<-tmp_mt_merge[,-1]
head(rownames(tmp_mt_merge))

gd_annotation<- gamma_delta@meta.data %>% rownames_to_column() %>% select(rowname, batch) %>% arrange(batch) %>% column_to_rownames('rowname')
gd_annotation

tmp_mt_merge <- tmp_mt_merge[,rownames(gd_annotation)]
aheatmap(tmp_mt_merge, scale = "row", annCol = gd_annotation, Rowv = NA, Colv = NA, height = 30, cexRow = 10)



