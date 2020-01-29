library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
library(gplots)
library(cowplot)

#Seurat-Integrating
#Setup control & met
data_dir <- './projects/00_cowork/jhpark_scRNA/control/outs/filtered_feature_bc_matrix'
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
immune.combined <- RunPCA(immune.combined, npcs = 100)

immune.combined <- JackStraw(immune.combined, dims = 100)
immune.combined <- ScoreJackStraw(immune.combined, dims = 100)
JackStrawPlot(immune.combined, dims = 1:100) 
ElbowPlot(immune.combined)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:100)
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:100)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:100)
immune.combined <- FindClusters(immune.combined, resolution = 0.7)

p1<-DimPlot(immune.combined, group.by = "batch", reduction = "umap")
p2<-DimPlot(immune.combined, group.by = "ident", reduction = "umap")
plot_grid(p1,p2)
DimPlot(immune.combined, reduction = "umap", split.by = "batch")
DimPlot(immune.combined, reduction = "tsne", group.by = "batch")
DimPlot(immune.combined, group.by = "ident", reduction = "umap", label = TRUE)

#Identify conserved cell type markers
DefaultAssay(immune.combined) <- "RNA"
marker.0 <-FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = "batch") #T
marker.1 <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "batch") #Microglia
marker.2 <- FindConservedMarkers(immune.combined, ident.1 = 2, grouping.var = "batch") #NK
marker.3 <- FindConservedMarkers(immune.combined, ident.1 = 3, grouping.var = "batch") #B
marker.4 <- FindConservedMarkers(immune.combined, ident.1 = 4, grouping.var = "batch")#
marker.5 <-FindConservedMarkers(immune.combined, ident.1 = 5, grouping.var = "batch") #Microglia
marker.6 <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "batch") #T
marker.7 <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "batch") #T
marker.8 <- FindConservedMarkers(immune.combined, ident.1 = 8, grouping.var = "batch") #T
marker.9 <- FindConservedMarkers(immune.combined, ident.1 = 9, grouping.var = "batch") #T
marker.10 <-FindConservedMarkers(immune.combined, ident.1 = 10, grouping.var = "batch") # ???
marker.11 <- FindConservedMarkers(immune.combined, ident.1 = 11, grouping.var = "batch") #Microglia
marker.12 <- FindConservedMarkers(immune.combined, ident.1 = 12, grouping.var = "batch") #T
marker.13 <-FindConservedMarkers(immune.combined, ident.1 = 13, grouping.var = "batch") #
marker.14 <- FindConservedMarkers(immune.combined, ident.1 = 14, grouping.var = "batch") #Microglia
marker.15 <- FindConservedMarkers(immune.combined, ident.1 = 15, grouping.var = "batch") #
marker.16 <- FindConservedMarkers(immune.combined, ident.1 = 16, grouping.var = "batch") #T
marker.17 <- FindConservedMarkers(immune.combined, ident.1 = 17, grouping.var = "batch") #
marker.18 <- FindConservedMarkers(immune.combined, ident.1 = 18, grouping.var = "batch") #T
marker.19 <- FindConservedMarkers(immune.combined, ident.1 = 19, grouping.var = "batch") #Fibroblast
marker.20 <- FindConservedMarkers(immune.combined, ident.1 = 20, grouping.var = "batch") #Plasmacytoid Dendritic/ pDC
marker.21 <- FindConservedMarkers(immune.combined, ident.1 = 21, grouping.var = "batch") #Neutrophil/PMN-MDSC
marker.22 <- FindConservedMarkers(immune.combined, ident.1 = 22, grouping.var = "batch") #Mast cell
marker.23 <- FindConservedMarkers(immune.combined, ident.1 = 23, grouping.var = "batch") #Erythroid precursor
marker.24 <- FindConservedMarkers(immune.combined, ident.1 = 24, grouping.var = "batch") #
marker.25 <- FindConservedMarkers(immune.combined, ident.1 = 25, grouping.var = "batch") #Basophil
marker.26 <- FindConservedMarkers(immune.combined, ident.1 = 26, grouping.var = "batch") #

immune.combined_ident <- RenameIdents(immune.combined, `0` = "T cell", `1` = "Microglia", `2` = "NK cell", 
                           `3` = "B cell", `4` = "Myeloid lineage agranulocyte", `5` = "Microglia", `6` = "T cell", `7` = "T cell", `8` = "T cell", `9` = "T cell", 
                           `10` = "NA", `11` = "Microglia", `12` = "T cell", `13` = "Myeloid lineage agranulocyte", `14` = "Microglia", `15` = "Myeloid lineage agranulocyte",
                           `16` = "T cell",`17` = "Myeloid lineage agranulocyte", `18` = "T cell", `19` = "Meningeal cell", `20` = "Myeloid lineage agranulocyte", `21`= "Neutrophil",
                           `22` = "Mast cell", `23` = "Erythroid precursor cell", `24` = "NA", `25` = "Basophil", `26` = "NA" )
DimPlot(immune.combined_ident, group.by = "ident", reduction = "umap", label = FALSE)

immune.freq_table <- prop.table(x = table(immune.combined_ident@active.ident, immune.combined_ident@meta.data[,"batch"]), margin = 2)
barplot(height = immune.freq_table)



#
FeaturePlot(immune.combined, features = c("Gpr34"))#,"Tmem119","P2ry12","Cx3cr1")) #Microglia = 1.5.11.14
FeaturePlot(immune.combined, features = c("Cd79b","Cd79a","Cd79b","Ly6d","Iglc2","Iglc3")) #B cell = 3
FeaturePlot(immune.combined, features = c("Igfbp6")) #,"Col5a2","Igfbp6","Col1a1","Col6a2")) #Fibroblast = 19
FeaturePlot(immune.combined, features = c("Retnlg")) #,"Cxcr2","Mmp8","S100a9","S100a8")) #Neutrophil = 21
FeaturePlot(immune.combined, features = c("Ace","Ear2","Adgre4","Treml4","Fabp4","Cd209a","Cd300e","Nxpe4","Cd16")) # = 
FeaturePlot(immune.combined, features = c("Cd8b1","Cd8a","Gzmk","Cxcr6","Klrc1")) # = cluster8
FeaturePlot(immune.combined, features = c("Pclaf","Ube2c","Birc5","Cenpf","Cdca8","Ccnb2")) # = cluster11
FeaturePlot(immune.combined, features = c("Il2ra","Ctla4","Foxp3","Izumo1r","Tnfrsf4")) #Treg = cluster9
FeaturePlot(immune.combined, features = c("Klra4"))#,"Ncr1","Klra8","Klrb1f")) #NK cell = 2
FeaturePlot(immune.combined, features = c("Kmo","Flt3","Kmo","Clec9a","Cd209a")) #Dendritic = cluster14 
FeaturePlot(immune.combined, features = c("Tcrg-C1","Trdc","Tcrg-V6","Trdv4"))
FeaturePlot(immune.combined, features = c("Cd3g"))#"Cd3d","Cd3g","Cd3e")) # T cells
FeaturePlot(immune.combined, features = c("Klra7","Klra6","Xcl1","Klra1"))
FeaturePlot(immune.combined, features = c("Tpsb2")) #Mast cell
FeaturePlot(immune.combined, features = c("Mcpt8")) #Basophil cell
FeaturePlot(immune.combined, features = c("Alas2")) #Erythroid precusor
FeaturePlot(immune.combined, features = c("Klk1")) #pDC
FeaturePlot(immune.combined, features = c("Tcrg-C1","Trdc","Tcrg-V6","Trdv4")) #cluster11 = gamma delta T 
FeaturePlot(immune.combined, features = c("Ccl22")) #17
FeaturePlot(immune.combined, features = c("Ace")) #13
FeaturePlot(immune.combined, features = c("Mrc1")) #15
FeaturePlot(immune.combined, features = c("Ifitm3")) #4
FeaturePlot(immune.combined, features = c("Il16")) #4
FeaturePlot(immune.combined, features = c("Arg2")) #PMN-MDSC = 21
FeaturePlot(immune.combined, features = c("Pax5","Fcrla")) #PMN-MDSC = 21
FeaturePlot(immune.combined, features = c("mt-Atp6")) #PMN-MDSC = 21



#cell 개수
nrow(immune.combined@meta.data) #Total = 11647
nrow(immune.combined@meta.data %>% filter(immune.combined@meta.data$seurat_clusters %in% c(0,5,10,13))) #Microglia = 2729
nrow(immune.combined_ident@meta.data %>% filter(immune.combined_ident@meta.data$seurat_clusters %in% c(21))) #T cell = 4649
nrow(t.combined@meta.data %>% filter(t.combined@meta.data$seurat_clusters %in% c(11))) #gamma delta T cell = 139 (control=75/ metformin=64)

immune.freq_table <- prop.table(x = table(immune.combined@active.ident, immune.combined@meta.data[,"batch"]))
barplot(height = immune.freq_table)


#T cell subset 찾기..
t.combined <- subset(immune.combined, idents = c("0","6","7","8","9","12","16","18"))
DimPlot(t.combined, reduction = "umap", group.by = "batch") %>% + ggtitle("T cell")
DimPlot(t.combined, reduction = "umap", group.by = "ident", label = TRUE) %>% + ggtitle("T cell")
DimPlot(t.combined, reduction = "umap", group.by = "ident", split.by = "batch" ,label = TRUE) %>% + ggtitle("T cell")

t.combined <- ScaleData(t.combined)
t.combined <- RunPCA(t.combined, npcs = 100)
t.combined <- JackStraw(t.combined, dims = 100)
t.combined <- ScoreJackStraw(t.combined, dims = 100)
JackStrawPlot(t.combined, dims = 1:100) 
ElbowPlot(t.combined)

t.combined <- RunUMAP(t.combined, reduction = "pca", dims = 1:100)
t.combined <- RunTSNE(t.combined, reduction = "pca", dims = 1:100)
t.combined <- FindNeighbors(t.combined, reduction = "pca", dims = 1:100)
t.combined <- FindClusters(t.combined, resolution = 0.7)

p_t1<-DimPlot(t.combined, group.by = "batch", reduction = "umap")
p_t2<-DimPlot(t.combined, group.by = "ident", reduction = "umap", label = TRUE)
plot_grid(p_t1,p_t2)
DimPlot(t.combined, reduction = "umap", split.by = "batch", label = TRUE)
DimPlot(t.combined, reduction = "tsne", group.by = "batch")
DimPlot(t.combined, group.by = "ident", reduction = "umap", label = TRUE)

t.marker.0 <-FindConservedMarkers(t.combined, ident.1 = 0, grouping.var = "batch") #CD4-/CD8+
t.marker.1 <-FindConservedMarkers(t.combined, ident.1 = 1, grouping.var = "batch") #CD3+/CD4-/CD8-/ 
t.marker.2 <-FindConservedMarkers(t.combined, ident.1 = 2, grouping.var = "batch") #CD4+
t.marker.3 <-FindConservedMarkers(t.combined, ident.1 = 3, grouping.var = "batch") #Treg
t.marker.4 <-FindConservedMarkers(t.combined, ident.1 = 4, grouping.var = "batch") #CD4-/Cd8+
t.marker.5 <-FindConservedMarkers(t.combined, ident.1 = 5, grouping.var = "batch") #CD4+
t.marker.6 <-FindConservedMarkers(t.combined, ident.1 = 6, grouping.var = "batch") #CD4+
t.marker.7 <-FindConservedMarkers(t.combined, ident.1 = 7, grouping.var = "batch") #CD4-/CD8+
t.marker.8 <-FindConservedMarkers(t.combined, ident.1 = 8, grouping.var = "batch") #CD4+
t.marker.9 <-FindConservedMarkers(t.combined, ident.1 = 9, grouping.var = "batch") #CD4-/CD8+
t.marker.10 <-FindConservedMarkers(t.combined, ident.1 = 10, grouping.var = "batch") #gamma delta
t.marker.11 <-FindConservedMarkers(t.combined, ident.1 = 11, grouping.var = "batch") #CD+ Bcell

FeaturePlot(t.combined, features = c("Izumo1r")) #"Foxp3","Izumo1r","Foxp3","Il2ra","Tnfrsf18")) # cluster3 = Treg ~PMC4384382
FeaturePlot(t.combined, features = c("Trdc")) #Tcrg-C1,"Trdc","Tcrg-V6","Trdv4")) #cluster11 = gamma delta T 
FeaturePlot(t.combined, features = c("Cd8a"))
FeaturePlot(t.combined, features = c("Cd4")) # 5.8.2
FeaturePlot(t.combined, features = c("Cd40lg","Tbc1d4")) #Tm
FeaturePlot(t.combined, features = c("Cd3d","Cd3g","Cd3e"))
FeaturePlot(t.combined, features = c("Fcrla"))#,"Pax5","Fcrla","Fcmr"))
FeaturePlot(immune.combined, features = c("Gja1"))
FeaturePlot(t.combined, features = c("Styk1","Cxcr4","Styk1","Il12rb2","Il12rb1","Tbx21"))

t.combined <- RenameIdents(t.combined, `0` = "CD8 T", `1` = "CD3+/CD4-/CD4-", `2` = "CD4 T", 
                                `3` = "CD4 Treg", `4` = "CD8 T", `5` = "CD4 T", `6` = "CD4 T", `7` = "CD8 T", `8` = "CD4 T", `9` = "CD8 T", 
                                `10` = "gamma delta T", `11` = "CD3+ B")

t.freq_table <- prop.table(x = table(t.combined@active.ident, t.combined@meta.data[,"batch"]), margin = 2)
barplot(height = t.freq_table)

#gamma delta T cell
gamma_delta <- subset(t.combined, idents = "10")
Idents(gamma_delta) <- "batch"
avg.gamma_delta <- log1p(AverageExpression(gamma_delta, verbose = FALSE)$RNA)
avg.gamma_delta$gene <- rownames(avg.gamma_delta)
avg.gamma_delta2 <- AverageExpression(gamma_delta, verbose = FALSE)$RNA
avg.gamma_delta2$gene <- rownames(avg.gamma_delta2)

gamma_delta.response <- FindMarkers(gamma_delta, ident.1 = "control", ident.2 = "metformin", verbose = TRUE)
gamma_delta.response <- gamma_delta.response %>% rownames_to_column()
gd_deg <- gamma_delta.response %>% filter(gamma_delta.response$p_val < 0.05)
head(gamma_delta.response)

tmp_mt <- gamma_delta@assays$RNA@data %>% as.matrix()
m_ids <- rownames(gamma_delta@meta.data)[gamma_delta@meta.data$batch == 'metformin']
c_ids <- rownames(gamma_delta@meta.data)[gamma_delta@meta.data$batch == 'control']
a <- tmp_mt['Gzmb',colnames(tmp_mt) %in% m_ids]
b <- tmp_mt['Gzmb',colnames(tmp_mt) %in% c_ids]
tmp_dt <- tmp_mt %>% as.data.frame() %>% rownames_to_column('gene')
tmp_dt <- tmp_dt %>% gather(-gene, key=sample, value=exp)
tmp_dt2 <- gamma_delta@meta.data %>% as.data.frame() %>% rownames_to_column('sample') %>% select(sample, batch)
tmp_dt <- left_join(tmp_dt, tmp_dt2)
tmp_dt <- tmp_dt %>% filter(gene == 'Gzmb')
colnames(tmp_dt)
ggplot(tmp_dt, aes(x=batch, y=log10(exp+0.01)))+
  geom_violin()+
  geom_jitter(width=0.2, height=0, alpha=0.3)

tmp_dt %>% filter(batch == 'metformin') %>% .$exp -> a
tmp_dt %>% filter(batch == 'control') %>% .$exp -> b
t.test(a,b)

tmp_dt <- tmp_dt %>% mutate(group = ifelse(exp == 0, 'nonExp','Exp'))

tmp_dt %>% group_by(batch, group) %>% count() %>% .$n -> a


gd_deg_gzmb_hilpda <- gd_deg %>% filter(gd_deg$rowname %in% c("Gzmb", "Hilpda"))

genes.to.label = c("Gzmb","Hilpda") #gamma delta
p_gd <-ggplot(avg.gamma_delta, aes(control, metformin)) + geom_point(size = 0.5) + ggtitle("gamma delta T cell")+ coord_cartesian(xlim = c(0:6), ylim = c(0:6))
p_gd <-LabelPoints(plot = p_gd, points = gd_deg$rowname, repel = TRUE)
 #p_gd <- LabelPoints(plot = p_gd, points = gd_deg_gzmb_hilpda$rowname , color= "red")
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


#특정 유전자. 
tmp_dt1 <-  gamma_delta@assays$RNA@data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column('gene_name') %>% filter(gene_name == 'Hilpda') %>% t() %>% as.data.frame() %>% rownames_to_column("cell_id")
tmp_dt1 <- tmp_dt1[-1,]
colnames(tmp_dt1) <- c("cell_id", "Hilpda")
tmp_dt2 <- gamma_delta@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% select(cell_id, batch)
tmp_dt <- left_join(tmp_dt1, tmp_dt2)
tmp_dt$Hilpda <- tmp_dt$Hilpda %>% as.numeric()
head(tmp_dt)

library(ggsignif)
ggplot(tmp_dt, aes(x=batch, y=Hilpda))+
  geom_boxplot(fatten = NULL) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  stat_summary(fun.y=mean) +
  geom_signif(xmin=start, xmax=end, y_position = 12, annotations = '*', tip_length = c(0,0))


RidgePlot(gamma_delta, features = c("Gzmb", "Hilpda"), ncol = 2)
RidgePlot(gamma_delta, features = c("Hilpda"), ncol = 2)



#cd3_b
cd3_b <- subset(t.combined, idents = "11")
Idents(cd3_b) <- "batch"
avg.cd3_b <- log1p(AverageExpression(cd3_b, verbose = FALSE)$RNA)
avg.cd3_b$gene <- rownames(avg.cd3_b)

cd3_b.response <- FindMarkers(cd3_b, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
cd3_b.response <- cd3_b.response %>% rownames_to_column()
cd3_b_deg <- cd3_b.response %>% filter(cd3_b.response$p_val < 0.05)
head(cd3_b.response)

genes.to.label = c("Lgals1") #gamma delta
p_cd3_b <-ggplot(avg.cd3_b, aes(control, metformin)) + geom_point() + ggtitle("Cd3+ B cell")+ coord_cartesian(xlim = c(0:6), ylim = c(0:6))
p_cd3_b <-LabelPoints(plot = p_cd3_b, points = gd_deg$rowname, repel = TRUE)
p_cd3_b <- LabelPoints(plot = p_cd3_b, points = genes.to.label, color= "red")
p_cd3_b

cd3_b_deg_ordr <- cd3_b_deg [ order(cd3_b_deg$avg_logFC), ] 
DoHeatmap(cd3_b, features = cd3_b_deg_ordr$rowname) 


library(NMF)

tmp_cd3_b <- as.matrix(cd3_b@assays$RNA@data) 
tmp_cd3_b <- tmp_cd3_b %>% as.data.frame %>% rownames_to_column() 
tmp_cd3_b_merge<- merge(tmp_cd3_b, cd3_b_deg_ordr, by.x = "rowname", by.y = "rowname")
tmp_cd3_b_merge <- tmp_cd3_b_merge [ order(tmp_cd3_b_merge$avg_logFC), ] 

tmp_cd3_b_merge<-tmp_cd3_b_merge[,-(113:117)]
rownames(tmp_cd3_b_merge)<-tmp_cd3_b_merge$rowname
tmp_cd3_b_merge<-tmp_cd3_b_merge[,-1]
head(rownames(tmp_cd3_b_merge))

cd3_b_annotation<- cd3_b@meta.data %>% rownames_to_column() %>% select(rowname, batch) %>% arrange(batch) %>% column_to_rownames('rowname')
head(cd3_b_annotation)

tmp_cd3_b_merge <- tmp_cd3_b_merge[,rownames(cd3_b_annotation)]
aheatmap(tmp_cd3_b_merge, scale = "row", annCol = cd3_b_annotation, Rowv = NA, Colv = NA, height = 30, cexRow = 10)


#cluster19
cluster19<-subset(immune.combined, idents = "19")
Idents(cluster19) <- "batch"
avg.cluster19 <- log1p(AverageExpression(cluster19, verbose = FALSE)$RNA)
avg.cluster19$gene <- rownames(avg.cluster19)

cluster19.response <- FindMarkers(cluster19, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
head(cluster19.response)

genes.to.label_c19 = c("Gh","Gzmb","Plac8")
p_c19 <- ggplot(avg.cluster19, aes(control, metformin)) + geom_point() + ggtitle("Fibroblast")
p_c19 <- LabelPoints(plot = p_c19, points = genes.to.label)
p_c19

#Treg
treg<- subset(t.combined, idents = c("3"))
Idents(treg) <- "batch"
avg.treg <- log1p(AverageExpression(treg, verbose = FALSE)$RNA)
avg.treg$gene <- rownames(avg.treg)

treg.response <- FindMarkers(treg, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
head(treg.response)

genes.to.label_treg = c()
p_treg <- ggplot(avg.treg, aes(control, metformin)) + geom_point() +ggtitle("Treg") + coord_cartesian(xlim = c(0:6), ylim = c(0:6))
p_treg <- LabelPoints(plot = treg, points = genes.to.label_treg)
p_treg

#Microglia
microglia <- subset(immune.combined, idents = c("1","5","11","14"))
Idents(microglia) <- "batch"
avg.microglia <- log1p(AverageExpression(microglia, verbose = FALSE)$RNA)
avg.microglia$gene <- rownames(avg.microglia)

microglia.response <- FindMarkers(microglia, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
head(microglia.response)

genes.to.label_microglia = c()
p_microglia <- ggplot(avg.microglia, aes(control, metformin)) + geom_point() +ggtitle("Microglia")
p_microglia <- LabelPoints(plot = microglia, points = genes.to.label_microglia)
p_microglia

#Neutrophil/PMN-MDSC
neutro_mdsc <- subset(immune.combined, idents = c("21"))
Idents(neutro_mdsc) <- "batch"
avg.neutro_mdsc <- log1p(AverageExpression(neutro_mdsc, verbose = FALSE)$RNA)
avg.neutro_mdsc$gene <- rownames(avg.neutro_mdsc)

neutro_mdsc.response <- FindMarkers(neutro_mdsc, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
neutro_mdsc_deg<- neutro_mdsc.response %>% filter(neutro_mdsc.response$p_val < 0.05)
head(neutro_mdsc.response)

genes.to.label_neutro_mdsc = c()
p_neutro_mdsc <- ggplot(avg.neutro_mdsc, aes(control, metformin)) + geom_point() +ggtitle("neutrophil")
p_neutro_mdsc <- LabelPoints(plot = neutro_mdsc, points = genes.to.label_neutro_mdsc)
p_neutro_mdsc

neutro_mdsc.response <- neutro_mdsc.response %>% rownames_to_column()
neutro_mdsc_deg<- neutro_mdsc.response %>% filter(neutro_mdsc.response$p_val < 0.05)
neutro_mdsc_ordr <- neutro_mdsc_deg [ order(neutro_mdsc_deg$avg_logFC), ] 
DoHeatmap(neutro_mdsc, features = neutro_mdsc_ordr$rowname) 

#basophil=25/ mast=22
basophil <- subset(immune.combined, idents = c("22"))
Idents(basophil) <- "batch"
avg.basophil <- log1p(AverageExpression(basophil, verbose = FALSE)$RNA)
avg.basophil$gene <- rownames(avg.basophil)

basophil.response <- FindMarkers(basophil, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
p_basophil <- ggplot(avg.basophil, aes(control, metformin)) + geom_point() +ggtitle("Mast cell")
p_basophil <- LabelPoints(plot = basophil, points = genes.to.label_basophil)
p_basophil

basophil.response <- basophil.response %>% rownames_to_column()
basophil_deg<- basophil.response %>% filter(basophil.response$p_val < 0.05)
basophil_ordr <- basophil_deg [ order(basophil_deg$avg_logFC), ] 
DoHeatmap(basophil, features = basophil_ordr$rowname) 

#NK
NK <- subset(immune.combined, idents = c("2"))
Idents(NK) <- "batch"
avg.NK <- log1p(AverageExpression(NK, verbose = FALSE)$RNA)
avg.NK$gene <- rownames(avg.NK)

NK.response <- FindMarkers(NK, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
NK.response <- NK.response %>% filter(NK.response$p_val < 0.05)
head(NK.response)

genes.to.label_NK = c()
p_NK <- ggplot(avg.NK, aes(control, metformin)) + geom_point() +ggtitle("NK")
p_NK <- LabelPoints(plot = NK, points = genes.to.label_NK)
p_NK

#B
B <- subset(immune.combined, idents = c("3"))
Idents(B) <- "batch"
avg.B <- log1p(AverageExpression(B, verbose = FALSE)$RNA)
avg.B$gene <- rownames(avg.B)

B.response <- FindMarkers(B, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
B.response <- B.response %>% filter(B.response$p_val < 0.05)
head(B.response)

genes.to.label_B = c()
p_B <- ggplot(avg.B, aes(control, metformin)) + geom_point() +ggtitle("B cell")
p_B <- LabelPoints(plot = B, points = genes.to.label_B)
p_B

#CD8+ cluster0,4,7
c0_4_7<- subset(t.combined, idents = c("0","4","7"))
Idents(c0_4_7) <- "batch"
avg.c0_4_7 <- log1p(AverageExpression(c0_4_7, verbose = FALSE)$RNA)
avg.c0_4_7$gene <- rownames(avg.c0_4_7)

c0_4_7.response <- FindMarkers(c0_4_7, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
head(c0_4_7.response)

genes.to.label_c0_4 = c()
p_c0_4_7 <- ggplot(avg.c0_4_7, aes(control, metformin)) + geom_point() +ggtitle("CD8+")  + coord_cartesian(xlim = c(0:6), ylim = c(0:6))
p_c0_4_7 <- LabelPoints(plot = c0_4_7, points = genes.to.label_c0_4_7)
p_c0_4_7

#CD4+ cluster2,5,8
c2_5_8<- subset(t.combined, idents = c("2","5","8","6"))
Idents(c2_5_8) <- "batch"
avg.c2_5_8 <- log1p(AverageExpression(c2_5_8, verbose = FALSE)$RNA)
avg.c2_5_8$gene <- rownames(avg.c2_5_8)

c2_5_8.response <- FindMarkers(c2_5_8, ident.1 = "control", ident.2 = "metformin", verbose = FALSE)
head(c2_5_8.response)

genes.to.label_c2_5_8 = c()
p_c2_5_8 <- ggplot(avg.c2_5_8, aes(control, metformin)) + geom_point() +ggtitle("CD4+")  + coord_cartesian(xlim = c(0:6), ylim = c(0:6))
p_c2_5_8 <- LabelPoints(plot = c2_5_8, points = genes.to.label_c2_5_8)
p_c2_5_8








> gamma_delta_barplot%>%as.data.frame() %>% write_csv("gamma_delta_readcount.csv")

ggboxplot(gd__bar_test, x = "batch", y = "B", color = "batch", palette = "npg", add = "jitter") + stat_compare_means() + facet_wrap(~rowname)




#
gamma_delta_barplot <- gamma_delta@assays$RNA@counts %>% as.data.frame()%>%rownames_to_column() 
gamma_delta_barplot <- merge(gd_deg_ordr, gamma_delta_barplot, by.x = "rowname", by.y = "rowname")
gamma_delta_barplot <- gamma_delta_barplot[,-(2:6)]
gd__bar_test<- gather(gamma_delta_barplot, "A","B",-rowname)
gd_batch<-gamma_delta@meta.data %>% as.data.frame() %>% rownames_to_column() %>% select(rowname, batch)
gd__bar_test <- merge(gd__bar_test, gd_batch, by.x = "A", by.y = "rowname")


gamma_delta_barplot_test <- merge(gd_deg_ordr, avg.gamma_delta2, by.x = "rowname", by.y = "gene")%>% select("rowname", "metformin", "control")%>%as.data.frame()
gamma_delta_barplot_test<- gather(gamma_delta_barplot, "A","B", -rowname)
gamma_delta_barplot_test$B<-as.numeric(gamma_delta_barplot_test$B)
ggplot(gamma_delta_barplot_test, aes(x = A, y = B)) + facet_wrap(~rowname) + geom_bar(aes(fill = A), stat = 'identity')
barchart(~B|rowname, group = A, data = gamma_delta_barplot_test, rectangles = TRUE, points = TRUE)
