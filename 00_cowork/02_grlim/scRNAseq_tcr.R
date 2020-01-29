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
obj_combined <- merge(obj_1, obj_2)
obj_combined <- merge(obj_combined, obj_3)
obj_combined <- merge(obj_combined, obj_4)
obj_combined <- merge(obj_combined, obj_5)
obj_combined <- merge(obj_combined, obj_6)

#####SCT
#object integration (batch balancing)
obj.list <- SplitObject(obj_combined, split.by = "batch")
for (i in names(obj.list)) {
  obj.list[[i]] <- SCTransform(obj.list[[i]], verbose = FALSE)
}
obj.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.features_re <- c(obj.features,"P2RY14", "ZBTB10", "PMCH", "CAMP", "ASGR2", "BST1", "CCR2", "CD1D", "BHLHE4", "DHX58", "CD180", "DPEP2", "HRH1", "ALOX15", "PTPRC")
length(obj.features_re)

obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = obj.features_re)
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = obj.features_re)
obj_combined <- IntegrateData(anchorset = obj.anchors_re, normalization.method = "SCT")

save.image("~/projects/00_cowork/02_krlim_scRNA/R.Data/scRNAseq_tcr_v4.RData")


obj_combined <- RunPCA(object = obj_combined, verbose = FALSE)
ElbowPlot(obj_combined, ndims = 50)
obj_combined <- FindNeighbors(obj_combined, dims = 1:30, verbose = FALSE)
obj_combined <- FindClusters(obj_combined, verbose = FALSE, resolution = 1.5)

obj_combined <- RunUMAP(object = obj_combined, dims = 1:30)
DimPlot(obj_combined, label = TRUE, split.by = "batch", ncol = 3) 
DimPlot(obj_combined, label = TRUE) 
plots <- DimPlot(obj_combined, group.by = c("batch"), combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
                                                                                                              byrow = TRUE, override.aes = list(size = 2.5))))
CombinePlots(plots)
DimPlot(obj_combined, group.by = "", split.by = "batch", ncol = 3)
DefaultAssay(obj) <- "RNA"

#Normalize RNA data for visualization purposes
obj_combined <- NormalizeData(obj_combined, verbose = FALSE)
all.markers <- FindAllMarkers(obj_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_top50 <- all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
FeaturePlot(obj_combined, c("CD3D","CD3E")) #Tcell ,"CD3G"
FeaturePlot(obj_combined, c("GNLY","NKG7","KLRC1")) #Tcell ,"CD3G"
FeaturePlot(obj_combined, c("CD79A","MS4A1")) #Bcell
FeaturePlot(obj_combined, c("AIF1","CD68")) #macrophage
FeaturePlot(obj_combined, c("LILRA4","PTCRA")) #plasmacytoid dendritic
FeaturePlot(obj_combined, c("IGHG1","IGHG3","IGLC2","IGKC"), ) #Bcell subtypes
FeaturePlot(obj_combined, c("CD8A","CD8B"))

FeaturePlot(obj_combined, c("LAMP3","MARCH3","CD209"))
FeaturePlot(obj_combined, c("CD80"))

FeaturePlot(obj_combined, c("CTLA4","FOXP3")) #plasmacytoid dendritic
FeaturePlot(obj_combined, c("LEF1","CCR7","CD4"))

#Cluster annotation
new.cluster.ids <- c("B cell", "Naive T cell", "Naive T cell", "CD8 T cell", "Naive B cell", "Regulatory T cell", "B cell", "Naive T cell", "B cell","Naive T cell","CD8 T cell","NK cell","Macrophage","Naive T cell","CD8 T cell","CD8 T cell","Regulatory T cell","Macrophage","CD8 T cell","Macrophage","CD8 T cell","CD8 T cell","NK cell","Naive T cell","Plasma B cell","Naive T cell", "Macrophage","Naive T cell","B cell","Regulatory T cell","Cancer","Plasmacytoid Dendritic cell","TCL1A+ B cell","CD8 T cell")
names(new.cluster.ids) <- levels(obj_combined)
plot <-FeaturePlot(obj_combined, c("CD79A")) 
plot
select_cluster <- CellSelector(plot = plot)
Idents(obj_combined, cells = select_cluster) <- "TCL1A+ B cell"

obj_combined_annotation <- RenameIdents(object = obj_combined_annotation,  'TCL1A+ B cell' = 'B cell')
obj_combined_annotation <- RenameIdents(object = obj_combined_annotation,  'Naive T cell' = 'CD4 T cell')

##obj_combined_annotation<-obj_combined
DimPlot(obj_combined, label = FALSE) 
DimPlot(obj_combined, split.by = "batch", ncol = 3) 

#Bar plot -> batch 마다 각 cell type frequency
freq_table <- prop.table(x = table(b_combined_ann@active.ident, b_combined_ann@meta.data[, "batch"]), margin = 2)
barplot(height = freq_table,col = brewer.pal(12,"Set3")[1:4]) #cluster 수에 맞게 color 수 도 설정
barplot(height = freq_table,col = brewer.pal(12,"Set3")[1:4], beside = TRUE) #beside option
legend('topright', legend = rownames(freq_table), fill = brewer.pal(12, "Set3")[1:4], cex = 1) #cluster 수에 맞게 color 수 도 설정

#Dot Plot
markers.to.plot <- c("LILRA4","ELF3","MZB1","AIF1","GNLY","FCER2","CD8A","CD3D", "FOXP3","MS4A1")
DotPlot(obj_combined, features = rev(markers.to.plot), cols = brewer.pal(12,"Accent")[1:], dot.scale = 8, 
        split.by = "batch") + RotatedAxis()


###############Load VDJ data (one csv per run)
tcr_1 <- read.csv("./projects/00_cowork/02_krlim_scRNA/00_data/scRNAseq_microgen/02_align/10xTCR_1/outs/filtered_contig_annotations.csv")
tcr_2 <- read.csv("./projects/00_cowork/02_krlim_scRNA/00_data/scRNAseq_microgen/02_align/10xTCR_2/outs/filtered_contig_annotations.csv")
tcr_3 <- read.csv("./projects/00_cowork/02_krlim_scRNA/00_data/scRNAseq_microgen/02_align/10xTCR_3/outs/filtered_contig_annotations.csv")
tcr_4 <- read.csv("./projects/00_cowork/02_krlim_scRNA/00_data/scRNAseq_microgen/02_align/10xTCR_4/outs/filtered_contig_annotations.csv")
tcr_5 <- read.csv("./projects/00_cowork/02_krlim_scRNA/00_data/scRNAseq_microgen/02_align/10xTCR_5/outs/filtered_contig_annotations.csv")
tcr_6 <- read.csv("./projects/00_cowork/02_krlim_scRNA/00_data/scRNAseq_microgen/02_align/10xTCR_6/outs/filtered_contig_annotations.csv")

# Create a function to trim unwanted "-1" and append sample information before barcodes (barcode에 "-1" 붙어있는것 떼줘야함. )
barcoder <- function(df, prefix,  trim="\\-1"){
  df$barcode <- gsub(trim, "", df$barcode)
  df$barcode <- paste0(prefix, df$barcode)
  df
}
tcr_1 <- barcoder(tcr_1, prefix = "pbmc_")
tcr_2 <- barcoder(tcr_2, prefix = "ovary_")
tcr_3 <- barcoder(tcr_3, prefix = "transverse.colon_")
tcr_4 <- barcoder(tcr_4, prefix = "bladder_")
tcr_5 <- barcoder(tcr_5, prefix = "omentum_")
tcr_6 <- barcoder(tcr_6, prefix = "rt.ex.ilian.lymphnode_")
all_tcr <- rbind(tcr_1, tcr_2, tcr_3, tcr_4, tcr_5, tcr_6)
data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "__")
}
suppressPackageStartupMessages(library(data.table))
all_tcr <- as.data.table(all_tcr)
all_tcr_collapsed <- all_tcr[, lapply(.SD, data_concater), by=barcode]
rownames(all_tcr_collapsed) <- all_tcr_collapsed$barcode

obj_combined <- AddMetaData(obj_combined, metadata = all_tcr_collapsed) # 메타데이터로 정보 넣어주기
DimPlot(obj_combined[,!(obj_combined@meta.data$chain%in%c(NA,"None"))], reduction = "umap", group.by = "chain")

################
FeaturePlot(object = obj_combined, features = "FOXP3", split.by = "batch")
FeaturePlot(object = obj_combined, features = "CTLA4", split.by = "batch")
FeaturePlot(object = obj_combined, features = "FCER2", split.by = "batch")
FeaturePlot(object = obj_combined, features = "RGS13", split.by = "batch")

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
b.cells <- subset(immune_combined, idents = "Plasma B cell")
Idents(b.cells) <- "pbmc"
avg.b.cells <- log1p(AverageExpression(b.cells, verbose = FALSE)$RNA)
avg.b.cells$gene <- rownames(avg.b.cells)

cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA)
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(pbmc, )) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2)

#Scatter plot
b_cell <- subset(obj_combined, idents = c("Plasmacytoid Dendritic cell"))
Idents(b_cell) <- "batch"
avg.b_cell <- log1p(AverageExpression(b_cell, verbose = FALSE)$RNA)
avg.b_cell$gene <- rownames(avg.b_cell)

b_cell.response <- FindMarkers(b_cell, ident.1 = "ovary", ident.2 = "transverse.colon", verbose = FALSE)
b_cell.response <- b_cell.response %>% rownames_to_column()
b_cell_deg <- b_cell.response %>% filter(b_cell.response$p_val_adj < 0.05)
b_cell_deg_ordr <- b_cell_deg [ order(b_cell_deg$avg_logFC), ] 

head(b_cell.response)

p1 <- ggplot(avg.b_cell, aes(ovary, transverse.colon)) + geom_point() + ggtitle("pDC")+ coord_cartesian(xlim = c(0:6), ylim = c(0:7))
p1 <- LabelPoints(plot = p1, repel = TRUE, points = c("NFKBIA"))
p1

b_cell_deg_ordr <- b_cell_deg [ order(b_cell_deg$avg_logFC), ] 
DoHeatmap(b_cell, features = b_cell_deg_ordr$rowname) 

#subset----------------------------------------------------------------------------------------------------------------

##T subset--------------------------------------------------------------------------------------
#t_combined은 backup
#t_combined_ann가 진짜!
t_combined <- subset(obj_combined, idents = c("CD4 T cell", "CD8 T cell", "Regulatory T cell"))
DefaultAssay(t_combined) <- "integrated"


#t_combined <- NormalizeData(obj_combined)
#t_combined <- FindVariableFeatures(t_combined)
#t_combined <- t_combined %>% ScaleData()
#dim(t_combined@assays$RNA@scale.data)

t_combined <- t_combined %>% RunPCA(npcs = 50)
ElbowPlot(t_combined,ndims = 50 )
t_combined <- t_combined %>%FindNeighbors(dims = 1:40)
t_combined <- t_combined %>% FindClusters(dims= 1:40, resolution = 1.5)
t_combined <- t_combined %>% RunUMAP(dims = 1:40)
DimPlot(t_combined, reduction = "umap", label = TRUE, group.by = "ident")
DimPlot(t_combined, reduction = "umap", split.by = "batch", label = FALSE)

FeaturePlot(object = t_combined, features = "GNLY", min.cutoff = 0, max.cutoff = 5)
FeaturePlot(object = t_combined, features = "TOX", min.cutoff = 0, max.cutoff = 5)
FeaturePlot(object = t_combined, features = c("CD8A","CD4"), min.cutoff = 2, max.cutoff = 5)
FeaturePlot(object = t_combined, features = "FOXP3")
FeaturePlot(object = t_combined, features = "CCL20", min.cutoff = 0, max.cutoff = 5) #memory activated
FeaturePlot(object = t_combined, features = c("IL21","CXCR5","LRMP","PDCD1"), min.cutoff = 0, max.cutoff = 5)#T follicular helper
FeaturePlot(object = t_combined, features = c("BACH2","TXK"), min.cutoff = 0, max.cutoff = 3) #CD4 naive
FeaturePlot(object = t_combined, features = c("PTGER2", "GNLY","ZFP36L2"), min.cutoff = 0, max.cutoff = 5) #CD4 T memory resting
FeaturePlot(object = t_combined, features = "TNFRSF4") #cd4 memory/ folicular helper/ treg
FeaturePlot(object = t_combined, features = "LTA", min.cutoff = 0, max.cutoff = 4)
FeaturePlot(object = t_combined, features = c("CTLA4","FOXP3"), min.cutoff = 0)
FeaturePlot(object = t_combined, features = c("TBX21", "GATA3"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(object = t_combined, features = c("CD4","CD8A"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(object = t_combined, features = c("SELL","CD8A","CCR7","GNLY","PDCD1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(object = t_combined, features = c("CD8A","CD4","FOXP3","GNLY","PDCD1","CCR7","CD69","IL21"), min.cutoff = 0, max.cutoff = 5, ncol = 4)

marker_9 <- FindMarkers(object = t_combined, ident.1 = "9", ident.2 = NULL)

t_combined_org<-t_combined_ann

##
plot <- DimPlot(t_combined)
plot <-FeaturePlot(t_combined, c("GNLY"), min.cutoff = 0, max.cutoff = 5) 
plot
select_cluster <- CellSelector(plot = plot)
Idents(t_combined, cells = select_cluster) <- "T cell"

#Cluster annotation
new.cluster.ids <- c("Naive CD4 T", "Follicular helper CD4 T", "Regulatory CD4 T", "Naive CD8 T", "Exhausted CD8 T", "Naive CD4 T", "Exhausted CD8 T", "Activated CD4 T", "Cytotoxic CD8 T","Activated CD8 T","Naive CD4 T","Naive CD4 T","Activated CD4 T","Naive CD4 T","Naive CD4 T","Activated CD8 T","Activated CD8 T","Naive CD8 T","Naive CD4 T","Activated CD8 T","Activated CD8 T","Activated CD8 T","Naive CD4 T","Naive CD4 T","24", "Regulatory CD4 T","Activated CD8 T","Regulatory CD4 T","Naive CD4 T","Naive CD4 T","Regulatory CD4 T","Exhausted CD8 T")
names(new.cluster.ids) <- levels(t_combined)
t_combined_ann <- RenameIdents(t_combined, new.cluster.ids)
DimPlot(t_combined_ann, label = TRUE)

plot <-FeaturePlot(t_combined_ann, c("PDCD1"), min.cutoff = 0, max.cutoff = 5) 
plot
select_cluster <- CellSelector(plot = plot)
Idents(t_combined_ann, cells = select_cluster) <- "Exhausted CD8 T"
t_combined_ann <- RenameIdents(object = t_combined_ann,  '24' = 'NK cell')
DimPlot(t_combined_ann, split.by = "batch", ncol = 3)

###1) Tcell subtyping for CD4(11889) & CD8
plot<-DimPlot(t_combined_ann, reduction = "umap", label = TRUE)
plot
select_cluster <- CellSelector(plot = plot)
Idents(t_combined_ann, cells = select_cluster) <- "CD4 T"

which(colnames(t(t_combined@assays$integrated@scale.data))=="CD8A")
rownames(t_combined@assays$integrated@scale.data)[227]
DimPlot(t_combined[,which(t_combined@assays$integrated@scale.data[227,]>0)], reduction = "umap", label = TRUE, group.by = "ident")
which(colnames(t(t_combined@assays$integrated@scale.data))=="CD4")
rownames(t_combined@assays$integrated@scale.data)[1937]
DimPlot(t_combined[,which(t_combined@assays$integrated@scale.data[1937,]>1)], reduction = "umap", label = TRUE, group.by = "ident")

DimPlot(t_combined[,which(t_combined@assays$integrated@scale.data[1937,]>0 ,t_combined@assays$integrated@scale.data[227,] < t_combined@assays$integrated@scale.data[1937,])], reduction = "umap", label = TRUE, group.by = "ident")
DimPlot(t_combined[,which(t_combined@assays$integrated@scale.data[227,]>0 ,t_combined@assays$integrated@scale.data[227,] > t_combined@assays$integrated@scale.data[1937,])], reduction = "umap", label = TRUE, group.by = "ident")

cd8_cells <- colnames(t_combined_ann[,which(t_combined_ann@assays$integrated@scale.data[227,]>0 ,t_combined_ann@assays$integrated@scale.data[227,] > t_combined@assays$integrated@scale.data[1937,])])
Idents(t_combined_ann, cells = cd8_cells) <- "CD8 T"

DimPlot(t_combined_ann, reduction = "umap", label = TRUE, split.by = "ident")

###2.1) Tcell subtyping in CD4------------------------------------------------------------------------------------------------------------------------------
cd4_combined <- subset(t_combined_ann, ident = "CD4 T")

FeaturePlot(cd4_combined, c("FOXP3", "IL2RA"), min.cutoff = 0, max.cutoff = 5)
which(colnames(t(cd4_combined@assays$integrated@scale.data))=="FOXP3")
rownames(cd4_combined@assays$integrated@scale.data)[245]
which(colnames(t(cd4_combined@assays$integrated@scale.data))=="IL2RA")
rownames(cd4_combined@assays$integrated@scale.data)[450]
DimPlot(cd4_combined[,which(cd4_combined@assays$integrated@scale.data[245,]>3|cd4_combined@assays$integrated@scale.data[450,]>3)], reduction = "umap", label = TRUE, group.by = "ident")
plot<-DimPlot(cd4_combined[,which(cd4_combined@assays$integrated@scale.data[245,]>3|cd4_combined@assays$integrated@scale.data[450,]>3)], reduction = "umap", label = TRUE, group.by = "ident")
plot
select_cluster <- CellSelector(plot = plot)
Idents(cd4_combined, cells = select_cluster) <- "Treg cells"
DimPlot(cd4_combined, reduction = "umap", label = TRUE)

FeaturePlot(cd4_combined, c("PDCD1", "TOX"), min.cutoff = 0, max.cutoff = 5)
which(colnames(t(cd4_combined@assays$integrated@scale.data))=="PDCD1")
rownames(cd4_combined@assays$integrated@scale.data)[601]
which(colnames(t(cd4_combined@assays$integrated@scale.data))=="TOX")
rownames(cd4_combined@assays$integrated@scale.data)[1082]
DimPlot(cd4_combined[,which((cd4_combined@assays$integrated@scale.data[601,]>3|cd4_combined@assays$integrated@scale.data[1082,]>3)&cd4_combined@active.ident != "Treg cells")], reduction = "umap", label = TRUE, group.by = "ident")
plot<-DimPlot(cd4_combined[,which((cd4_combined@assays$integrated@scale.data[601,]>3|cd4_combined@assays$integrated@scale.data[1082,]>3)&cd4_combined@active.ident != "Treg cells")], reduction = "umap", label = TRUE, group.by = "ident")
plot
select_cluster <- CellSelector(plot = plot)
Idents(cd4_combined, cells = select_cluster) <- "Exhausted CD4 T cells"
DimPlot(cd4_combined, reduction = "umap", label = TRUE)FeaturePlot(cd4_combined, c("PDCD1", "TOX"), min.cutoff = 0, max.cutoff = 5)

FeaturePlot(cd4_combined, c("CD69","CCL20", "IL2RA", "IFNG","GZMB"), min.cutoff = 0, max.cutoff = 1)
which(colnames(t(cd4_combined@assays$integrated@scale.data))=="CD69")
rownames(cd4_combined@assays$integrated@scale.data)[331]
which(colnames(t(cd4_combined@assays$integrated@scale.data))=="CCL20")
rownames(cd4_combined@assays$integrated@scale.data)[69]
which(colnames(t(cd4_combined@assays$integrated@scale.data))=="IFNG")
rownames(cd4_combined@assays$integrated@scale.data)[164]
which(colnames(t(cd4_combined@assays$integrated@scale.data))=="GZMB")
rownames(cd4_combined@assays$integrated@scale.data)[93]
DimPlot(cd4_combined[,which(cd4_combined@assays$integrated@scale.data[331,]>0 & cd4_combined@active.ident == "CD4 T")], reduction = "umap", label = TRUE, group.by = "ident")
plot<- DimPlot(cd4_combined[,which(cd4_combined@assays$integrated@scale.data[93,]>1 & cd4_combined@active.ident =="CD4 T")], reduction = "umap", label = TRUE, group.by = "ident")
plot
select_cluster <- CellSelector(plot = plot)
Idents(cd4_combined, cells = select_cluster) <- "activated CD4 Tme cells"
DimPlot(cd4_combined, reduction = "umap", label = TRUE)

FeaturePlot(cd4_combined, c("CCR7", "SELL", "PTGER2", "ZFP36L2","CD44"), min.cutoff = 0, max.cutoff = 5)
which(colnames(t(cd4_combined@assays$integrated@scale.data))=="CCR7")
rownames(cd4_combined@assays$integrated@scale.data)[320]
which(colnames(t(cd4_combined@assays$integrated@scale.data))=="SELL")
rownames(cd4_combined@assays$integrated@scale.data)[710]
DimPlot(cd4_combined[,which((cd4_combined@assays$integrated@scale.data[320,]>2|cd4_combined@assays$integrated@scale.data[710,]>2)&cd4_combined@active.ident == "CD4 T")], reduction = "umap", label = TRUE, group.by = "ident")
plot<-DimPlot(cd4_combined[,which((cd4_combined@assays$integrated@scale.data[320,]>1|cd4_combined@assays$integrated@scale.data[710,]>1)&cd4_combined@active.ident == "CD4 T")], reduction = "umap", label = TRUE, group.by = "ident")
plot
select_cluster <- CellSelector(plot = plot)
Idents(cd4_combined, cells = select_cluster) <- "Naive CD4 T cells"
DimPlot(cd4_combined, reduction = "umap", label = TRUE)

plot<-DimPlot(cd4_combined[,which(cd4_combined@active.ident == "Naive CD4 T cells" & cd4_combined@assays$integrated@scale.data[1082,]>0)], reduction = "umap", label = TRUE, group.by = "ident")
plot<-DimPlot(cd4_combined[,cd4_combined@active.ident== "Activated CD4 Tm cells" & cd4_combined@assays$integrated@scale.data[320,]>1 & cd4_combined@assays$integrated@scale.data[710,]>1])
plot
select_cluster <- CellSelector(plot = plot)
Idents(cd4_combined, cells = select_cluster) <- "Naive CD4 T cells"
DimPlot(cd4_combined, reduction = "umap", split.by = "batch")
cd4_combined <- RenameIdents(object = cd4_combined, "Naive CD4 T cells" = "CD4_T_naive", "Treg cells"= "CD4_Treg", "Exhausted CD4 Tm cells"= "CD4_Tem_exhausted", "Activated CD4 Tm cells"= "CD4_Tem_activated")

Idents(t_combined_ann, cells = colnames(cd4_combined[,cd4_combined@active.ident=="CD4_T_naive"])) <- "CD4_T_naive"
Idents(t_combined_ann, cells = colnames(cd4_combined[,cd4_combined@active.ident=="CD4_Treg"])) <- "CD4_Treg"
Idents(t_combined_ann, cells = colnames(cd4_combined[,cd4_combined@active.ident=="CD4_Tem_exhausted"])) <- "CD4_Tem_exhausted"
Idents(t_combined_ann, cells = colnames(cd4_combined[,cd4_combined@active.ident=="CD4_Tem_activated"])) <- "CD4_Tem_activated"
DimPlot(t_combined_ann)

##2)Tcell subtyping by CD8A(15719)------------------------------------------------------------------------------------------------------------------------------- 
cd8_combined <- t_combined[,t_combined_ann@active.ident == "CD8 T"]
DimPlot(cd8_combined)

FeaturePlot(cd8_combined, c("CCR7", "SELL"), min.cutoff = 0, max.cutoff = 5)
which(colnames(t(cd8_combined@assays$integrated@scale.data))=="CCR7")
rownames(cd8_combined@assays$integrated@scale.data)[320]
which(colnames(t(cd8_combined@assays$integrated@scale.data))=="SELL")
rownames(cd8_combined@assays$integrated@scale.data)[710]
DimPlot(cd8_combined[,which((cd8_combined@assays$integrated@scale.data[320,]>0|cd8_combined@assays$integrated@scale.data[710,]>0))], reduction = "umap", label = TRUE, group.by = "ident")
plot<- DimPlot(cd8_combined[,which((cd8_combined@assays$integrated@scale.data[320,]>1|cd8_combined@assays$integrated@scale.data[710,]>1))], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(cd8_combined, cells = select_cluster) <- "Naive CD8 T cells"
DimPlot(cd8_combined, reduction = "umap", label = TRUE)

FeaturePlot(cd8_combined, c("TOX","PDCD1"), min.cutoff = 0, max.cutoff = 5)
which(colnames(t(cd8_combined@assays$integrated@scale.data))=="TOX")
rownames(cd8_combined@assays$integrated@scale.data)[1082]
which(colnames(t(cd8_combined@assays$integrated@scale.data))=="PDCD1")
rownames(cd8_combined@assays$integrated@scale.data)[601]
DimPlot(cd8_combined[,which((cd8_combined@assays$integrated@scale.data[1082,]>1|cd8_combined@assays$integrated@scale.data[601,]>1))], reduction = "umap", label = TRUE, group.by = "ident")
plot<- DimPlot(cd8_combined[,which((cd8_combined@assays$integrated@scale.data[1082,]>1|cd8_combined@assays$integrated@scale.data[601,]>1))], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(cd8_combined, cells = select_cluster) <- "Exhausted CD8 Tem cells"
DimPlot(cd8_combined, reduction = "umap", label = TRUE)

which(colnames(t(cd8_combined@assays$integrated@scale.data))=="FOXP3")
rownames(cd8_combined@assays$integrated@scale.data)[245]
DimPlot(cd8_combined[,which(cd8_combined@assays$integrated@scale.data[245,]>3)], reduction = "umap", label = TRUE, group.by = "ident")
plot<- DimPlot(cd8_combined[,which(cd8_combined@assays$integrated@scale.data[245,]>3)], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(cd8_combined, cells = select_cluster) <- "Treg cells"
DimPlot(cd8_combined, reduction = "umap", label = TRUE)

which(colnames(t(cd8_combined@assays$integrated@scale.data))=="GZMB")
rownames(cd8_combined@assays$integrated@scale.data)[93]
which(colnames(t(cd8_combined@assays$integrated@scale.data))=="PRF1")
rownames(cd8_combined@assays$integrated@scale.data)[396]
DimPlot(cd8_combined[,which((cd8_combined@assays$integrated@scale.data[93,]>2 | cd8_combined@assays$integrated@scale.data[396,]>2) & cd8_combined@active.ident != "Exhausted CD8 Tem cells")], reduction = "umap", label = TRUE, group.by = "ident")
plot<- DimPlot(cd8_combined[,which((cd8_combined@assays$integrated@scale.data[93,]>2 | cd8_combined@assays$integrated@scale.data[396,]>2) & cd8_combined@active.ident != "Exhausted CD8 Tem cells")], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(cd8_combined, cells = select_cluster) <- "CD8 Tc cells"
DimPlot(cd8_combined, reduction = "umap")

FeaturePlot(cd8_combined, c("CD69", "IL2RA"), min.cutoff = 0, max.cutoff = 1)
which(colnames(t(cd8_combined@assays$integrated@scale.data))=="IL2RA")
rownames(cd8_combined@assays$integrated@scale.data)[450]
DimPlot(cd8_combined[,which(cd8_combined@assays$integrated@scale.data[331,]>2 | cd8_combined@assays$integrated@scale.data[450,]>2)], reduction = "umap", label = TRUE, group.by = "ident")
DimPlot(cd8_combined)
plot<- DimPlot(cd8_combined[,which(cd8_combined@assays$integrated@scale.data[331,]>2 | cd8_combined@assays$integrated@scale.data[450,]>2)], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(cd8_combined, cells = select_cluster) <- "Activated CD8 Tem cells"
DimPlot(cd8_combined, reduction = "umap", split.by = "ident")

DimPlot(cd8_combined[,!cd8_combined@active.ident %in% c("Activated CD8 Tem cells","CD8 Tc cells","Exhausted CD8 Tem cells","Treg cells", "Naive CD8 T cells")], label = TRUE)

cd8_combined <- RenameIdents(object = cd8_combined,  '0' = 'Naive CD8 T cells', "1" = "Treg cells", "2"= "Treg cells", "3" = "3", "4" = "Exhausted CD8 Tem cells","6"= "Exhausted CD8 Tem cells","7"= "Activated CD8 Tem cells","8"= "CD8 Tc cells", "10"= "Naive CD9 T cells")
cd8_combined <- RenameIdents(object = cd8_combined,  'Exhausted CD9 Tem cells' = 'Exhausted CD8 Tem cells', "5"= "Activated CD8 Tem cells", "3" = "Naive CD8 T cells", "9" = "CD8 Tc cells", "11" = "Activated CD8 Tem cells", "14"= "Naive CD8 T cells", "15"= "Activated CD8 Tem cells", "18" = "Activated CD8 Tem cells", "22"= "Naive CD8 T cells", "26"= "Activated CD8 Tem cells", "28" = "Naive CD8 T cells")
cd8_combined <- RenameIdents(object = cd8_combined, "Exhausted CD8 Tem cells" = "CD8_Tem_exhausted", "Activated CD8 Tem cells"= "CD8_Tem_activated", "Naive CD8 T cells" = "CD8_T_naive", "Treg cells"= "CD4_Treg")
cd8_combined <- RenameIdents(object = cd8_combined, "CD8 Tc cells" = "CD8_Tc")
DimPlot(cd8_combined, reduction = "umap")

Idents(t_combined_ann, cells = colnames(cd8_combined[,cd8_combined@active.ident=="CD8_Tem_exhausted"])) <- "CD8_Tem_exhausted"
Idents(t_combined_ann, cells = colnames(cd8_combined[,cd8_combined@active.ident=="CD8_Tem_activated"])) <- "CD8_Tem_activated"
Idents(t_combined_ann, cells = colnames(cd8_combined[,cd8_combined@active.ident=="CD8_T_naive"])) <- "CD8_T_naive"
Idents(t_combined_ann, cells = colnames(cd8_combined[,cd8_combined@active.ident=="CD4_Treg"])) <- "CD4_Treg"
Idents(t_combined_ann, cells = colnames(cd8_combined[,cd8_combined@active.ident=="CD8_Tc"])) <- "CD8_Tc"
DimPlot(t_combined_ann)

DimPlot(cd8_combined[,cd8_combined@active.ident !="CD4_Treg"], split.by = "ident")
DimPlot(t_combined_ann[,t_combined_ann@active.ident %in% c("CD4_Treg","CD4_Tem_exhausted","CD4_Tem_activated","CD4_T_naive")])

FeaturePlot(cd8_combined[,cd8_combined@active.ident !="CD4_Treg"], "TOX",min.cutoff = 0, max.cutoff = 5, split.by = "batch")

##B subset--------------------------------------------------------------------------------------
b_combined <- subset(obj_combined, idents = c("Plasma B cell", "Naive B cell", "B cell"))
DefaultAssay(b_combined) <- "integrated"
b_combined <- b_combined %>% RunPCA(npcs = 50)
ElbowPlot(b_combined,ndims = 50 )
b_combined <- b_combined %>%FindNeighbors(dims = 1:30)
b_combined <- b_combined %>% FindClusters(dims= 1:30, resolution = 1)
b_combined <- b_combined %>% RunUMAP(dims = 1:30)
DimPlot(b_combined, reduction = "umap", label = TRUE)
DimPlot(b_combined, reduction = "umap", split.by = "batch", label = TRUE)

FeaturePlot(object = b_combined, features = c("MZB1", "FCER2","P2RY14","CD72","TCL1A"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(object = b_combined, features = c("MZB1"), min.cutoff = 0, max.cutoff = 3) #plasma
FeaturePlot(object = b_combined, features = c("TCL1A","FCER2"), min.cutoff = 0, max.cutoff = 3) #naive
FeaturePlot(object = b_combined, features = c("CD27",'IGHD',"IGHE","IGHG1","CD19","FCER2","TCL1A","IL2RA"), min.cutoff = 0, max.cutoff = 5, ncol = 2) #total "IGHM"
FeaturePlot(object = b_combined, features = c("SELL","LY86","IRF8","IL4R"), min.cutoff = 0, max.cutoff = 3) #plasma
FeaturePlot(object = b_combined, features = c("TLR4"), min.cutoff = 0, max.cutoff = 3) #plasma

FeaturePlot(object = b_combined, features = c("TCL1A","FCER2","TRIB2"), min.cutoff = 0, max.cutoff = 3) #Naive
which(colnames(t(b_combined@assays$integrated@scale.data))=="TCL1A")
rownames(b_combined@assays$integrated@scale.data)[239]
which(colnames(t(b_combined@assays$integrated@scale.data))=="FCER2")
rownames(b_combined@assays$integrated@scale.data)[552]
which(colnames(t(b_combined@assays$integrated@scale.data))=="TRIB2")
rownames(b_combined@assays$integrated@scale.data)[2199]
DimPlot(b_combined[,which(b_combined@assays$integrated@scale.data[239,]>0 & b_combined@assays$integrated@scale.data[552,]>0 & b_combined@assays$integrated@scale.data[2199,]>0)], reduction = "umap", label = TRUE, group.by = "ident")
DimPlot(b_combined[,which(b_combined@assays$integrated@scale.data[239,]>2)], reduction = "umap", label = TRUE, group.by = "ident")
plot<- DimPlot(b_combined[,which(b_combined@assays$integrated@scale.data[2199,]>3)], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(b_combined, cells = select_cluster) <- "B_naive"
DimPlot(b_combined, reduction = "umap")

FeaturePlot(object = b_combined, features = c("AIM2","GPR183","CD27"), min.cutoff = 0, max.cutoff = 1) #Memory
which(colnames(t(b_combined@assays$integrated@scale.data))=="AIM2")
rownames(b_combined@assays$integrated@scale.data)[2269]
which(colnames(t(b_combined@assays$integrated@scale.data))=="GPR183")
rownames(b_combined@assays$integrated@scale.data)[430]
which(colnames(t(b_combined@assays$integrated@scale.data))=="CD27")
rownames(b_combined@assays$integrated@scale.data)[761]
DimPlot(b_combined[,which(b_combined@assays$integrated@scale.data[2269,]>0 & b_combined@assays$integrated@scale.data[430,]>0 & b_combined@assays$integrated@scale.data[761,]>0)], reduction = "umap", label = TRUE, group.by = "ident")
DimPlot(b_combined[,which(b_combined@assays$integrated@scale.data[2269,]>2 & b_combined@active.ident != "B_naive" )], reduction = "umap", label = TRUE, group.by = "ident")
plot<- DimPlot(b_combined[,which(b_combined@assays$integrated@scale.data[761,]>1 & b_combined@active.ident != "B_naive" )], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(b_combined, cells = select_cluster) <- "B_memory"
DimPlot(b_combined, reduction = "umap")

FeaturePlot(object = b_combined, features = c("MZB1","TNFRSF17"), min.cutoff = 0, max.cutoff = 5) #plasma
which(colnames(t(b_combined@assays$integrated@scale.data))=="MZB1")
rownames(b_combined@assays$integrated@scale.data)[82]
which(colnames(t(b_combined@assays$integrated@scale.data))=="TNFRSF17")
rownames(b_combined@assays$integrated@scale.data)[670]
DimPlot(b_combined[,which(b_combined@assays$integrated@scale.data[82,]>2 & b_combined@assays$integrated@scale.data[670,]>2)], reduction = "umap", label = TRUE, group.by = "ident")
DimPlot(b_combined[,which(b_combined@assays$integrated@scale.data[82,]>1)], reduction = "umap", label = TRUE, group.by = "ident")
plot<- DimPlot(b_combined[,which(b_combined@assays$integrated@scale.data[82,]>1)], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(b_combined, cells = select_cluster) <- "Plasma Cells"
DimPlot(b_combined, reduction = "umap")

DimPlot(b_combined, reduction = "umap", split.by = "ident")
DimPlot(b_combined, reduction = "umap", split.by = "batch")
DimPlot(b_combined[,!b_combined@active.ident %in% c("B_memory","Plasma Cells","B_naive")], reduction = "umap", label = TRUE)

b_combined <- RenameIdents(object = b_combined, "0"= "B_memory", "1"= "B_naive", "2"= "B_memory", "3"= "B_memory", "4"= "B_memory", "5"= "B_memory", "6"= "B_memory", "7"= "B_naive", "8"= "B_memory", "9"= "B_naive",
                           "10"= "Plasma Cells", "11"= "B_memory", "12"= "B_memory", "13"= "B_memory", "14"= "B_memory", "15"= "Plasma Cells", "16"= "Plasma Cells", "17"= "B_naive", "18"= "B_memory", "19"= "Plasma Cells")


FeaturePlot(b_combined, c("TCL1A","FCER2","AIM2","CD27","MZB1"), min.cutoff = 0, max.cutoff = 3, split.by = "batch")
FeaturePlot(b_combined, c("TCL1A"), min.cutoff = 0, max.cutoff = 3, split.by = "batch", pt.size = 0.3)
FeaturePlot(b_combined, c("FCER2"), min.cutoff = 0, max.cutoff = 3, split.by = "batch", pt.size = 0.3)
FeaturePlot(b_combined, c("AIM2"), min.cutoff = 0, max.cutoff = 3, split.by = "batch", pt.size = 0.3)
FeaturePlot(b_combined, c("CD27"), min.cutoff = 0, max.cutoff = 3, split.by = "batch", pt.size = 0.3)
FeaturePlot(b_combined, c("MZB1"), min.cutoff = 0, max.cutoff = 3, split.by = "batch", pt.size = 0.3)

DimPlot(b_combined, pt.size = 1)
DimPlot(t_combined_ann, pt.size = 0.6)

FeaturePlot(t_combined_ann[,t_combined_ann@active.ident == "CD8"], "CD69", min.cutoff = 0, max.cutoff = 5)
plot<-FeaturePlot(t_combined_ann[,t_combined_ann@active.ident == "CD8"], "CD69", min.cutoff = 0, max.cutoff = 5)
select_cluster <- CellSelector(plot = plot)
Idents(t_combined_ann, cells = select_cluster) <- "CD8_Tem_activated"
DimPlot(t_combined_ann, reduction = "umap")
t_combined_ann <- RenameIdents(object = t_combined_ann,  'CD8_Tem_Exhausted' = 'CD8_Tem_exhausted')

#M subset--------------------------------------------------------------------------------------
m_combined <- subset(obj_combined_annotation, idents = c("Macrophage"))
DefaultAssay(m_combined) <- "integrated"
m_combined <- m_combined %>% RunPCA(npcs = 50)
ElbowPlot(m_combined,ndims = 50 )
m_combined <- m_combined %>%FindNeighbors(dims = 1:20)
m_combined <- m_combined %>% FindClusters(dims= 1:20, resolution = 1)
m_combined <- m_combined %>% RunUMAP(dims = 1:20)
DimPlot(m_combined, reduction = "umap", label = TRUE)
DimPlot(m_combined, reduction = "umap", split.by = "batch", label = TRUE)

#Signature
FeaturePlot(object = m_combined, features = c("MARCO","CXCL3"), min.cutoff = 0, max.cutoff = 5) #M0
FeaturePlot(object = m_combined, features = c("CD86","SOCS3","TLR2","TNFAIP6"), min.cutoff = 0, max.cutoff = 5, ncol = 4) #M1
FeaturePlot(object = m_combined, features = c("CD163","PTGER2","MS4A6A","CLEC10A"), min.cutoff = 0, max.cutoff = 5) #M2

#Transcriptional Profiling of the Human Monocyte-to-Macrophage Differentiation and Polarization: New Molecules and Patterns of Gene Expression
FeaturePlot(object = m_combined, features = c("CCR7","BCL2A1","TNF","OASL"), min.cutoff = 0, max.cutoff = 5) #M1 (BCL2A1:M0, M1 모두발현)
FeaturePlot(object = m_combined, features = c("MAF","FGL2","MS4A4A","CD36","EGR2","MS4A6A"), min.cutoff = 0, max.cutoff = 5) #M2 (EGR2: M0도 발현)

#Final
FeaturePlot(object = m_combined, features = c("MARCO","CXCL3"), min.cutoff = 0, max.cutoff = 5) #M0
FeaturePlot(object = m_combined, features = c("CD163","MS4A6A"), min.cutoff = 0, max.cutoff = 5) #M2
FeaturePlot(object = m_combined, features = c("CCR7","BCL2A1","TNF","CD86"), min.cutoff = 0, max.cutoff = 2) #M1 (BCL2A1:M0, M1 모두발현)

#
FeaturePlot(object = m_combined, features = c("MARCO","CXCL3"), min.cutoff = 0, max.cutoff = 3) #M0
which(colnames(t(m_combined@assays$integrated@scale.data))=="MARCO")
rownames(m_combined@assays$integrated@scale.data)[720]
which(colnames(t(m_combined@assays$integrated@scale.data))=="CXCL3")
rownames(m_combined@assays$integrated@scale.data)[123]
which(colnames(t(m_combined@assays$integrated@scale.data))=="BCL2A1")
rownames(m_combined@assays$integrated@scale.data)[286]
DimPlot(m_combined[,which(m_combined@assays$integrated@scale.data[720,]>1 & m_combined@assays$integrated@scale.data[123,]>1 & m_combined@assays$integrated@scale.data[286,]>1),], reduction = "umap", label = TRUE, group.by = "ident")
plot<- DimPlot(m_combined[,which(m_combined@assays$integrated@scale.data[720,]>1 & m_combined@assays$integrated@scale.data[123,]>1 & m_combined@assays$integrated@scale.data[286,]>1),], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(m_combined, cells = select_cluster) <- "M0"
DimPlot(m_combined, reduction = "umap")

FeaturePlot(object = m_combined, features = c("CD163","MS4A6A"), min.cutoff = 0, max.cutoff = 3) #M2
which(colnames(t(m_combined@assays$integrated@scale.data))=="CD163")
rownames(m_combined@assays$integrated@scale.data)[487]
which(colnames(t(m_combined@assays$integrated@scale.data))=="MS4A6A")
rownames(m_combined@assays$integrated@scale.data)[226]
DimPlot(m_combined[,which(m_combined@assays$integrated@scale.data[487,]>1 & m_combined@assays$integrated@scale.data[226,]>1 & m_combined@active.ident != "M0")], reduction = "umap", label = TRUE, group.by = "ident")
plot<- DimPlot(m_combined[,which(m_combined@assays$integrated@scale.data[487,]>1 & m_combined@assays$integrated@scale.data[226,]>1 & m_combined@active.ident != "M0")], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(m_combined, cells = select_cluster) <- "M2"
DimPlot(m_combined, reduction = "umap")

FeaturePlot(object = m_combined, features = c("CCR7","BCL2A1","TNF"), min.cutoff = 0, max.cutoff = 3) #M
which(colnames(t(m_combined@assays$integrated@scale.data))=="CCR7")
rownames(m_combined@assays$integrated@scale.data)[320]
which(colnames(t(m_combined@assays$integrated@scale.data))=="BCL2A1")
rownames(m_combined@assays$integrated@scale.data)[286]
which(colnames(t(m_combined@assays$integrated@scale.data))=="TNF")
rownames(m_combined@assays$integrated@scale.data)[220]

DimPlot(m_combined[,which(m_combined@assays$integrated@scale.data[320,]>2)], reduction = "umap", label = TRUE, group.by = "ident")
plot<-DimPlot(m_combined[,which(m_combined@assays$integrated@scale.data[320,]>2)], reduction = "umap", label = TRUE, group.by = "ident")
select_cluster <- CellSelector(plot = plot)
Idents(m_combined, cells = select_cluster) <- "M1"
DimPlot(m_combined, reduction = "umap")

DimPlot(m_combined, reduction = "umap", label = TRUE)
m_combined <- RenameIdents(object = m_combined, '0' = 'M2', '2' = 'M0', '8' = 'M0', '11' = 'M0', '12' = 'M1')
DimPlot(m_combined, reduction = "umap", label = TRUE, split.by = "batch")
DimPlot(m_combined, reduction = "umap", label = TRUE, split.by = "ident")

#NK subset--------------------------------------------------------------------------------------
nk_combined <- subset(obj_combined_annotation, idents = c("NK cell"))
DefaultAssay(nk_combined) <- "integrated"
nk_combined <- nk_combined %>% RunPCA(npcs = 50)
ElbowPlot(nk_combined,ndims = 50 )
nk_combined <- nk_combined %>%FindNeighbors(dims = 1:20)
nk_combined <- nk_combined %>% FindClusters(dims= 1:20, resolution = 1)
nk_combined <- nk_combined %>% RunUMAP(dims = 1:20)
DimPlot(nk_combined, reduction = "umap", label = TRUE)
DimPlot(nk_combined, reduction = "umap", split.by = "batch", label = TRUE)

FeaturePlot(object = nk_combined, features = c("EGR2"), min.cutoff = 0, max.cutoff = 5) #M0
FeaturePlot(object = nk_combined, features = c("CD163","MS4A6A"), min.cutoff = 0, max.cutoff = 5) #M2