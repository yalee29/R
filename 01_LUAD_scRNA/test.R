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

plot <-FeaturePlot(prox_epi_seurat, reduction = "umap", c("ENSMUSG00000064057")) #
plot
select_cluster <- CellSelector(plot = plot)
Idents(prox_epi_seurat, cells = select_cluster) <- "Clara"

nrow(prox_pnec_basal@meta.data[prox_pnec_basal@active.ident =="Basal",])
nrow(prox_epi_seurat@meta.data[prox_epi_seurat@active.ident =="PNEC",])
nrow(prox_pnec_basal@meta.data[prox_pnec_basal@active.ident =="Ciliated",])
nrow(prox_epi_seurat@meta.data[prox_epi_seurat@active.ident =="Clara",])
DimPlot(prox_epi_seurat, reduction = "umap", group.by = "ident", label = TRUE)
DimPlot(prox_pnec_basal, reduction = "umap", group.by = "ident", label = TRUE)
test1 <-E_seurat1_ad_epi_bb
test2 <- subset(prox_pnec_basal, idents = c("PNEC", "Basal")) #PNEC= 36/ Basal = 228
test4 <- subset(prox_epi_seurat, idents = c("PNEC", "Basal")) #PNEC = 70/ Basal = 239

#Clara, ciliated를 추가해서 batch balancing 해볼 것. (일종의 anchor.?)
test6 <- subset(prox_pnec_basal, idents = c("PNEC", "Basal","Clara","Ciliated")) #PNEC = 70/ Basal = 239/ Ciliated = 215/ Clara =323
test7 <- subset(prox_epi_seurat, idents = c("PNEC", "Basal","Clara","Ciliated")) #PNEC = 70/ Basal = 239/ Ciliated = 431/ Clara = 285

##Preprocessing before merging mouse traj + PNEC & Basal(gene set) -> 11528 genes
gene_set <- intersect(intersect(rownames(test1), rownames(test2)), rownames(test4)) %>% as.data.frame() #11528
colnames(gene_set) <- "gene set"
dim(gene_set)
test1 <- E_seurat1_ad_epi_bb@assays$RNA@counts[rownames(E_seurat1_ad_epi_bb) %in% gene_set$`gene set`,] %>% CreateSeuratObject()
test2 <- test2@assays$RNA@counts[rownames(test2) %in% gene_set$`gene set`,] %>% CreateSeuratObject()
test4 <- test4@assays$RNA@counts[rownames(test4) %in% gene_set$`gene set`,] %>% CreateSeuratObject()

test6 <- test6@assays$RNA@counts[rownames(test6) %in% gene_set$`gene set`,] %>% CreateSeuratObject()
test7 <- test7@assays$RNA@counts[rownames(test7) %in% gene_set$`gene set`,] %>% CreateSeuratObject()

dim(test1)
dim(test2)
dim(test4)
dim(test6)
dim(test7)

#test5 = test1 + test2 + test4 (batch balancing with basal, PNEC)
batch_vector = stringr::str_replace(colnames(test1), "_.*","")
test1$batch = batch_vector
test2$batch <- "P6to8"
test4$batch <- "P6to12"
test5 <- merge(test1, merge(test2, test4))
dim(test5)
unique(test5$batch)

#test8 = test1 +test6 + test7 (batch balancing with basal, PNEC, ciliated, clara)
batch_vector = stringr::str_replace(colnames(test1), "_.*","")
test1$batch = batch_vector
test6$batch <- "P6to8"
test7$batch <- "P6to12"
test8 <- merge(test1, merge(test6, test7))
dim(test8)
unique(test8$batch)

test5 <-NormalizeData(test5)
test5 <-FindVariableFeatures(test5)
test5 <-ScaleData(test5, features = rownames(test5@assays$RNA@counts))
test5 <- RunPCA(test5, npcs = 50)
test5 <- JackStraw(test5, dims = 50)
test5 <- ScoreJackStraw(test5, dims = 1:50)
JackStrawPlot(test5, dims = 1:50)
ElbowPlot(test5, ndims = 50)
test5 <- FindNeighbors(test5, dims = 1:20)
test5 <- FindClusters(test5, resolution = 1) #resolution 변경한 것임. 
test5 <- RunUMAP(test5, dims = 1:20)
DimPlot(test5, reduction = "umap")

pdf("test5.pdf", 7,5)
for(k in c(4:6)){
  for(trim in c(25,30,35,40,45)){
    print(
      test5 %>% RunBBKNN(dims.use = 1:20,
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

test5 <- test5 %>% RunBBKNN(dims.use = 1:20,
                            neighbors_within_batch = 6,
                            trim = 50,
                            batch.key = "batch",
                            python.path = "/home/users/yunah1029/anaconda3/bin/python")
DimPlot(test5, reduction = "bbknn", group.by = "batch") + ggtitle(paste0("bbknn1; neighbors=",6,",trim=",50,",add_PNEC&Basal"))

##Clustering
test5 <- RenameIdents(object = test5, "0" = "AT2", "1" = "AT2", "2" = "AT2", "3" = "AT1", "4" = "AT1", "5" = "Clara", "6" = "Basal", "7" = "Ciliated", "8" = "Basal", "10" = "AT2", "11" = "AT2", "12" = "AT2", "13" = "PNEC", "15" = "AT2" ,"17" ="AT1")
test5 <- RenameIdents(object = test5, "14" = "Fetal lung")

FeaturePlot(test5, c("ENSMUSG00000023951"), reduction = "bbknn") #AT1 "ENSMUSG00000059325","ENSMUSG00000028583","ENSMUSG00000023039","ENSMUSG00000023951"
FeaturePlot(test5, c("ENSMUSG00000041247"), reduction = "bbknn") #AT2 "ENSMUSG00000029375","ENSMUSG00000037071","ENSMUSG00000042784"
FeaturePlot(test5, c("ENSMUSG00000038791","ENSMUSG00000031722"), reduction = "bbknn", max.cutoff = 8) #Clara "ENSMUSG00000022595"
FeaturePlot(test5, c("ENSMUSG00000027676","ENSMUSG00000040703"), reduction = "bbknn") #Ciliated "ENSMUSG00000040703", "ENSMUSG00000024653" ,"ENSMUSG00000027676"
FeaturePlot(test5, reduction = "bbknn", c("ENSMUSG00000021194","ENSMUSG00000030669", "ENSMUSG00000020052"), max.cutoff = 5) #PNEC
FeaturePlot(test5, reduction = "bbknn", c("ENSMUSG00000028435","ENSMUSG00000022676"), max.cutoff = 5) #BASAL

plot <-DimPlot(test5[,test5$batch %in% c("E8.25","E9.5to11.5","E14.5.2","E16.5","E18.5","LungMap.E18.5")], reduction = "bbknn", group.by = "ident")
plot
select_cluster <- CellSelector(plot = plot)
Idents(test5, cells = select_cluster) <- "Fetal lung"

##Gene set 결정
test5_markers <- FindAllMarkers(test5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
table(test5_markers$cluster)
test5_markers <- test5_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) #700
test5_markers_rm_nonepi<-setdiff(test5_markers$gene, nonepi_marker_dt$gene_name) #658

DoHeatmap(test5, features = test5_markers$gene)
DoHeatmap(test5, features = test5_markers_rm_nonepi)

###epi marker + lung subtype markers = final_genes ->738
test5_final <- union(test5_markers_rm_nonepi, cur_markers$gene_name ) 
test5_final<-setdiff(test5_final, luad_genes) #luad_gene들은 제거/ 그대로 유지.
length(test5_final)

###+final gene 수정; 738 (그대로 유지)
### cell cycle related gene 21개 + luad related gene 1개 제거
test5_final <- setdiff(test5_final, c("ENSG00000118971","ENSG00000039068","ENSG00000198668","ENSG00000135446","ENSG00000105976","ENSG00000171848","ENSG00000073350","ENSG00000100714","ENSG00000150753","ENSG00000111057","ENSG00000148834","ENSG00000144381","ENSG00000155368","ENSG00000173207","ENSG00000010278","ENSG00000075131","ENSG00000129757","ENSG00000105968","ENSG00000188486","ENSG00000160752","ENSG00000110092","ENSG00000147883"))
DoHeatmap(test5, features = test5_final$mouse_ensembl)

###Final_gene to human ensembl; 736
test5_final_df <- test5_final%>%as.data.frame()
colnames(test5_final_df) <- "epi_mouse"
test5_final<-merge(test5_final_df, ensembl_mouse_human, by.x = "epi_mouse", by.y = "mouse_ensembl")
table(duplicated(test5_final$human_ensembl))
test5_final<- test5_final[!duplicated(test5_final[2]),]
head(test5_final)
names(test5_final)[1]<-"mouse_ensembl"
dim(test5_final) #736gene


test5_final_markers <- test5_markers_rm_nonepi %>% as.data.frame()
colnames(test5_final_markers) <- "epi_mouse"
test5_final_markers<-merge(test5_final_markers, ensembl_mouse_human, by.x = "epi_mouse", by.y = "mouse_ensembl")
table(duplicated(test5_final_markers$human_ensembl))
test5_final_markers<- test5_final_markers[!duplicated(test5_final_markers[2]),]
head(test5_final_markers)
names(test5_final_markers)[1]<-"mouse_ensembl"
dim(test5_final_markers) #736gene


##single cell data processing -> with 700 genes (expression data)   #gene set; test5_final -> test5_final_markers(539genes)
mouse_epi_traj<-test5@assays$RNA@data %>%as.data.frame() %>%rownames_to_column() # cell; 6476 / gene; 11528
names(mouse_epi_traj)[1] <- c("ensembl")
mouse_epi_700 <- merge(mouse_epi_traj, test5_final, by.x = "ensembl", by.y = "mouse_ensembl") #734
mouse_epi_700 <- mouse_epi_700[,-1]
head(mouse_epi_700[,1:10])
dim(mouse_epi_700)

##mouse traj knn pooling (sparcity 보정)
mouse_epi_700 <- mouse_epi_700%>%column_to_rownames(var = "human_ensembl")
mouse_epi_700_knn <- knn(mouse_epi_700, 5)

##Plot
mouse_epi_traj<-test5
group1 <- c('Fetal lung')
group2 <- c('AT1')
group3 <- c('AT2')
group4 <- c('Ciliated')
group5 <- c('Clara')
group6 <- c('Basal')
group7 <- c('PNEC')
cell1 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group1])
cell2 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group2])
cell3 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group3])
cell4 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group4])
cell5 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group5])
cell6 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group6])
cell7 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group7])
bbknn_dt <- mouse_epi_traj@reductions$bbknn@cell.embeddings %>% as.data.frame() %>% rownames_to_column('cell_id')
bbknn_dt <- bbknn_dt %>% mutate(group = ifelse(cell_id %in% cell1, paste(group1,collapse=','), 
                                               ifelse(cell_id %in% cell2, paste(group2, collapse=','),
                                                      ifelse(cell_id %in% cell3, paste(group3, collapse=','),
                                                             ifelse(cell_id %in% cell4, paste(group4, collapse=','),
                                                                    ifelse(cell_id %in% cell5, paste(group5, collapse=','),
                                                                           ifelse(cell_id %in% cell6, paste(group6, collapse=','),
                                                                                  "PNEC")))))))
ggplot(bbknn_dt, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), size = 0.3, alpha = 0.7)+
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#87d4f5","#7dab78", "#f57c16",
                                 "#F0E442", "#cf99f7","#bf046a","#a8850f","#0e5e05","#9e5416","#cf99f7")) +
  theme_classic()

#4.Human normal과 비교 (# 부분 조절해서 AT1~ciliated까지 가능)-----------------------------------------------------------------
mouse_epi_700_knn <- mouse_epi_700_knn %>% as.data.frame() %>% rownames_to_column(var = "human_ensembl")
hu_AT1 <- mouse_epi_700_knn$human_ensembl %>% as.data.frame()
hu_AT1$id  <- 1:nrow(hu_AT1)
colnames(hu_AT1) <- c("gene","id")
hu_AT1 <- merge(hu_AT1, donor_cili, by.x = "gene", by.y = "gene", all.x = TRUE, all.y = FALSE) #여기 바꾸면서 조절
sum(is.na(hu_AT1))
hu_AT1[is.na(hu_AT1)] <- 0
sum(is.na(hu_AT1))
hu_AT1 <- hu_AT1[order(hu_AT1$id),]
rownames(hu_AT1) <-NULL
hu_AT1 <- hu_AT1 %>% as.data.frame() %>% column_to_rownames(var = "gene")
hu_AT1 <- hu_AT1[,-1]
head(rownames(hu_AT1))
head(colnames(hu_AT1))
dim(hu_AT1)
mouse_epi_700_knn <- mouse_epi_700_knn%>%column_to_rownames(var = "human_ensembl")

###human lung cluster마다 mouse traj와 cor 계산
m_h_AT1 <- cbind(mouse_epi_700_knn, hu_AT1)
cor_m_h_AT1 <- cor(m_h_AT1, method = "pearson") #method = "spearman"/ pearson(default)
cor_m_h_AT1 <- cor_m_h_AT1 %>% as.data.frame()
cor_m_h_AT1_sort <- cor_m_h_AT1[1:6476,6477:nrow(cor_m_h_AT1)]
#cor_m_h_AT1_max <- apply(cor_m_h_AT1_sort,2, max) #correlation cutoff 설정 
cor_m_h_AT1_sort <- cor_m_h_AT1_sort[,which(cor_m_h_AT1_max > 0.45)]
cor_m_h_AT1_sort_rank <- apply(cor_m_h_AT1_sort,2,function(x){order(-x)})[1:7,]%>%t()
cor_m_h_AT1_sort_rank <- cor_m_h_AT1_sort_rank %>% as.data.frame()

###plot(highly correlated human cell positioning on mouse traj)
a<-Embeddings(test5,reduction = "bbknn") %>% as.data.frame()
cor_m_h_AT1_sort_rank <- cor_m_h_AT1_sort_rank[,1:3] %>% t()
matched <- data.frame()

for(n in 1:ncol(cor_m_h_AT1_sort_rank)){
  tmp_matched <- a[cor_m_h_AT1_sort_rank[,n],]
  tmp_matched$group <- "h_Ciliated" #여기도 cluster 바뀔때마다 조절
  matched<- rbind(matched, tmp_matched)
}
matched <- matched %>% rownames_to_column(var = "cell_id")

bbknn_dt_plot <- rbind(bbknn_dt, matched)
ggplot(bbknn_dt_plot, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), size = .5, alpha = 1)+
  #scale_shape_manual(values = c( 16, 16, 16,16, 16, 2, 17))+
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#87d4f5","#7dab78", "#f57c16",
                                 "#F0E442", "#0018c9","#cf99f7","#bf046a","#a8850f","#0e5e05","#9e5416","#cf99f7"))+
  theme_classic()

#5.cancer와 비교--------------------------------------------------------------------------------------------------------------------------------
###LUAD_bulk_RNA/ LUSC_bulk/ SCLC_bulk (expression data_FPKM) ;  dir, file_list와 plot에 표기되는 부분만 바꿔주면 ok
dir <- ("./projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUSC_551_FPKM//")
file_list<-list.files("/home/users/yunah1029/projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUAD_551_FPKM",, recursive=TRUE, pattern="*.FPKM.txt.gz$")
file_list<-list.files("/home/users/yunah1029/projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUSC_551_FPKM",, recursive=TRUE, pattern="*.FPKM.txt.gz$")
####sclc의 preprocessing과정은 sclc_plot.R 참고

mouse_epi_700_knn <- mouse_epi_700_knn %>% as.data.frame() %>% rownames_to_column(var = "human_ensembl")
tcga <- mouse_epi_700_knn$human_ensembl %>% as.data.frame()
tcga$id  <- 1:nrow(tcga)
colnames(tcga) <- "gene"
dim(tcga)

for (file in file_list) {
  temp <- read.csv(paste0(dir, file), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = FALSE)
  colnames(temp) <- c("gene", strsplit(file, "/")[[1]][1])
  temp$gene <- sapply(strsplit(temp$gene, "[.]"), function(x) x[1])
  tcga <- merge(tcga, temp, by.x = "gene", by.y = "gene", all.x = TRUE, all.y = FALSE)
}
sum(is.na(tcga))
tcga[is.na(tcga)] <- 0
sum(is.na(tcga))

###
tcga<-tcga[order(tcga$"NA"), ] #mouse_epi_600의 gene 순서와 동일하게 정렬 -> 나중에 cbind
rownames(tcga) <- NULL
tcga<- tcga %>% as.data.frame() %>% column_to_rownames(var = "gene")
tcga<-tcga[,-1]
head(rownames(tcga))
head(colnames(tcga))
head(tcga[,1:10])
dim(tcga)

####purity cutoff = 0.7 -> sample 156개
library(readr)
luad_purity_0.7 <- read_csv("lusc_sum_purity.csv") #luad는 luad_purity_0.7.csv
head(luad_purity_0.7)
tcga_tmp <- tcga[,colnames(tcga) %in% luad_purity_0.7$File_id]
dim(tcga_tmp)
head(tcga_tmp[,1:10])
##tcga<-tcga_tmp
##rm(tcga_tmp)
mouse_epi_700_knn <- mouse_epi_700_knn%>%column_to_rownames(var = "human_ensembl")

###tcga sample 별 mouse_epi_700과 cor -> 가장 가까운 cell 찾기. (cutoff 0.4 -> sample 80개)
single_bulk_luad <- cbind(mouse_epi_700_knn, tcga_tmp)
cor_single_bulk_luad <- cor(single_bulk_luad, method = "pearson") #method = "spearman"/ pearson(default)
colnames(cor_single_bulk_luad)[6476]
cor_single_bulk_luad_sort <- cor_single_bulk_luad[1:6476,6477:nrow(cor_single_bulk_luad)]
cor_max_luad <- apply(cor_single_bulk_luad_sort,2, max) #correlation cutoff 설정 -> 0.4
cor_single_bulk_luad_sort <- cor_single_bulk_luad_sort[,which(cor_max_luad > 0.4)]
cor_single_bulk_luad_sort_rank <- apply(cor_single_bulk_luad_sort,2,function(x){order(-x)})[1:7,]%>%t()

# original gene set으로 correlation 볼 때 각 tcga sample 별 Median
original_matched <- data.frame()
original_avg <- data.frame()

for (x in 1:nrow(cor_single_bulk_luad_sort_rank)) {
  tmp_original_matched <- a[cor_single_bulk_luad_sort_rank[x,],]
  tmp_original_matched$group <- "LUSC:Origianl genes"
  tmp_original_matched$cancer_sample <- rownames(cor_single_bulk_luad_sort_rank)[x]
  tmp_original_matched <- tmp_original_matched %>% rownames_to_column(var = "cell_id") #추가
  tmp_original_matched$chull <- x
  
  tmp_original_avg <- apply(tmp_original_matched[,2:3],2,median) %>% t() %>% as.data.frame()
  tmp_original_avg$group <- "LUSC:Origianl genes_Median"
  tmp_original_avg$cell_id <- rownames(cor_single_bulk_luad_sort_rank)[x]
  tmp_original_avg$cancer_sample <- rownames(cor_single_bulk_luad_sort_rank)[x]
  tmp_original_avg <- tmp_original_avg[,c(4,1,2,3,5)]
  tmp_original_avg$chull_group <- x
  
  original_matched <- rbind(original_matched, tmp_original_matched)
  original_avg <- rbind(original_avg, tmp_original_avg)
}

dim(original_avg)
head(original_avg)

##random sampling + convex hull 그리기
single_bulk_luad_cut<-cbind(mouse_epi_700_knn,single_bulk_luad[,colnames(single_bulk_luad) %in% colnames(cor_single_bulk_luad_sort)])

i = 1
matched_final <- data.frame()
avg_final <- data.frame()

###cutoff 하려면 single bulk luad 조절 ㄱ
#single_bulk_luad_cutoff <- cbind(mouse_epi_700_knn, tcga[,which(cor_max_luad > 0.45)])

repeat{
  random_700<-sample(rownames(single_bulk_luad_cut), 734, replace = TRUE)
  single_bulk_random <- single_bulk_luad_cut %>% filter(rownames(single_bulk_luad_cut) %in% random_700) 
  cor_single_bulk_random <- cor(single_bulk_random, method = "pearson") %>% as.data.frame() #method = "spearman"/ pearson(default)
  cor_single_bulk_sort_random <- cor_single_bulk_random[1:6476,6477:nrow(cor_single_bulk_random)]
  cor_single_bulk_sort_rank_random <- apply(cor_single_bulk_sort_random,2,function(x){order(-x)})[1:7,]%>%t() #1등부터 7등까지
  cor_single_bulk_sort_rank_random <- cor_single_bulk_sort_rank_random %>% as.data.frame()
  cor_single_bulk_sort_rank_random <- cor_single_bulk_sort_rank_random %>% t()
  
  matched <- data.frame()
  avg <- data.frame()
  
  for (n in 1:ncol(cor_single_bulk_sort_rank_random)){
    tmp_matched <- a[cor_single_bulk_sort_rank_random[,n],]
    tmp_matched$group <- "LUSC:Highly correlated"
    tmp_matched$cancer_sample <-colnames(cor_single_bulk_sort_rank_random)[n]
    tmp_matched <- tmp_matched %>% rownames_to_column(var = "cell_id")
    tmp_matched$chull_group <- paste0(i,"_",n)
    
    tmp_avg <- apply(tmp_matched[,2:3],2,median) %>% t() %>% as.data.frame()
    tmp_avg$group <- "LUSC:Median"
    tmp_avg$cell_id <- colnames(cor_single_bulk_sort_rank_random)[n]
    tmp_avg$cancer_sample <-colnames(cor_single_bulk_sort_rank_random)[n]
    tmp_avg <- tmp_avg[,c(4,1,2,3,5)]
    tmp_avg$chull_group <- "median"
    
    tmp_matched_dist <- tmp_matched
    for ( n in 1:nrow(tmp_matched_dist)) {
      tmp_matched_dist[n,7] <- sqrt((tmp_matched_dist[n,2] - tmp_avg$BBKNN_1)^2 + (tmp_matched_dist[n,3] - tmp_avg$BBKNN_2)^2)
    }
    
    tmp_matched_final <- tmp_matched[tmp_matched$cell_id %in% tmp_matched_dist[which(rank(tmp_matched_dist$V7) < nrow(tmp_matched_dist)*0.9),1],]
    
    matched <- rbind(matched, tmp_matched_final) #하나의 gene set에 대해
    avg <- rbind(avg, tmp_avg)
  }
  
  matched_final <- rbind(matched_final, matched)
  avg_final <- rbind(avg_final, avg)
  
  if(i == 100) break
  i = i +1
}

bbknn_dt_tmp <- bbknn_dt
bbknn_dt_tmp$cancer_sample <- "."
bbknn_dt_tmp$chull_group <- "."

# Median spot으로 convex hull/ sampled gene set Med <-> original gene set Med => convex hull 그릴 준비. 
# avg_final = 100번 random sampling 했을 때 각 gene set에 대한 median(top 7) // original_avg = original gene set에서 cor을 계산했을 때 top 7에 대한 Median
med_dist_final <- data.frame()
for (y in 1:nrow(original_avg)) {
  original_avg_dist <- avg_final %>% filter(cell_id %in% original_avg[y,1])
  original_avg_dist$dist <- sqrt((original_avg[y,2] - original_avg_dist$BBKNN_1)^2 + (original_avg[y,3] - original_avg_dist$BBKNN_2)^2)
  med_dist_tmp <- original_avg_dist[which(rank(original_avg_dist$dist) < nrow(original_avg_dist)*0.9),]
  med_dist_final <- rbind(med_dist_final,med_dist_tmp)
}

med_dist_final <- med_dist_final[,-7]
head(med_dist_final)
hull2 <- med_dist_final %>% group_by(cancer_sample) %>% slice(chull(BBKNN_1, BBKNN_2)) # 

pdf("perturbation(lusc)_sampled_genes3_v3_chull_medians_0.7_0.4.pdf",8,5)
for (k in 1:length(unique(original_avg$cancer_sample))) {
  bbknn_dt_tmp2 <- rbind(bbknn_dt_tmp,med_dist_final%>% filter(cancer_sample == unique(original_avg$cancer_sample)[k]),original_avg%>%filter(cancer_sample == unique(original_avg$cancer_sample)[k])) #med_dist_final%>% filter(cancer_sample == unique(original_avg$cancer_sample)[k]
  hull3<- hull2 %>% filter(cancer_sample == unique(original_avg$cancer_sample)[k])
  tcga_sample_name <- unique(original_avg$cancer_sample)[k]
  print(
    bbknn_dt_tmp2 %>%
      ggplot(aes(x = BBKNN_1, y = BBKNN_2)) +
      ggtitle(tcga_sample_name)+
      geom_point(aes(color=group), size = 0.8, alpha = 0.7)+
      scale_colour_manual(values = c("#CC79A7", "#E69F00", "#87d4f5","#7dab78", "#f57c16",
                                     "#F0E442", "#4674e8","#000000","#cf99f7","#a8850f","#0e5e05","#9e5416","#cf99f7"))+ # 뒤에서 두번째 "#0066CC" 
      geom_polygon(data = hull3, size = 0.5, color = "#0066CC", alpha = 0.5, fill = NA) +
      theme_classic()
  )
}
dev.off()

#6.origin 별 cancer DEG--------------------------------------------------------------------------------------------------------------------------------
library( "DESeq2" )
library(ggplot2)
##tcga read count data with all genes(60483)
dir <- ("./projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga(HTseq)_LUAD_551/")
file_list<-list.files("/home/users/yunah1029/projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga(HTseq)_LUAD_551",, recursive=TRUE, pattern="*.htseq.counts$")
head(file_list)
tcga_count <- read.csv(paste0(dir, file_list[1]), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = FALSE)
colnames(tcga_count) <- c("gene", strsplit(file, "/")[[1]][1])
tcga_count$gene <- sapply(strsplit(tcga_count$gene, "[.]"), function(x) x[1])
for (file in file_list) {
  temp <- read.csv(paste0(dir, file), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = FALSE)
  colnames(temp) <- c("gene", strsplit(file, "/")[[1]][1])
  temp$gene <- sapply(strsplit(temp$gene, "[.]"), function(x) x[1])
  tcga_count <- merge(tcga_count, temp, by.x = "gene", by.y = "gene", all.x = TRUE, all.y = FALSE)
}
sum(is.na(tcga_count))

head(tcga_count[,1:10])
tcga_count <- tcga_count[-c(1:5),]
head(tcga_count[,1:10])
which(colnames(tcga_count) == "ffdf38f8-864c-4354-bcb3-83617bbbcff5.y")
names(tcga_count)[553] <- c("ffdf38f8-864c-4354-bcb3-83617bbbcff5")
rownames(tcga_count) <- NULL
tcga_count<- tcga_count %>% as.data.frame() %>% column_to_rownames(var = "gene")
tcga_count<-tcga_count[,-1]
head(rownames(tcga_count))
head(colnames(tcga_count))

## Meta data
library(readxl)
luad_summary<- read_csv("luad_sum_purity.csv")
head(luad_summary)
dim(luad_summary)

##Before DESeq2-------------------------------------------------------------------------------------------
###classifing tcga samples 
tcga_filter <- tcga_count %>% select(which(colnames(tcga_count) %in% luad_summary$File_id_count))
tcga_filter_gene <- merge(tcga_filter%>%rownames_to_column(var = "gene"), human_ensmebl_genename, by.x = "gene", by.y = "Gene stable ID", all.x = TRUE, all.y = FALSE)
tcga_filter_gene <- tcga_filter_gene[,-1]
tcga_filter_gene<-na.omit(tcga_filter_gene)
tcga_filter_gene<- tcga_filter_gene[!duplicated(tcga_filter_gene[ncol(tcga_filter_gene)]),]
rownames(tcga_filter_gene) <- NULL
tcga_filter_gene <- tcga_filter_gene %>% column_to_rownames(var = "Gene name")
dim(tcga_filter_gene)

###tcga_filtere with 734 genes(->539 gene ->534(ensembl to gene name))
tcga_filter_700 <- merge(mouse_epi_700_knn%>%rownames_to_column(var = "gene")%>%select(gene),tcga_filter%>%rownames_to_column(var = "gene"), by.x = "gene", by.y = "gene", all.x = "TRUE", all.y = "FALSE")
sum(is.na(tcga_filter_700))
tcga_filter_700[is.na(tcga_filter_700)] <- 0
sum(is.na(tcga_filter_700))
tcga_filter_700 <- tcga_filter_700 %>% column_to_rownames(var = "gene")
dim(tcga_filter_700)
head(tcga_filter_700[,1:10])

tcga_filter_700_gene <- merge(tcga_filter_700%>%rownames_to_column(var = "gene"), human_ensmebl_genename, by.x = "gene", by.y = "Gene stable ID", all.x = TRUE, all.y = FALSE)
tcga_filter_700_gene <- tcga_filter_700_gene[,-1]
tcga_filter_700_gene[is.na(tcga_filter_700_gene)] <- 0
tcga_filter_700_gene<- tcga_filter_700_gene[!duplicated(tcga_filter_700_gene[ncol(tcga_filter_700_gene)]),]
rownames(tcga_filter_700_gene) <- NULL
tcga_filter_700_gene <- tcga_filter_700_gene %>% column_to_rownames(var = "Gene name")
dim(tcga_filter_700_gene)
head(tcga_filter_700_gene[,1:10])

###at2/fetal type tcga sample sorting
at2_fetal_filename <-luad_summary%>%filter(Type %in% c("AT2","Fetal")) #Basal
at2_fetal_filename <- at2_fetal_filename$File_id_count
at2_fetal_700<-tcga_filter_gene[,at2_fetal_filename]
dim(at2_fetal_700)

luad_summary_filter <- luad_summary %>% filter(luad_summary$File_id_count %in% colnames(at2_fetal_700))
dim(luad_summary_filter)

##DESeq2
tcga_deseq2 <- DESeqDataSetFromMatrix(
  countData = at2_fetal_700,
  colData = luad_summary_filter,
  design = ~Type
)
dim(tcga_deseq2)

keep <- rowSums(counts(tcga_deseq2)) >= 10
tcga_deseq2 <- tcga_deseq2[keep,]
dim(tcga_deseq2)

tcga_deseq2 <- DESeq(tcga_deseq2)
res <- results(tcga_deseq2, contrast = c("Type", "AT2", "Fetal")) #비교할 대상 바꿔줄 수 있음

#
metadata(res)$alpha
metadata(res)$filterThreshold
plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
resNoFilt <- results(tcga_deseq2, independentFiltering=FALSE)
addmargins(table(filtering=(res$padj < .1),
                 noFiltering=(resNoFilt$padj < .1)))
#
mcols(res, use.names = TRUE)
head(results(tcga_deseq2, tidy=TRUE))
sum( res$pvalue < 0.01, na.rm=TRUE)
table(is.na(res$pvalue))
sum( res$padj < 0.05, na.rm=TRUE ) #0.1
resSig <- res[ which(res$padj < 0.05 ), ] #0.1
head( resSig[ order( resSig$log2FoldChange ), ] ) #significant genes with the strongest down-regulation
tail( resSig[ order( resSig$log2FoldChange ), ] ) #significant genes with the strongest up-regulation

basal_fetal_deg_pos <- resSig %>% subset(log2FoldChange > 2)
basal_fetal_deg_neg <- resSig %>% subset(log2FoldChange < -2)
rownames(basal_fetal_deg_pos) %>% as.data.frame %>% write_csv("basal_fetal_deg_pos.csv")
rownames(basal_fetal_deg_neg) %>% as.data.frame %>% write_csv("basal_fetal_deg_neg.csv")

#par(mfrow=c(1,1))
plotCounts(tcga_deseq2, gene="SFTPC",intgroup="Type") #
plotCounts(tcga_deseq2, gene="TMEM213", intgroup="Type")
plotCounts(tcga_deseq2, gene="LNCAROD", intgroup="Type")
plotCounts(tcga_deseq2, gene="TAC1", intgroup="Type") #AT2
plotCounts(tcga_deseq2, gene="SLCO1B3", intgroup="Type")
plotCounts(tcga_deseq2, gene="NOL4", intgroup="Type") #fetal
plotCounts(tcga_deseq2, gene="ECT2", intgroup="Type") #fetal
plotCounts(tcga_deseq2, gene="EGFR", intgroup="Type") #fetal
plotCounts(tcga_deseq2, gene="TP53", intgroup="Type") #fetal

#dev.off()

##Volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20,  , main="Volcano plot", xlim=c(-25,25)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=1.0)
abline(v=2, col="black", lty=4, lwd=1.0)
abline(h=-log10(max(res$padj[res$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=1.0)
with(subset(res, padj<.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(res), res$padj<0.05 & abs(res$log2FoldChange)>2), cex=0.8, pos=3))

##Diagnostic plot
plotMA( res, ylim = c(-5, 5))
plotDispEsts( tcga_deseq2, ylim = c(1e-2, 1e1) )
hist( res$pvalue, breaks=20, col="grey" )

##rlog(regularized-logarithm) transform => Sample distances
rld <- vst(tcga_deseq2, blind = FALSE)
plotPCA(rld, intgroup=c("Type"))
head(assay(rld))
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Type,
                                     rld$`Sample ID`, sep="_" )
colnames(sampleDistMatrix) <- NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

####Gene clustering
library( "genefilter")
library("RColorBrewer")
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ),100)
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( AT2="gray", Fetal="orange" )[
             colData(rld)$Type ] )

##tcga read count data with all genes(LUSC)---------------------------------------------------------------------------------------------------------------------------------------------------------------
dir <- ("./projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUSC_551_count/")
file_list<-list.files("/home/users/yunah1029/projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUSC_551_count",, recursive=TRUE, pattern="*.htseq.counts.gz$")
head(file_list)
length(file_list)
tcga_count_lusc <- read.csv(paste0(dir, file_list[1]), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = FALSE)
colnames(tcga_count_lusc) <- c("gene", strsplit(file, "/")[[1]][1])
tcga_count_lusc$gene <- sapply(strsplit(tcga_count_lusc$gene, "[.]"), function(x) x[1])
for (file in file_list) {
  temp <- read.csv(paste0(dir, file), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = FALSE)
  colnames(temp) <- c("gene", strsplit(file, "/")[[1]][1])
  temp$gene <- sapply(strsplit(temp$gene, "[.]"), function(x) x[1])
  tcga_count_lusc <- merge(tcga_count_lusc, temp, by.x = "gene", by.y = "gene", all.x = TRUE, all.y = FALSE)
}
sum(is.na(tcga_count_lusc))

head(tcga_count_lusc[,1:10])
tcga_count_lusc <- tcga_count_lusc[-c(1:5),]

rownames(tcga_count_lusc) <- NULL
tcga_count_lusc<- tcga_count_lusc %>% as.data.frame() %>% column_to_rownames(var = "gene")
tcga_count_lusc<-tcga_count_lusc[,-1]
head(rownames(tcga_count_lusc))
head(colnames(tcga_count_lusc))
dim(tcga_count_lusc)

## Meta data
library(readxl)
lusc_summary<- read_csv("lusc_sum_purity.csv")
head(lusc_summary)
dim(lusc_summary)

###classifing tcga samples 
tcga_filter_lusc <- tcga_count_lusc [,lusc_summary$File_id_count]
tcga_filter_gene_lusc <- merge(tcga_filter_lusc%>%rownames_to_column(var = "gene"), human_ensmebl_genename, by.x = "gene", by.y = "Gene stable ID", all.x = TRUE, all.y = FALSE)
tcga_filter_gene_lusc <- tcga_filter_gene_lusc[,-1]
tcga_filter_gene_lusc<-na.omit(tcga_filter_gene_lusc)
tcga_filter_gene_lusc<- tcga_filter_gene_lusc[!duplicated(tcga_filter_gene_lusc[ncol(tcga_filter_gene_lusc)]),]
rownames(tcga_filter_gene_lusc) <- NULL
tcga_filter_gene_lusc <- tcga_filter_gene_lusc %>% column_to_rownames(var = "Gene name")
dim(tcga_filter_gene_lusc)

##DESeq2
tcga_deseq2_lusc <- DESeqDataSetFromMatrix(
  countData = tcga_filter_gene_lusc,
  colData = lusc_summary,
  design = ~Type
)
dim(tcga_deseq2_lusc)

keep_lusc <- rowSums(counts(tcga_deseq2_lusc)) >= 10
tcga_deseq2_lusc <- tcga_deseq2_lusc[keep_lusc,]
dim(tcga_deseq2_lusc)

tcga_deseq2_lusc <- DESeq(tcga_deseq2_lusc)
res_lusc <- results(tcga_deseq2_lusc, contrast = c("Type", "Basal", "Fetal")) #비교할 대상 바꿔줄 수 있음

mcols(res_lusc, use.names = TRUE)
head(results(tcga_deseq2_lusc, tidy=TRUE))
sum( res_lusc$pvalue < 0.01, na.rm=TRUE)
table(is.na(res_lusc$pvalue))
sum( res_lusc$padj < 0.05, na.rm=TRUE ) #0.1
resSig_lusc <- res_lusc[ which(res_lusc$padj < 0.05 ), ] #0.1
head( resSig_lusc[ order( resSig_lusc$log2FoldChange ), ] ) #significant genes with the strongest down-regulation
tail( resSig_lusc[ order( resSig_lusc$log2FoldChange ), ] ) #significant genes with the strongest up-regulation

basal_fetal_deg_pos_lusc <- resSig_lusc %>% subset(log2FoldChange > 2)
basal_fetal_deg_neg_lusc <- resSig_lusc %>% subset(log2FoldChange < -2)
rownames(basal_fetal_deg_pos_lusc) %>% as.data.frame %>% write_csv("basal_fetal_deg_pos_lusc")
rownames(basal_fetal_deg_neg_lusc) %>% as.data.frame %>% write_csv("basal_fetal_deg_neg_lusc")

##Volcano plot
with(res_lusc, plot(log2FoldChange, -log10(padj), pch=20,  , main="Volcano plot", xlim=c(-13,13)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res_lusc, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res_lusc, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=1.0)
abline(v=2, col="black", lty=4, lwd=1.0)
abline(h=-log10(max(res_lusc$padj[res_lusc$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=1.0)
with(subset(res_lusc, padj<.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(res_lusc), res_lusc$padj<0.05 & abs(res_lusc$log2FoldChange)>2), cex=0.8, pos=3))

#par(mfrow=c(1,1))
plotCounts(tcga_deseq2_lusc, gene="PIK3CA",intgroup="Type") #
plotCounts(tcga_deseq2_lusc, gene="NFE2L2",intgroup="Type") #
plotCounts(tcga_deseq2_lusc, gene="PTEN",intgroup="Type") #
#KEAP1, NFE2L2 /TP53, PI3K, RB1 and NFE2L2/KEAP1

##Diagnostic plot
plotMA( res_lusc, ylim = c(-5, 5))
plotDispEsts( tcga_deseq2_lusc, ylim = c(1e-2, 1e1) )
hist( res_lusc$pvalue, breaks=20, col="grey" )

##rlog(regularized-logarithm) transform => Sample distances
rld_lusc <- vst(tcga_deseq2_lusc, blind = FALSE)
plotPCA(rld_lusc, intgroup=c("Type"))
head(assay(rld_lusc))
sampleDists_lusc <- dist(t(assay(rld_lusc)))
sampleDists_lusc
sampleDistMatrix_lusc <- as.matrix( sampleDists_lusc )
rownames(sampleDistMatrix_lusc) <- paste( rld_lusc$Type,
                                          rld_lusc$`Sample ID`, sep="_" )
colnames(sampleDistMatrix_lusc) <- NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix_lusc, trace="none", col=colours)

##gene clustering
library( "genefilter")
library("RColorBrewer")
topVarGenes <- head( order( rowVars( assay(rld_lusc) ), decreasing=TRUE ),100)
heatmap.2( assay(rld_lusc)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( AT2="gray", Fetal="orange" )[
             colData(rld_lusc)$Type ] )

## clinic/exposure data 
library(readr)
library(dplyr)
library(tidyverse)

luad_meta<- read_csv("luad_sum_purity.csv")
head(luad_meta)
luad_meta <- luad_meta[,1:18]
luad_meta <- luad_meta[1:32,]
head(luad_meta)

lusc_meta <- read_csv("lusc_sum_purity.csv")
head(lusc_meta)

ggplot(lusc_meta%>%select(Type,sex),aes(x=factor(Type),fill=factor(sex)))+
  geom_bar(position = "fill")
geom_text(aes(label=..count..),stat="count",position=position_stack(0.5))
