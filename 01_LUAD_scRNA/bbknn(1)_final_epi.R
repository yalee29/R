#resolution 높여서 다시 보기. 
E_seurat1_ad <- FindClusters(E_seurat1_ad, resolution = 2.5)
E_seurat1_ad_bb<- E_seurat1_ad %>% RunBBKNN(dims.use = 1:40,
                                            neighbors_within_batch =  5,
                                            trim=45,
                                            batch.key = "batch",
                                            python.path = "/home/users/yunah1029/anaconda3/bin/python") 
DimPlot(E_seurat1_ad_bb, reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; dims=40, k=5, trim=45, all")
DimPlot(E_seurat1_ad_bb, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=40, k=5, trim=45, all")


#Fetal lung 세분화
DimPlot(E_seurat1_ad_bb[,(E_seurat1_ad_bb@active.ident %in% c('48','19','27','33','31','45','42','43','41','21') & E_seurat1_ad_bb$batch %in% c("E8.25","E14.5.1","E9.5to11.5","LungMap.E18.5","E14.5.2","E16.5","E18.5"))], reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=40, k=5, trim=45, fetal")
DimPlot(E_seurat1_ad_bb[,(E_seurat1_ad_bb@active.ident %in% c('48','19','27','33','31','45','42','43','41','21') & E_seurat1_ad_bb$batch %in% c("E8.25","E14.5.1","E9.5to11.5","LungMap.E18.5","E14.5.2","E16.5","E18.5"))], reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; dims=40, k=5, trim=45, fetal")

marker.48 <- FindMarkers(E_seurat1_ad_bb, ident.1 = 48)
marker.19 <- FindMarkers(E_seurat1_ad_bb, ident.1 = 19)
marker.27 <- FindMarkers(E_seurat1_ad_bb, ident.1 = 27)
marker.33 <- FindMarkers(E_seurat1_ad_bb, ident.1 = 33)
marker.31 <- FindMarkers(E_seurat1_ad_bb, ident.1 = 31)
marker.45 <- FindMarkers(E_seurat1_ad_bb, ident.1 = 45)
marker.42 <- FindMarkers(E_seurat1_ad_bb, ident.1 = 42)
marker.43 <- FindMarkers(E_seurat1_ad_bb, ident.1 = 43)
marker.41 <- FindMarkers(E_seurat1_ad_bb, ident.1 = 41)
marker.21 <- FindMarkers(E_seurat1_ad_bb, ident.1 = 21)
markers.all <- FindAllMarkers(E_seurat1_ad_bb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

FeaturePlot(E_seurat1_ad_bb[,(E_seurat1_ad_bb@active.ident %in% c('48','19','27','33','31','45','42','43','41','21') & E_seurat1_ad_bb$batch %in% c("E8.25","E14.5.1","E9.5to11.5","LungMap.E18.5","E14.5.2","E16.5","E18.5"))], c("ENSMUSG00000016494"), reduction = "bbknn")
FeaturePlot(E_seurat1_ad_bb, c("ENSMUSG00000023868"), reduction = "bbknn")

#
DimPlot(E_seurat1_ad_bb[,E_seurat1_ad_bb@active.ident == c('31','45','43','23','25','0','4','35','5','8','44','32','11')], reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; dims=40, k=5, trim=45, all")
DimPlot(E_seurat1_ad_bb[,E_seurat1_ad_bb@active.ident == c('31','45','43','23','25','0','4','35','5','8','44','32','11')|E_seurat1_ad_bb$batch %in% c("E8.25","E14.5.2","E16.5","E18.5")], reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; dims=40, k=5, trim=45, all")

# 나중에! (annotation)
new.cluster.ids <- c("Epithelial", "Mesenchymal", "Mesenchymal", "Epithelial", "Immune", "Immune", "Immune", "Immune", "Immune","9", "Epithelial", "Epithelial", "Mesenchymal","Immune", "Mesenchymal","Immune","Mesenchymal","Mesenchymal","18","19","Mesenchymal","21","Epithelial","23","Epithelial","Immune","Immune","Mesenchymal") 
names(new.cluster.ids) <- levels(E_seurat1_ad_bb)
E_seurat1_ad_bb <- RenameIdents(E_seurat1_ad_bb, new.cluster.ids)
DimPlot(E_seurat1_ad_bb, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=40, k=5, trim=45, all")

#bbknn1_Epithelial cells
E_seurat1_ad_epi <- E_seurat1_ad_bb[,(E_seurat1_ad_bb@active.ident %in% c('31','45','23','25','0','4','35','5','8','44','32','11')|E_seurat1_ad_bb$batch %in% c("E8.25","E14.5.2","E16.5","E18.5"))]
DimPlot(E_seurat1_ad_epi, reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; dims=40, k=5, trim=45, epi")
unique(E_seurat1_ad_epi$batch)
ElbowPlot(E_seurat1_ad_epi)
E_seurat1_ad_epi <- JackStraw(E_seurat1_ad_epi,dims = 30)
E_seurat1_ad_epi <- ScoreJackStraw(E_seurat1_ad_epi, dims = 1:30)
JackStrawPlot(E_seurat1_ad_epi, dims = 1:30)

E_seurat1_ad_epi <- FindNeighbors(E_seurat1_ad_epi, dims = 1:20)
E_seurat1_ad_epi <- FindClusters(E_seurat1_ad_epi, resolution = 1.5)

pdf("bbknn1_epi2_3(PND1.2).pdf")
for(k in c(3:6)){
  for(trim in c(20,25,30,35,40,45)){
    print(
      E_seurat1_ad_epi %>% RunBBKNN(dims.use = 1:20,
                                    neighbors_within_batch = k,
                                    trim = trim,
                                    batch.key = "batch",
                                    python.path = "/home/users/yunah1029/anaconda3/bin/python") %>%  
        DimPlot(reduction = "bbknn", group.by = "batch") %>% 
        LabelClusters(id="batch") + ggtitle(paste0("bbknn1; neighbors=",k,",trim=",trim,",epi"))
    )
  }
}
dev.off()

E_seurat1_ad_epi_bb<- E_seurat1_ad_epi %>% RunBBKNN(dims.use = 1:20,
                                                    neighbors_within_batch =5,
                                                    trim=35,
                                                    batch.key = "batch",
                                                    python.path = "/home/users/yunah1029/anaconda3/bin/python") 
DimPlot(E_seurat1_ad_epi_bb, reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; dims=20, k=5, trim=35, epi")
DimPlot(E_seurat1_ad_epi_bb, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=20, k=5, trim=35, epi")

#Find clusters
FeaturePlot(E_seurat1_ad_epi_bb, c("ENSMUSG00000042784"), reduction = "bbknn")
FeaturePlot(E_seurat1_ad_epi_bb, c("ENSMUSG00000079465"), reduction = "bbknn") #AT1 "ENSMUSG00000059325","ENSMUSG00000028583","ENSMUSG00000023039","ENSMUSG00000023951"
FeaturePlot(E_seurat1_ad_epi_bb, c("ENSMUSG00000041247"), reduction = "bbknn") #"ENSMUSG00000029375","ENSMUSG00000037071","ENSMUSG00000042784"), reduction = "bbknn") #AT2
FeaturePlot(E_seurat1_ad_epi_bb, c("ENSMUSG00000038791"), reduction = "bbknn") #Clara #,"ENSMUSG00000022595")
FeaturePlot(E_seurat1_ad_epi_bb, c("ENSMUSG00000027676"), reduction = "bbknn") #Ciliated "ENSMUSG00000040703", "ENSMUSG00000024653" ,"ENSMUSG00000027676"

new.cluster.ids <- c("AT2", "AT2", "AT2", "AT2","AT2", "AT1","AT1","AT1","Ciliated", "Clara","Fetal lung","AT1","AT1","AT2","Clara","Fetal lung","Fetal lung","Fetal lung","Fetal lung") 
names(new.cluster.ids) <- levels(E_seurat1_ad_epi_bb)
E_seurat1_ad_epi_bb <- RenameIdents(E_seurat1_ad_epi_bb, new.cluster.ids)
DimPlot(E_seurat1_ad_epi_bb, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=20, k=5, trim=35, epi")
DimPlot(E_seurat1_ad_epi_bb[,E_seurat1_ad_epi_bb@active.ident=="Fetal lung"], reduction = "bbknn")

##Epi_fetal lung specification (fetal lung 부분 더 clustering 해볼것 -> fetal_AT2/fetal_AT1/fetal_Clara/fetal_Cili)
E_seurat1_ad_epi_bb_fetal <- E_seurat1_ad_epi_bb[,E_seurat1_ad_epi_bb$batch %in% c("E8.25","E9.5to11.5","E14.5.2","E16.5","E18.5","LungMap.E18.5")]
E_seurat1_ad_epi_bb_fetal <- FindNeighbors(E_seurat1_ad_epi_bb_fetal, dims = 1:40)
E_seurat1_ad_epi_bb_fetal <- FindClusters(E_seurat1_ad_epi_bb_fetal, resolution = 5) 
DimPlot(E_seurat1_ad_epi_bb_fetal, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id="ident")

FeaturePlot(E_seurat1_ad_epi_bb_fetal, c("ENSMUSG00000041247","ENSMUSG00000056370","ENSMUSG00000029375","ENSMUSG00000056370"), reduction = "bbknn") #AT2 cell
FeaturePlot(E_seurat1_ad_epi_bb_fetal, c("ENSMUSG00000037846","ENSMUSG00000000216","ENSMUSG00000067158","ENSMUSG00000021057"), reduction = "bbknn") #AT1 cell
FeaturePlot(E_seurat1_ad_epi_bb_fetal, c("ENSMUSG00000038791","ENSMUSG00000053279"), reduction = "bbknn") #Clara cell
FeaturePlot(E_seurat1_ad_epi_bb_fetal, c("ENSMUSG00000049811","ENSMUSG00000027676"), reduction = "bbknn") #Ciliated cell

plot <- FeaturePlot(E_seurat1_ad_epi_bb_fetal, "ENSMUSG00000056370", reduction = "bbknn")
plot
select.AT2 <- CellSelector(plot = plot)
head(select.AT2)
Idents(E_seurat1_ad_epi_bb, cells = select.AT2) <- "Fetal_AT2"

plot <- FeaturePlot(E_seurat1_ad_epi_bb_fetal, "ENSMUSG00000056370", reduction = "bbknn")
plot
select.AT1 <- CellSelector(plot = plot)
head(select.AT1)
Idents(E_seurat1_ad_epi_bb, cells = select.AT1) <- "Fetal_AT1"

plot <- FeaturePlot(E_seurat1_ad_epi_bb_fetal, "ENSMUSG00000038791", reduction = "bbknn")
plot
select.Clara <- CellSelector(plot = plot)
head(select.Clara)
Idents(E_seurat1_ad_epi_bb, cells = select.Clara) <- "Fetal_Clara"

plot <- FeaturePlot(E_seurat1_ad_epi_bb_fetal, "ENSMUSG00000049811", reduction = "bbknn")
plot
select.Cili <- CellSelector(plot = plot)
head(select.Cili)
Idents(E_seurat1_ad_epi_bb, cells = select.Cili) <- "Fetal_Cili"

DimPlot(E_seurat1_ad_epi_bb, reduction = "bbknn", group.by = "ident")

#epi_gene
##epi<->nonepi에서 markergene 뽑기
test<-E_seurat1_ad_bb
tmp_epithelial <- test$seurat_clusters %in% c('31','45','23','25','0','4','35','5','8','44','32','11')
tmp_fetalepithelial <- test$batch %in% c("E8.25","E14.5.2","E16.5","E18.5")
table(tmp_epithelial, tmp_fetalepithelial)
test$second_ident <- ifelse(tmp_epithelial | tmp_fetalepithelial, "epi", "nonepi")

test$second_ident<-test$second_ident%>%as.factor()
test@active.ident <- test$second_ident

group1 <- c('epi')
group2 <- NULL  
marker_dt <- FindMarkers(test, ident.1 = "epi", ident.2 = "nonepi") %>% as.data.frame %>% rownames_to_column('gene_name') %>% as.tibble() #1447
marker_dt <- marker_dt %>% filter(marker_dt$p_val <0.001) #1376
marker_dt_ordr <- marker_dt[ order(marker_dt$avg_logFC), ] 
marker_dt_ordr <- marker_dt_ordr %>% filter(abs(marker_dt_ordr$avg_logFC)>0.4) #623 (절대값 avg_logFC >0.4)
cur_markers <- marker_dt_ordr$gene_name

cell1 <- names(Idents(test)[Idents(test) %in% group1])
cell2 <- names(Idents(test)[Idents(test) %in% group2])
bbknn_dt <- test@reductions$bbknn@cell.embeddings %>% as.data.frame() %>% rownames_to_column('cell_id')
bbknn_dt <- bbknn_dt %>% mutate(group = ifelse(cell_id %in% cell1, paste(group1,collapse=','), ifelse(cell_id %in% cell2, paste(group2, collapse=','),'nonepi')))
ggplot(bbknn_dt, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), alpha=0.3)

DoHeatmap(test,features = final_genes)

cur_markers <- marker_dt_ordr %>% filter(marker_dt_ordr $avg_logFC > 0.4) %>% select("gene_name") #376 only epi high



#AT1/AT2/Ciliated/Club/Fetal lung marker gene뽑기 & 그 중에서 nonepi marker gene은 제거. 
#방법1
epi.markers <- FindAllMarkers(E_seurat1_ad_epi_bb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top100 <- epi.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) #500
DoHeatmap(E_seurat1_ad_epi_bb, features = final$mouse_ensembl) 
nonepi_marker_dt <- marker_dt%>% filter(marker_dt$avg_logFC < 0) %>% select("gene_name")
top100_rm_nonepi<-setdiff(top100$gene, nonepi_marker_dt$gene_name) #409

#방법2
marker_dt <- FindMarkers(E_seurat1_ad_epi_bb, ident.1 = c("AT2","Fetal_AT2")) %>% as.data.frame %>% rownames_to_column('gene_name') %>% as.tibble() 
marker_dt <- marker_dt %>% filter(marker_dt$p_val <0.001) #1376
marker_dt_ordr <- marker_dt[ order(marker_dt$avg_logFC), ] 
marker_dt_ordr <- marker_dt_ordr %>% filter(marker_dt_ordr$avg_logFC>0.4) #623 (절대값 avg_logFC >0.4)
at2_markers <- marker_dt_ordr$gene_name

DoHeatmap(E_seurat1_ad_epi_bb,features = final_genes$mouse_ensembl)

#epi marker + lung subtype markers = final_genes #623
fianl_genes <- union(top100_rm_nonepi, cur_markers$gene_name )

##luad related gene 들은 final_genes에서 제거. 
m_ensembl_h_gene<-read_tsv("./projects/ensembl/m_ensembl_h_gene.tsv") 
m_ensembl_h_gene<-m_ensembl_h_gene[!duplicated(m_ensembl_h_gene[2]),] #remove duplicated gene name

luad_genes<-c("APEX1", "AXIN2","CHRNA5", "CHRNA3", "CXCR2","CYP1B1", "CYP2E1","EPHX1","GSTM1","GSTP1","MTHFR","NAT2","NBN","SOD2","TP53","XRCC1","ATM","CYP1A1","ERCC1","ERCC2","GSTM1","OGG1","XRCC1", "HYKK", "PON1", "REV3L", "ATM", "CD3EAP", "CYP2A6", "HIF1A", "PDCD5", "PROM1", "TP53", "TP63", "WWOX", "XRCC1", "AGER", "BCL2", "CHRNA3", "CHRNA5", "CLPTM1L", "CYP1A1", "CYP1B1", "CYP2A6", "CYP2E1", "ELANE", "ERCC1", "ERCC2", "ERCC5", "ERCC6", "FGFR4", "GSTM1", "GSTP1", "GSTT1", "HRAS1", "IL10", "MAPKAPK2", "MDM2", "MIR146A", "MMP2", "MTRR", "NOD2", "SFTPB", "SOD2", "TERT", "UGT1A6", "XRCC1")
luad_genes<-luad_genes%>%as.data.frame()
luad_genes<-luad_genes[!duplicated(luad_genes[1]),] %>% as.data.frame()
colnames(luad_genes)<-"gene_name"
luad_genes<-merge(luad_genes, m_ensembl_h_gene, by.x = "gene_name", by.y = "Gene name") #luad_genes= 42개. 
luad_genes<-luad_genes[2]
luad_genes<-luad_genes%>%t()%>%as.vector()
final_genes<-setdiff(fianl_genes, luad_genes) #luad_gene들은 제거/ 그대로 유지.

#Final_gene to human ensembl; 
final_genes_df <- final_genes%>%as.data.frame()
colnames(final_genes_df) <- "epi_mouse"
names(ensembl_mouse_human)<-c("mouse_ensembl","human_ensembl") 
final<-merge(final_genes_df, ensembl_mouse_human, by.x = "epi_mouse", by.y = "mouse_ensembl")
#final<- final[!duplicated(final[1]),]
final<- final[!duplicated(final[2]),]
names(final)[1]<-"mouse_ensembl"

# +final gene 수정.
# cell cycle related gene 21개 + luad related gene 1개 제거
final_test <- final %>% filter(human_ensembl %in% setdiff(final$human_ensembl, c("ENSG00000118971","ENSG00000039068","ENSG00000198668","ENSG00000135446","ENSG00000105976","ENSG00000171848","ENSG00000073350","ENSG00000100714","ENSG00000150753","ENSG00000111057","ENSG00000148834","ENSG00000144381","ENSG00000155368","ENSG00000173207","ENSG00000010278","ENSG00000075131","ENSG00000129757","ENSG00000105968","ENSG00000188486","ENSG00000160752","ENSG00000110092","ENSG00000147883")))
final<- final_test

#single cell data processing -> with 600 genes (expression data)
mouse_epi_traj<-E_seurat1_ad_epi_bb@assays$RNA@data # cell; 5903 / gene; 11854
mouse_epi_traj<-mouse_epi_traj%>%as.data.frame()
mouse_epi_traj<-mouse_epi_traj%>%rownames_to_column()
names(mouse_epi_traj)[1] <- c("ensembl")
mouse_epi_600 <- merge(mouse_epi_traj, final, by.x = "ensembl", by.y = "mouse_ensembl")
mouse_epi_600 <- mouse_epi_600[,-1]

## + E8.25 제거 
#mouse_epi_traj<- E_seurat1_ad_epi_bb[,!(E_seurat1_ad_epi_bb$batch %in% c("E8.25"))]
#mouse_epi_traj<- mouse_epi_traj@assays$RNA@data # cell; 5718 / gene; 11854
#mouse_epi_traj<-mouse_epi_traj%>%as.data.frame()
#mouse_epi_traj<-mouse_epi_traj%>%rownames_to_column()
#names(mouse_epi_traj)[1] <- c("ensembl")
#mouse_epi_600 <- merge(mouse_epi_traj, final, by.x = "ensembl", by.y = "mouse_ensembl")
#mouse_epi_600 <- mouse_epi_600[,-1]

#LUAD_bulk_RNAseq with 600 genes (expression data_FPKM)
dir <- ("./projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUAD_551_FPKM/")
file_list<-list.files("/home/users/yunah1029/projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUAD_551_FPKM",, recursive=TRUE, pattern="*.FPKM.txt.gz$")

tcga <- mouse_epi_600_knn$human_ensembl %>% as.data.frame()
tcga$id  <- 1:nrow(tcga)
colnames(tcga) <- "gene"

for (file in file_list) {
  temp <- read.csv(paste0(dir, file), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = FALSE)
  colnames(temp) <- c("gene", strsplit(file, "/")[[1]][1])
  temp$gene <- sapply(strsplit(temp$gene, "[.]"), function(x) x[1])
  tcga <- merge(tcga, temp, by.x = "gene", by.y = "gene", all.x = TRUE, all.y = FALSE)
}
tcga[is.na(tcga)] <- 0
sum(is.na(tcga))

#tcga sample 별 mouse_epi_600과 cor -> 가장 가까운 cell 찾기. 
tcga<-tcga[order(tcga$"NA"), ] #mouse_epi_600의 gene 순서와 동일하게 정렬 -> 나중에 cbind
rownames(tcga) <- NULL
tcga<- tcga %>% as.data.frame() %>% column_to_rownames(var = "gene")
tcga<-tcga[,-1]
head(rownames(tcga))
head(colnames(tcga))

mouse_epi_600 <- mouse_epi_600%>%column_to_rownames(var = "human_ensembl")

single_bulk <- cbind(mouse_epi_600, tcga)
cor_single_bulk <- cor(single_bulk, method = "spearman") #method = "spearman"/ pearson(default)
cor_single_bulk <- cor_single_bulk %>% as.data.frame()
cor_single_bulk_sort <- cor_single_bulk[1:5903,5904:nrow(cor_single_bulk)]
cor_max <- apply(cor_single_bulk_sort,2, max) #correlation cutoff 설정 
#cor_single_bulk_sort<- cor_single_bulk_sort[,which(cor_max > 0.5)]
cor_single_bulk_sort_rank <- apply(cor_single_bulk_sort,2,function(x){order(-x)})[1:10,]%>%t()
cor_single_bulk_sort_rank <- cor_single_bulk_sort_rank %>% as.data.frame()

#Plot
mouse_epi_traj<-E_seurat1_ad_epi_bb
group1 <- c('Fetal lung')
group2 <- c('AT1')
group3 <- c('AT2')
group4 <- c('Ciliated')
group5 <- c('Clara')
group6 <- c('Fetal_AT1')
group7 <- c('Fetal_AT2')
group8 <- c('Fetal_Clara')
group9 <- c('Fetal_Cili')
cell1 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group1])
cell2 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group2])
cell3 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group3])
cell4 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group4])
cell5 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group5])
cell6 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group6])
cell7 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group7])
cell8 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group8])
cell9 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group9])
bbknn_dt <- mouse_epi_traj@reductions$bbknn@cell.embeddings %>% as.data.frame() %>% rownames_to_column('cell_id')
bbknn_dt <- bbknn_dt %>% mutate(group = ifelse(cell_id %in% cell1, paste(group1,collapse=','), 
                                               ifelse(cell_id %in% cell2, paste(group2, collapse=','),
                                                      ifelse(cell_id %in% cell3, paste(group3, collapse=','),
                                                             ifelse(cell_id %in% cell4, paste(group4, collapse=','),
                                                                    ifelse(cell_id %in% cell5, paste(group5, collapse=','),
                                                                           ifelse(cell_id %in% cell6, paste(group6, collapse=','),
                                                                                  ifelse(cell_id %in% cell7, paste(group7, collapse=','),
                                                                                         ifelse(cell_id %in% cell8, paste(group8, collapse=','),
                                                                                                "Fetal_Cili")))))))))
ggplot(bbknn_dt, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), alpha=0.5)

ggplot(bbknn_dt, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), size = 0.8, alpha = 0.7)+
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#7dab78", "#f57c16",
                                 "#F0E442", "#bf046a","#fa2a00","#0e5e05","#9e5416")) #0072B2

##
#a<-Embeddings(mouse_epi_traj,reduction = "bbknn") %>% as.data.frame()
cor_single_bulk_sort_rank <- cor_single_bulk_sort_rank[,1:3] %>% t()
matched <- data.frame()
matched_avg <- data.frame()

for(n in 1:ncol(cor_single_bulk_sort_rank)){
  tmp_matched <- a[cor_single_bulk_sort_rank[,n],]
  tmp_matched$group <- "LUAD:Highly correlated"
  tmp_matched <- tmp_matched %>% rownames_to_column(var = "cell_id")
  
  tmp_avg <- colMeans(tmp_matched[2:3]) %>% t() %>% as.data.frame()
  tmp_avg$group <- "LUAD:Average"
  tmp_avg$cell_id <- colnames(cor_single_bulk_sort_rank_3)[n]
  tmp_avg <- tmp_avg[,c(4,1,2,3)]
  
  matched<- rbind(matched, tmp_matched)
  matched_avg <- rbind(matched_avg, tmp_avg)
}

DimPlot(E_seurat1_ad_epi_bb[,E_seurat1_ad_epi_bb@active.ident %in% c("Fetal lung","Fetal_AT1","Fetal_AT2","Fetal_Clara","Fetal_Cili")], reduction = "bbknn", group.by = "ident")
DimPlot(E_seurat1_ad_epi_bb, reduction = "bbknn", group.by = "ident")

bbknn_dt_plot <- rbind(bbknn_dt, matched)
unique(bbknn_dt_plot$group)
ggplot(bbknn_dt_plot, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), size = 0.8, alpha = 0.7)+
  #scale_shape_manual(values = c( 16, 16, 16,16, 16, 2, 17))+
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#7dab78", "#f57c16",
                                 "#F0E442", "#bf046a","#fa2a00","#0e5e05","#9e5416","#3e0ffa")) #0072B2


#Hightly correlated cells 정리 -> expression pattern 확인
high_cor_cells<-table(cor_single_bulk_sort_rank[2,]) #1/2/3등 조절 가능
high_cor_cells%<>% as.data.frame()
high_cor_cells$cell_name <- colnames(cor_single_bulk)[high_cor_cells$Var1 %>% as.character() %>% as.numeric()]
#high_cor_cells<-high_cor_cells[-1]
high_cor_cells<-high_cor_cells[,c(3,1,2)] #column 순서 조절
head(high_cor_cells)
#high_cor_cells %>% write_csv("high_cor_cells3.csv") #1/2/3조절

high_cor_cells_sort <- single_bulk[,c(4987,5638,5654,5722,5728,5729,5878)]
head(high_cor_cells_sort)
high_cor_cells_sort$tcga_avg <- rowMeans(single_bulk[,5904:ncol(single_bulk)])

high_cor_cells_sort %>% write.csv("high_cor_cells_tcga_expression.csv")

## group 별 600 gene expression pattern; group 별 avg.exp with 600 genes
mouse_epi_traj<- E_seurat1_ad_epi_bb
group1 <- c("Fetal lung")
cell1 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group1])
fetal_600 <- mouse_epi_traj@assays$RNA@data %>% as.data.frame() 
fetal_600 <- fetal_600[,which(colnames(fetal_600) %in% cell1)]
fetal_600 <- fetal_600 %>% rownames_to_column( "ensembl")
fetal_600 <- merge(fetal_600, final, by.x = "ensembl", by.y = "mouse_ensembl")
fetal_600 <- fetal_600[,-1] #rm mouse_ensembl (human ensembl left)
fetal_600 <- fetal_600[!(duplicated(fetal_600$human_ensembl)),]
fetal_600 <- fetal_600 %>% mutate(fetal_avg.exp = rowMeans(fetal_600[,-(ncol(fetal_600))]))

head(fetal_600[,c(ncol(fetal_600), ncol(fetal_600)-1)])
tmp<-fetal_600[,c(ncol(fetal_600), ncol(fetal_600)-1)]
head(tmp)
#group_avg.exp_test <- merge(, tmp, all.x = TRUE)
#head(group_avg.exp_test)
group_avg.exp <- fetal_600 %>% select(human_ensembl, Clara_avg.exp) #처음 시작할 때만!
#group_avg.exp <- group_avg.exp_test
rm(group_avg.exp_test)
#ggplot(data = group_avg.exp) + 
#geom_bar(mapping = aes(x = human_ensembl ), stat = "identity", position = "dodge")
group_avg.exp %>% write_csv("group_avg_exp.csv")

#600gene expression heatmap 
library(NMF)
mouse_epi_600_hm <- mouse_epi_600_knn 
mouse_epi_annotation <- E_seurat1_ad_epi_bb@active.ident %>% as.data.frame() %>% rownames_to_column()
colnames(mouse_epi_annotation) <- c("cell", "ident")
mouse_epi_annotation <- mouse_epi_annotation %>% arrange(ident)
mouse_epi_annotation %<>% column_to_rownames('cell')
mouse_epi_600_hm <- mouse_epi_600_hm[,rownames(mouse_epi_annotation)] #ident 순서대로 정렬되도록
aheatmap(mouse_epi_600_hm, scale = "col", annCol = mouse_epi_annotation, Rowv = NA, Colv = NA, height = 30, cexRow = 10)

#LUSC cor
#LUSC_bulk_RNAseq with 600 genes (expression data_FPKM)
dir <- ("./projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUSC_551_FPKM/")
file_list<-list.files("/home/users/yunah1029/projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUSC_551_FPKM/",, recursive=TRUE, pattern="*.FPKM.txt.gz$")

mouse_epi_600_knn %<>% rownames_to_column(var = "human_ensembl")
tcga_lusc <- mouse_epi_600_knn$human_ensembl %>% as.data.frame()
tcga_lusc$id  <- 1:nrow(tcga_lusc)
colnames(tcga_lusc) <- "gene"
for (file in file_list) {
  temp_lusc <- read.csv(paste0(dir, file), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = FALSE)
  colnames(temp_lusc) <- c("gene", strsplit(file, "/")[[1]][1])
  temp_lusc$gene <- sapply(strsplit(temp_lusc$gene, "[.]"), function(x) x[1])
  tcga_lusc <- merge(tcga_lusc, temp_lusc, by.x = "gene", by.y = "gene", all.x = TRUE, all.y = FALSE)
}
tcga_lusc[is.na(tcga_lusc)] <- 0
sum(is.na(tcga_lusc))

#tcga_lusc sample 별 mouse_epi_600과 cor -> 가장 가까운 cell 찾기. 
tcga_lusc<-tcga_lusc[order(tcga_lusc$"NA"), ] #mouse_epi_600의 gene 순서와 동일하게 정렬 -> 나중에 cbind
rownames(tcga_lusc) <- NULL
tcga_lusc<- tcga_lusc %>% as.data.frame() %>% column_to_rownames(var = "gene")
tcga_lusc<-tcga_lusc[,-1]
head(rownames(tcga_lusc))
head(colnames(tcga_lusc))
mouse_epi_600_knn <- mouse_epi_600_knn%>%column_to_rownames(var = "human_ensembl")

single_bulk_lusc <- cbind(mouse_epi_600_knn, tcga_lusc)
cor_single_bulk_lusc <- cor(single_bulk_lusc, method = "pearson") #method = "spearman"/ pearson(default)
cor_single_bulk_lusc_sort <- cor_single_bulk_lusc[1:5903,5904:nrow(cor_single_bulk_lusc)]
cor_max_lusc <- apply(cor_single_bulk_lusc_sort,2, max) #correlation cutoff 설정 
#cor_single_bulk_lusc_sort_0.3 <- cor_single_bulk_lusc_sort[,which(cor_max_lusc > 0.45)]
cor_single_bulk_lusc_sort_rank <- apply(cor_single_bulk_lusc_sort,2,function(x){order(-x)})[1:7,]%>%t()
cor_single_bulk_lusc_sort_rank <- cor_single_bulk_lusc_sort_rank %>% as.data.frame()


#Plot
group1 <- c('Fetal lung')
group2 <- c('AT1')
group3 <- c('AT2')
group4 <- c('Ciliated')
group5 <- c('Clara')
cell1 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group1])
cell2 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group2])
cell3 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group3])
cell4 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group4])
cell5 <- names(Idents(mouse_epi_traj)[Idents(mouse_epi_traj) %in% group5])

bbknn_dt <- mouse_epi_traj@reductions$bbknn@cell.embeddings %>% as.data.frame() %>% rownames_to_column('cell_id')
bbknn_dt <- bbknn_dt %>% mutate(group = ifelse(cell_id %in% cell1, paste(group1,collapse=','), 
                                               ifelse(cell_id %in% cell2, paste(group2, collapse=','),
                                                      ifelse(cell_id %in% cell3, paste(group3, collapse=','),
                                                             ifelse(cell_id %in% cell4, paste(group4, collapse=','),
                                                                    "Clara")))))
ggplot(bbknn_dt, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), alpha=0.5)

##
a<-Embeddings(mouse_epi_traj,reduction = "bbknn") %>% as.data.frame()
cor_single_bulk_lusc_sort_rank <- cor_single_bulk_lusc_sort_rank[,1:3] %>% t()
matched_lusc <- data.frame()
matched_avg_lusc <- data.frame()

for(n in 1:ncol(cor_single_bulk_lusc_sort_rank)){
  tmp_matched_lusc <- a[cor_single_bulk_lusc_sort_rank[,n],]
  tmp_matched_lusc$group <- "LUSC:Highly correlated"
  tmp_matched_lusc <- tmp_matched_lusc %>% rownames_to_column(var = "cell_id")
  
  #tmp_avg_lusc <- colMeans(tmp_matched_lusc[2:3]) %>% t() %>% as.data.frame()
  #tmp_avg_lusc$group <- "LUSC:Average"
  #tmp_avg_lusc$cell_id <- colnames(cor_single_bulk_lusc_sort_rank_3)[n]
  #tmp_avg_lusc <- tmp_avg_lusc[,c(4,1,2,3)]
  
  matched_lusc<- rbind(matched_lusc, tmp_matched_lusc)
  #matched_avg_lusc <- rbind(matched_avg_lusc, tmp_avg_lusc)
}

DoHeatmap(bbknn_dt)

bbknn_dt_lusc_plot <- rbind(bbknn_dt, matched_lusc)
ggplot(bbknn_dt_lusc_plot, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), size = 1, alpha = 0.3)+
  #scale_shape_manual(values = c( 16, 16, 16,16, 16, 2, 17))+
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#52854C",  "#D16103",
                                 "#F0E442", "#0072B2", "#D55E00"))


#human noraml scRNAseq to mouse epi traj
##Find AT1/AT2/Ciliated/Clara cluster from human normal scRNAseq data
donor1<-Read10X_h5("./projects/01_LUAD_scRNA/donor_single_cell/GSM3489185_Donor_02_filtered_gene_bc_matrices_h5.h5")
rownames(donor1)<-make.unique(rownames(donor1))
colnames(donor1)=paste("D1",colnames(donor1))
donor1 %<>% as.data.frame()  %<>% rownames_to_column()

human_ensmebl_genename<-read_tsv("./projects/01_LUAD_scRNA/donor_single_cell/human_ENSEMBL_genename.tsv")
donor1 <- merge(donor1, human_ensmebl_genename, by.x = "rowname", by.y = "Gene name")
table(duplicated(donor1[ncol(donor1)]))
donor1 <- donor1[,-1]
donor1 %<>% column_to_rownames(var = "Gene stable ID")

donor1<-CreateSeuratObject(count=donor1)
donor1<-NormalizeData(object=donor1)
donor1<-FindVariableFeatures(object=donor1,x.low.cutoff=0.0125,x.high.cutoff=3,y.cutoff=0.5)
donor1<-ScaleData(object=donor1,do.par=T,num.cores=2)
donor1<-RunPCA(object=donor1, npcs = 50)
ElbowPlot(donor1, ndims = 50)
donor1<-FindNeighbors(object=donor1, dims = 1:30)
donor1<-FindClusters(object=donor1,dims=1:30,resolution=0.7,save.SNN=T,force.recalc=T)
donor1<-RunUMAP(object = donor1, dims = 1:30)
donor1<-RunTSNE(object = donor1, dims = 1:30)
DimPlot(object=donor1,reduction = "umap") %>% LabelClusters(id = "ident")
DimPlot(object=donor1,reduction = "tsne") %>% LabelClusters(id = "ident")

markers.all.donor1 <- FindAllMarkers(donor1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker.6 <- FindMarkers(donor1, ident.1 = c("6"))
marker.13 <- FindMarkers(donor1, ident.1 = c("13"))
marker.16 <- FindMarkers(donor1, ident.1 = c("16"))
marker.10 <- FindMarkers(donor1, ident.1 = c("10"))
marker.11 <- FindMarkers(donor1, ident.1 = c("11"))
marker.12 <- FindMarkers(donor1, ident.1 = c("12"))
marker.15 <- FindMarkers(donor1, ident.1 = c("15"))
marker.17 <- FindMarkers(donor1, ident.1 = c("17"))

FeaturePlot(donor1, c("ENSG00000204305","ENSG00000182010", "ENSG00000012171","ENSG00000112782"), reduction = "umap") #AT1/ cluster 6
FeaturePlot(donor1, c("ENSG00000164265","ENSG00000149021","ENSG00000145113","ENSG00000196188"), reduction = "umap") #Clara / cluster 16 
FeaturePlot(donor1, c("ENSG00000129654","ENSG00000117477","ENSG00000159588","ENSG00000248712"), reduction = "umap") #Cili / cluster 13
FeaturePlot(donor1, c("ENSG00000096088", "ENSG00000131400"), reduction = "umap") #AT2 / cluster 1,2,3,4,5,8,9,14
FeaturePlot(donor1, c("ENSG00000019169","ENSG00000173391"), reduction = "umap") #AM / cluster 0,7
FeaturePlot(donor1, c("ENSG00000081237"), reduction = "umap") #immune(cd45+)
FeaturePlot(donor1, c("ENSG00000169413","ENSG00000132514"), reduction = "umap") #Dendritic cell/ cluster 12
FeaturePlot(donor1, c("ENSG00000085265","ENSG00000038427"), reduction = "umap") #Monocyte/ cluster 11
FeaturePlot(donor1, c("ENSG00000198851","ENSG00000271503"), reduction = "umap") #T cell/ cluster 15
FeaturePlot(donor1, c("ENSG00000011465","ENSG00000139329"), reduction = "umap") #Fibroblast / cluster 17
FeaturePlot(donor1, c("ENSG00000011465","ENSG00000139329"), reduction = "umap") #Endothelial cell / cluster 10

new.cluster.ids <- c("h_Alveolar macrophage", "h_AT2", "h_AT2","h_AT2","h_AT2","h_AT2","h_AT1","h_Alveolar macrophage","h_AT2","h_AT2","h_Endothelial","h_Monocyte","h_Dendritic","h_Ciliated","h_AT2","h_T","h_Clara","h_Fibroblast") 
names(new.cluster.ids) <- levels(donor1)
donor1 <- RenameIdents(donor1, new.cluster.ids)
DimPlot(donor1, reduction = "umap", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("Human lung")

donor_AT1 <- donor1@assays$RNA@data[,donor1@active.ident == "h_AT1"] %>% as.data.frame() %>% rownames_to_column(var = "gene") #166
donor_AT2 <- donor1@assays$RNA@data[,donor1@active.ident == "h_AT2"] %>% as.data.frame() %>% rownames_to_column(var = "gene") #2050
donor_cili <- donor1@assays$RNA@data[,donor1@active.ident == "h_Ciliated"] %>% as.data.frame() %>% rownames_to_column(var = "gene") #79
donor_clara <- donor1@assays$RNA@data[,donor1@active.ident == "h_Clara"] %>% as.data.frame() %>% rownames_to_column(var = "gene") #46

##human lung cluster마다 mouse traj에 붙여보기 (#표시한곳 조건 바꿔가면서 cluster별 그릴 수 있음)
mouse_epi_traj<-E_seurat1_ad_epi_bb@assays$RNA@data # cell; 5903 / gene; 11854
mouse_epi_traj<-mouse_epi_traj%>%as.data.frame()
mouse_epi_traj<-mouse_epi_traj%>%rownames_to_column()
names(mouse_epi_traj)[1] <- c("ensembl")
mouse_epi_600 <- merge(mouse_epi_traj, final, by.x = "ensembl", by.y = "mouse_ensembl")
mouse_epi_600 <- mouse_epi_600[,-1]

###-### bbknn(1)_final_troubleshooting 과 연계 -> plot
mouse_epi_600_knn <- mouse_epi_600_knn %>% as.data.frame() %>% rownames_to_column(var = "human_ensembl")
hu_AT1 <- mouse_epi_600_knn$human_ensembl %>% as.data.frame()
hu_AT1$id  <- 1:nrow(hu_AT1)
colnames(hu_AT1) <- c("gene","id")
hu_AT1 <- merge(hu_AT1, donor_AT1, by.x = "gene", by.y = "gene", all.x = TRUE, all.y = FALSE) #여기 바꾸면서 조절
hu_AT1[is.na(hu_AT1)] <- 0
sum(is.na(hu_AT1))
hu_AT1 <- hu_AT1[order(hu_AT1$id),]
rownames(hu_AT1) <-NULL
hu_AT1 <- hu_AT1 %>% as.data.frame() %>% column_to_rownames(var = "gene")
hu_AT1 <- hu_AT1[,-1]
head(rownames(hu_AT1))
head(colnames(hu_AT1))
dim(hu_AT1)
mouse_epi_600_knn <- mouse_epi_600_knn%>%column_to_rownames(var = "human_ensembl")
###

m_h_AT1 <- cbind(mouse_epi_600, hu_AT1)
cor_m_h_AT1 <- cor(m_h_AT1, method = "pearson") #method = "spearman"/ pearson(default)
cor_m_h_AT1 <- cor_m_h_AT1 %>% as.data.frame()
cor_m_h_AT1_sort <- cor_m_h_AT1[1:5903,5904:nrow(cor_m_h_AT1)]
cor_m_h_AT1_max <- apply(cor_m_h_AT1_sort,2, max) #correlation cutoff 설정 
#cor_single_bulk_sort_0.3 <- cor_m_h_AT1_sort[,which(cor_max > 0.45)]
cor_m_h_AT1_sort_rank <- apply(cor_m_h_AT1_sort,2,function(x){order(-x)})[1:10,]%>%t()
cor_m_h_AT1_sort_rank <- cor_m_h_AT1_sort_rank %>% as.data.frame()

##
#a<-Embeddings(mouse_epi_traj,reduction = "bbknn") %>% as.data.frame()
cor_m_h_AT1_sort_rank <- cor_m_h_AT1_sort_rank[,1:3] %>% t()
matched <- data.frame()

for(n in 1:ncol(cor_m_h_AT1_sort_rank)){
  tmp_matched <- a[cor_m_h_AT1_sort_rank[,n],]
  tmp_matched$group <- "h_AT2" #여기도 cluster 바뀔때마다 조절
  tmp_matched <- tmp_matched %>% rownames_to_column(var = "cell_id")
  matched<- rbind(matched, tmp_matched)
}

bbknn_dt_plot <- rbind(bbknn_dt, matched)
ggplot(bbknn_dt_plot, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), size = 1, alpha = 0.3)+
  #scale_shape_manual(values = c( 16, 16, 16,16, 16, 2, 17))+
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#52854C",  "#D16103",
                                 "#F0E442", "#0072B2", "#D55E00"))



