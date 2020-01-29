#Mouse embryo data loading & processing (extract foregut/ gene name>ENSEMBL)
#E8.25_foregut & E9.5_11.5 & E14.5 & P6~P10(1/2/3) & collected(E14.5/E16.5/E18.5/PND107)
E8.25<-read_tsv("./Mouse_embryo/E8.25_E-MTAB-6153.tsv")
E9.5_11.5<-read_csv("./Mouse_embryo/GSE87038_E9.5_E10.5_E11.5_counts_matrix.csv")
E14.5<-read.csv("./Mouse_embryo/GSM2906422_FetalLung_dge.txt", sep = ' ', stringsAsFactors = FALSE, na.strings = "")
P6_10_1<-read.csv("./Mouse_embryo/GSM2906429_Lung1_dge.txt", sep = ' ', stringsAsFactors = FALSE, na.strings = "")
P6_10_2<-read.csv("./Mouse_embryo/GSM2906430_Lung2_dge.txt", sep = ' ', stringsAsFactors = FALSE, na.strings = "")
P6_10_3<-read.csv("./Mouse_embryo/GSM2906431_Lung3_dge.txt", sep = ' ', stringsAsFactors = FALSE, na.strings = "")
E8.25_foregut<-E8.25%>%select(gene_id, starts_with('foregut')) #ncol(E8.25_foregut)=185
E8.25_foregut_ensembl<-E8.25%>%select(gene_id, starts_with('foregut')) #ncol(E8.25_foregut)=185 
P1_1<-read_tsv("./projects/LUAD_scRNA/mouse_embryo/yan/GSM3464101_Batch1_out_gene_exon_tagged.dge.txt")
P1_2<-read_tsv("./projects/LUAD_scRNA/mouse_embryo/yan/GSM3464102_Batch2_out_gene_exon_tagged.dge.txt")
collected<-read_tsv("./projects/LUAD_scRNA/mouse_embryo/treulein/GSE52583/collected.txt")
LungMap_E18.5 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_E18.5.xlsx")
LungMap_P1 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_P01.xlsx")
LungMap_P3 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_P03.xlsx")
LungMap_P7 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_P07.xlsx")
LungMap_P10 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_P10.xlsx")
LungMap_P14 <- read_excel("projects/LUAD_scRNA/mouse_embryo/lungmap/LungMap_P14.xlsx")
LungMap_E18.5 <- LungMap_E18.5[-1,] #LungMap data는 1행에 필요없는 내용 있음. 
LungMap_P1 <- LungMap_P1[-1,]
LungMap_P3 <- LungMap_P3[-1,]
LungMap_P7 <- LungMap_P7[-1,]
LungMap_P10 <- LungMap_P10[-1,]
LungMap_P14 <- LungMap_P14[-1,]
LungMap_E18.5 <- rename(LungMap_E18.5, gene_name = ...1) 
LungMap_P1 <- rename(LungMap_P1, gene_name = ...1) 
LungMap_P3 <- rename(LungMap_P3, gene_name = ...1) 
LungMap_P7 <- rename(LungMap_P7, gene_name = ...1) 
LungMap_P10 <- rename(LungMap_P10, gene_name = ...1)
LungMap_P14 <- rename(LungMap_P14, gene_name = "Type") 

#mouse gene_name ↔ gene ensembl table (reference)
gene_name_ensembl_mouse<-read_tsv("./Mouse_embryo/gene_name_ensembl_mouse.tsv")
gene_name_ensembl_mouse<- gene_name_ensembl_mouse[!duplicated(gene_name_ensembl_mouse[1]),] #remove duplicated gene name

#gene name > gene_ensembl
E9.5_11.5_ensembl<-merge(E9.5_11.5,gene_name_ensembl_mouse,by.x = "Gene",by.y = "Gene name")
E9.5_11.5_ensembl<-E9.5_11.5_ensembl[,-1] # remove mouse gene name # Gene stable ID in [,225]
E14.5<-rownames_to_column(E14.5,"gene_name")
E14.5_ensembl<-E14.5_ensembl[,-1]
P6_10_1<-P6_10_1%>%rownames_to_column("gene_name")
P6_10_1_ensembl<-merge(P6_10_1,gene_name_ensembl_mouse,by.x = "gene_name",by.y = "Gene name")
P6_10_1_ensembl<-P6_10_1_ensembl[,-1]
P6_10_2<-P6_10_2%>%rownames_to_column("gene_name")
P6_10_2_ensembl<-merge(P6_10_2,gene_name_ensembl_mouse,by.x = "gene_name",by.y = "Gene name")
P6_10_2_ensembl<-P6_10_2_ensembl[,-1]
P6_10_3<-P6_10_3%>%rownames_to_column("gene_name")
P6_10_3_ensembl<-merge(P6_10_3,gene_name_ensembl_mouse,by.x = "gene_name",by.y = "Gene name")
P6_10_3_ensembl<-P6_10_3_ensembl[,-1]
P1_1_ensembl<-merge(P1_1, gene_name_ensembl_mouse, by.x = "GENE", by.y = "Gene name")
P1_1_ensembl<-P1_1_ensmebl[,-1]
names(P1_1_ensembl)[ncol(P1_1_ensembl)] <- c("Gene stable ID")
P1_2_ensembl<-merge(P1_2, gene_name_ensembl_mouse, by.x = "GENE", by.y = "Gene name")
P1_2_ensembl<-P1_2_ensembl[,-1]
names(P1_2_ensembl)[ncol(P1_2_ensembl)] <- c("Gene stable ID")
collected_ensembl<-merge(collected, gene_name_ensembl_mouse, by.x = "gene_symbol", by.y = "Gene name")
collected_ensembl<-collected_ensembl[,-1]
names(collected_ensembl)[ncol(collected_ensembl)] <- c("Gene stable ID")
LungMap_E18.5<-merge(LungMap_E18.5, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_E18.5_ensembl<-LungMap_E18.5[,-1]
LungMap_P1<-merge(LungMap_P1, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_P1_ensembl<-LungMap_P1[,-1]
LungMap_P3<-merge(LungMap_P3, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_P3_ensembl<-LungMap_P3[,-1]
LungMap_P7<-merge(LungMap_P7, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_P7_ensembl<-LungMap_P7[,-1]
LungMap_P10<-merge(LungMap_P10, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_P10_ensembl<-LungMap_P10[,-1]
LungMap_P14<-merge(LungMap_P14, gene_name_ensembl_mouse, by.x = "gene_name", by.y = "Gene name")
LungMap_P14_ensembl<-LungMap_P14[,-1]

# Embryo datas 공통된 genes로 합치기(cbind) & batch vector 부여
union_gene_vector = intersect(E8.25_foregut_ensembl$gene_id, E9.5_11.5_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,E14.5_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,P6_10_2_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,P6_10_1_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,P6_10_3_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,P1_1_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,P1_2_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,collected_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_E18.5_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_P1_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_P3_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_P7_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_P10_ensembl$`Gene stable ID`)
union_gene_vector = intersect(union_gene_vector,LungMap_P14_ensembl$`Gene stable ID`)

A = E8.25_foregut_ensembl[match(union_gene_vector, E8.25_foregut_ensembl$gene_id),-1]
colnames(A) = paste0("E8.25_", colnames(A))
B = E9.5_11.5_ensembl[match(union_gene_vector, E9.5_11.5_ensembl$`Gene stable ID`),-ncol(E9.5_11.5_ensembl)]
colnames(B) = paste0("E9.5to11.5_", colnames(B))
C = E14.5_ensembl[match(union_gene_vector, E14.5_ensembl$'Gene stable ID'), -ncol(E14.5_ensembl)]
colnames(C) = paste0("E14.5.1_", colnames(C))
D = P6_10_1_ensembl[match(union_gene_vector,P6_10_1_ensembl$`Gene stable ID`),-ncol(P6_10_1_ensembl)]
colnames(D) = paste0("P6to10.1_", colnames(D))
G = P6_10_2_ensembl[match(union_gene_vector,P6_10_2_ensembl$`Gene stable ID`),-ncol(P6_10_2_ensembl)]
colnames(G) = paste0("P6to10.2_", colnames(G))
H = P6_10_3_ensembl[match(union_gene_vector,P6_10_3_ensembl$`Gene stable ID`),-ncol(P6_10_3_ensembl)]
colnames(H) = paste0("P6to10.3_", colnames(H))
I = P1_1_ensembl[match(union_gene_vector,P1_1_ensembl$`Gene stable ID`),-ncol(P1_1_ensembl)]
colnames(I) = paste0("PND1.1_", colnames(I))
J = P1_2_ensembl[match(union_gene_vector,P1_2_ensembl$`Gene stable ID`),-ncol(P1_2_ensembl)]
colnames(J) = paste0("PND1.2_", colnames(J))
L = LungMap_E18.5_ensembl[match(union_gene_vector,LungMap_E18.5_ensembl$`Gene stable ID`),-ncol(LungMap_E18.5_ensembl)]
colnames(L) = paste0("LungMap.E18.5_", colnames(L))
M = LungMap_P1_ensembl[match(union_gene_vector,LungMap_P1_ensembl$`Gene stable ID`),-ncol(LungMap_P1_ensembl)]
colnames(M) = paste0("LungMap.PND1_", colnames(M))
N = LungMap_P3_ensembl[match(union_gene_vector,LungMap_P3_ensembl$`Gene stable ID`),-ncol(LungMap_P3_ensembl)]
colnames(N) = paste0("LungMap.PND3_", colnames(N))
O = LungMap_P7_ensembl[match(union_gene_vector,LungMap_P7_ensembl$`Gene stable ID`),-ncol(LungMap_P7_ensembl)]
colnames(O) = paste0("LungMap.PND7_", colnames(O))
P = LungMap_P10_ensembl[match(union_gene_vector,LungMap_P10_ensembl$`Gene stable ID`),-ncol(LungMap_P10_ensembl)]
colnames(P) = paste0("LungMap.PND10_", colnames(P))
Q = LungMap_P14_ensembl[match(union_gene_vector,LungMap_P14_ensembl$`Gene stable ID`),-ncol(LungMap_P14_ensembl)]
colnames(Q) = paste0("LungMap.PND14_", colnames(Q))

K = collected_ensembl[match(union_gene_vector,collected_ensembl$`Gene stable ID`),-ncol(collected_ensembl)]
K1<-K[1:80]
K2<-K[81:125]
K3<-K[126:171]
K4<-K[172:198]
colnames(K1) = paste0("E18.5_", colnames(K1))
colnames(K2) = paste0("E14.5.2_", colnames(K2))
colnames(K3) = paste0("P15_", colnames(K3))
colnames(K4) = paste0("E16.5_", colnames(K4))
K = cbind(K1, K2, K3, K4)

# downsampling → E14.5.1 data만 1k cell로 downsampling
C_1k<-knn(C,10)
sample<-sample(ncol(C_1k),1000)
C_1k<-C_1k[,c(sample)]

#bbknn1 pnd1.2 추가. (pnd1.1 제거)/ total cell 
E1 = cbind(A,B,C_1k,D,G,H,J,K,L,M,N,O,P,Q)
rownames(E1) = union_gene_vector
table(duplicated(t(E1)))
E1 = E1[,!duplicated(t(E1))]
batch_vector = stringr::str_replace(colnames(E1), "_.*","")

E_seurat1_ad = CreateSeuratObject(E1)
E_seurat1_ad$batch = batch_vector
E_seurat1_ad <-NormalizeData(E_seurat1_ad)
E_seurat1_ad <-FindVariableFeatures(E_seurat1_ad)
E_seurat1_ad <-ScaleData(E_seurat1_ad)
E_seurat1_ad <- RunPCA(E_seurat1_ad, npcs = 50)
E_seurat1_ad <- JackStraw(E_seurat1_ad,dims = 50)
E_seurat1_ad <- ScoreJackStraw(E_seurat1_ad, dims = 1:50)
JackStrawPlot(E_seurat1_ad, dims = 1:50)
ElbowPlot(E_seurat1_ad, ndims = 50)
E_seurat1_ad <- FindNeighbors(E_seurat1_ad, dims = 1:40)
E_seurat1_ad <- FindClusters(E_seurat1_ad, resolution = 2) #resolution 변경한 것임. 
E_seurat1_ad <- RunTSNE(E_seurat1_ad, dims = 1:40)
E_seurat1_ad <- RunUMAP(E_seurat1_ad, dims = 1:40)

#bbknn1_total
pdf("bbknn1_total(PND1.2)_35.pdf", 6,5)
for(k in c(3:12)){
  for(trim in c(20,25,30,35,40,45)){
    print(
      E_seurat1_ad %>% RunBBKNN(dims.use = 1:35,
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

E_seurat1_ad_bb<- E_seurat1_ad %>% RunBBKNN(dims.use = 1:40,
                                            neighbors_within_batch =  5,
                                            trim=45,
                                            batch.key = "batch",
                                            python.path = "/home/users/yunah1029/anaconda3/bin/python") 
DimPlot(E_seurat1_ad_bb, reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; dims=40, k=5, trim=45, all")


DimPlot(E_seurat1_ad_bb, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=40, k=5, trim=45, all")




#Find cluster_total 
group1 <- c('10')
group2 <- NULL  
marker_dt <- FindMarkers(E_seurat1_ad_bb, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
marker_10<-marker_dt

FeaturePlot(E_seurat1_ad_bb, "ENSMUSG00000056973", reduction = "bbknn")

new.cluster.ids <- c("Epithelial", "Mesenchymal", "Mesenchymal", "Epithelial", "Immune", "Immune", "Immune", "Immune", "Immune","Fetal lung", "Epithelial", "Epithelial", "Mesenchymal","Immune", "Mesenchymal","Immune","Mesenchymal","Mesenchymal","Fetal lung","Fetal lung","Mesenchymal","Fetal lung","Epithelial","Fetal lung","Epithelial","Immune","Immune","Mesenchymal") 
names(new.cluster.ids) <- levels(E_seurat1_ad_bb)
E_seurat1_ad_bb <- RenameIdents(E_seurat1_ad_bb, new.cluster.ids)
DimPlot(E_seurat1_ad_bb, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=40, k=5, trim=45, all")

#bbknn1_Epithelial cells
E_seurat1_ad_epi <- E_seurat1_ad_bb[,(E_seurat1_ad_bb@active.ident %in% c('Epithelial')|E_seurat1_ad_bb$batch %in% c("E8.25","E9.5to11.5","E14.5.2","E16.5"))]
DimPlot(E_seurat1_ad_epi, reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; dims=40, k=5, trim=45, epi")


unique(E_seurat1_ad_epi$batch)
ElbowPlot(E_seurat1_ad_epi)

E_seurat1_ad_epi <- FindNeighbors(E_seurat1_ad_epi, dims = 1:20)
E_seurat1_ad_epi <- FindClusters(E_seurat1_ad_epi, resolution = 0.7)
E_seurat1_ad_epi <- RunTSNE(E_seurat1_ad_epi, dims = 1:20)
E_seurat1_ad_epi <- RunUMAP(E_seurat1_ad_epi, dims = 1:20)

pdf("bbknn1_epi(PND1.2).pdf", 6,5)
for(k in c(3:12)){
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
                                                    neighbors_within_batch =  4,
                                                    trim=45,
                                                    batch.key = "batch",
                                                    python.path = "/home/users/yunah1029/anaconda3/bin/python") 
DimPlot(E_seurat1_ad_epi_bb, reduction = "bbknn", group.by = "batch") %>% LabelClusters(id = "batch") + ggtitle("bbknn1; dims=20, k=4, trim=45, epi")
DimPlot(E_seurat1_ad_epi_bb, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=20, k=4, trim=45, epi")

#Find cluster_total 
group1 <- c('7','8','12')
group2 <- NULL  
marker_dt <- FindMarkers(E_seurat1_ad_epi_bb, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
View(marker_dt)

FeaturePlot(E_seurat1_ad_epi_bb, c("ENSMUSG00000028583"), reduction = "bbknn")
FeaturePlot(E_seurat1_ad_epi_bb, c("ENSMUSG00000042784","ENSMUSG00000028583"), reduction = "bbknn")

new.cluster.ids <- c("AT2", "AT2", "AT2", "AT1", "AT1", "Clara","Ciliated","Fetal lung","Fetal lung","AT1","AT2","Fetal lung","Fetal lung") 
names(new.cluster.ids) <- levels(E_seurat1_ad_epi_bb)
E_seurat1_ad_epi_bb <- RenameIdents(E_seurat1_ad_epi_bb, new.cluster.ids)
DimPlot(E_seurat1_ad_epi_bb, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=20, k=4, trim=45, epi")
