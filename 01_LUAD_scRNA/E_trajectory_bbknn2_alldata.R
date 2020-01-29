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

# downsampling → 500 cell
C_500<-knn(C,10)
D_500<-knn(D,10)
G_500<-knn(G,10)
H_500<-knn(H,10)
I_500<-knn(I,10)
J_500<-knn(J,10)

sample<-sample(ncol(C_500),500)
C_500<-C_500[,c(sample)]
sample<-sample(ncol(D_500),500)
D_500<-D_500[,c(sample)]
sample<-sample(ncol(G_500),500)
G_500<-G_500[,c(sample)]
sample<-sample(ncol(H_500),500)
H_500<-H_500[,c(sample)]
sample<-sample(ncol(I_500),500)
I_500<-I_500[,c(sample)]
sample<-sample(ncol(J_500),500)
J_500<-J_500[,c(sample)]

#collected_ensembl 자료는 E8.5~P15까지 data가 다 섞여있기 때문에 추가적인 과정 필요. 
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

E = cbind(A,B,C_500,D_500,G_500,H_500,I_500,J_500,K,L,M,N,O,P,Q)
rownames(E) = union_gene_vector
table(duplicated(t(E)))
E = E[,!duplicated(t(E))]
table(duplicated(t(E)))
batch_vector = stringr::str_replace(colnames(E), "_.*","")

#Seurat (+bbknn2)
#function(FindBatchBalanceNeighbor)등록 해야함. script → bbknn2
E_seurat2 <- CreateSeuratObject(E)
E_seurat2$batch = batch_vector
E_seurat2 <- E_seurat2%>% NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()
E_seurat2_bb <- E_seurat2 %>% FindBatchBalancedNeighbors(batch.key = "batch",
                                                         neighbors_within_batch = 9,
                                                         prune.SNN = 1/15) 
E_seurat2_bb <- E_seurat2_bb %>%  FindClusters(graph.name = "bbknn_snn") 
E_seurat2_bb <- E_seurat2_bb %>%  RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_")

DimPlot(E_seurat2_bb, reduction = "bbknn", group.by = "batch")%>%LabelClusters(id = "batch") + ggtitle("k=6, prune= 1/15, all cells")
DimPlot(E_seurat2_bb, reduction = "bbknn", group.by = "ident")%>%LabelClusters(id = "ident")

#bbknn 가장 적절한 조건 찾기(neighbors, prune 조절)
pdf("bbknn2_test.pdf", 6,5)
for(k in c(5:15)){
  for(prune in c(1,2)){
    print(
      E_seurat2 %>% FindBatchBalancedNeighbors(batch.key = "batch",
                                               neighbors_within_batch = k,
                                               prune.SNN = prune/15) %>%  
        FindClusters(graph.name = "bbknn_snn") %>%
        RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_") %>% 
        DimPlot(reduction = "bbknn", group.by = "batch") %>% 
        LabelClusters(id="batch") + ggtitle(paste0("neighbors=",k,",prune=",prune,"/15",",every cells"))
    )
  }
}
dev.off()

#원하는 조건으로 다시 그리기
E_seurat2_bb <- E_seurat2 %>% FindBatchBalancedNeighbors(batch.key = "batch",
                                                         neighbors_within_batch = 9,
                                                         prune.SNN = 1/15) 
E_seurat2_bb <- E_seurat2_bb %>%  FindClusters(graph.name = "bbknn_snn", resolution = 1.2) 
E_seurat2_bb <- E_seurat2_bb %>%  RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_")

DimPlot(E_seurat2_bb, reduction = "bbknn", group.by = "batch")%>%LabelClusters(id = "batch") + ggtitle("bbknn2; k=9, prune= 1/15, resolution=1.2, all cells")
DimPlot(E_seurat2_bb, reduction = "bbknn", group.by = "ident")%>%LabelClusters(id = "ident")+ ggtitle("bbknn2; k=9, prune= 1/15, resolution=1.2, all cells")

#Cluster3 = immune
DimPlot(E_seurat2_bb, reduction = "bbknn", group.by = "batch", cells = WhichCells(E_seurat2_bb, idents = c(0,1,2,4))) %>% LabelClusters(id = "batch") + ggtitle("bbknn2; k=9, prune= 1/15, resolution=1.2, all cells")
DimPlot(E_seurat2_bb, reduction = "bbknn", group.by = "ident", cells = WhichCells(E_seurat2_bb, idents = c(0,1,2,4))) %>% LabelClusters(id = "ident") + ggtitle("bbknn2; k=9, prune= 1/15, resolution=1.2, all cells")

#Find Marker gene_1 #그룹간 비교
#group1 cluster와 나머지 cluster 비교
group1 <- c('2')
group2 <- NULL  
marker_dt <- FindMarkers(E_seurat2_bb, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
marker_2<-marker_dt

#Immune cell 제거 
test <- subset(E_seurat2_bb, idents = c("0","1","2","4"))
test_bb <- test %>% FindBatchBalancedNeighbors(batch.key = "batch",
                                               neighbors_within_batch = 10,
                                               prune.SNN = 1/15) 
test_bb <- test_bb %>%  FindClusters(graph.name = "bbknn_snn", resolution = 1.2) 
test_bb <- test_bb %>%  RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_")

DimPlot(test_bb, reduction = "bbknn", group.by = "batch")%>%LabelClusters(id = "batch") + ggtitle("bbknn2; k=10, prune= 1/15, resolution=1.2, no im")
DimPlot(test_bb, reduction = "bbknn", group.by = "ident")%>%LabelClusters(id = "ident")+ ggtitle("bbknn2; k=10, prune= 1/15, resolution=1.8, no im")

E_seurat2_no_im<-test
E_seurat2_no_im_bb<-test_bb
rm(test)
rm(test_bb)

group1 <- c('3')
group2 <- NULL  
marker_dt <- FindMarkers(E_seurat2_no_im_bb, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
marker_3<-marker_dt

#Epithelial cell만 추출
test <- E_seurat2_no_im[,(E_seurat2_no_im@meta.data$seurat_clusters %in% c(0)|E_seurat2_no_im$batch %in% c("E8.25","E9.5to11.5","E14.5.2"))]
unique(test$batch)

#bbknn2 original version
pdf("bbknn2_epi.pdf", 8,6)
for(k in c(7:12)){
  for(prune in c(1,2)){
    print(
      test %>% FindBatchBalancedNeighbors(batch.key = "batch",
                                               neighbors_within_batch = k,
                                               prune.SNN = prune/15) %>%  
        FindClusters(graph.name = "bbknn_snn") %>%
        RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_") %>% 
        DimPlot(reduction = "bbknn", group.by = "batch") %>% 
        LabelClusters(id="batch") + ggtitle(paste0("neighbors=",k,",prune=",prune,"/15",",epi"))
    )
  }
}
dev.off()


test_bb<-test%>%FindBatchBalancedNeighbors(batch.key = "batch",
                                           neighbors_within_batch = 5,
                                           prune.SNN = 1/15)   
test_bb<- test_bb %>% FindClusters(graph.name = "bbknn_snn", resolution = 1.2) 
test_bb<- test_bb %>% RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_") 

DimPlot(test_bb, reduction = "bbknn", group.by = "batch")%>%LabelClusters(id = "batch") + ggtitle("bbknn2; k=5, prune= 1/15, resolution=1.2, epi")

group1 <- c('2')
group2 <- NULL  
marker_dt <- FindMarkers(test_bb, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
marker_2<-marker_dt

#bbknn2 upgrade version
test <- E_seurat2_no_im[,(E_seurat2_no_im@meta.data$seurat_clusters %in% c(0)|E_seurat2_no_im$batch %in% c("E8.25","E9.5to11.5","E14.5.2"))]
unique(test$batch)

pdf("bbknn2_up_epi.pdf", 8,6)
for(k in c(5:13)){
  for(kb in c((k-1)%/%2:(k-1))){
    for(prune in c(1,2)){
      print(
        test %>% FindBatchBalancedNeighbors(batch.key = "batch",
                                          neighbors_within_own = k,
                                          neighbors_within_batch = kb,
                                          prune = prune/15) %>% 
            FindClusters(graph.name = "bbknn_snn") %>%
            RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_") %>% 
            DimPlot(reduction = "bbknn", group.by = "batch") %>% 
            LabelClusters(id="batch") + ggtitle(paste0("bbknn2_up; neighbors=",k,"(",kb,")",",prune=",prune,"/15",",all"))
      )
    }
  }
}
dev.off()

#적당한 조건으로 그려보기 
test_bb_up<-test%>%FindBatchBalancedNeighbors(batch.key = "batch",
                                              neighbors_within_own = 8,
                                              neighbors_within_batch = 4,
                                              prune = 1/15)
test_bb_up<-test_bb_up%>%FindClusters(graph.name = "bbknn_snn", resolution = 0.8)
test_bb_up<-test_bb_up%>%RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_")
DimPlot(test_bb_up, reduction = "bbknn", group.by = "batch") %>%
  LabelClusters(id="batch") + ggtitle("bbknn2_up; neighbors = 8(4), prune = 1/15, resolution = 1.2")
DimPlot(test_bb_up, reduction = "bbknn", group.by = "ident") %>%
  LabelClusters(id="ident") + ggtitle("bbknn2_up; neighbors = 8(4), prune = 1/15, resolution = 0.8")

group1 <- c('6','7')
group2 <- NULL  
marker_dt <- FindMarkers(test_bb_up, ident.1 = group1, ident.2 = group2) %>% as.data.frame %>% rownames_to_column('gene_id') %>% as_tibble()
marker_67<-marker_dt

#bbknn2 upgrade version; re-down sampling & reducing batches
# downsampling → E14.5.1 data만 1k cell로 downsampling
C_1k<-knn(C,10)
sample<-sample(ncol(C_1k),1000)
C_1k<-C_1k[,c(sample)]

C_500<-knn(C,10)
sample<-sampel(ncol(C_500), 500)
C_500<-C_500[,c(sample)]

E2 = cbind(A,B,C_1k,D,G,H,I,K,L,M,N,O,P,Q)
rownames(E2) = union_gene_vector
table(duplicated(t(E2)))
E2 = E2[,!duplicated(t(E2))]
table(duplicated(t(E2)))
batch_vector = stringr::str_replace(colnames(E2), "_.*","")

E_seurat2<- CreateSeuratObject(E2)
E_seurat2$batch = batch_vector
E_seurat2<- E_seurat2%>% NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()
  
E_seurat2_up_bb<- E_seurat2 %>% FindBatchBalancedNeighbors(batch.key = "batch",
                                                           neighbors_within_own = 10,
                                                           neighbors_within_batch = 7,
                                                           prune = 1/15) 
E_seurat2_up_bb <- E_seurat2_up_bb %>%  FindClusters(graph.name = "bbknn_snn", resolution = 1.2) 
E_seurat2_up_bb <- E_seurat2_up_bb %>%  RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_")

DimPlot(E_seurat2_up_bb, reduction = "bbknn", group.by = "batch")%>%LabelClusters(id = "batch") + ggtitle("bbknn2_up;k=10(7), prune= 1/15, resolution = 1.2,all cells")
DimPlot(E_seurat2_up_bb, reduction = "bbknn", group.by = "ident")%>%LabelClusters(id = "ident") + ggtitle("bbknn2_up;k=10(7), prune= 1/15, resolution = 1.2,all cells")

#bbknn2 version
E_seurat2_bb<-E_seurat2%>%FindBatchBalancedNeighbors(batch.key = "batch",
                                           neighbors_within_batch = 5,
                                           prune.SNN = 1/15)   
E_seurat2_bb<-E_seurat2_bb %>% FindClusters(graph.name = "bbknn_snn", resolution = 1.2) 
E_seurat2_bb<-E_seurat2_bb %>% RunUMAP(graph = "bbknn_snn", reduction.name = "bbknn", reduction.key = "BBKNN_") 

DimPlot(E_seurat2_bb, reduction = "bbknn", group.by = "batch")%>%LabelClusters(id = "batch") + ggtitle("bbknn2; k=5, prune= 1/15, resolution=1.2, all")
DimPlot(E_seurat2_bb, reduction = "bbknn", group.by = "ident")%>%LabelClusters(id = "ident") + ggtitle("bbknn2; k=5, prune= 1/15, resolution=1.2, all")
