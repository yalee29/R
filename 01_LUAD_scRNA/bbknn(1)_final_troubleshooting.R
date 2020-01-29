#single cell data processing -> with 600 genes (expression data)
mouse_epi_traj<-E_seurat1_ad_epi_bb@assays$RNA@data # cell; 5903 / gene; 11854
mouse_epi_traj<-mouse_epi_traj%>%as.data.frame()
mouse_epi_traj<-mouse_epi_traj%>%rownames_to_column()
names(mouse_epi_traj)[1] <- c("ensembl")
mouse_epi_600 <- merge(mouse_epi_traj, final, by.x = "ensembl", by.y = "mouse_ensembl")
mouse_epi_600 <- mouse_epi_600[,-1]

#LUAD_bulk_RNAseq with 600 genes (expression data_FPKM)
dir <- ("./projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUAD_551_FPKM/")
file_list<-list.files("/home/users/yunah1029/projects/01_LUAD_scRNA/LUAD/bulkRNAseq/tcga_LUAD_551_FPKM",, recursive=TRUE, pattern="*.FPKM.txt.gz$")

tcga <- mouse_epi_600$human_ensembl %>% as.data.frame()
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


#1. mouse_epi_600 knn pooling
##knn function
knn <- function(A,k){
  i = nrow(A)
  j = ncol(A)
  B = matrix(0,nrow = i, ncol = j)
  D = cor(A)
  for (m in 1:j){
    r = rank(-D[,m],ties.method = "first") <= k+1
    B[,m] = rowSums(A[,r,drop=F])
  }
  colnames(B) <- colnames(A)
  rownames(B) <- rownames(A)
  B
}
##mouse_epi_600 knn
mouse_epi_600_knn <- knn(mouse_epi_600, 5)

##
single_bulk_knn <- cbind(mouse_epi_600_knn, tcga)
cor_single_bulk_knn <- cor(single_bulk_knn, method = "pearson") #method = "spearman"/ pearson(default)
cor_single_bulk_knn <- cor_single_bulk_knn %>% as.data.frame()
cor_single_bulk_sort_knn <- cor_single_bulk_knn[1:5903,5904:nrow(cor_single_bulk_knn)]
cor_max_knn <- apply(cor_single_bulk_sort_knn,2, max) #correlation cutoff 설정 
#cor_single_bulk_sort_knn<- cor_single_bulk_sort_knn[,which(cor_max_knn > 0.45)]
cor_single_bulk_sort_rank_knn <- apply(cor_single_bulk_sort_knn,2,function(x){order(-x)})[1:3,]%>%t()
cor_single_bulk_sort_rank_knn <- cor_single_bulk_sort_rank_knn %>% as.data.frame()

##Plot
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
  geom_point(aes(color=group), size = 0.8, alpha = 0.7)+
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#7dab78", "#f57c16",
                                 "#F0E442", "#bf046a","#fa2a00","#0e5e05","#9e5416")) +
  theme_classic()

##
#a<-Embeddings(mouse_epi_traj,reduction = "bbknn") %>% as.data.frame()
cor_single_bulk_sort_rank_knn <- cor_single_bulk_sort_rank_knn[,1:3] %>% t()
matched <- data.frame()
matched_avg <- data.frame()

for(n in 1:ncol(cor_single_bulk_sort_rank_knn)){
  tmp_matched <- a[cor_single_bulk_sort_rank_knn[,n],]
  tmp_matched$group <- "LUAD:Highly correlated"
  tmp_matched <- tmp_matched %>% rownames_to_column(var = "cell_id")
  
  #tmp_avg <- colMeans(tmp_matched[2:3]) %>% t() %>% as.data.frame()
  #tmp_avg$group <- "LUAD:Average"
  #tmp_avg$cell_id <- colnames(cor_single_bulk_sort_rank_3)[n]
  #tmp_avg <- tmp_avg[,c(4,1,2,3)]
  
  matched<- rbind(matched, tmp_matched)
  #matched_avg <- rbind(matched_avg, tmp_avg)
}

bbknn_dt_plot <- rbind(bbknn_dt, matched)
unique(bbknn_dt_plot$group)
ggplot(bbknn_dt_plot, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), size = 0.8, alpha = 0.7)+
  #scale_shape_manual(values = c( 16, 16, 16,16, 16, 2, 17))+
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#7dab78", "#f57c16",
                                 "#F0E442", "#bf046a","#fa2a00","#0e5e05","#9e5416","#3e0ffa")) #0072B2


##Hightly correlated cells 정리 -> expression pattern 확인
high_cor_cells_knn<-table(cor_single_bulk_sort_rank_knn[1,]) #1/2/3등 조절 가능
high_cor_cells_knn%<>% as.data.frame()
high_cor_cells_knn$cell_name <- colnames(cor_single_bulk_knn)[high_cor_cells_knn$Var1 %>% as.character() %>% as.numeric()]
#high_cor_cells<-high_cor_cells[-1]
high_cor_cells_knn<-high_cor_cells_knn[,c(3,1,2)] #column 순서 조절
head(high_cor_cells_knn)
#high_cor_cells %>% write_csv("high_cor_cells3.csv") #1/2/3조절

high_cor_cells_sort <- single_bulk[,c(4987,5638,5654,5722,5728,5729,5878)]
head(high_cor_cells_sort)
high_cor_cells_sort$tcga_avg <- rowMeans(single_bulk[,5904:ncol(single_bulk)])

high_cor_cells_sort %>% write.csv("high_cor_cells_tcga_expression.csv")

##human lung cluster마다 mouse traj에 붙여보기 (#표시한곳 조건 바꿔가면서 cluster별 그릴 수 있음)  #전과정은 bbknn(1)_final_epi.R 에 정리돼있음. 여기는 knn pooling한 data와 human noraml lung을 비교하기 위함
m_h_AT1 <- cbind(mouse_epi_600_knn, hu_AT1)
cor_m_h_AT1 <- cor(m_h_AT1, method = "pearson") #method = "spearman"/ pearson(default)
cor_m_h_AT1 <- cor_m_h_AT1 %>% as.data.frame()
cor_m_h_AT1_sort <- cor_m_h_AT1[1:5903,5904:nrow(cor_m_h_AT1)]
cor_m_h_AT1_max <- apply(cor_m_h_AT1_sort,2, max) #correlation cutoff 설정 
#cor_single_bulk_sort_0.3 <- cor_m_h_AT1_sort[,which(cor_max > 0.45)]
cor_m_h_AT1_sort_rank <- apply(cor_m_h_AT1_sort,2,function(x){order(-x)})[1:8,]%>%t()
cor_m_h_AT1_sort_rank <- cor_m_h_AT1_sort_rank %>% as.data.frame()

##
#a<-Embeddings(mouse_epi_traj,reduction = "bbknn") %>% as.data.frame()
cor_m_h_AT1_sort_rank <- cor_m_h_AT1_sort_rank[,1:3] %>% t()
matched <- data.frame()

for(n in 1:ncol(cor_m_h_AT1_sort_rank)){
  tmp_matched <- a[cor_m_h_AT1_sort_rank[,n],]
  tmp_matched$group <- "h_Cili" #여기도 cluster 바뀔때마다 조절
  matched<- rbind(matched, tmp_matched)
}

bbknn_dt_plot <- rbind(bbknn_dt, matched)
ggplot(bbknn_dt_plot, aes(BBKNN_1, BBKNN_2))+
  geom_point(aes(color=group), size = 1, alpha = 0.3)+
  #scale_shape_manual(values = c( 16, 16, 16,16, 16, 2, 17))+
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#7dab78", "#f57c16",
                                 "#F0E442", "#bf046a","#fa2a00","#0e5e05","#9e5416","#3e0ffa"))

#2. Perturbation
##600gene random sampling
random_600<-sample(rownames(single_bulk_knn), 613, replace = TRUE)
length(unique(random_600))
single_bulk_random <- single_bulk_knn %>% filter(rownames(single_bulk_knn) %in% random_600) 

##option1; cor with reduced(sampled) gene sets
cor_single_bulk_random <- cor(single_bulk_random, method = "pearson") #method = "spearman"/ pearson(default)
cor_single_bulk_random <- cor_single_bulk_random %>% as.data.frame()
cor_single_bulk_sort_random <- cor_single_bulk_random[1:5903,5904:nrow(cor_single_bulk_random)]
cor_max_random <- apply(cor_single_bulk_sort_random,2, max) #for correlation cutoff 설정 
#cor_single_bulk_sort_knn<- cor_single_bulk_sort_knn[,which(cor_max_knn > 0.45)]
cor_single_bulk_sort_rank_random <- apply(cor_single_bulk_sort_random,2,function(x){order(-x)})[1:7,]%>%t() #1등부터 8등까지
cor_single_bulk_sort_rank_random <- cor_single_bulk_sort_rank_random %>% as.data.frame()

##option2; original gene sets
single_bulk_knn <- cbind(mouse_epi_600_knn, tcga)
cor_single_bulk_knn <- cor(single_bulk_knn, method = "pearson") #method = "spearman"/ pearson(default)
cor_single_bulk_knn <- cor_single_bulk_knn %>% as.data.frame()
cor_single_bulk_sort_knn <- cor_single_bulk_knn[1:5903,5904:nrow(cor_single_bulk_knn)]
cor_max_knn <- apply(cor_single_bulk_sort_knn,2, max) #for correlation cutoff 설정 
#cor_single_bulk_sort_knn<- cor_single_bulk_sort_knn[,which(cor_max_knn > 0.45)]
cor_single_bulk_sort_rank_knn <- apply(cor_single_bulk_sort_knn,2,function(x){order(-x)})[1:7,]%>%t()
cor_single_bulk_sort_rank_knn <- cor_single_bulk_sort_rank_knn %>% as.data.frame()


##convex hull 그리기
#option1 일 경우 아래 그림그릴 때 cor_single_bulk_sort_rank_random 사용
cor_single_bulk_sort_rank_random <- cor_single_bulk_sort_rank_random %>% t()
#option2 일 경우같은 방법으로 cor_single_bulk_sort_rank 사용
cor_single_bulk_sort_rank_knn <- cor_single_bulk_sort_rank_knn%>% t()

pdf("perturbation_sampled_genes.pdf",8,5)
for(n in 1:ncol(cor_single_bulk_sort_rank_random)){
  tmp_matched <- a[cor_single_bulk_sort_rank_random[,n],]
  tmp_matched$group <- "LUAD:Highly correlated"
  tmp_matched <- tmp_matched %>% rownames_to_column(var = "cell_id")
  tcga_sample_name <-colnames(cor_single_bulk_sort_rank_random)[n]
  
  tmp_avg <- apply(tmp_matched[,2:3],2,median) %>% t() %>% as.data.frame()
  tmp_avg$group <- "LUAD:Median"
  tmp_avg$cell_id <- colnames(cor_single_bulk_sort_rank_random)[n]
  tmp_avg <- tmp_avg[,c(4,1,2,3)]
  
  tmp_matched_dist <- tmp_matched
  for ( n in 1:nrow(tmp_matched_dist)) {
    tmp_matched_dist[n,5] <- sqrt((tmp_matched_dist[n,2] - tmp_avg$BBKNN_1)^2 + (tmp_matched_dist[n,3] - tmp_avg$BBKNN_2)^2)
  }
  
  tmp_matched_final <- tmp_matched[tmp_matched$cell_id %in% tmp_matched_dist[which(rank(tmp_matched_dist$V5) < nrow(tmp_matched_dist)*0.9),1],]
  tmp_matched_chull <- chull(tmp_matched_final[,2:3], NULL) #highly correlated 10 cell로 convex hull의 outline 찾는과정
  
  bbknn_dt_tmp <- rbind(bbknn_dt, tmp_matched, tmp_avg)
  print(
    bbknn_dt_tmp %>%
      ggplot(aes(x = BBKNN_1, y = BBKNN_2)) +
      ggtitle(tcga_sample_name)+
      geom_point(aes(color=group), size = 0.8, alpha = 0.7)+
      scale_colour_manual(values = c("#CC79A7", "#E69F00", "#7dab78", "#f57c16",
                                     "#F0E442", "#bf046a","#fa2a00","#0e5e05","#9e5416","#0066CC","#000000")) +
      geom_polygon(data = tmp_matched_final[tmp_matched_chull,], size = 0.5, color = "black", alpha = 1, fill = NA)
  )
}
dev.off()

