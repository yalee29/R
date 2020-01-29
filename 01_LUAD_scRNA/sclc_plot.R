
#Preprocessing ###########################################################################################################################################################
sclc_n <- read_tsv("./projects/01_LUAD_scRNA/LUAD/bulkRNAseq/GSE60052_SCLC/GSE60052_7normal.normalized.log2.data.Rda.tsv")
sclc_t <- read_tsv("./projects/01_LUAD_scRNA/LUAD/bulkRNAseq/GSE60052_SCLC/GSE60052_79tumor.normalized.log2.data.Rda 2.tsv")
sclc_n <- sclc_t

#sclc_n
 #human_ensmebl_genename<-read_tsv("./projects/01_LUAD_scRNA/donor_single_cell/human_ENSEMBL_genename.tsv")
sclc_n <- merge(sclc_n, human_ensmebl_genename, by.x = "gene", by.y = "Gene name")
table(duplicated(sclc_n[ncol(sclc_n)]))
sclc_n<- sclc_n[,-1]

mouse_epi_700_knn <-mouse_epi_700_knn %>% rownames_to_column(var = "human_ensembl")
sclc <- mouse_epi_700_knn$human_ensembl %>% as.data.frame()
sclc$id  <- 1:nrow(sclc)
colnames(sclc) <- c("gene","id")
sclc_n<- merge(sclc, sclc_n, by.x = "gene", by.y = "Gene stable ID", all.x = TRUE, all.y = FALSE)
sclc_n<-sclc_n[order(sclc_n$id), ]
rownames(sclc_n) <- NULL
sclc_n <- sclc_n %>% column_to_rownames(var = "gene")
sclc_n <- sclc_n[,-1]
mouse_epi_700_knn <-mouse_epi_700_knn %>% column_to_rownames("human_ensembl")

#test.R 에 사용하기전 -> 바로 "###tcga sample 별 mouse_epi_700과 cor -> 가장 가까운 cell 찾기." 부터 시작ㄱㄱ
tcga<-sclc_n
sum(is.na(tcga))
tcga[is.na(tcga)] <- 0
sum(is.na(tcga))
############################################################################################################################################################################

#sclc_n cor
single_bulk_sclc_n <- cbind(mouse_epi_700_knn, sclc_n)
single_bulk_sclc_n[is.na(single_bulk_sclc_n)] <- 0
cor_single_bulk_sclc_n <- cor(single_bulk_sclc_n, method = "pearson")
cor_single_bulk_sort_sclc_n <- cor_single_bulk_sclc_n[1:5903, 5904:ncol(cor_single_bulk_sclc_n)]
cor_max_sclc_n <- apply(cor_single_bulk_sort_sclc_n, 2, max)
cor_single_bulk_sort_sclc_n_rank <- apply(cor_single_bulk_sort_sclc_n,2,function(x){order(-x)})[1:7,]%>%t()

# original gene set으로 correlation 볼 때 각 sclc_n sample 별 Median
original_matched <- data.frame()
original_avg <- data.frame()

for (x in 1:nrow(cor_single_bulk_sort_sclc_n_rank)) {
  tmp_original_matched <- a[cor_single_bulk_sort_sclc_n_rank[x,],]
  tmp_original_matched$group <- "SCLC_t:Origianl genes_Median"
  tmp_original_matched$cancer_sample <- rownames(cor_single_bulk_sort_sclc_n_rank)[x]
  tmp_original_matched$chull <- x
  
  tmp_original_avg <- apply(tmp_original_matched[,2:3],2,median) %>% t() %>% as.data.frame()
  tmp_original_avg$group <- "SCLC_t:Origianl genes_Median"
  tmp_original_avg$cell_id <- rownames(cor_single_bulk_sort_sclc_n_rank)[x]
  tmp_original_avg$cancer_sample <- rownames(cor_single_bulk_sort_sclc_n_rank)[x]
  tmp_original_avg <- tmp_original_avg[,c(4,1,2,3,5)]
  tmp_original_avg$chull_group <- x
  
  original_matched <- rbind(original_matched, tmp_original_matched)
  original_avg <- rbind(original_avg, tmp_original_avg)
}

dim(original_avg)
head(original_avg)

##100번
i = 1
matched_final <- data.frame()
avg_final <- data.frame()

repeat{
  random_600<-sample(rownames(single_bulk_sclc_n), 591, replace = TRUE)
  single_bulk_random <- single_bulk_sclc_n %>% filter(rownames(single_bulk_sclc_n) %in% random_600) 
  cor_single_bulk_random <- cor(single_bulk_random, method = "pearson") %>% as.data.frame() #method = "spearman"/ pearson(default)
  cor_single_bulk_sort_random <- cor_single_bulk_random[1:5903,5904:nrow(cor_single_bulk_random)]
  #cor_max_random <- apply(cor_single_bulk_sort_random,2, max) #for correlation cutoff 설정 
  #cor_single_bulk_sort_random <- cor_single_bulk_sort_random[,which(cor_max_random > median(cor_max_random))]
  #cor_single_bulk_sort_random<- cor_single_bulk_sort_random[,which(cor_max_knn > 0.4)] #original gene set에서 cor_max가 0.4(median) 보다 큰 sample만 -> 223/551
  cor_single_bulk_sort_rank_random <- apply(cor_single_bulk_sort_random,2,function(x){order(-x)})[1:7,]%>%t() #1등부터 7등까지
  cor_single_bulk_sort_rank_random %<>% as.data.frame()
  cor_single_bulk_sort_rank_random <- cor_single_bulk_sort_rank_random %>% t()
  
  matched <- data.frame()
  avg <- data.frame()
  
  for (n in 1:ncol(cor_single_bulk_sort_rank_random)){
    tmp_matched <- a[cor_single_bulk_sort_rank_random[,n],]
    tmp_matched$group <- "SCLC:Highly correlated"
    tmp_matched$cancer_sample <-colnames(cor_single_bulk_sort_rank_random)[n]
    tmp_matched$chull_group <- paste0(i,"_",n)
    
    tmp_avg <- apply(tmp_matched[,2:3],2,median) %>% t() %>% as.data.frame()
    tmp_avg$group <- "SCLC:Median"
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



## Median spot으로 convex hull/ sampled gene set Med <-> original gene set Med => convex hull 그릴 준비. 
# avg_final = 100번 random sampling 했을 때 각 gene set에 대한 median(top 7)/ original_avg = original gene set에서 cor을 계산했을 때 top 7에 대한 Median
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

pdf("perturbation(SCLC_t)_sampled_genes2_v3_chull_medians.pdf",8,5)
for (k in 1:length(unique(original_avg$cancer_sample))) {
  bbknn_dt_tmp2 <- rbind(bbknn_dt_tmp,med_dist_final%>% filter(cancer_sample == unique(original_avg$cancer_sample)[k]),original_avg%>%filter(cancer_sample == unique(original_avg$cancer_sample)[k])) #med_dist_final%>% filter(cancer_sample == unique(original_avg$cancer_sample)[k]
  hull3<- hull2 %>% filter(cancer_sample == unique(original_avg$cancer_sample)[k])
  tcga_sample_name <- unique(original_avg$cancer_sample)[k]
  print(
    bbknn_dt_tmp2 %>%
      ggplot(aes(x = BBKNN_1, y = BBKNN_2)) +
      ggtitle(tcga_sample_name)+
      geom_point(aes(color=group), size = 0.8, alpha = 0.7)+
      scale_colour_manual(values = c("#CC79A7", "#E69F00", "#7dab78", "#f57c16",
                                     "#F0E442", "#bf046a","#fa2a00","#0e5e05","#9e5416","#000000","#0066CC"))+ # 뒤에서 두번째 "#0066CC" 
      geom_polygon(data = hull3, size = 0.5, color = "#0066CC", alpha = 0.5, fill = NA)
  )
}
dev.off()
