gse69405_bulk <- read_csv("./projects/01_LUAD_scRNA/LUAD/GSE69405_bulk_single/GSE69405_processed_tpm_bulks.csv")
gse69405_bulk <- gse69405_bulk[,-c(1,3)]

#bulk
human_ensmebl_genename<-read_tsv("./projects/01_LUAD_scRNA/donor_single_cell/human_ENSEMBL_genename.tsv")
gse69405_bulk<-gse69405_bulk[!duplicated(gse69405_bulk$gene_name),]
gse69405_bulk <- merge(gse69405_bulk, human_ensmebl_genename, by.x = "gene_name", by.y = "Gene name")
gse69405_bulk<- gse69405_bulk[,-1]

mouse_epi_600_knn %<>% rownames_to_column(var = "human_ensembl")
gse69405 <- mouse_epi_600_knn$human_ensembl %>% as.data.frame()
gse69405$id  <- 1:nrow(gse69405)
colnames(gse69405) <- c("gene","id")
gse69405_bulk<- merge(gse69405, gse69405_bulk, by.x = "gene", by.y = "Gene stable ID", all.x = TRUE, all.y = FALSE)
gse69405_bulk<-gse69405_bulk[order(gse69405_bulk$id), ]
rownames(gse69405_bulk) <- NULL
gse69405_bulk %<>% column_to_rownames(var = "gene")
gse69405_bulk <- gse69405_bulk[,-1]
mouse_epi_600_knn %<>% column_to_rownames("human_ensembl")

#bulk corsingle_bulk_gse <- cbind(mouse_epi_600_knn, gse69405_bulk)
single_bulk_gse[is.na(single_bulk_gse)] <- 0
cor_single_bulk_gse <- cor(single_bulk_gse, method = "pearson")
cor_single_bulk_sort_gse <- cor_single_bulk_gse[1:5903, 5904:ncol(cor_single_bulk_gse)]
cor_max_gse <- apply(cor_single_bulk_sort_gse, 2, max)
cor_single_bulk_sort_gse_rank <- apply(cor_single_bulk_sort_gse,2,function(x){order(-x)})[1:7,]%>%t()

# original gene set으로 correlation 볼 때 각 sclc_n bulk 별 Median
original_matched <- data.frame()
original_avg <- data.frame()

for (x in 1:nrow(cor_single_bulk_sort_gse_rank)) {
  tmp_original_matched <- a[cor_single_bulk_sort_gse_rank[x,],]
  tmp_original_matched$group <- "scLCMBT15:Origianl genes_Median"
  tmp_original_matched$cancer_sample <- rownames(cor_single_bulk_sort_gse_rank)[x]
  tmp_original_matched$chull <- x
  
  tmp_original_avg <- apply(tmp_original_matched[,2:3],2,median) %>% t() %>% as.data.frame()
  tmp_original_avg$group <- "scLCMBT15:Origianl genes_Median"
  tmp_original_avg$cell_id <- rownames(cor_single_bulk_sort_gse_rank)[x]
  tmp_original_avg$cancer_sample <- rownames(cor_single_bulk_sort_gse_rank)[x]
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
  random_600<-sample(rownames(single_bulk_gse), 591, replace = TRUE)
  single_bulk_random <- single_bulk_gse %>% filter(rownames(single_bulk_gse) %in% random_600) 
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
    tmp_matched$group <- "scLCMBT15:Highly correlated"
    tmp_matched$cancer_sample <-colnames(cor_single_bulk_sort_rank_random)[n]
    tmp_matched$chull_group <- paste0(i,"_",n)
    
    tmp_avg <- apply(tmp_matched[,2:3],2,median) %>% t() %>% as.data.frame()
    tmp_avg$group <- "scLCMBT15:Median"
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

pdf("perturbation(gse69405_scLCMBT15)_sampled_genes2_v3_chull_medians.pdf",8,5)
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
                                     "#F0E442", "#bf046a","#fa2a00","#0e5e05","#9e5416","#0066CC","#000000"))+ # 뒤에서 두번째 "#0066CC" 
      geom_polygon(data = hull3, size = 0.5, color = "#0066CC", alpha = 0.5, fill = NA)
  )
}
dev.off()








