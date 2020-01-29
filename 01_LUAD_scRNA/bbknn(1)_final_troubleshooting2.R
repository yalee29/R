#
single_bulk_knn <- cbind(mouse_epi_600_knn, tcga)
cor_single_bulk_knn <- cor(single_bulk_knn, method = "pearson") #method = "spearman"/ pearson(default)
cor_single_bulk_knn <- cor_single_bulk_knn %>% as.data.frame()
cor_single_bulk_sort_knn <- cor_single_bulk_knn[1:5903,5904:nrow(cor_single_bulk_knn)]
cor_max_knn <- apply(cor_single_bulk_sort_knn,2, max) #correlation cutoff 설정 
#cor_single_bulk_sort_knn <- cor_single_bulk_sort_knn[,which(cor_max_knn > 0.45)]
cor_single_bulk_sort_rank_knn <- apply(cor_single_bulk_sort_knn,2,function(x){order(-x)})[1:7,]%>%t()
#cor_single_bulk_sort_rank_knn <- cor_single_bulk_sort_rank_knn %>% as.data.frame()

# original gene set으로 correlation 볼 때 각 tcga sample 별 Median
original_matched <- data.frame()
original_avg <- data.frame()

for (x in 1:nrow(cor_single_bulk_sort_rank_knn)) {
  tmp_original_matched <- a[cor_single_bulk_sort_rank_knn[x,],]
  tmp_original_matched$group <- "LUAD:Origianl genes_Median"
  tmp_original_matched$cancer_sample <- rownames(cor_single_bulk_sort_rank_knn)[x]
  tmp_original_matched$chull <- x
  
  tmp_original_avg <- apply(tmp_original_matched[,2:3],2,median) %>% t() %>% as.data.frame()
  tmp_original_avg$group <- "LUAD:Origianl genes_Median"
  tmp_original_avg$cell_id <- rownames(cor_single_bulk_sort_rank_knn)[x]
  tmp_original_avg$cancer_sample <- rownames(cor_single_bulk_sort_rank_knn)[x]
  tmp_original_avg <- tmp_original_avg[,c(4,1,2,3,5)]
  tmp_original_avg$chull_group <- x
  
  original_matched <- rbind(original_matched, tmp_original_matched)
  original_avg <- rbind(original_avg, tmp_original_avg)
}

dim(original_matched)
dim(original_avg)

##random sampling + convex hull 그리기
i = 1
matched_final <- data.frame()
avg_final <- data.frame()

repeat{
  random_600<-sample(rownames(single_bulk_knn), 591, replace = TRUE)
  single_bulk_random <- single_bulk_knn %>% filter(rownames(single_bulk_knn) %in% random_600) 
  cor_single_bulk_random <- cor(single_bulk_random, method = "pearson") %>% as.data.frame() #method = "spearman"/ pearson(default)
  cor_single_bulk_sort_random <- cor_single_bulk_random[1:5903,5904:nrow(cor_single_bulk_random)]
  #cor_max_random <- apply(cor_single_bulk_sort_random,2, max) #for correlation cutoff 설정 
  #cor_single_bulk_sort_random <- cor_single_bulk_sort_random[,which(cor_max_random > median(cor_max_random))]
  #cor_single_bulk_sort_random<- cor_single_bulk_sort_random[,which(cor_max_knn > 0.4)] #original gene set에서 cor_max가 0.4(median) 보다 큰 sample만 -> 223/551
  cor_single_bulk_sort_rank_random <- apply(cor_single_bulk_sort_random,2,function(x){order(-x)})[1:11,]%>%t() #1등부터 7등까지
  cor_single_bulk_sort_rank_random %<>% as.data.frame()
  cor_single_bulk_sort_rank_random <- cor_single_bulk_sort_rank_random %>% t()
  
  matched <- data.frame()
  avg <- data.frame()
  
  for (n in 1:ncol(cor_single_bulk_sort_rank_random)){
    tmp_matched <- a[cor_single_bulk_sort_rank_random[,n],]
    tmp_matched$group <- "LUAD:Highly correlated"
    #tmp_matched <- tmp_matched %>% rownames_to_column(var = "cell_id")
    tmp_matched$cancer_sample <-colnames(cor_single_bulk_sort_rank_random)[n]
    tmp_matched$chull_group <- paste0(i,"_",n)
    
    tmp_avg <- apply(tmp_matched[,2:3],2,median) %>% t() %>% as.data.frame()
    tmp_avg$group <- "LUAD:Median"
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

# Highly correlated cell 로 convex hull 
hull <- matched_final %>% group_by(chull_group) %>% slice(chull(BBKNN_1, BBKNN_2)) #highly correlated cell 로 convex hull 그릴 때 

pdf("perturbation_sampled_genes2_v1_dots_0.45.pdf",8,5)

for (k in 1:length(unique(matched_final$cancer_sample))) {
  bbknn_dt_tmp2 <- rbind(bbknn_dt_tmp,avg_final%>%filter(cancer_sample == unique(matched_final$cancer_sample)[k])) #matched_final %>% filter(cancer_sample == unique(matched_final$cancer_sample)[k])
  hull2<- hull %>% filter(cancer_sample == unique(matched_final$cancer_sample)[k])
  tcga_sample_name <- unique(matched_final$cancer_sample)[k]
  print(
    bbknn_dt_tmp2 %>%
      ggplot(aes(x = BBKNN_1, y = BBKNN_2)) +
      ggtitle(tcga_sample_name)+
      geom_point(aes(color=group), size = 0.8, alpha = 0.7)+
      scale_colour_manual(values = c("#CC79A7", "#E69F00", "#7dab78", "#f57c16",
                                     "#F0E442", "#bf046a","#fa2a00","#0e5e05","#9e5416","#000000")) # 뒤에서 두번째 "#0066CC" 
    #    geom_polygon(data = hull2, size = 0.5, color = "#0066CC", alpha = 0.5, fill = NA)
  )
}
dev.off()


# Median spot으로 convex hull/ sampled gene set Med <-> original gene set Med => convex hull 그릴 준비. 
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

pdf("perturbation_sampled_genes2_v3_chull_medians_0.45.pdf",8,5)
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


# bronchioalveolar stem cells (BASCs) 확인
FeaturePlot(E_seurat1_ad_epi_bb, c("ENSMUSG00000056370","ENSMUSG00000038791","ENSMUSG00000020044"), reduction = "bbknn")

#Doublet check on mouse traj
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
library(DoubletFinder)
sweep.res.E_seurat1_ad_epi_bb <- paramSweep_v3(E_seurat1_ad_epi_bb, PCs = 1:20, sct = FALSE)
sweep.stats_E_seurat1_ad_epi_bb <- summarizeSweep(sweep.res.E_seurat1_ad_epi_bb, GT = FALSE)
bcmvn_E_seurat1_ad_epi_bb <- find.pK(sweep.stats_E_seurat1_ad_epi_bb)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotation <- E_seurat1_ad_epi_bb@active.ident
homotypic.prop <- modelHomotypic(annotation)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*length(colnames(E_seurat1_ad_epi_bb)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
E_seurat1_ad_epi_bb <- doubletFinder_v3(E_seurat1_ad_epi_bb, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
E_seurat1_ad_epi_bb <- doubletFinder_v3(E_seurat1_ad_epi_bb, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(E_seurat1_ad_epi_bb, reduction = "bbknn", group.by = "DF.classifications_0.25_0.09_443")

E_seurat1_ad_epi_bb_single <- E_seurat1_ad_epi_bb[,E_seurat1_ad_epi_bb@meta.data$DF.classifications_0.25_0.09_443 == "Singlet"]

E_seurat1_ad_epi_bb_single<- E_seurat1_ad_epi_bb_single %>% RunBBKNN(dims.use = 1:20,
                                                                     neighbors_within_batch =  4,
                                                                     trim=45,
                                                                     batch.key = "batch",
                                                                     python.path = "/home/users/yunah1029/anaconda3/bin/python") 
DimPlot(E_seurat1_ad_epi_bb_single, reduction = "bbknn", group.by = "batch") %>% LabelC
lusters(id = "batch") + ggtitle("bbknn1; dims=20, k=4, trim=45, epi")
DimPlot(E_seurat1_ad_epi_bb_single, reduction = "bbknn", group.by = "ident") %>% LabelClusters(id = "ident") + ggtitle("bbknn1; dims=20, k=4, trim=45, epi")
