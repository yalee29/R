library(readxl)

#load purity data (sequenza & vaf purity info)
LungCancerIO_summary <- read_excel("projects/04_lungCa_IO_SNUH/LungCancerIO_summary.xlsx")
LungCancerIO_summary <- LungCancerIO_summary %>% select(sample_id, Seqz_purity, vaf_peak, vaf_cellularity)
LungCancerIO_summary$vaf_peak[22] <- 0.28
LungCancerIO_summary$vaf_peak[23] <- 0.07
LungCancerIO_summary$vaf_cellularity[23] <- 0.14
LungCancerIO_summary$vaf_cellularity[22] <- 0.56
LungCancerIO_summary$vaf_peak<- gsub("NA", 0, LungCancerIO_summary$vaf_peak)
LungCancerIO_summary$vaf_cellularity<- gsub("NA", 0, LungCancerIO_summary$vaf_cellularity)
LungCancerIO_summary$Seqz_purity<- LungCancerIO_summary$Seqz_purity %>% as.numeric()
LungCancerIO_summary$vaf_cellularity<- LungCancerIO_summary$vaf_cellularity %>% as.numeric()
LungCancerIO_summary$diff_seqz_vaf <- LungCancerIO_summary$Seqz_purity - LungCancerIO_summary$vaf_cellularity
LungCancerIO_summary <- LungCancerIO_summary [-11,] #sample11 is removed

#sequenza data
dir <- ("/home/users/team_projects/LungCancer_IO_WGS/01_WGS/08_sequenza/result/")
file_list <- list.files("/home/users/team_projects/LungCancer_IO_WGS/01_WGS/08_sequenza/result/",recursive=TRUE, pattern="*segments.txt$")
file_list
file_list <- setdiff(file_list, file_list[11:13]) #sample11 related files are removed

#Plotting diploid_proportion <-> diff_purity(seqz, vaf)
purity_ploidy <- purity_ploidy %>% as.data.frame()
for (n in 1:length(file_list)) {
  temp <- read.csv(paste0(dir, file_list[n]), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = TRUE)
  temp$diploid <- paste0(temp$CNt,"_",temp$A,"_",temp$B)
  temp$seg_length <- temp$end.pos - temp$start.pos
  
  temp_purity_ploidy <- temp_purity_ploidy %>% as.data.frame()
  colnames(temp_purity_ploidy) <- c("sample_id","diploid_proportion","diff_seqz_vaf")
  temp_purity_ploidy[1,1] <- LungCancerIO_summary$sample_id[n]
  temp_purity_ploidy[1,2] <- sum(temp%>%filter(diploid == "2_1_1")%>%select(seg_length)) / sum(temp$seg_length)
  temp_purity_ploidy[1,3] <- LungCancerIO_summary$diff_seqz_vaf[n]
  
  purity_ploidy <- rbind(purity_ploidy, temp_purity_ploidy)
}

purity_ploidy <- transform(purity_ploidy, sample = substr(sample_id, 11, 12))
ggplot(purity_ploidy, aes(x=diploid_proportion,y=diff_seqz_vaf, color = sample, label = purity_ploidy$sample)) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 3) +
  geom_text(hjust = 0.02, nudge_x = 0.01) +
  theme_minimal() 
  #geom_label(check_overlap = TRUE)

#data export
purity_ploidy_sum <- cbind(LungCancerIO_summary, purity_ploidy[,2:3])
purity_ploidy_sum <- purity_ploidy_sum[,-ncol(purity_ploidy_sum)]
write.table(purity_ploidy_sum, file = "/home/users/yunah1029/projects/04_lungCa_IO_SNUH/purity_ploidy_info.tsv ", row.names=FALSE, sep="\t")
