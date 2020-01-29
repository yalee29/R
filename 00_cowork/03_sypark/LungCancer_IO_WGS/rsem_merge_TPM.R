library(readr)
library(tidyverse)
library(dplyr)

#RNAseq sample별 TPM merge-------------------------------------------------------------------
dir <- ("/home/users/team_projects/LungCancer_IO_WGS/02_RNAseq/HN00103479/RSEM_result/")
file_list<-list.files("/home/users/team_projects/LungCancer_IO_WGS/02_RNAseq/HN00103479/RSEM_result",, recursive=TRUE, pattern="*.genes.results$")
tpm_matrix <- read.csv(paste0(dir, file_list[1]), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = TRUE) %>% select(gene_id, TPM)
colnames(tpm_matrix) <- c("gene_id","ref")
for (file in file_list) {
  temp <- read.csv(paste0(dir, file), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = TRUE)
  temp <- temp %>% select(gene_id,TPM)
  colnames(temp) <- c("gene_id", paste0(strsplit(file[1], "_")[[1]][1],"_",strsplit(file[1], "_")[[1]][2]))
  #temp$gene <- sapply(strsplit(temp$gene, "[.]"), function(x) x[1])
  tpm_matrix <- merge(tpm_matrix, temp, by.x = "gene_id", by.y = "gene_id", all.x = TRUE, all.y = FALSE)
}
dim(tpm_matrix)
sum(is.na(tpm_matrix))
head(tpm_matrix)
tpm_matrix <- tpm_matrix[,-2] #ref 제거

#RNAseq sample 별 count merge-------------------------------------------------------------------
count_matrix <- read.csv(paste0(dir, file_list[1]), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = TRUE) %>% select(gene_id, expected_count)
colnames(count_matrix) <- c("gene_id","ref")
for (file in file_list) {
  temp <- read.csv(paste0(dir, file), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = TRUE)
  temp <- temp %>% select(gene_id,expected_count)
  colnames(temp) <- c("gene_id", paste0(strsplit(file[1], "_")[[1]][1],"_",strsplit(file[1], "_")[[1]][2]))
  #temp$gene <- sapply(strsplit(temp$gene, "[.]"), function(x) x[1])
  count_matrix <- merge(count_matrix, temp, by.x = "gene_id", by.y = "gene_id", all.x = TRUE, all.y = FALSE)
}
dim(count_matrix)
sum(is.na(count_matrix))
head(count_matrix)
count_matrix <- count_matrix[,-2] #ref 제거

#Extract variable genes--------------------------------------------------------------------------
