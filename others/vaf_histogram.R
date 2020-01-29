library(cowplot)
library(ggplot2)
library(dplyr)

#SNV data
dir <- ("/home/users/team_projects/LungCancer_IO_WGS/01_WGS/07_pointMT/01_snv/")
file_list<-list.files("/home/users/team_projects/LungCancer_IO_WGS/01_WGS/07_pointMT/01_snv/",recursive=TRUE, pattern="*.seqzcn$")

#INDEL data
dir_2 <- ("/home/users/team_projects/LungCancer_IO_WGS/01_WGS/07_pointMT/02_indel/")
file_list_2 <- list.files("/home/users/yunah1029/projects/04_lungCa_IO_SNUH/07_pointMT/02_indel/",recursive=TRUE, pattern="*.seqzcn$")

#Plot vaf histogram (luad + lusc combined)
par(mfrow = c(4,6))
par("mar")
par(mar=c(2,2,2,2))

for (n in setdiff(1:length(file_list),11)) {
  print(file_list[n])
  temp_snv <- read.csv(paste0(dir, file_list[n]), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = TRUE)
  temp_snv$CN <- paste0(temp_snv$totCN,"_",temp_snv$majCN,"_",temp_snv$minCN)
  temp_snv <- temp_snv %>% filter(temp_snv$CN == "2_1_1")
  temp_snv$ref_cor<-lapply(strsplit(temp_snv$repeat_unit.ref_repeat_count.lt_seq_size.rt_seq_size.t_ref.t_var.t_unknown.n_other.n_consen.t_ref_cor.t_var_cor.t_ukn_cor,"[;]"),function(x) x[10]) %>% as.numeric()
  temp_snv$var_cor<-lapply(strsplit(temp_snv$repeat_unit.ref_repeat_count.lt_seq_size.rt_seq_size.t_ref.t_var.t_unknown.n_other.n_consen.t_ref_cor.t_var_cor.t_ukn_cor,"[;]"),function(x) x[11]) %>% as.numeric()
  temp_snv$vaf <- temp_snv$var_cor/(temp_snv$ref_cor + temp_snv$var_cor)
  
  temp_indel <- read.csv(paste0(dir_2, file_list_2[n]), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = TRUE)
  temp_indel$CN <- paste0(temp_indel$totCN,"_",temp_indel$majCN,"_",temp_indel$minCN)
  temp_indel <- temp_indel %>% filter(temp_indel$CN == "2_1_1")
  temp_indel$ref_cor<-lapply(strsplit(temp_indel$repeat_unit.ref_repeat_count.lt_seq_size.rt_seq_size.t_ref.t_var.t_unknown.n_other.n_consen.t_ref_cor.t_var_cor.t_ukn_cor,"[;]"),function(x) x[10]) %>% as.numeric()
  temp_indel$var_cor<-lapply(strsplit(temp_indel$repeat_unit.ref_repeat_count.lt_seq_size.rt_seq_size.t_ref.t_var.t_unknown.n_other.n_consen.t_ref_cor.t_var_cor.t_ukn_cor,"[;]"),function(x) x[11]) %>% as.numeric()
  temp_indel$vaf <- temp_indel$var_cor/(temp_indel$ref_cor + temp_indel$var_cor)
  
  sample_name <- strsplit(file_list[n],"[.]")[[1]][1]
  
  temp_snv_indel_vaf <- rbind(temp_indel$vaf%>%as.data.frame(), temp_snv$vaf%>%as.data.frame())
  colnames(temp_snv_indel_vaf) <- "vaf"
  
  hist(temp_snv_indel_vaf$vaf, breaks=seq(0,1,by=0.02), col=rgb(0,0,1,1/4), xlim=c(0,1), xaxt="n",xlab = "Vaf", ylab = "Frequency", main = sample_name)
  axis(1, seq(0,1,0.04))
  
  #p1 <- hist(temp_snv$vaf, breaks=seq(0,1,by=0.02)) 
  #p2 <- hist(temp_indel$vaf, breaks=seq(0,1,by=0.02))                     
  #with(temp_snv, plot(p1, col=rgb(0,0,1,1/4), xlim=c(0,1), xlab = "Vaf", ylab = "Frequency", main = sample_name))
  #par(new = T)
  #with(temp_indel, plot( p2, col=rgb(1,0,0,1/4), axes=F, xlab=NA, ylab=NA, cex=1.2, main = NULL))
  #axis(side = 4)
  #legend("topright", legend=c("SNV","INDEL"), pch = c(15,15), col = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

}

  