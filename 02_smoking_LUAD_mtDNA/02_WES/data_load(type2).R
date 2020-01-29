#WES variant calling 결과 Varscan은 모든 sample에 대해 variant를 찾지만 Mutect2는 아님 -> column 수가 안맞음. 

dir="/home/users/yunah1029/projects/smoking_LUAD_mtDNA/02_WES_MT/04_merged/01_annotated/"
threshold <- 0.1

#df_merged <- data.frame(ncol=length(default_col_names), nrow=0)
#colnames(df_merged) <- default_col_names
#print(colnames(df_merged))
df_merged <- setNames(data.frame(matrix(ncol = 84, nrow = 0)), default_col_names)

for(f in list.files(dir, pattern="_anv.tsv$")) {
  
  df_temp <- read.csv(paste0(dir, f), sep='\t', stringsAsFactors = FALSE)

  if(ncol(df_temp)==83) {
    
      df_temp$sample <- strsplit(f, "[_]")[[1]][1]
      df_temp <- df_temp[, default_col_names]
      df_merged <- rbind(df_merged, df_temp)
      
  } else {
    
    df_temp$sample <- strsplit(f, "[_]")[[1]][1]
    
    diff <- setdiff(default_col_names, colnames(df_temp))
    for(d in diff) df_temp[,d] <- NA
    
    df_temp <- df_temp[, default_col_names]
    df_merged <- rbind(df_merged, df_temp)
    
  }
  
}

#information data
LungWES_info <- read_excel("projects/smoking_LUAD_mtDNA/02_WES_MT/LungWES_info.xlsx")
names(LungWES_info)<-LungWES_info[1,]
LungWES_info<-LungWES_info[-1,]

Sample<-t(data.frame(strsplit(LungWES_info$`Patient ID`, split = "[-]")))[,3]
LungWES_info_test<-cbind(Sample, LungWES_info)
LungWES_info_test<-LungWES_info_test[,-2]
LungWES_info_test$purity<-sample_purity$cellularity
LungWES_info_test<-LungWES_info_test[,-15]
rownames(LungWES_info_test)<-c()
LungWES_info<-LungWES_info_test

#variant calling data + information data 정리 -> wes_sorted
wes_sorted<-df_merged%>%select(sample,POS,REF,ALT,NORMAL_AF_Mutect2,TUMOR_AF_Mutect2,NORMAL_AD_VarScan2,NORMAL_RD_VarScan2,TUMOR_AD_VarScan2,TUMOR_RD_VarScan2,NORMAL_DP_Mutect2,TUMOR_DP_Mutect2,NORMAL_DP_VarScan2,TUMOR_DP_VarScan2,SS_VarScan2,Func_ensGene, Gene_ensGene, GeneDetail_ensGene, ExonicFunc_ensGene, AAChange_ensGene,NORMAL_FREQ_VarScan2,TUMOR_FREQ_VarScan2)
wes_sorted <- wes_sorted %>% mutate(NORMAL_AF_VarScan2 = wes_sorted$NORMAL_AD_VarScan2/(wes_sorted$NORMAL_AD_VarScan2+wes_sorted$NORMAL_RD_VarScan2))
wes_sorted <- wes_sorted %>% mutate(TUMOR_AF_VarScan2 = wes_sorted$TUMOR_AD_VarScan2/(wes_sorted$TUMOR_AD_VarScan2+wes_sorted$TUMOR_RD_VarScan2))
wes_sorted <- wes_sorted%>%select(sample, POS, REF, ALT, NORMAL_AF_Mutect2, TUMOR_AF_Mutect2, NORMAL_AF_VarScan2, TUMOR_AF_VarScan2,NORMAL_DP_Mutect2,TUMOR_DP_Mutect2,NORMAL_DP_VarScan2,TUMOR_DP_VarScan2,SS_VarScan2, Func_ensGene, Gene_ensGene, ExonicFunc_ensGene)
wes_sorted_test <- merge(LungWES_info, wes_sorted, by.x = "Sample" , by.y = "sample", all.x = FALSE, all.y = TRUE)
wes_sorted<-wes_sorted_test%>%select(Sample, POS, REF, ALT, NORMAL_AF_Mutect2, TUMOR_AF_Mutect2, NORMAL_AF_VarScan2, TUMOR_AF_VarScan2, NORMAL_DP_Mutect2,TUMOR_DP_Mutect2,NORMAL_DP_VarScan2,TUMOR_DP_VarScan2, SS_VarScan2, purity ,Func_ensGene, Gene_ensGene, ExonicFunc_ensGene, Gender, `Age at Surgery`,`Smoking Status`, `Tumor Stage`, `Lymph node or distant Metastasis`)










