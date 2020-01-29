dir="/home/users/yunah1029/projects/smoking_LUAD_mtDNA/02_WES_MT/04_merged/01_annotated/"
threshold <- 0.1

df_merged <- data.frame(POS=NA, REF=NA, ALT=NA)

for(f in list.files(dir, pattern="_anv.tsv$")) {
  
  df_temp <- read.csv(paste0(dir, f), sep='\t', stringsAsFactors = FALSE)
  df_temp$VAF <- df_temp$TUMOR_AD_VarScan2 / (df_temp$TUMOR_AD_VarScan2 + df_temp$TUMOR_RD_VarScan2) 
  
  # filtering option...
  df_temp <- subset(df_temp, VAF >= threshold)
  
  
  df_merged <- merge(df_merged, df_temp[,c("POS", "REF", "ALT", "VAF")], by=c("POS", "REF", "ALT"), all=TRUE)
  colnames(df_merged)[ncol(df_merged)] <-  strsplit(f, "[_]")[[1]][1]
  
  
}

df_merged[is.na(df_merged)] <- 0
