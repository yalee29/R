dir="/home/users/yunah1029/projects/smoking_LUAD_mtDNA/03_WES/02_sequenza/output/"

sample_purity <- data.frame()

for (f in list.files(dir, pattern = "*_confints_CP.txt")){
  temp1<-read.csv(paste0(dir,f),sep='\t', stringsAsFactors = FALSE, na.strings = "")
  temp1$sample<-strsplit(f,"[_]")[[1]][1]
  temp1<-temp1[1,]
  sample_purity<-rbind(sample_purity,temp1)
}
  
}

