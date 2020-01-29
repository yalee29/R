#Data load & make data frame(sample & variant number)
l <- list()

for(f in list.files("./projects/smoking_LUAD_mtDNA/04_Merged/")) {
  df_temp <- read.csv(paste0("./projects/smoking_LUAD_mtDNA/04_Merged/",f), sep='\t', stringsAsFactors = FALSE, na.strings = "")
  l[[length(l) + 1]] <- list(sample=f, count=nrow(df_temp))
}

library(plyr)
df_temp<-ldply(l,data.frame)

#Edit data(sample name 마다 붙은 .tsv를 없애기)
df_temp$sample<-as.character(df_temp$sample)
for(i in 1:nrow(df_temp)){
  df_temp[i,1]<-strsplit(df_temp[i,1],"[.]")[[1]][1]
}

#Data(information) load & Edit data (sample name에서 .merged 없애기)
somatic_mut<-read.csv("./projects/smoking_LUAD_mtDNA/LungWGSinfo.csv", stringsAsFactors = FALSE, na.strings = "")

for(i in 1:nrow(somatic_mut)){
   somatic_mut[i,1]<-strsplit(somatic_mut[i,1],"[.]")[[1]][1]
}

#information data + count 
test<-merge(df_temp,somatic_mut,by.x = "sample",by.y = "ID",all.x = TRUE, all.y = FALSE)


#.tsv data 모두 merge & 필요한 column만 
library(dplyr)
dir<-("../smoking_LUAD_mtDNA/04_Merged/")
file_list<-list.files(dir)
mt_var_rbind<-data.frame()

for (file in file_list){
  temp<-read.csv(paste0(dir,file),sep='\t', stringsAsFactors = FALSE, na.strings = "")
  temp$sample<-strsplit(file,"[.]")[[1]][1]
  mt_var_rbind<-rbind(mt_var_rbind,temp)
}

mt_var_rbind_sort<-mt_var_rbind[,c(79,1,2,3,4,26,62,24,60,53,27,63,30,66)]

#filtering
#rm frequent false-positive variant
mt_var_filt<- mt_var_rbind_sort 
mt_var_filt$NORMAL_FREQ_VarScan2 <- mt_var_filt$NORMAL_FREQ_VarScan2 %>% gsub('%', '', .) %>% as.numeric()
mt_var_filt$TUMOR_FREQ_VarScan2 <- mt_var_filt$TUMOR_FREQ_VarScan2 %>% gsub('%', '', .) %>% as.numeric()

mt_var_filt_test<-mt_var_filt%>% filter(!(POS %in% c(302:315, 513:525, 3105:3109)))
mt_var_filt_test<-mt_var_filt_test[!grepl(",",mt_var_filt_test$ALT),]
mt_var_filt_test$TUMOR_AF_Mutect2<-mt_var_filt_test$TUMOR_AF_Mutect2%>%as.numeric()
mt_var_filt_test$NORMAL_AF_Mutect2<-mt_var_filt_test$NORMAL_AF_Mutect2%>%as.numeric()
mt_var_filt_test$NORMAL_AF_Mutect2<-mt_var_filt_test$NORMAL_AF_Mutect2*100
mt_var_filt_test$TUMOR_AF_Mutect2<-mt_var_filt_test$TUMOR_AF_Mutect2*100

mt_var_filt_test<-mt_var_filt_test%>%filter(mt_var_filt_test$SS_VarScan2==2 | is.na(mt_var_filt_test$TUMOR_FREQ_VarScan2)==TRUE) #varscan에서 somatic mutation만 

mt_var_filt_test<-mt_var_filt_test[mt_var_filt_test$TUMOR_AF_Mutect2>=1 | is.na(mt_var_filt_test$TUMOR_AF_Mutect2)==TRUE,]
mt_var_filt_test<-mt_var_filt_test[mt_var_filt_test$TUMOR_FREQ_VarScan2>=1| is.na(mt_var_filt_test$TUMOR_FREQ_VarScan2)==TRUE,]

mt_var_both<-filter(mt_var_filt_test,is.na(mt_var_filt_test$NORMAL_AF_Mutect2)==FALSE,is.na(mt_var_filt_test$TUMOR_AF_Mutect2)==FALSE, is.na(mt_var_filt_test$NORMAL_FREQ_VarScan2)==FALSE, is.na(mt_var_filt_test$TUMOR_FREQ_VarScan2)==FALSE)
mt_var_filt_test<-filter(mt_var_filt_test,is.na(mt_var_filt_test$NORMAL_AF_Mutect2)==TRUE | is.na(mt_var_filt_test$TUMOR_AF_Mutect2)==TRUE | is.na(mt_var_filt_test$NORMAL_FREQ_VarScan2)==TRUE | is.na(mt_var_filt_test$TUMOR_FREQ_VarScan2)==TRUE)
mt_var_both<-mt_var_both %>% filter(mt_var_both$TUMOR_FREQ_VarScan2 %/% mt_var_both$NORMAL_FREQ_VarScan2 > 2)

mt_var_filt_test<-mt_var_filt_test %>% filter(mt_var_filt_test$TUMOR_FREQ_VarScan2 %/% mt_var_filt_test$NORMAL_FREQ_VarScan2 > 2 | is.na(mt_var_filt_test$TUMOR_FREQ_VarScan2)==TRUE)
mt_var_filt_test<-mt_var_filt_test %>% filter(mt_var_filt_test$TUMOR_AF_Mutect2 %/% mt_var_filt_test$NORMAL_AF_Mutect2 >3 | is.na(mt_var_filt_test$TUMOR_AF_Mutect2)==TRUE)

mt_var_filt_test<-rbind(mt_var_filt_test,mt_var_both)

# filtering 후 sample 별 variants 개수 정리  
sample_count<-data.frame(table(mt_var_filt_test$sample))
names(sample_count)[1] <- c("sample")
names(sample_count)[2] <- c("count")

#information data + count 
test2<-merge(sample_count,somatic_mut,by.x = "sample",by.y = "ID",all.x = TRUE, all.y = FALSE)

#비교 $ t-test
ggplot(data=test2) + geom_point(aes(x=count, y=revised_SNV + revised_indel), color="red", alpha=3, size=5) +
theme(axis.title.x=element_text(size=15)) + geom_smooth(method='lm', aes(x=count, y=revised_SNV+revised_indel))

ggplot(data=test2[!(is.na(test2$Smoking) | test2$Smoking=="NA"),]) + geom_bar(aes(x=1, y=count, fill=Smoking), position="dodge", stat="summary", fun.y="mean") 
+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())

ggplot(data=test2) + geom_point(aes(x=count, y=Age)) + geom_smooth(method="lm", aes(x=count, y=Age)) 

ggplot(data=test2[!(is.na(test2$Smoking) | test2$Smoking=="NA"),]) + geom_bar(aes(x=1, y=revised_SNV + revised_indel, fill=Smoking), position="dodge", stat="summary", fun.y="mean") 
+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())

ggplot(data=test2[!(is.na(test2$Smoking) | test2$Smoking=="NA"),]) + geom_bar(aes(x=1, y=ctSig4, fill=Smoking), position="dodge", stat="summary", fun.y="mean") 
+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())

t.test(test[!(is.na(test2$Smoking) | test2$Smoking=="NA"),]$count~test2[!(is.na(test2$Smoking) | test2$Smoking=="NA"),]$Smoking)

