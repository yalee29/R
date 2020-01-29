library(tidyverse)
library(ggplot2)

#Data(+information) load & Edit data (sample name에서 .merged 없애기)
re_filtering_total<-read.csv("projects/smoking_LUAD_mtDNA/01_WGS/re_filtering.csv", stringsAsFactors = FALSE, na.strings = "")
re_filtering_info<-read.csv("projects/smoking_LUAD_mtDNA/01_WGS/LungWGSinfo.csv", stringsAsFactors = FALSE, na.strings = "")

#information data와 merge하기 전 sample name, ID 일치시켜줘야함
for(i in 1:nrow(re_filtering_info)){
  re_filtering_info[i,1]<-strsplit(re_filtering_info[i,1],"[.]")[[1]][1]
}

re_filtering_merged<-merge(re_filtering_total, re_filtering_info, by.x = "sample", by.y = "ID", all.x = TRUE, all.y = FALSE)
re_filtering_merged<-re_filtering_merged%>%select("sample","POS","REF","ALT","NORMAL_AF_Mutect2","TUMOR_AF_Mutect2","NORMAL_FREQ_VarScan2","TUMOR_FREQ_VarScan2","Func_ensGene","Gene_ensGene","ExonicFunc_ensGene","Gender.x","Age.x","Smoking.x","ctSig4.x","revised_purity")
re_filtering<-re_filtering_merged%>%rename(c("Age.x"="Age","Gender.x"="Gender","ctSig4.x"="ctSig4", "revised_purity"="purity","Smoking.x"="Smoking"))

#vaf 값 지정 & vaf/purity 값 지정
re_filtering$vaf<-apply(re_filtering[,c("TUMOR_AF_Mutect2","TUMOR_FREQ_VarScan2")],1, max)
for (f in which(is.na(re_filtering$vaf)==TRUE)){
  if (is.na(re_filtering[f,6])==TRUE){
    re_filtering[f,17]<-re_filtering[f,8]
  } else {re_filtering[f,17] <- re_filtering[f,6]}
} 

re_filtering$purity <- re_filtering$purity*100  
re_filtering <- re_filtering %>% mutate(v_d_p = vaf/purity)

#sample별로 가장 높은 vaf만 추출
df_test <- do.call(rbind, by(re_filtering, re_filtering$sample, function(x) {
  return(x[which.max(x$v_d_p),])
}))
re_filtering_top<-df_test

#Smoker <-> Nonsmoker group 지정  
smoker<-re_filtering_top%>%filter(Smoking =="Yes")
nonsmoker<-re_filtering_top%>%filter(Smoking =="No")

#비교 분석.
library(ggplot2)
library(ggpubr)

smoker_sig4<-smoker%>%filter(ctSig4>=15000)
nonsmoker_sig4<-nonsmoker%>%filter(ctSig4<200)
smoker_sig4$ID <- "smoker"
nonsmoker_sig4$ID <- "nonsmoker"
test<-rbind(smoker_sig4, nonsmoker_sig4)
ggboxplot(test, x="ID", y="v_d_p", color = "ID", palette = "jco", add = "jitter", width = 0.5) + stat_compare_means(method = t.test) + ggtitle("smoker, nonsmoker")


#
geom_jitter() + ggtitle("smoker 30000 , nonsmoker 50")
ggplot(test, aes(ID, Age)) + geom_boxplot() + ggtitle("smoker, nonsmoker")
ggplot(test, aes(ID, v_d_p)) + geom_violin() + ggtitle("smoker 10000 , nonsmoker 1000")

#비교분석. age


ggplot(re_filtering_top, aes(Smoking, Age)) + geom_violin() + ggtitle("Age")




t.test(test[test$ID==smoker,18]~ test[test$ID==nonsmoker,18])
t.test(test$ID ~ test$v_d_p)

t.test(test[!(is.na(test2$Smoking) | test2$Smoking=="NA"),]$count~test2[!(is.na(test2$Smoking) | test2$Smoking=="NA"),]$Smoking)




