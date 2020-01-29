#_anv.tsv data 모두 load, merge & 필요한 column만 
library(dplyr)
dir<-("./projects/smoking_LUAD_mtDNA/04_Merged/01_Annotated/")
file_list<-list.files(dir)
annotation_rbind<-data.frame()

for (file in file_list){
  temp<-read.csv(paste0(dir,file),sep='\t', stringsAsFactors = FALSE, na.strings = "")
  temp$sample<-strsplit(file,"[_]")[[1]][1]
  annotation_rbind<-rbind(annotation_rbind,temp)
}

annotation_sort<-annotation_rbind%>%select(sample,CHROM, POS, REF, ALT, NORMAL_AF_Mutect2, TUMOR_AF_Mutect2, NORMAL_FREQ_VarScan2, TUMOR_FREQ_VarScan2, SS_VarScan2, Func_ensGene, Gene_ensGene, GeneDetail_ensGene, ExonicFunc_ensGene, AAChange_ensGene)

#annotation_sort filtering -> annotation_filt_test 
annotation_filt<- annotation_sort 
annotation_filt$NORMAL_FREQ_VarScan2 <- annotation_filt$NORMAL_FREQ_VarScan2 %>% gsub('%', '', .) %>% as.numeric()
annotation_filt$TUMOR_FREQ_VarScan2 <- annotation_filt$TUMOR_FREQ_VarScan2 %>% gsub('%', '', .) %>% as.numeric()

annotation_filt_test<-annotation_filt%>% filter(!(POS %in% c(302:315, 513:525, 3105:3109))) #특정 위치 제거
annotation_filt_test<-annotation_filt_test[!grepl(",",annotation_filt_test$ALT),]
annotation_filt_test$TUMOR_AF_Mutect2<-annotation_filt_test$TUMOR_AF_Mutect2%>%as.numeric()
annotation_filt_test$NORMAL_AF_Mutect2<-annotation_filt_test$NORMAL_AF_Mutect2%>%as.numeric()
annotation_filt_test$NORMAL_AF_Mutect2<-annotation_filt_test$NORMAL_AF_Mutect2*100
annotation_filt_test$TUMOR_AF_Mutect2<-annotation_filt_test$TUMOR_AF_Mutect2*100

annotation_filt_test<-annotation_filt_test%>%filter(annotation_filt_test$SS_VarScan2==2 | is.na(annotation_filt_test$TUMOR_FREQ_VarScan2)==TRUE) #varscan에서 somatic mutation만 

#MUtect, Varscan 모두 vaf 1%이상 만 가져오기. 
annotation_filt_test<-annotation_filt_test[annotation_filt_test$TUMOR_AF_Mutect2>=1 | is.na(annotation_filt_test$TUMOR_AF_Mutect2)==TRUE,]
annotation_filt_test<-annotation_filt_test[annotation_filt_test$TUMOR_FREQ_VarScan2>=1| is.na(annotation_filt_test$TUMOR_FREQ_VarScan2)==TRUE,]

#varscan, mutect이 모두 detect한 것들은 따로 뽑아서 정리.
#둘 다 찾은 것들은 varscan기준으로 tumor가 normal vaf보다 최소 2배 이상 큰 것들만. 
annotation_both<-filter(annotation_filt_test,is.na(annotation_filt_test$NORMAL_AF_Mutect2)==FALSE,is.na(annotation_filt_test$TUMOR_AF_Mutect2)==FALSE, is.na(annotation_filt_test$NORMAL_FREQ_VarScan2)==FALSE, is.na(annotation_filt_test$TUMOR_FREQ_VarScan2)==FALSE)
annotation_filt_test<-filter(annotation_filt_test,is.na(annotation_filt_test$NORMAL_AF_Mutect2)==TRUE | is.na(annotation_filt_test$TUMOR_AF_Mutect2)==TRUE | is.na(annotation_filt_test$NORMAL_FREQ_VarScan2)==TRUE | is.na(annotation_filt_test$TUMOR_FREQ_VarScan2)==TRUE)
annotation_both<-annotation_both %>% filter(annotation_both$TUMOR_FREQ_VarScan2 %/% annotation_both$NORMAL_FREQ_VarScan2 > 2)

#mutect은 tumor vaf/ normal vaf >2, varscan은 tumor vaf/ normal vaf >3 인 것들만. 
#annotation_both + annotation_filt_test
annotation_filt_test<-annotation_filt_test %>% filter(annotation_filt_test$TUMOR_FREQ_VarScan2 %/% annotation_filt_test$NORMAL_FREQ_VarScan2 > 2 | is.na(annotation_filt_test$TUMOR_FREQ_VarScan2)==TRUE)
annotation_filt_test<-annotation_filt_test %>% filter(annotation_filt_test$TUMOR_AF_Mutect2 %/% annotation_filt_test$NORMAL_AF_Mutect2 >3 | is.na(annotation_filt_test$TUMOR_AF_Mutect2)==TRUE)

annotation_filt_test<-rbind(annotation_filt_test,annotation_both)

#annotation info + smoking info(=somatic_mut)
sample_annotation<-annotation_filt_test%>%select(sample, Func_ensGene, ExonicFunc_ensGene, Gene_ensGene)
test3<-merge(sample_annotation,somatic_mut,by.x = "sample",by.y = "ID",all.x = TRUE, all.y = FALSE)
snv_smoking<-test3
rm(test3)

# filtering 후 sample 별 variants 개수 정리  
sample_count<-data.frame(table(annotation_filt_test$sample))
names(sample_count)[1] <- c("sample")
names(sample_count)[2] <- c("count")

#information data + count 
test2<-merge(sample_count,somatic_mut,by.x = "sample",by.y = "ID",all.x = TRUE, all.y = FALSE)
sample_count_info<-test2
rm(test2)

#re_filtering data 불러오기
re_filtering<-read.csv("./projects/smoking_LUAD_mtDNA/re_feiltering.csv")

sample_count2<-data.frame(table(re_filtering$sample)) #refilt_count data; 각 sample 별 count 
names(sample_count2)[1] <-c("sample")
names(sample_count2)[2] <- c("count_refilt")
cf_count<-merge(sample_count, sample_count2, by.x = "sample", by.y = "sample") #refeiltering 이후 count 수 변화. 

#ggplot input으로 넣어줄 data (refilt_count + info)
test3<-merge(sample_count2,somatic_mut,by.x = "sample",by.y = "ID",all.x = TRUE, all.y = FALSE)

#1) Smoking==NA 바꿔주기 방법 1
for(k in c(which(test3$Smoking=="NA"))){
  if (test3$ctSig4[k] >4000){
    test3$Smoking[k]<-gsub("NA","Yes",test3$Smoking[k])
    }
  else {
    test3$Smoking[k]<-gsub("NA","No",test3$Smoking[k])
    }
} 

#1) Smoking==NA 바꿔주기 방법 2
test3$Smoking[test3$Smoking=="NA"] <- ifelse(test3$ctSig4[test3$Smoking=="NA"]>4000,"Yes","No")

sample_count2_info<-test3
rm(test3)

#ggplot
library(ggplot2)

#sig4 <-> count
ggplot(data=sample_count2_info) + geom_point(aes(x=count_refilt, y=ctSig4), color="red") +
  theme(axis.title.x=element_text(size=15)) + geom_smooth(method='lm', aes(x=count_refilt, y=ctSig4))
t.test(sample_count2_info$count_refilt~sample_count2_info$ctSig4)

#Smoking <-> count
ggplot(data=sample_count2_info) + geom_bar(aes(x=1, y=count_refilt, fill=Smoking), position="dodge", stat="summary", fun.y="mean") + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())
t.test(count_refilt ~ Smoking, data=sample_count2_info, var.equal = TRUE, conf.level = 0.95)

#StrandBias
#StrandBias script 함수등록
re_filtering$REF<-as.character(re_filtering$REF)
re_filtering$ALT<-as.character(re_filtering$ALT)
StrandBias(re_filtering)

nonsmoker <-re_filtering%>%filter(Smoking=="No")
smoker <- re_filtering%>%filter(Smoking=="Yes")
StrandBias(smoker)
StrandBias(nonsmoker)

test4<-merge(annotation_filt_test, somatic_mut, by.x = "sample", by.y = "ID", all.x = TRUE, all.y = FALSE)
smoker<-test4%>%filter(Smoking=="Yes")
nonsmoker <-test4%>%filter(Smoking=="No")
StrandBias(smoker)
StrandBias(nonsmoker)





#
ggplot(sample_count2_info, aes(x=Func_ensGene,1, fill=Smoking))+geom_bar(stat='identity')
t.test(sample_count2_info$Smoking~re_filtering$Func_ensGene=="exonic")

ggplot(test3[!(is.na(test3$Smoking)|test3$Smoking=="NA"),], aes(Gene_ensGene, 1, fill=Smoking))+geom_bar(stat='identity')+theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggplot(test3[!(is.na(test3$Smoking)|test3$Smoking=="NA"),], aes(Func_ensGene,n(Func_ensGene), fill=Smoking)) + 
  geom_bar(aes(Func_ensGene, fill = as.factor(Smoking)), 
           position = "dodge", stat = "summary", fun.y = "mean")

ggplot(test3[!(is.na(test3$Smoking)|test3$Smoking=="NA"),], aes(Func_ensGene, ?mean, fill=Smoking))+geom_bar(stat='identity', position = 'dodge')+stat_summary(fun.y = 'mean')
ggplot(test3[!(is.na(test3$Smoking)|test3$Smoking=="NA"),], aes(Gene_ensGene, 1, fill=Smoking))+geom_bar(stat='identity', position = 'dodge')+theme(axis.text.x = element_text(angle = 60, hjust = 1))+stat_summary(fun.y = mean)

ggplot(test3[!(is.na(test3$Smoking)|test3$Smoking=="NA"),], aes(x=Func_ensGene, fill=Smoking))+geom_bar(aes(x=Func_ensGene, y=, fill=Smoking), position="dodge", stat="summary", fun.y="mean") 

library(wesanderson)
test_rk$Func_ensGene[grepl("upstream", test_rk$Func_ensGene)] <- "Intergenic"
ggplot(data=test_rk[!(is.na(test_rk$Smoking) | test_rk$Smoking=="NA"),], aes(x=Smoking)) + geom_bar(position="fill", aes(y=(..count..)/sum(..count..), fill=Func_ensGene), stat="count") +scale_fill_manual(values=wes_palette(n=3, name="Moonrise1"))  


#StrandBias - Smoker/Nonsmoker
test4<-merge(annotation_filt_test, somatic_mut, by.x = "sample", by.y = "ID", all.x = TRUE, all.y = FALSE)
smoker<-test4%>%filter(Smoking=="Yes")
nonsmoker <-test4%>%filter(Smoking=="No")
StrandBias(smoker)
StrandBias(nonsmoker)
