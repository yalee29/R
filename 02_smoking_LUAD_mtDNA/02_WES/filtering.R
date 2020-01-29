#filtering
wes_filtering<-wes_sorted
wes_filtering<-wes_filtering%>% filter(!(POS %in% c(302:315, 513:525, 3105:3109)))
wes_filtering<-wes_filtering%>% filter(wes_filtering$SS_VarScan2==2) #varscan에서 somatic mutation만 

#MUtect, Varscan 모두 vaf 1%이상 만 가져오기. (2707 -> 164)
wes_filtering<-wes_filtering[wes_filtering$TUMOR_AF_Mutect2>=0.01 | is.na(wes_filtering$TUMOR_AF_Mutect2)==TRUE,]
wes_filtering<-wes_filtering[wes_filtering$TUMOR_AF_VarScan2>=0.01| is.na(wes_filtering$TUMOR_AF_VarScan2)==TRUE,]

#filtering 1 -> candidate1 ; (both calling은 tumor depth 높은것 기준)  depth 100미만은 vaf 10%이상, depth 100이상은 vaf 4%이상 / (164 -> 124)
wes_filtering$Sample<-as.character(wes_filtering$Sample) #Sample class가 factor -> wes_filtering에서 하나씩 가져오면 sample 이름이 바뀜. 
class(wes_filtering$Sample)

colname_wes <-colnames(wes_filtering)
candidate1<-data.frame()
candidate1<- setNames(data.frame(matrix(ncol = 22, nrow = 0)), colname_wes)

for (i in c(1:nrow(wes_filtering))){
  if (sum(is.na(wes_filtering[i,c(6,8)]))==0 & wes_filtering[i,10] > wes_filtering[i,12]){
    if (wes_filtering[i,10] <100){
      if (wes_filtering[i,6]> 0.1){
        candidate1<-rbind(candidate1,wes_filtering[i,])
      }
    }
    else if (wes_filtering[i,10] >100){
      if (wes_filtering[i,6]>0.04){
        candidate1<-rbind(candidate1,wes_filtering[i,])
      }
    }
  }
  else if (sum(is.na(wes_filtering[i,c(6,8)]))==0 & wes_filtering[i,10] < wes_filtering[i,12]){
    if (wes_filtering[i,12] <100){
      if (wes_filtering[i,8]> 0.1){
        candidate1<-rbind(candidate1,wes_filtering[i,])
      }
    }
    else if (wes_filtering[i,12] >100){
      if (wes_filtering[i,8]>0.04){
        candidate1<-rbind(candidate1,wes_filtering[i,])
      }
    }
  }
  else if (sum(is.na(wes_filtering[i,c(6,8)]))==1) {
    if (wes_filtering[i,12] <100){
      if (wes_filtering[i,8]> 0.1){
        candidate1<-rbind(candidate1,wes_filtering[i,])
      }
    }
    else if (wes_filtering[i,12] >100){
      if (wes_filtering[i,8]>0.04){
        candidate1<-rbind(candidate1,wes_filtering[i,])
      }
    }
  }  
}

#filtering2 ; Panel of normal(2707) ; variant 중 pon에 2개 이상 있으면 제외/(124 -> 106)
pon<-wes_sorted%>%select(Sample, POS, REF, ALT, NORMAL_AF_Mutect2, NORMAL_AF_VarScan2, NORMAL_DP_Mutect2, NORMAL_DP_VarScan2)
candidate2<-data.frame()
candidate2<-setNames(data.frame(matrix(ncol = 22, nrow = 0)), colname_wes)

for (i in c(1:nrow(candidate1))){
  if (sum(pon$POS %in% candidate1[i,2])<=2){
    candidate2<-rbind(candidate2, candidate1[i,])
  }
}

#vaf 값 지정 & vaf/purity 값 지정
candidate2[which(is.na(candidate2$TUMOR_AF_Mutect2)==TRUE),c(5,6)]<-0
candidate2$vaf<-apply(candidate2[,c("TUMOR_AF_Mutect2","TUMOR_AF_VarScan2")],1, max)
candidate2$vaf<-candidate2$vaf%>%as.numeric()
candidate2 <- candidate2 %>% mutate(v_d_p = vaf/purity)

#variant calling 되지 않은 sample 추가 (총 19개 sample 추가; nrow 106 ->125)
candidate2_test<-merge(candidate2, LungWES_info, all.y = TRUE)
candidate2_test<-candidate2_test%>%select(colnames(candidate2))
candidate2<-candidate2_test
rm(candidate2_test)
candidate2[which(is.na(candidate2$vaf)==TRUE),c(23,24)]<-0

#sample별로 가장 높은 vaf만 추출-> can2_top_vaf로 저장
can2_top <- do.call(rbind, by(candidate2, candidate2$Sample, function(x) {
  return(x[which.max(x$v_d_p),])
}))


#Smoker <-> Nonsmoker group 지정
#Smoking=NA인 sample(RS24)은 나중에 sig4보고 추가해줄 것. 
smoker<-can2_top%>%filter(can2_top$`Smoking Status` =="Yes") #29
nonsmoker<-can2_top%>%filter(can2_top$`Smoking Status` =="No") #39

#비교 분석.
library(ggplot2)
library(ggpubr)
smoker$ID <- "smoker"
nonsmoker$ID <- "nonsmoker"
test<-rbind(smoker, nonsmoker)
ggboxplot(test, x="ID", y="v_d_p", color = "ID", palette = "jco", add = "jitter", width = 0.5) + stat_compare_means() + ggtitle("smoker, nonsmoker")


#top_vaf 말고 calling mutation 개수로 비교. 



#sig4로 비교

smoker_sig4<-smoker%>%filter(ctSig4>=15000)
nonsmoker_sig4<-nonsmoker%>%filter(ctSig4<200)










