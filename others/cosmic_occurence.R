#cosmic data 이용 -> filtering 전 sample 별 driver mutation 확인

library(readr)
#cosmic_occurence function 
cosmic_occurence <- function(x){
  if (x == '.'){
    return(0)
  } else{
    y <- x %>% strsplit(., 'OCCURENCE=') %>% .[[1]] %>% .[2] %>%
      gsub('_','',.) %>%
      gsub('\\([A-z0-9]*\\)', '', .) %>%
      gsub('\\([A-z0-9]*\\)', '', .) %>%
      strsplit(',') %>% .[[1]] %>%
      as.numeric() %>%
      sum()
    return(y)
  }
}

#sample마다 cosmic_occurence 구한 후 내림차순 정렬 ->바로 저장. 
for(f in list.files("/home/users/team_projects/LungCancer_IO_WGS/01_WGS/07_pointMT/02_indel",pattern = "*.anv$")) {
  temp <- read_tsv(paste0("/home/users/team_projects/LungCancer_IO_WGS/01_WGS/07_pointMT/02_indel/",f))
  for (i in 1:nrow(temp)){
    temp$cosmic_coding_occ[i] <- cosmic_occurence(temp$cosmic86_coding[i])
  }
  temp <- temp[order(temp$cosmic_coding_occ, decreasing = TRUE),]
  temp %>% write_tsv(paste0("/home/users/team_projects/LungCancer_IO_WGS/01_WGS/07_pointMT/02_indel/",f,".occ"))
}

##test
##t1_2 <- setdiff(t2_3, t1_1)
#for (i in 1:nrow(test)){
#  test$cosmic_coding_occ[i] <- cosmic_occurence(test$cosmic86_coding[i])
#}
#test <- test[order(test$cosmic_coding_occ, decreasing = TRUE),]