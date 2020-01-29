df_result <- tibble(break_point = character())


for(f in list.files("/home/users/yunah1029/notebook/output", patter="*.tsv")) {
  df_temp <- read_tsv(paste0("/home/users/yunah1029/notebook/output/", f))
  colnames(df_temp)[2] <- f
  
  df_result %<>% full_join(df_temp)
}

write_tsv(df_result, "/home/users/yunah1029/notebook/output/test.tsv")

result_filtering <-df_result %>% filter(break_point %in% c("120_16241","16118_174","16128_165","16134_155","16189_65","16192_66","16224_67","16332_151","16365_42","16449_104","16483_171","45_16243","59_16243","54_16243","56_16243","64_16190","65_16193")) %>%
  select(1, which(colSums(is.na(.[2:ncol(.)]))!=nrow(.)) + 1)
  
write_tsv(result_filtering, "/home/users/yunah1029/notebook/output/result_filtering.tsv")

re_filtering_t <- result_filtering %>% t() %>% as.data.frame()
class(re_filtering_t)
head(re_filtering_t)
re_filtering_t<-re_filtering_t[,which(colSums(is.na(re_filtering_t))<100)]
head(re_filtering_t)
re_filtering_t <- re_filtering_t %>% rownames_to_column()
write_tsv(re_filtering_t, "/home/users/yunah1029/projects/03_long_read_mt/re_filtering_t.tsv")

re_filtering_t <- read_tsv("./projects/03_long_read_mt/re_filtering_t.tsv")
re_filtering_t[is.na(re_filtering_t)] <-0

ggplot(re_filtering_t)+
  geom_histogram(aes(break_point), binwidth=0.1, stat = "count")





