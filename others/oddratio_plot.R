library(cowplot)
library(ggplot2)
library(dplyr)

dir <- ("/home/users/yunah1029/projects/04_lungCa_IO_SNUH/07_pointMT/01_snv/")
file_list<-list.files("/home/users/yunah1029/projects/04_lungCa_IO_SNUH/07_pointMT/01_snv",recursive=TRUE, pattern="*.oddp$")

#pdf("/home/users/yunah1029/projects/LungCancer_IO_or_3.pdf",8,5)
p2 <- list()
for (n in setdiff(1:length(file_list),11)) {
  print(file_list[n])
  temp <- read.csv(paste0(dir, file_list[n]), sep = '\t', stringsAsFactors = FALSE, na.strings = "", header = TRUE)
  temp$mutation <- paste0(temp$REF,">",temp$ALT)
  #for (f in 1:nrow(temp)){
  #  temp$mutation[f] <- paste0(strsplit(temp$REF, "")[[f]][1],">",strsplit(temp$ALT, "")[[f]][1])
  #}
  temp$mutation <- gsub('G>T','C>A', temp$mutation)
  temp$mutation <- gsub('G>C','C>G', temp$mutation)
  temp$mutation <- gsub('G>A','C>T', temp$mutation)
  temp$mutation <- gsub('A>T','T>A', temp$mutation)
  temp$mutation <- gsub('A>G','T>C', temp$mutation)
  temp$mutation <- gsub('A>C','T>G', temp$mutation)
  temp<- temp[is.na(temp$log2OR)==0,]
  sample_name <- strsplit(file_list[n],"[.]")[[1]][1]

  p2[[n]] <-ggplot(sample_n(temp,2000), aes(var_readN, log2OR)) + #color = mutation
      ggtitle(sample_name)+
      theme_bw()+
      coord_cartesian(xlim=c(0,20), y=c(-12,0))+
      #scale_color_manual(values=c('C>A'='#3bb2ed','C>G'='#000000', 'C>T'='#ff0000','T>A'="#ada6a6",'T>C'="#7dcf19",'T>G'="#f29bcb")) +
      scale_y_continuous(breaks=seq(-12,0,1))+
      scale_x_continuous(breaks=seq(0,20,2))+
      geom_jitter(height=0, width=0.5, size=3, alpha=0.3)
  
}

plot_grid(plotlist=p2, ncol=6)

#dev.off()