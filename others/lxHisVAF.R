
lxHistVAF <- function(POS, ALT, dir="/home/users/yunah1029/projects/smoking_LUAD_mtDNA/02_WES_MT/05_pileup/", MT="MT", show.x.label=FALSE, label.size=10) {
  
  out <- system(sprintf("grep %s[[:blank:]]%s[[:blank:]] %s/*.mp", MT, POS, dir), intern = TRUE)
  out <- strsplit(out, '\t')
  
  REF <- out[[1]][3]
  
  if(startsWith(ALT, "-") || startsWith(ALT, "+")) {
    
    VAFs <- sapply(out, function(x) {
      
      s <- gsub("[\\{\\}]", "", x[10])
      s <- gsub(":", "=", s)
      s <- eval(parse(text=paste("list(", s, ")")))
      
      upper <- if(!is.null(s[[toupper(ALT)]])) s[[toupper(ALT)]] else 0
      lower <- if(!is.null(s[[tolower(ALT)]])) s[[tolower(ALT)]] else 0
      
      return((upper + lower)/as.numeric(x[4]) * 100)
    })
    
    if(startsWith(ALT, "-")) Variant <- paste(POS,
                                              paste0(paste0(REF, substr(ALT, 2, nchar(ALT))),
                                                     ">",
                                                     REF))
    if(startsWith(ALT, "+")) Variant <- paste(POS,
                                              paste0(REF,
                                                     ">",
                                                     paste0(REF, substr(ALT, 2, nchar(ALT)))))
  

  } else {
    
    switch(ALT,
           "A"={ column = 5 },
           "G"={ column = 6 },
           "C"={ column = 7 },
           "T"={ column = 8 })
    
    VAFs    <- sapply(out, function(x) { as.numeric(x[column])/as.numeric(x[4]) * 100 })
    
    Variant <- paste(POS, paste0(REF,">",ALT))
  }
  
  df <- data.frame(sample=sapply(out, function(x) { gsub(paste0(".mp:", MT), "", basename(x[1])) }), VAF=VAFs)
  df$sample <- factor(df$sample, levels=df$sample[order(df$VAF, decreasing = TRUE)])
  
  print(df)
  
  g <- ggplot(data=df, aes(x=sample, y=VAF)) + geom_bar(stat="identity") + theme_minimal() +
    theme(axis.text.x=if(show.x.label) element_text(size=label.size, angle=45, hjust=1.0, vjust=1.0) else element_blank(),
          legend.title=element_blank(),
          plot.margin = margin(0,0,0,1.5,"cm")
    )
  
  print(g)
  
}
