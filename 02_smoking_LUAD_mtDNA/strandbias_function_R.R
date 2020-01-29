# Strand bias

StrandBias <- function(df_varmat) {
  
  df <- subset(df_varmat, nchar(REF) == nchar(ALT))
  
  AG = 0; AC = 0; AT = 0; GA = 0; GC = 0; GT = 0; CA = 0; CG = 0; CT = 0; TA = 0; TG = 0; TC = 0
  
  for(i in 1:nrow(df)) {
    
    ref = strsplit(df$REF[i], "")[[1]]
    alt = strsplit(df$ALT[i], "")[[1]]
    
    for(j in 1:length(ref)) {
      
      switch(ref[j],
             "A"={switch(alt[j],
                         "G"={ AG = AG + 1 },
                         "C"={ AC = AC + 1 },
                         "T"={ AT = AT + 1 })},
             "G"={switch(alt[j],
                         "A"={ GA = GA + 1 },
                         "C"={ GC = GC + 1 },
                         "T"={ GT = GT + 1 })},
             "C"={switch(alt[j],
                         "A"={ CA = CA + 1 },
                         "G"={ CG = CG + 1 },
                         "T"={ CT = CT + 1 })},
             "T"={switch(alt[j],
                         "A"={ TA = TA + 1 },
                         "G"={ TG = TG + 1 },
                         "C"={ TC = TC + 1 })}
      )
    }
  }
  
  df <- data.frame(context=c("AG", "AC", "AT", "GA", "GC", "GT", "CA", "CG", "CT", "TA", "TC", "TG"), count=NA, strand=NA)
  
  df$count[df$context=="AG"] <- AG
  df$count[df$context=="AC"] <- AC
  df$count[df$context=="AT"] <- AT
  df$count[df$context=="GA"] <- GA
  df$count[df$context=="GC"] <- GC
  df$count[df$context=="GT"] <- GT
  df$count[df$context=="CA"] <- CA
  df$count[df$context=="CG"] <- CG
  df$count[df$context=="CT"] <- CT
  df$count[df$context=="TA"] <- TA
  df$count[df$context=="TC"] <- TC
  df$count[df$context=="TG"] <- TG
  
  df$strand <- replace(df$strand, df$context %in% c("CA", "CG", "CT", "TA", "TC", "TG"), "L") # L-strand가 reference
  df$strand <- replace(df$strand, df$context %in% c("GT", "GC", "GA", "AT", "AG", "AC"), "H")   
  
  df$context <- replace(df$context, df$context == "GT", "CA")
  df$context <- replace(df$context, df$context == "GC", "CG")
  df$context <- replace(df$context, df$context == "GA", "CT")
  df$context <- replace(df$context, df$context == "AT", "TA")
  df$context <- replace(df$context, df$context == "AG", "TC")
  df$context <- replace(df$context, df$context == "AC", "TG")
  
  df$context <- as.character(df$context)
  
  l <- list()
  for(context in unique(df$context)) {
    l[[length(l) + 1]] <- list(context = context, count = sum(df$count[(df$context == context) & (df$strand=="H")]), strand="H")
    l[[length(l) + 1]] <- list(context = context, count = sum(df$count[(df$context == context) & (df$strand=="L")]), strand="L")
  }
  
  df <- ldply(l, data.frame)
  
  g <- ggplot(data=df, aes(x=context, y=count, fill=strand)) + geom_bar(position=position_dodge(), stat="identity", width=0.7) +
    ylab("Counts") +
    scale_x_discrete(labels=c("CA"="C>A", "CG"="C>G", "CT"="C>T", "TA"="T>A", "TC"="T>C", "TG"="T>G")) +
    theme_minimal() + 
    theme(
      axis.text.x  = element_text(size=15),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=15),
      axis.text.y  = element_text(size=15)
    )
  
  print(g)
  
}
