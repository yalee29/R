library(sequenza)

for(f in list.files("/home/users/yunah1029/R/x86_64-redhat-linux-gnu-library/3.6/sequenza/extdata/lungCa_IO_SNUH/",pattern = "*.small.seqz.gz$")) {
  temp <- sequenza.extract(paste0("/home/users/yunah1029/R/x86_64-redhat-linux-gnu-library/3.6/sequenza/extdata/lungCa_IO_SNUH/",f), verbose = TRUE)
  cp_temp <- sequenza.fit(temp)
  ID=strsplit(f,split = "[.]")[[1]][1]
  sequenza.results(sequenza.extract = temp, cp.table = CP_temp, sample.id = ID, out.dir = "/home/users/team_projects/LungCancer_IO_WGS/01_WGS/08_sequenza/result")
}

rm(temp)
rm(cp_temp)
rm(ID)