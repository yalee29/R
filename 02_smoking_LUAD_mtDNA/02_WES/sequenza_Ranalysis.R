library(sequenza)

for(f in list.files("/home/users/yunah1029/R/x86_64-redhat-linux-gnu-library/3.5/sequenza/smoking_wes", pattern="*.small.seqz.gz$")) {
  temp <- sequenza.extract(paste0("/home/users/yunah1029/R/x86_64-redhat-linux-gnu-library/3.5/sequenza/smoking_wes/",f), verbose = TRUE)
  cp_temp <- sequenza.fit(temp)
  ID=strsplit(f, split="[.]")[[1]][1]
  sequenza.results(sequenza.extract = temp, cp.table = CP_temp, sample.id = ID, out.dir = "/home/users/yunah1029/projects/smoking_LUAD_mtDNA/03_WES/02_sequenza/output")
}

rm(temp)
rm(cp_temp)
rm(ID)



