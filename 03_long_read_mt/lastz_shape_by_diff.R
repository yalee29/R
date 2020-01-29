#CHM13 16k에서 v/^ form 찾기(미분사용)
dir = "/home/users/mspark/long_read/data/CHM13/test/16k/rdp/"
pdf("/home/users/yunah1029/projects/images/long_read_mt/CHM13_16k_Vform.pdf")

for (f in list.files(dir, pattern = ".rdp")) {
  temp_rdp <- read.table(paste0(dir, f), header = TRUE)
  
  #dy/dx
  dy = diff(temp_rdp[,2])
  dx = diff(temp_rdp[,1])
  
  dydx = dy/dx
  dydx <- dydx[!is.na(dydx)]
  
  if(length(unique(dydx))>=2) {
    print(sprintf("%s: %d", f, length(unique(dydx))))
    plot(temp_rdp,type = "l", main = "CHM13_10to15k", asp = 1, xaxs = 'i', yaxs = 'i')
  }
}

dev.off()