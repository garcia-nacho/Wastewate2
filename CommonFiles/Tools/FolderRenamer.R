library(readxl)


files<-list.files(pattern = "xlsx")
if(length(files)==1){
  df<-read_xlsx(files)
  if(length(which(df$Status=="FAIL"))){
    df<-df[-which(df$Status=="FAIL"),]
  }
  all.files<-list.files(recursive = TRUE)
  
  for (i in 1:nrow(df)) {
    if(length(grep(df$Barcode[i],all.files))>0){
      dir.create(paste(df$`NGS Oppsett`[i],".", df$`Geo-Navn`[i],sep = ""))
      file.rename(all.files[grep(df$Barcode[i],all.files)], gsub(df$Barcode[i],paste(df$`NGS Oppsett`[i],".", df$`Geo-Navn`[i],sep = ""),all.files[grep(df$Barcode[i],all.files)]))
      file.remove(df$Barcode[i])
    }
  }
  
}



