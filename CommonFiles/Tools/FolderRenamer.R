library(readxl)


files<-list.files(pattern = "xlsx")
if(length(files)==1){
  df<-read_xlsx(files)
  colnames(df)<-df[1,]
  df<-df[-1,] 
  if(length(which(is.na(df$Barcode)))>0) df<-df[-which(is.na(df$Barcode)),]
  all.files<-list.files(recursive = TRUE)
  for (i in 1:nrow(df)) {
    if(length(grep(df$Barcode[i],all.files))>0){
      dir.create(df$SequenceID[i])
      file.rename(all.files[grep(df$Barcode[i],all.files)], gsub(df$Barcode[i],df$SequenceID[i],all.files[grep(df$Barcode[i],all.files)]))
      file.remove(df$Barcode[i])
    }
  }
}