sam<-list.files(pattern ="*voc.sam" )
all<-read.csv(sam, sep = "\t", header = FALSE, skip = 7)
if(length(which(all$V3==""))>0) all<-all[-which(all$V3==""),]
results<-as.data.frame(table(all$V3))
write.csv(results, gsub("voc.sam","voc_depth.csv",sam),row.names = FALSE  )
