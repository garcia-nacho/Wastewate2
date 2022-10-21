library(data.table)
library(ggplot2)
library("doParallel")
library("parallel")
library("foreach")
library(doSNOW)
library(progress)
library(digest)
library(writexl)
library(seqinr)
path=commandArgs(TRUE)
co.n<-0.15
pathf<-path[1]
co.n<-as.numeric(path[2])
n.start<-as.numeric(path[3])
n.end<- as.numeric(path[4])

noise.path<-list.files(pathf, pattern = ".*noise.tsv",full.names = TRUE, recursive = TRUE)

for (i in 1:length(noise.path)) {
  noise.table.d<-read.csv(noise.path[i], sep = "\t", header = FALSE)  
  if(!exists("noise.table")){
    noise.table<-noise.table.d
  }else{
    noise.table<-rbind(noise.table,noise.table.d)
  }
}


noise.table.agg<-aggregate(V2~V1, noise.table, sum)

poi.co<-  (n.end-n.start)*0.02
poit<-noise.table[which(noise.table$V2>co.n & noise.table$V1>n.start & noise.table$V1<n.end),]
poi<-poit$V1[order(poit$V2, decreasing=TRUE)][1:round(poi.co)]
write.csv(poi, "poi.temp",row.names = FALSE)
