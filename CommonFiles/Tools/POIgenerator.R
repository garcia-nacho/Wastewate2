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


noise.table.agg<-aggregate(V2~V1, noise.table, mean)
depth.table.agg<-aggregate(V3~V1, noise.table, mean)

poi.co<-  length(depth.table.agg$V3>20)*0.025

poit<-noise.table.agg[which(noise.table.agg$V2>co.n & noise.table.agg$V1>n.start & noise.table.agg$V1 < n.end & depth.table.agg$V3>20),]
poi<-poit$V1[order(poit$V2, decreasing=TRUE)][1:min(round(poi.co),nrow(poit))]
write.csv(poi, "poi.temp",row.names = FALSE)
