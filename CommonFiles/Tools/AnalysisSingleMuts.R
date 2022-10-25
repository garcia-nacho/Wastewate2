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
library(reshape2)

path=commandArgs(TRUE)
co.n<-0.15
co.n<-as.numeric(path[2])
n.start<-1
n.end<- 3822
path<-path[1]


noise.path<-list.files(path, pattern = ".*noise.tsv",full.names = TRUE, recursive = TRUE)

for (i in 1:length(noise.path)) {
  noise.table.d<-read.csv(noise.path[i], sep = "\t", header = FALSE)  
  if(!exists("noise.table")){
    noise.table<-noise.table.d
  }else{
    noise.table<-rbind(noise.table,noise.table.d)
  }
}

#co.n<-co.n/length(noise.path)

noise.table.agg<-aggregate(V2~V1, noise.table, mean)
depth.table.agg<-aggregate(V3~V1, noise.table, mean)



poi.co<-  length(depth.table.agg$V3>20)*0.025
poit<-noise.table.agg[which(noise.table.agg$V2>co.n & noise.table.agg$V1>n.start & noise.table.agg$V1 < n.end & depth.table.agg$V3>20),]
poi<-poit$V1[order(poit$V2, decreasing=TRUE)][1:min(round(poi.co),nrow(poit))]

bam.path.total<-list.files(path, pattern = ".*bam$",full.names = TRUE, recursive = TRUE)

cores<-as.numeric(detectCores())-2

allpos<-TRUE

for (b in 1:length(bam.path.total)) {
  bam.path<-bam.path.total[b]

if(length(poi)>0){
  
  samples.to.analyze<-poi
  
  pb <- progress_bar$new(
    format = "Position: :samp.pb [:bar] :elapsed | eta: :eta",
    total = length(samples.to.analyze),    # 100 
    width = 60)
  
  samp <- samples.to.analyze
  
  progress <- function(n){
    pb$tick(tokens = list(samp.pb = samp[n]))
  } 
  
  opts <- list(progress = progress)
  
  ###
  gc()
  cluster.cores<-makeCluster(cores)
  registerDoSNOW(cluster.cores)
  
  out.par<-foreach(k=1:length(poi), .verbose=FALSE, .options.snow = opts) %dopar%{
    system(paste(path,"bbasereader -p ",poi[k]  ,
                 " -f ", bam.path, "> ", gsub(".*/",path,gsub("\\.sorted.bam","_P",bam.path)), poi[k],".tsv",sep = ""))
  }
  
  stopCluster(cluster.cores)
  



}
}
position.files<-list.files(path, full.names = TRUE, pattern = "_P.*\\.tsv")

basediff<-as.data.frame(matrix(NA, ncol = 6, nrow=length(position.files)))
colnames(basediff)<-c("Base", "A","T","C","G", "Count")
for (i in 1:length(position.files)) {
  dummy<-fread(position.files[i],sep = "\t", header = FALSE)
  dummy<-as.data.frame(dummy[,-1])
  colnames(dummy)<-c("Position","Base")
  
  if(length(which(dummy$Base=="D"))>0) dummy<-dummy[-which(dummy$Base=="D"),]
  if(length(which(dummy$Base=="I"))>0) dummy<-dummy[-which(dummy$Base=="I"),]
  
  if(nrow(dummy)>0){
    dummy<-as.data.frame(dummy)
    basediff$Base[i]<-unique(dummy$Position)[1]
    dummy<-as.data.frame(table(dummy$Base))
    
    basediff$Count[i]<-sum(dummy$Freq)
    
    if(length(which(dummy$Var1=="A"))>0){
      basediff$A[i]<-dummy$Freq[which(dummy$Var1=="A")]/sum(dummy$Freq)
    }else{
      basediff$A[i]<-0
    }
    
    if(length(which(dummy$Var1=="T"))>0){
      basediff$T[i]<-dummy$Freq[which(dummy$Var1=="T")]/sum(dummy$Freq)
    }else{
      basediff$T[i]<-0
    }
    
    if(length(which(dummy$Var1=="C"))>0){
      basediff$C[i]<-dummy$Freq[which(dummy$Var1=="C")]/sum(dummy$Freq)
    }else{
      basediff$C[i]<-0
    }
    
    if(length(which(dummy$Var1=="G"))>0){
      basediff$G[i]<-dummy$Freq[which(dummy$Var1=="G")]/sum(dummy$Freq)
    }else{
      basediff$G[i]<-0
    }
    basediff$File[i]<-gsub("_P.*","",gsub(".*/","",position.files[i]))
  }
  
  
}
basediff.melted<-melt(basediff, id=c("Base","Count","File"))
s<-read.fasta("/home/docker/CommonFiles/reference/SpikeRef.fa")
s<-as.character(s$Spike)

basediff.melted$CI<-NA
basediff.melted$AA<-NA

basediff.melted$CI<-qt(0.975,df=basediff.melted$Count-1)*(sd(basediff.melted$value,na.rm = TRUE))/sqrt(basediff.melted$Count)    

basediff.melted<-basediff.melted[which(basediff.melted$CI<basediff.melted$value),]

bases.to.check<-unique(paste(basediff.melted$Base,basediff.melted$variable,sep = "-"))
original<-translate(s)

for (i in 1:length(bases.to.check)){
mutated<-s
mutated[as.numeric(gsub("-.*","",bases.to.check[i]))]<-tolower(gsub(".*-","",bases.to.check[i]))
mutated<-translate(mutated)
if(length(which(original!=mutated))){
  basediff.melted$AA[which(paste(basediff.melted$Base,basediff.melted$variable,sep = "-") == bases.to.check[i])]<- 
    paste("S:",original[which(original!=mutated)], which(original!=mutated), mutated[which(original!=mutated)],sep = "")
}else{
  basediff.melted$AA[which(paste(basediff.melted$Base,basediff.melted$variable,sep = "-") == bases.to.check[i])]<-"SynMut"
}

}

basediff.melted.bk<-basediff.melted
  
if(length(which(basediff.melted$AA=="SynMut"))>0) basediff.melted<-basediff.melted[-which(basediff.melted$AA=="SynMut"),]
if(length(grep("\\*", basediff.melted$AA))>0) basediff.melted<-basediff.melted[-grep("\\*", basediff.melted$AA),]

basediff.melted$AApos<-gsub(".$","",basediff.melted$AA)
basediff.melted$AAnew<-gsub("S:.[0-9]","",basediff.melted$AA)
basediff.melted$AAnew<-gsub("^[0-9]","",basediff.melted$AAnew)
basediff.melted$AAnew<-gsub("^[0-9]","",basediff.melted$AAnew)
basediff.melted$AAnew<-gsub("^[0-9]","",basediff.melted$AAnew)


 

basediff.melted$AApos<-factor( basediff.melted$AApos, levels =unique(basediff.melted$AApos[order(as.numeric(gsub(".$","",gsub("S:.","",basediff.melted$AApos))))]) )


colnames(basediff.melted)[9]<-"Mutation"
if(length(unique(basediff.melted$File))<10){
ggplot(basediff.melted)+
  geom_bar(aes(AApos,value, fill=Mutation), stat = "identity")+
  theme_minimal()+
  facet_wrap(~File,ncol = 1)+
    scale_fill_manual(values = rainbow(length(unique(basediff.melted$Mutation))))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Position Spike")+
    ylab("Ratio mutant reads")
  ggsave(paste(path,"SingleMutationsBarplots.pdf",sep = "") ,width =  8.27, height =  11.69)
}

if(length(unique(basediff.melted$File))>10){
  files.to.plot<-unique(basediff.melted$File)
printing<-TRUE
counter<-0
start<-1
stop<-10
  while(printing){
    counter<-counter+1
    if(stop>=length(files.to.plot)){
      stop<-length(files.to.plot)
      printing<-FALSE
    }
    ggplot(basediff.melted[which(basediff.melted$File %in% files.to.plot[start:stop]),])+
      geom_bar(aes(AApos,value, fill=Mutation), stat = "identity")+
      theme_minimal()+
      facet_wrap(~File,ncol = 1)+
      scale_fill_manual(values = rainbow(length(unique(basediff.melted$Mutation))))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      xlab("Position Spike")+
      ylab("Ratio mutant reads")
    
    ggsave(paste(path,"SingleMutationsBarplots_part",counter, ".pdf",sep = "") ,width =  8.27, height =  11.69)
    start<-stop+1
    stop<-start+9
    if(start>length(files.to.plot)) printing<-FALSE
}
  
}

file.remove(position.files)
file.remove(paste(path,"bbasereader",sep = ""))
write_xlsx(basediff.melted.bk, paste(path,"SingleMutationsResults.xlsx",sep = ""))     