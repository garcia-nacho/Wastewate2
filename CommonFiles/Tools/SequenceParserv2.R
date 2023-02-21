library(seqinr)
library(digest)

df<-read.csv("/media/nacho/Data/wastewater/AFMode/refs/references.csv",sep = ";")
sqs<-read.fasta("/media/nacho/Data/wastewater/AFMode/refs/references.aligned.fasta")

df<-df[-which(duplicated(names(sqs))),]
sqs<-sqs[-which(duplicated(names(sqs)))]

df<-df[which(df$totalMissing<300),]
sqs<-sqs[which(names(sqs) %in% gsub(" .*","",df$seqName))]

sqs<-lapply(sqs, function(x) x[c((21563+1250):(21563+2250))] )

hashes<-lapply(sqs, function(x) digest(paste(x,collapse = ""),algo="md5"))

df$hash<-unlist(hashes)

uniqueregions<-unique(hashes)
df$QC<-NA
df$Corrected<-NA
unique.index<-vector()

for (i in 1:length(uniqueregions)) {
  pangos<- df$Nextclade_pango[ which(df$hash== uniqueregions[i])]
  if(length(unique(pangos))==1){
    df$QC[which(df$hash==uniqueregions[i])]<-"OK"
    df$Corrected[which(df$hash==uniqueregions[i])]<-df$Nextclade_pango[which(df$hash==uniqueregions[i])]
    unique.index<-c(unique.index,which(df$hash==uniqueregions[i])[1])
  }else{
    df$QC[which(df$hash==uniqueregions[i])]<-"Warning"
    
    lineages<- df$Nextclade_pango[which(df$hash==uniqueregions[i])]
    lineages<-strsplit(lineages,"\\.")
    noe<-max(unlist(lapply(lineages, length)))
    
    running<-TRUE
    while(running){
      new.lineages<- unlist(lapply(lineages, function(x) paste(x[(1:min(length(x),noe))], collapse = ".") ))  
      if(length(unique(new.lineages))==1){
        new.lin<-unique(new.lineages)
        unique.index<-c(unique.index,which(df$hash==uniqueregions[i])[1])
        running<-FALSE
        
        
      }else{
        new.lin<-NA
        if(noe==1) running<-FALSE
      }
        noe<-noe-1
    }
    if(!is.na(new.lin))df$Corrected[which(df$hash==uniqueregions[i])]<-paste(new.lin,".X",sep = "")
    
  }
}
df<-df[unique.index,]

sqs<-sqs[which(names(sqs) %in% gsub(" .*","",df$seqName))]

total.lineages<-as.data.frame(table(df$Corrected))

names(sqs)<-paste("SID","_",as.character(df$Corrected),sep = "")
write.fasta(sqs, "/media/nacho/Data/wastewater/AFMode/references.1Kregion.fasta",names = names(sqs))
