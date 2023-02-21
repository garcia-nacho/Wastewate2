library(data.table)
library("doParallel")
library("parallel")
library("foreach")
library(doSNOW)
library(progress)
library(seqinr)
library(umap)
library(ggplot2)

#src/bam2msa /media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/SpikeRef.fa ../fullwaste18/results/bam/waste018.43_Oslo.sorted.bam Spike:1260-2000 bam2msa


# Parallel generation of fastas -------------------------------------------

file.copy("/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/Tools/bam2msa","bam2msa")
file.copy("/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/SpikeRef.fa","SpikeRef.fa")
bam.files<-list.files(pattern = ".bam$",full.names = TRUE, recursive = TRUE)
dir.create("MSAFastas")
samples.to.analyze<-c(1:length(bam.files))

pb <- progress_bar$new(
  format = "File: :samp.pb [:bar] :elapsed | eta: :eta",
  total = length(samples.to.analyze),    # 100 
  width = 60)
samp <- samples.to.analyze

progress <- function(n){
  pb$tick(tokens = list(samp.pb = samp[n]))
} 

opts <- list(progress = progress)

gc()
cores.n<-detectCores()
cores<- cores.n -2

cluster.cores<-makeCluster(cores)
registerDoSNOW(cluster.cores)

out.par<-foreach(i=1:length(samples.to.analyze), .verbose=FALSE, .options.snow = opts) %dopar%{
system(paste("./bam2msa ", "./SpikeRef.fa ",bam.files[i]," Spike:1250-2250",
             " \\| cut -f 1-2 > ", gsub(".*/","MSAFastas/",gsub(".bam","_b2f.tsv",bam.files[i])), sep = ""))
}
stopCluster(cluster.cores)















df<-fread("/media/nacho/Data/wastewater/bam2msa/test.al.fasta.tsv", sep = "\t")
poi<-c(130,  81,  86, 105, 206, 207, 219, 223, 292, 332, 419)
#poi<-c(10:600)
pmatrix<-read.csv("/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/ProbMatrix.csv")
reference<-read.fasta("/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/SpikeRef.fa")
reference<-unlist(reference)






samples.to.analyze<-c(1:nrow(df))

pb <- progress_bar$new(
  format = "Read: :samp.pb [:bar] :elapsed | eta: :eta",
  total = length(samples.to.analyze),    # 100 
  width = 60)
samp <- samples.to.analyze

progress <- function(n){
  pb$tick(tokens = list(samp.pb = samp[n]))
} 

opts <- list(progress = progress)


###
gc()
cores.n<-detectCores()
cores<- cores.n -2

cluster.cores<-makeCluster(cores)
registerDoSNOW(cluster.cores)

out.par<-foreach(i=1:nrow(df), .verbose=FALSE, .options.snow = opts) %dopar%{
  ins<-which(unlist(strsplit(df$ref_msa[i], ""))=="-")
  if(length(ins)>0 ){
    dum.out<-unlist(strsplit(df$`#query_msa`[i], ""))[-ins][poi]
  }else{
    dum.out<-unlist(strsplit(df$`#query_msa`[i], ""))[poi]
  }
  out.ohe<-rep(0,(length(poi)*5))
  out.ohe[match(dum.out, c("A","T","C","G","-"))+seq(from=0, to=(length(poi)-1)*5, by=5)]<-1
  return(list(dum.out,out.ohe))
}
stopCluster(cluster.cores)

lineages<-lapply(out.par, function(x)x[[1]])
lineages.ohe<-lapply(out.par, function(x)x[[2]])

lineages<-do.call( rbind, lineages)
lineages.ohe<-do.call( rbind, lineages.ohe)
lineages<-as.data.frame(lineages)
lineages.ohe<-as.data.frame(lineages.ohe)

colnames(lineages)<-paste("P",poi,sep = "")
colnames(lineages.ohe)<-paste(rep(poi+1250,5)[order(rep(poi+1250,5))],c("A","T","C","G","D"),sep = "")
#write.csv(lineages,"Lineages.csv",row.names = FALSE)








# Pango Lineage probability calculation -----------------------------------------

rownames(pmatrix)<-pmatrix$X
pmatrix$X<-NULL
pmatrix<-pmatrix[which(row.names(pmatrix) %in% colnames(lineages.ohe)),]

mut.aa<-vector()
for (i in 1:nrow(pmatrix)) {
  if(length(which( gsub("^[0-9]{1,}","",rownames(pmatrix)[i])%in% c("A","T","C","G")))>0){
    r<-translate(reference)
    rm<-reference
    rm[as.numeric(gsub(".$","",rownames(pmatrix)[i]))]<-tolower(gsub("^[0-9]{1,}","",rownames(pmatrix)[i])) 
    m<-translate(rm)
    mut.site<-which(r!=m)
    if(length(mut.site)>0){
     mut.aa<-c(mut.aa, paste("S:",r[mut.site],mut.site,m[mut.site],sep = ""))   
    }else{
      mut.aa<-c(mut.aa,NA)
    }
  }else{
    mut.aa<-c(mut.aa,NA)
  }
}


pmatrix[which(is.na(pmatrix), arr.ind = TRUE)]<-0

plinmatrix<-as.matrix(lineages.ohe) %*% as.matrix(pmatrix)

lin.index<-apply(plinmatrix, 1, function(x)which(x==max(x,na.rm = TRUE)))
lin.max<-apply(plinmatrix, 1, function(x)max(x,na.rm = TRUE))
lin.tab<-unlist(lapply(lin.index, function(x) paste(colnames(plinmatrix)[x],collapse = "/")))

# Compacted Lineage calculator --------------------------------------------

lineages.compact<- as.data.frame(apply(lineages,1, function(x)paste(x, collapse = "")))
colnames(lineages.compact)<-"Seq"
#lineages.compact
lineages.compact$PangoLineage<-lin.tab
lineages.compact$Count<-0
lineages.compact$Score<-lin.max
lineages.compact.agg<-aggregate(Count ~Seq+PangoLineage+Score+ScoreDif, lineages.compact, length)
lineages.compact.agg$Positions<-paste(gsub("P","",colnames(lineages)),collapse = "/")

index<-vector()
for (i in 1:nrow(lineages.compact.agg)) {
   index<-c(index,which(lineages.compact$Seq==lineages.compact.agg$Seq[i])[1])
}

plin.umap<-umap(plinmatrix[index,])
lineages.compact.agg$X<-plin.umap$layout[,1]
lineages.compact.agg$Y<-plin.umap$layout[,2]

# Single Mutation Calculation ----------------------------------------------

colnames(lineages.ohe)<-mut.aa
lineages.ohe<-lineages.ohe[,-which(is.na(mut.aa))]

sin.mut<-as.data.frame(colnames(lineages.ohe)) 
colnames(sin.mut)<-"Mutation"
sin.mut$Count<-NA
sin.mut$Ratio<-NA

for (i in 1:nrow(sin.mut)) {
  sin.mut$Count[i]<-sum(lineages.ohe[,which(colnames(lineages.ohe)==sin.mut$Mutation[i])])
  sin.mut$Ratio[i]<-sum(lineages.ohe[,which(colnames(lineages.ohe)==sin.mut$Mutation[i])])/nrow(lineages.ohe)
}


# Plotting ---------------------------------



ggplot(lineages.compact.agg)+
  geom_point(aes(X,Y,col=PangoLineage),alpha=0.2)+
  scale_colour_manual(values = rainbow(length(unique(lineages.compact.agg$PangoLineage))))+
  theme_minimal()
