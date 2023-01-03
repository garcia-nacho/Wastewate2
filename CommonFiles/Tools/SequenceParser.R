library(seqinr)
library(digest)

fastas<-seqinr::read.fasta("Seqs_aligned.fa")
metadata<-read.csv("Seqs_aligned_NextClade.csv",sep = ";", stringsAsFactors = FALSE)

metadata<-metadata[which(metadata$totalMissing<200),]

fastas<-fastas[which(names(fastas) %in% metadata$seqName)]

fastas.trimmed<-lapply(fastas, function(x) x[c(21563:25384)] )

for (i in 1:length(fastas.trimmed)) {
  if(length(which(fastas.trimmed[[i]]=="-"))>0) fastas.trimmed[[i]] <- fastas.trimmed[[i]][-which(fastas.trimmed[[i]]=="-")] 
  names(fastas.trimmed)[i]<-paste("SID",i,"_",as.character(metadata$Nextclade_pango[which(metadata$seqName==names(fastas.trimmed)[i])]),sep = "")
  
}

#CleanDups
hashes<-lapply(fastas.trimmed, function(x) digest(paste(x,collapse = ""),algo="md5"))
hashes<-unlist(hashes)
dups<-names(hashes)[which( duplicated(hashes))]
fastas.trimmed<-fastas.trimmed[-which(names(fastas.trimmed) %in% dups )]

write.fasta(fastas.trimmed,"/media/nacho/Data/wastewater/LastkSeqs.fasta", names = names(fastas.trimmed))
