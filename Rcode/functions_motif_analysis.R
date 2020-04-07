## motif analysis using Bioconductor package
## Created by ZheZhen Wang
## Modified by Holly Yang # add teh parameter path
## Last updated: March 10, 2020


require("rtracklayer")
library("ggplot2")
library("PWMEnrich")

mywriteFasta <- function(enhancers, BSgenome, fasta.file, width.cut=5)
{
  if(class(BSgenome) != "BSgenome" )
    stop("the object 'BSgenome.Hsapiens.UCSC.hg19' requires a 'BSgenome' object ")
  
  if(class(enhancers)== "GRanges") Sequences=enhancers else {  
    if(class(enhancers) %in% c("matrix","data.frame")) {
      if(!all(c("start","end","peakNames") %in% colnames(enhancers))) 
        stop("please ensure that enhancers has the columns 'start','end','peakNames'")
      rgBED <- IRanges(start=enhancers$start,end=enhancers$end,names=enhancers$peakNames)
      Sequences <-  GenomicData(rgBED, chrom= enhancers$chr)
    } 
  }
  if (is(Sequences, "GRanges") & is.null(genome)) {
    stop("You have specified a RangedData object but no genome is specified")
  }
  if (is(Sequences, "GRanges")) {
    spSeq <- as.vector(chrom(Sequences))
    stSeq <- start(Sequences)
    edSeq <- end(Sequences)
  }  
  FastaXstring <- getSeq(BSgenome, 
                         names=spSeq, start = stSeq, end = edSeq, 
                         as.character = FALSE)
  names(FastaXstring) <- rownames(Sequences)            
  ## check the width
  if(any(width(FastaXstring) < width.cut))
    FastaXstring <- FastaXstring[-which(width(FastaXstring)< width.cut)]          
  writeXStringSet(FastaXstring, 
                  file = fasta.file, 
                  format="fasta", width=80)
}

wirte_table = function(res,top=0.10,saveN, path=NULL){
  myreport <- groupReport(res, top=top)
  tmp = as.data.frame(myreport)
  x <- which(tmp$top.motif.prop>0)
  tmp$fdr<- NA
  tmp$fdr[x] = p.adjust(tmp$p.value[x],'fdr')
  write.table(tmp,file = paste0(path, saveN,'.txt'),sep = '\t',row.names = F)
  return(tmp)
}

plotMotif = function(tmp,saveN,myreport,q = 0.05,score = 0.1,top = 0.1,save = TRUE, path=NULL){
  idx = row.names(tmp[!is.na(tmp$fdr) &tmp$fdr < q & tmp$top.motif.prop>top & tmp$raw.score>=score,])
  tmp2 = tmp[idx,]
  tmp3 = row.names(tmp2[!duplicated(tmp2$target),])
  len = ifelse(length(tmp3)>30,30,length(tmp3))
  plot(myreport[as.numeric(tmp3)][1:len], fontsize=7, id.fontsize=6)
  if(length(tmp3)>30) warning(paste0('the plot only contains top 30 motifs, ',
                                     length(tmp3)-30,' not shown'))
  if(save) dev.copy2pdf(file = paste0(path, saveN,'.pdf'))
  return(tmp3)
}

getplotMotif = function(FastaXstring,saveN,motif,q = 0.05,score = 0.1,top = 0.1,save = TRUE, path=NULL){
  res <-  motifEnrichment(FastaXstring, motif)
  save(res,file = paste0(path,'res+',saveN,'.rData'))
  tmp = wirte_table(res,top,saveN=saveN, path=path)
  #myreport <- groupReport(res)
  myreport <- groupReport(res, top=top)
  tmp2 = plotMotif(tmp,saveN=saveN,myreport,q,score,top,save = save,path=path)
  return(myreport)
}

motif_analysis = function(gr,BSgenome,saveN,
               motif = motifs,q = 0.05,score = 0.1,top = 0.1,
               save = TRUE, path=NULL){
  mywriteFasta(gr,BSgenome,
               fasta.file = paste0(path, saveN,'.fasta'))
  FastaXstring = readDNAStringSet(file=paste0(path, saveN,".fasta"),
                                  format="fasta")
  old = getplotMotif(FastaXstring, 
         saveN=paste0('motif_',saveN),motif,q,score,save = save, path=path)
  return(old)
}

# created on 5/11/2018
get_just_motif = function(gr,BSgenome,saveN,motif = motifs,path=NULL){
  mywriteFasta(gr,BSgenome,
               fasta.file = paste0(path, saveN,'.fasta'))
  FastaXstring = readDNAStringSet(file=paste0(path, saveN,".fasta"),
                                  format="fasta")
  res <-  motifEnrichment(FastaXstring, motif)
  return(res)
}

# created on 5/15/2018
get_instance_byq = function(res,q){
  idx = c()
  for(i in 1:length(FastaXstring)){
    sigmotif = sequenceReport(res,1)[sequenceReport(res,i)$p.value<q]$target
    if(target %in% sigmotif) idx = c(idx,i)
  }
  return(idx)
}
  
