library(Matrix)
library(reshape2)
library(gplots)
library(DESeq2)
library(ggplot2)
library(monocle)
library(DESeq)
library(RColorBrewer)
#library(Vennerable)
library(edgeR)
library(EDASeq)
#library(vioplot)
library(scatterplot3d)
library(scde)
library(tsne)
library(statmod)
library(genefilter)
library(caret)
library(Biobase)
library(apcluster)
library(clValid)
library(ggdendro)
library(heatmap3)
library(sva)
#library(jackstraw)
library(fastICA)
library(monocle)
library(pheatmap)
library(psych)
library(pvclust)
#library(SIMLR)
#library(diffusionMap)
#setwd('~/rosetting')

std.stages.col.list=list('ring'='green','late.ring'='green4',
                         'early.trophozoite'='blue',
                         'late.trophozoite'='blue4',
                         'early.schizont'='deeppink3',
                         'schizont'='red')
abbrev.std.stages.col.list=list('R'='cyan','LR'='green4','ET'='blue4',
                                'T'='orange','ES'='purple','S'='firebrick')
abbrev.std.stages.col.vec=c('R'='cyan','LR'='green4','ET'='blue4',
                            'T'='orange','ES'='purple','S'='firebrick')
std.stages.col.legend.str=c('ring','late.ring','early.trophozoite',
                            'late.trophozoite','early.schizont','schizont')
std.stages.col.abbrev.legend.str=c('R','LR','ET','T','ES','S')
std.stages.ordered.str=c('R','LR','ET','T','ES','S')
std.stages.col.legend.fill=c('cyan','green4','blue4','orange','purple','firebrick')
#markers.peaking.stages.list=list("Early ring transcripts"='R,LR,S','Cytoplasmic Translation machinery'='R,LR,ET',
#"Transcription machinery"='R,LR,ET',"Actin myosin motors"='S',"Ribonucleotide synthesis"='R,LR,ET', 
#"Glycolytic pathway"='R,LR,ET',"Merozoite Invasion"='S',"Proteasome"='T,ES',"DNA replication"='T,ES', 
#"Organellar Translation machinery"='S',"Mitochondrial"='T,ES',"TCA cycle"='T,ES', "Deoxynucleotide synthesis"='T,ES')
markers.peaking.stages.list=list("Early ring transcripts"='R,S',
                                 'Cytoplasmic Translation machinery'='LR,ET',
                                 "Transcription machinery"='R,LR,ET',
                                 "Actin myosin motors"='R,ES,S',
                                 "Ribonucleotide synthesis"='R,LR,ET', 
                                 "Glycolytic pathway"='R, LR, ET',
                                 "Merozoite Invasion"='R,ES,S',
                                 "Proteasome"='T,ES',"DNA replication"='T,ES', 
                                 "Organellar Translation machinery"='S',
                                 "Mitochondrial"='T,ES',"TCA cycle"='T,ES', 
                                 "Deoxynucleotide synthesis"='T')
surface.markers.col.list=c('rifin'='red','var'='blue','surfin'="yellow",
                           'stevor'="green")
schizont.stage.sub.group.samples.vec=c("sample_40h_4_one","sample_40h_11_one",
                                       "sample_40h_15_one" ,"sample_40h_18_one",
                                       "sample_40h_20_one","sample_40h_26_one",
                                       "sample_40h_27_one","sample_40h_28_one",
                                       "sample_40h_29_one","sample_40h_sc_sp_1_two",
                                       "sample_40h_sc_sp_2_two","sample_40h_sc_sp_3_two",
                                       "sample_40h_sc_sp_6_two","sample_40h_sc_sp_7_two",
                                       "sample_40h_sc_sp_8_two","sample_40h_10_two",
                                       "sample_40h_18_two","sample_40h_20_two",
                                       "sample_40h_21_two","sample_40h_22_two",
                                       "sample_40h_23_two","sample_40h_24_two",
                                       "sample_40h_26_two","sample_40h_28_two",
                                       "sample_40h_31_two","sample_40h_32_two",
                                       "sample_40h_33_two")
#sub.populations.grp.names=c('SC.grp.1','SC.grp.2','SC.grp.3','SC.grp.4','SC.grp.5','SC.grp.6','SC.grp.7','SC.grp.8')
sub.populations.grp.names=c('SP.1','SP.2','SP.3','SP.4','SP.5','SP.6','SP.7','SP.8')
sub.populations.grp.order.vec=c('8'='SP.1','1'='SP.2','7'='SP.3','5'='SP.4',
                                '2'='SP.5','6'='SP.6','4'='SP.7','3'='SP.8')
sub.populations.grp.order.names.vec=c('SP.1'='8','SP.2'='1','SP.3'='7','SP.4'='4','SP.5'='2','SP.6'='6',
                                      'SP.7'='5','SP.8'='3')
bulk.sub.grps.order.vec=c('SP:1'='3', 'SP:2'='5', 'SP:3'='6', 'SP:4'='7', 
                          'SP:5'='1' ,'SP:6'='4' ,'SP:7'='2')
abbrev.std.stages.timepoints.vec=c('R' = '8-12' ,'LR' = '12-16', 'ET'= '20-24', 
                                   'T' = '30-34', 'ES'= '36-40', 'S' = '44-48')
ordered.sp.cols=c("SP.1"="#FFFFB3", "SP.2"="#80B1D3", "SP.3"="#FCCDE5", 
                  "SP.4"="#FB8072", "SP.5"="#B3DE69", "SP.6"="#FDB462", 
                  "SP.7"="#BEBADA", "SP.8"="#8DD3C7")
sp.outlier.limits=c('Actin myosin motors'=2500,'Cytoplasmic Translation machinery'=400,
                     'Early ring transcripts'=450,'Merozoite Invasion'=750)
sp.outlier.limits=c('Actin myosin motors'=2500,'Cytoplasmic Translation machinery'=1500,
                    'Early ring transcripts'=4500,'Merozoite Invasion'=1500)
gene.cat.names=c('Actin myosin motors' ,  'Cytoplasmic Translation machinery' ,  'Deoxynucleotide synthesis' , 
  'DNA replication' ,  'Early ring transcripts' ,  'Glycolytic pathway' ,  'Merozoite Invasion' ,  
  'Mitochondrial' ,  'none' ,  'Organellar Translation machinery' ,  'Proteasome' ,  'Ribonucleotide synthesis' , 
  'TCA cycle' ,  'Transcription machinery')
gene.cat.col.names=c("#54db15" ,  "#ba5407" ,  "#49ddc0" ,  "#8cd9f7" ,"#cda9f2" ,  "#34f9aa" ,  "#ffffba" , 
  "#b25701","#e0f7a5" ,  "#c6f98b" ,  "#8ffcbd" ,  "#00efb3","#66c6c9" ,  "#7cef23")
gene.cat.col.names=c( "#079933",  "#ecfc5f" ,  "#d16729"  ,  'firebrick' ,"#e76ff2"  ,  "#e0597d" ,  "#2a3fb7" , 
                     "#7eb7e5","gray" ,  "#7be0b4" ,  "#edbe95" ,  "#4291c9","#5497ce" ,  "#efd5a5")
gene.cat.col.names=c( "#079933",  "#ecfc5f" ,  "gray28"  ,  'firebrick' ,"#e76ff2"  ,  "#e0597d" ,  "#2a3fb7" , 
                      "coral","gray" ,  "#7be0b4" ,  "#edbe95" ,  "#4291c9","#5497ce" ,  "purple")
names(gene.cat.col.names)=gene.cat.names


#Mapping proportions
plot.samples.mapping.proportions.barplot=function(df,cut.off=0,title.str='Test'){
  samples=rownames(df)
  input.reads.vec=as.numeric(as.character(df[,'Number_of_input_reads_']))
  uniquely.mapped.reads.vec=as.numeric(as.character(df[,'Uniquely_mapped_reads_number_']))
  multimaps.vec=as.numeric(as.character(df[,'Number_of_reads_mapped_to_multiple_loci_']))
  multimaps.too.many.loci.vec=as.numeric(as.character(df[,'Number_of_reads_mapped_to_too_many_loci_']))
  total.multimaps=multimaps.vec+multimaps.too.many.loci.vec
  unmapped.reads.vec=input.reads.vec-(uniquely.mapped.reads.vec+total.multimaps)
  mapping.df=data.frame(input.reads=input.reads.vec,unique.reads=uniquely.mapped.reads.vec,multimaps=total.multimaps)
  rownames(mapping.df)=rownames(df)
  mapping.df$unmapped.reads=unmapped.reads.vec
  unmapped.too.short.reads.vec=as.numeric(as.character(df[,'X._of_reads_unmapped._too_short_']))
  unmapped.too.many.mismatches.reads.vec=as.numeric(as.character(df[,'X._of_reads_unmapped._too_many_mismatches_']))
  unmapped.other.reasons.vec=as.numeric(as.character(df[,'X._of_reads_unmapped._other_']))
  unmapped.df=data.frame(unmapped.too.short.reads=unmapped.too.short.reads.vec,unmapped.too.many.mismatches.reads=unmapped.too.many.mismatches.reads.vec,unmapped.other.reasons=unmapped.other.reasons.vec)
  unmapped.df=(unmapped.df/(apply(unmapped.df,1,sum)))*100
  rownames(unmapped.df)=rownames(mapping.df)
  temp.mat=as.matrix(mapping.df)
  filt.mat =temp.mat[which(temp.mat[,'input.reads']>=cut.off),]
  #temp.mat=filt.mat
  read.mapping.mat=as.matrix(data.frame(temp.mat[,2:4]))
  names.label=gsub(pattern = '^sample_','',rownames(read.mapping.mat))
  #kit.str=as.character(df[rownames(temp.mat),'kit'])
  #kit.col.factor=get.col.factor(col.factor = kit.str)
  barplot.input.reads=as.numeric(temp.mat[,1])
  barplot.unique.reads=ifelse(as.numeric(temp.mat[,2])==0,1,as.numeric(temp.mat[,2]))
  #temp.barplot=barplot2(barplot.input.reads,main='No. of reads',col=c('lightblue'),ylab='',xlab='',beside=F,log='y',las=2, cex.names = .6,names.arg = names.label)
  title.str=convert.to.title.case(title.str)
  #temp.barplot=barplot2(barplot.input.reads,main=title.str,ylab='',xlab='',beside=F,log='y',las=2, cex.names = .2,names.arg = names.label,col=kit.col.factor$col.str)
  #legend('topright',legend=kit.col.factor$legend.str,fill=kit.col.factor$legend.col,cex=.5)
  #temp.barplot=barplot2(barplot.unique.reads,main='No. of unique reads',col=c('lightblue'),ylab='',xlab='',beside=F,log='y',las=2, cex.names = .6,names.arg = names.label)
  #temp.barplot=barplot2(barplot.unique.reads,main=title.str,ylab='',xlab='',beside=F,log='y',las=2, cex.names = .2,names.arg = names.label,col=kit.col.factor$col.str)
  #legend('topright',legend=kit.col.factor$legend.str,fill=kit.col.factor$legend.col,cex=.5)
  input_reads=as.numeric(as.character(temp.mat[,1]))
  temp.mat=round(temp.mat/input_reads*100,2)
  temp.mat=data.frame(temp.mat[,2:4])
  temp.mat=as.matrix(temp.mat)
  #temp.barplot.2=barplot2(t(temp.mat),main=title.str,col=c('red','green','blue'),ylab='',xlab='',beside=F,las=2, cex.names = .2,names.arg = names.label)
  temp.barplot.2=barplot2(t(temp.mat),main=title.str,col=c('red','green','blue'),ylab='',xlab='',beside=F,las=2, 
                          cex.names = .2,las=2,border = NA,xaxt='none')
  #text(temp.barplot.2, par("usr")[3], labels = samples.names, srt = 90, adj = c(1.1,1.1), xpd = T, cex=.3)
  #plot.new()
  #legend('center',legend=c('Unique','Multimap','Unmapped'),fill=c('red','green','blue'),cex=1.5,box.lty = 0)
  #legend('center',legend=c('','',''),fill=c('red','green','blue'),cex=1.5,box.lty = 0)
  #legend('center',legend=c('','','',''),fill=c('lightblue','red','green','blue'),cex=1.5,box.lty = 0)
  #Unmapped proportions plot
  #temp.barplot.3=barplot(t(as.matrix(unmapped.df)),main='Unmapped proportions',col=c('red','green','blue'),ylab='Proportion (%)',xlab='',xaxt='none',beside=F,ylim=c(0,115))
  #text(temp.barplot.2, par("usr")[3], labels = rownames(unmapped.df), srt = 90, adj = c(1.1,1.1), xpd = T, cex=.3)
  #legend('topcenter',legend=c('Too short reads','Too many mismatches','Other reasons'),fill=c('red','green','blue'),cex=.8)
  header.str=c('Group','Input.reads','%unique.reads', '%multimaps', '%unmapped.reads', '%unmapped_short', '%unmapped_mismatches')
  header.str=paste(header.str,collapse=' & ')
  header.str=paste(header.str,'\\ \\hline')
  #hist.breaks=seq(0,max(gene.counts)+200,by = 200)
  #hist(x = gene.counts,xlab = 'Counts', main ='Gene detected',breaks = hist.breaks)
  #plot.new()
  inhouse.samples=samples[grepl('_inhouse_',samples)]
  len.inhouse.samples=length(inhouse.samples)
  x.points=c()
  y.points=c()
  labels.str=c()
#   for(m in 1:len.inhouse.samples){
#     
#     in.house.sample=inhouse.samples[m]
#     
#     in.house.sample.reads=as.numeric(df[in.house.sample,'Uniquely_mapped_reads_number_'])
#     
#     nextera.sample=gsub(pattern = '_inhouse',replacement = '',in.house.sample)
#     
#     nextera.sample.reads=as.numeric(df[nextera.sample,'Uniquely_mapped_reads_number_'])
#     
#     x.points=append(x.points,in.house.sample.reads,length(x.points))
#     
#     y.points=append(y.points,nextera.sample.reads,length(y.points))
#     
#     labels.str=append(labels.str,paste(in.house.sample,nextera.sample,sep=':'),length(labels.str))
#     
#     
#   }
#   
#   #plot(x=log10(x.points),y=log10(y.points),xlab='Log10 (Tn5 reads)',ylab='Log10 (Nextera reads)',main=title.str,pch=19,cex=1.0)
# 
#   batch.five.reads.df=data.frame(tn5=x.points,nextera=y.points)
#   
#   rownames(batch.five.reads.df)=labels.str
#   
#   barplot.col.list=get.col.factor(colnames(batch.five.reads.df))
#   
#  
#   batch.five.reads.prop.df=batch.five.reads.df/as.numeric(apply(batch.five.reads.df,1,sum))
#   
#   show(head(batch.five.reads.prop.df))
  #barplot2(t(batch.five.reads.prop.df),beside = F,las=2,cex.names =.2,log = 'y',main=title.str,col = barplot.col.list$col.str)
  #legend('topright',legend=barplot.col.list$legend.str,fill=barplot.col.list$legend.col,cex=.5)
  out.list=list(mapping=mapping.df,unmapped=unmapped.df)
  return(out.list)
}

#Compares the batch eight samples libraries and mapping profiles
#No. of detected genes
analyze.batch.eight.samples=function(batch.eight.rpkm.df,batch.eight.counts.df,batch.eight.meta.df){
  samples=intersect(colnames(batch.eight.rpkm.df),colnames(batch.eight.counts.df))
  rpkm.df=batch.eight.rpkm.df[,samples]
  count.df=batch.eight.counts.df[,samples]
  meta.df=batch.eight.meta.df[samples,]
  meta.df$staining=ifelse(grepl('_Vy_',samples),'vy','mito')
  meta.df$kit=ifelse(grepl('_inhouse_',samples),'tn5','nextera')
  meta.list=split(meta.df,f=meta.df$development.stage)
  stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  stages=intersect(ordered.development.stages,stages)
  len.stages=length(stages)
  for(m in 1:len.stages){
    stage=stages[m]
    temp.meta.df=meta.list[[stage]]
    temp.samples=rownames(temp.meta.df)
    temp.meta.list=split(temp.meta.df,f=temp.meta.df$staining)
    stains=names(temp.meta.list)
    len.stains=length(stains)
    for(n in 1:len.stains){
      stain=stains[n]
      temp.stain.meta.df=temp.meta.list[[stain]]
      title.str=convert.to.title.case(paste(stage,stain,sep='_'))
      plot.samples.mapping.proportions.barplot(df = temp.stain.meta.df,title.str = title.str)
      samples=rownames(temp.stain.meta.df)
      inhouse.samples=samples[grepl('_inhouse_',samples)]
      len.inhouse.samples=length(inhouse.samples)
      for(l in 1:len.inhouse.samples){
        in.house.sample=inhouse.samples[l]
        nextera.sample=gsub(pattern = '_inhouse',replacement = '',in.house.sample)
        x.points=as.numeric(batch.eight.rpkm.df[,in.house.sample])
        y.points=as.numeric(batch.eight.rpkm.df[,nextera.sample])
        cor.score=round(cor(x =x.points ,y= y.points,method = 'spearman'),2)
        x.points=log10(x.points+1)
        y.points=log10(y.points+1)
        temp.title.str=paste(title.str,'(Corr: ',cor.score,')',collapse = '')
        plot(x=x.points,y=y.points,ylab=nextera.sample,xlab=in.house.sample,main=temp.title.str,pch=19,cex=.5)
      }
    }
  }
}


plot.samples.unique.mapping.hist=function(df,cut.off=0,title.str='Test'){
  input.reads.vec=as.numeric(as.character(df[,'Number_of_input_reads_']))
  uniquely.mapped.reads.vec=as.numeric(as.character(df[,'Uniquely_mapped_reads_number_']))
  multimaps.vec=as.numeric(as.character(df[,'Number_of_reads_mapped_to_multiple_loci_']))
  multimaps.too.many.loci.vec=as.numeric(as.character(df[,'Number_of_reads_mapped_to_too_many_loci_']))
  total.multimaps=multimaps.vec+multimaps.too.many.loci.vec
  unmapped.reads.vec=input.reads.vec-(uniquely.mapped.reads.vec+total.multimaps)
  mapping.df=data.frame(input.reads=input.reads.vec,unique.reads=uniquely.mapped.reads.vec,multimaps=total.multimaps)
  rownames(mapping.df)=rownames(df)
  mapping.df$unmapped.reads=unmapped.reads.vec
  temp.mat=as.matrix(mapping.df)
  filt.mat =temp.mat[which(temp.mat[,'input.reads']>=cut.off),]
  read.mapping.mat=as.matrix(data.frame(temp.mat[,2:4]))
  names.label=gsub(pattern = '^sample_','',rownames(read.mapping.mat))
  barplot.input.reads=as.numeric(temp.mat[,1])
  barplot.unique.reads=ifelse(as.numeric(temp.mat[,2])==0,1,as.numeric(temp.mat[,2]))
  input_reads=as.numeric(as.character(temp.mat[,1]))
  total.mapped.reads.counts.vec=as.numeric(as.character(temp.mat[,1]))+as.numeric(as.character(temp.mat[,2]))
  read.counts.breaks=seq(from = 0,to = max(total.mapped.reads.counts.vec)+100000,by = 100000)
  hist(total.mapped.reads.counts.vec,main = '',xlab='',ylab='',read.counts.breaks,las=2)
  temp.mat=round(temp.mat/input_reads*100,2)
  temp.mat=data.frame(temp.mat[,2:4])
  temp.mat=as.matrix(temp.mat)
  total.mapped.reads.vec=as.numeric(as.character(temp.mat[,1]))+as.numeric(as.character(temp.mat[,2]))
  temp.barplot.2=barplot2(t(temp.mat),main='',col=c('red','green','blue'),ylab='',xlab='',beside=F,las=2, cex.names = .6,names.arg = rep('',times=length(names.label)))
  plot.new()
  legend('center',legend=c(rep('',times = 3)),fill=c('red','green','blue'),cex = 1.5,box.lty = 0)
}

plot.samples.chrom.mapping.proportions.barplot.for.subset.samples=function(meta.df,chrom.mapping.df,title.str='test',
                                                                           breaks.interval=100){
  no.samples=dim(meta.df)[1]
  breaks=seq(1,to =no.samples,by = breaks.interval)
  no.breaks=length(breaks)
  for(m in 1:no.breaks){
    start=breaks[m]
    end=start+(breaks.interval-1)
    if(end<no.samples){
      temp.meta.df=meta.df[c(start:end),]
      temp.sc.pfa.hg.mapped.reads.sum.df=chrom.mapping.df[rownames(temp.meta.df),]
      barplot2(height = t(log10(as.matrix(temp.sc.pfa.hg.mapped.reads.sum.df))),
               las=2,col=c('green','red'),border = F,xaxt='none', beside = T,main=title.str)
    }
    else{
      temp.meta.df=meta.df[c(start:no.samples),]
      temp.sc.pfa.hg.mapped.reads.sum.df=chrom.mapping.df[rownames(temp.meta.df),]
      barplot2(height = t(log10(as.matrix(temp.sc.pfa.hg.mapped.reads.sum.df))),
               las=2,col=c('green','red'),border = F,xaxt='none', beside = T,main=title.str)
    }
  }
}
plot.samples.mapping.proportions.barplot.for.subset.samples=function(df,cut.off,title.str,breaks.interval=100){
  no.samples=dim(df)[1]
  breaks=seq(1,to =no.samples,by = breaks.interval)
  no.breaks=length(breaks)
  for(m in 1:no.breaks){
    start=breaks[m]
    end=start+(breaks.interval-1)
    if(end<no.samples){
      temp.meta.df=df[c(start:end),]
      temp.out.list=plot.samples.mapping.proportions.barplot(df =temp.meta.df,cut.off = cut.off,title.str =  title.str)
    }
    else{
      temp.meta.df=df[c(start:no.samples),]
      temp.out.list=plot.samples.mapping.proportions.barplot(df =temp.meta.df,cut.off = cut.off,title.str =  title.str)
    }
  }
  plot.new()
  #legend('center',legend=c('Unique','Multimap','Unmapped'),fill=c('red','green','blue'),cex=1.5,box.lty = 0)
  legend('center',legend=c('','',''),fill=c('red','green','blue'),cex=2.0,box.lty = 0,border=NA)
  unique.reads.vec=as.numeric(as.character(df$Uniquely_mapped_reads_number_))
  total.mapped.vec=as.numeric(as.character(df$Uniquely_mapped_reads_number_))+as.numeric(as.character(df$Number_of_reads_mapped_to_multiple_loci_))+as.numeric(as.character(df$Number_of_reads_mapped_to_too_many_loci_))
  unique.intervals = seq(0,max(unique.reads.vec+10000),10000)
  mapped.intervals = seq(0,max(total.mapped.vec+10000),10000)
  #par(mfrow=c(1,1))
  #hist(x = total.mapped.vec,main = 'Total mapped reads distribution',xlab='Read counts',breaks = mapped.intervals)
  #hist(x = unique.reads.vec,main = 'Unique reads distribution',breaks = unique.intervals,xlab='Read counts')
}


plot.samples.reads.hist=function(df,cut.off,title.str,breaks.interval=100){
  no.samples=dim(df)[1]
  breaks=seq(1,to =no.samples,by = breaks.interval)
  no.breaks=length(breaks)
  for(m in 1:no.breaks){
    start=breaks[m]
    end=start+(breaks.interval-1)
    if(end<no.samples){
      temp.meta.df=df[c(start:end),]
      temp.out.list=plot.samples.mapping.proportions.barplot(df =temp.meta.df,cut.off = cut.off,title.str =  title.str)
    }
    else{
      temp.meta.df=df[c(start:no.samples),]
      temp.out.list=plot.samples.mapping.proportions.barplot(df =temp.meta.df,cut.off = cut.off,title.str =  title.str)
    }
  }
  unique.reads.vec=as.numeric(as.character(df$Uniquely_mapped_reads_number_))
  total.mapped.vec=as.numeric(as.character(df$Uniquely_mapped_reads_number_))+as.numeric(as.character(df$Number_of_reads_mapped_to_multiple_loci_))+as.numeric(as.character(df$Number_of_reads_mapped_to_too_many_loci_))
  total.mapped.vec =total.mapped.vec[total.mapped.vec<=2000000]
  log.total.mapped.vec=log10(total.mapped.vec)
  unique.reads.vec=unique.reads.vec[unique.reads.vec<2000000]
  log.unique.reads.vec=log10(unique.reads.vec)
  unique.intervals = seq(0,max(unique.reads.vec+10000),10000)
  mapped.intervals = seq(0,max(total.mapped.vec+10000),10000)
  hist(x = total.mapped.vec,main = 'Total mapped reads distribution',xlab='Read counts',breaks = mapped.intervals)
  hist(x = log.total.mapped.vec,main = '',xlab='',ylab='')
  #hist(x = unique.reads.vec,main = 'Unique reads distribution',breaks = unique.intervals,xlab='Read counts')
  hist(x = log.unique.reads.vec,main = '',xlab='',ylab='')
}


plot.samples.reads.hist=function(counts.df){
  no.detected.genes=as.numeric(apply(counts.df,2,nnzero))
  hist(x = log10(no.detected.genes),main = 'Detected genes',xlab='No. of genes counts',ylab='No. of samples')
  hist(x = log10(no.detected.genes),main = '',xlab='',ylab='')
  #hist(x = log.total.mapped.vec,main = '',xlab='',ylab='')
  #hist(x = unique.reads.vec,main = 'Unique reads distribution',breaks = unique.intervals,xlab='Read counts')
  #hist(x = log.unique.reads.vec,main = '',xlab='',ylab='')
}


plot.pca=function(in.rpkm.df,trans=F,meta.data,title.str='Test pca plot',var.cut.off=0,first.pc='PC1',
                  second.pc='PC2',log.pca.scores=F,plot.pairs=T,three.dim.scatter=F,log.rpkm=F){
  rpkm.df=filter.none.expressed.samples(in.rpkm.df)
  pseudo.value=min(rpkm.df[rpkm.df!=0])/2
  samples=as.character(colnames(rpkm.df))
  meta.data=meta.data[samples,]
  meta.data=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.data)
  meta.data$development.stage=factor(meta.data$development.stage,
                                     levels = factor(c('R','LR','ET','T','ES','S')))
  meta.data$no.of.detected.genes=as.numeric(as.character(apply(rpkm.df,2,nnzero)))
  rpkm.df=if(trans){t(rpkm.df)} else{rpkm.df}
  #var.cut.off=as.numeric(quantile(as.numeric(apply(rpkm.df,1,sd)),.3))
  #rpkm.df=filter.non.variable.rows(rpkm.df,cut.off=var.cut.off)
  rpkm.df=filter.genes.with.zero.variance(df = rpkm.df)
  no.samples=dim(rpkm.df)[2]
  no.genes=dim(rpkm.df)[1]
  out.results.list=list()
  if(log.rpkm){
    log.rpkm.df=rpkm.df
    #log.rpkm.df=log2(log.rpkm.df+pseudo.value)
    log.rpkm.df=log2(log.rpkm.df+1)
    rpkm.df=log.rpkm.df
  }
  samples.pca=prcomp(t(rpkm.df),retx=T,center=T,scale.=T)
  gene.rotations.df=subset(samples.pca$rotation,select = c(PC1,PC2))
  #samples.pca=prcomp(t(log.rpkm.df),retx=T,center=T)
  samples.pca.summary=summary(samples.pca)$importance
  out.results.list[['pca.summary']]=samples.pca.summary
  ylim.vec=c(0,max(100*samples.pca.summary[2,])+10)
  barplot(100*samples.pca.summary[2,],main='Components explained variance',ylab='%',las=2,cex.names = .5,ylim = ylim.vec)
  abline(h=100/length(colnames(samples.pca.summary)),b=1)
  first.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,first.pc]))*100,2)
  second.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,second.pc]))*100,2)
  out.results.list[['pca.obj']]= samples.pca
  samples.pca.scores=samples.pca$x
  samples=rownames(samples.pca.scores)
  pc.scores.meta.df=data.frame(samples.pca.scores,meta.data)
  x.points=as.numeric(as.character(pc.scores.meta.df[,first.pc]))
  y.points=as.numeric(as.character(pc.scores.meta.df[,second.pc]))
  if(log.pca.scores){
    pseudo.x=min(x.points[x.points!=0])/2
    pseudo.y=min(y.points[y.points!=0])/2
    x.points=log2(x.points+pseudo.x)
    y.points=log2(y.points+pseudo.y)
  }
  #abbrev.std.stages.col=c('R'='green','LR'='green4','ET'='blue','LT'='blue4','ES'='deeppink3','S'='red')
  abbrev.std.stages.col=c('R'='cyan','LR'='green4','ET'='blue4','T'='orange','ES'='purple','S'='firebrick')
  pca.col.str=get.abbrev.std.stage.cols(in.stages.vec = as.character(pc.scores.meta.df$development.stage))
  # pca.col.vec=get.color.list.for.pheatmap(in.vec  = 
  #                                           as.character(pc.scores.meta.df$markers.cluster.groups))
  # 
  #pca.col.str=pca.col.vec[as.character(pc.scores.meta.df$markers.cluster.groups)]
  #mod.col.str=ifelse(as.character(pc.scores.meta.df$mRFP1.expr)=='Y','deeppink4','darkcyan')
  #pca.col.str=mod.col.str
  #biplot(x = samples.pca,col=pca.col.str,main=title.str,cex=.2,arrow.len = 0)
  first.pc.explained.var=paste('(Explained var: ',paste(first.pc.explained.var,'%)',sep=''),sep='')
  second.pc.explained.var=paste('(Explained var: ',paste(second.pc.explained.var,'%)',sep=''),sep='')
  first.pc.lab=paste(first.pc,first.pc.explained.var,'')
  second.pc.lab=paste(second.pc,second.pc.explained.var,'')
  x.lim=c(min(x.points,y.points),max(x.points,y.points))
  y.lim=c(min(x.points,y.points),max(x.points,y.points))
  temp.pc.scores.meta.df=pc.scores.meta.df
  points.size=as.numeric(lapply(strsplit(x=
                                           as.character(temp.pc.scores.meta.df[,'no.of.detected.genes']),
                                         split=':'),function(counts.category){
    counts.category=as.character(counts.category)
    return(counts.category[1])
  }))
  points.size= as.numeric(as.character(temp.pc.scores.meta.df$no.of.detected.genes))
  point.size.fact=points.size/sum(points.size)
  point.size.fact=point.size.fact*(2.5/max(point.size.fact))
  #pch.fact.vec=ifelse(as.character(temp.pc.scores.meta.df$protocol)=='tn5',19,17)
  #protocol.name.vec=as.character(temp.pc.scores.meta.df$protocol)
  plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,
       ylab=second.pc.lab,main=title.str,pch=19,cex=1.5,frame.plot=F)
  #text(x=x.points,y=y.points,labels =  protocol.name.vec,cex = .5)
  abline(h=0,b=1,col = "lightgray")
  abline(v=0,b=1,col = "lightgray")
  plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,
       ylab=second.pc.lab,main=title.str,pch=19,cex=3.0,frame.plot=F)
  plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,ylab=second.pc.lab,
       main=title.str,pch=19,cex=3.0,frame.plot=F)
  abline(h=0,b=1,col = "lightgray")
  abline(v=0,b=1,col = "lightgray")
  plot.new()
  #legend('center',legend=rep('',2),fill=unique(pca.col.str),cex=1.5,box.lty = 0,border = NA)
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,cex=1.0,box.lty = 0,border = NA)
  #x.points=as.numeric(as.character(pc.scores.meta.df[,first.pc]))
  #y.points=as.numeric(as.character(pc.scores.meta.df[,second.pc]))
  temp.pc.scores.meta.df$genes.detected.prop=point.size.fact
  #temp.pc.scores.meta.df$Protocol=ifelse(grepl(pattern = 'inhouse',x = rownames(temp.pc.scores.meta.df)),'Tn5','Nextera')
  #pca.ggplot=ggplot(data = temp.pc.scores.meta.df,aes(x = PC1,y=PC2,colour=development.stage))+geom_point(aes(size=genes.detected.prop))+ggtitle(label = title.str)+xlab(first.pc.lab)+ylab(second.pc.lab)+scale_colour_manual(name='Development stage',values = abbrev.std.stages.col)
  #pca.ggplot=ggplot(data = temp.pc.scores.meta.df,aes(x = PC1,y=PC2,colour=development.stage,shape=protocol))+geom_point(aes(size=1.5))+ggtitle(label = title.str)+xlab(first.pc.lab)+ylab(second.pc.lab)+scale_colour_manual(name='Development stage',values = abbrev.std.stages.col)+scale_size_manual(name='Proportion',values = genes.detected.prop)
  #pca.ggplot+xlab(first.pc.lab)
  #pca.ggplot+ylab(second.pc.lab)
  #print(pca.ggplot)
  #plot.new()
  #legend('center',legend=rep('',length(pca.col.factor[['legend.str']])),fill=pca.col.factor[['legend.cols']],cex = 1.5,box.lty = 0,border = NA)
  #legend('top',legend=pca.col.factor$legend.str,fill=pca.col.factor$legend.cols,cex = 1.5,box.lty = 0,border = NA)
  total.no.pc=dim(samples.pca.scores)[2]
  pca.pair.mat=samples.pca.scores
  if(total.no.pc>=6){
    pca.pair.mat=samples.pca.scores[,1:5]
  }
  if(plot.pairs){
    pairs(pca.pair.mat,cex=1,pch=19,col=pca.col.str,main=title.str)
    #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill)
  }
  if(three.dim.scatter){
    scatterplot3d(x=temp.pc.scores.meta.df[,'PC1'],y=temp.pc.scores.meta.df[,'PC2'],z=temp.pc.scores.meta.df[,'PC3'],xlab='PC1',ylab='PC2',zlab='PC3',pch=19,color=pca.col.str)
  }
  return(out.results.list)
}


plot.pca.biplot=function(in.rpkm.df,trans=F,meta.data,title.str='Test pca plot',gene.annotation.df,var.cut.off=0,first.pc='PC1',second.pc='PC2',
                         log.pca.scores=F,plot.pairs=T,log.rpkm=F,marker='PF3D7_1317200'){
  rpkm.df=filter.none.expressed.samples(in.rpkm.df)
  #rpkm.df=in.rpkm.df
  pseudo.value=min(rpkm.df[rpkm.df!=0])/2
  samples=as.character(colnames(rpkm.df))
  meta.data=meta.data[samples,]
  meta.data=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.data)
  #meta.data$development.stage=factor(meta.data$development.stage,levels = factor(c('R','LR','ET','T','ES','S')))
  meta.data$no.of.detected.genes=as.numeric(as.character(apply(rpkm.df,2,nnzero)))
  rpkm.df=if(trans){t(rpkm.df)} else{rpkm.df}
  #var.cut.off=as.numeric(quantile(as.numeric(apply(rpkm.df,1,sd)),.5))
  rpkm.df=filter.non.variable.rows(rpkm.df,cut.off=var.cut.off)
  no.samples=dim(rpkm.df)[2]
  no.genes=dim(rpkm.df)[1]
  out.results.list=list()
  if(log.rpkm){
    log.rpkm.df=rpkm.df
    log.rpkm.df=log2(log.rpkm.df+pseudo.value)
    rpkm.df=log.rpkm.df
  }
  samples.pca=prcomp(t(rpkm.df),retx=T,center=T,scale.=T)
  gene.rotations.df=data.frame(subset(samples.pca$rotation,select = c(PC1,PC2,PC3,PC4)))
  #gene.rotations.df=subset(gene.rotations.df,PC2>0&PC3>0)
  sex.markers.vec=rownames(gene.annotation.df)
  all.genes.col=ifelse(rownames(gene.rotations.df) %in%  sex.markers.vec, 'red','blue')
  all.pcs.vec=colnames(gene.rotations.df)
  all.pcs.vec.len=length(all.pcs.vec)
  for(m in 1:all.pcs.vec.len){
    temp.pc=all.pcs.vec[m]
    temp.loading.vec=as.numeric(gene.rotations.df[,temp.pc])
    temp.df=data.frame(loading=temp.loading.vec,names = rownames(gene.rotations.df))
    top.genes.df=subset(temp.df,loading>0)
    filt.top.genes.df=top.genes.df[with(top.genes.df,order(loading,decreasing = T)),]
    filt.top.genes.df=head(filt.top.genes.df,200)
    top.col.str=ifelse(filt.top.genes.df$names %in% sex.markers.vec,'red', 'lightgray')
    barplot(height = as.numeric(filt.top.genes.df$loading),ylab='Loading',main = temp.pc,col=top.col.str,border=NA,names.arg = as.character(filt.top.genes.df$names),cex.names = .2,las=2)
    down.genes.df=subset(temp.df,loading<0)
    filt.down.genes.df=down.genes.df[with(down.genes.df,order(loading)),]
    filt.down.genes.df=head(filt.down.genes.df,200)
    #barplot(height = as.numeric(filt.down.genes.df$loading),ylab='Loading',main = temp.pc,col=top.col.str,border=NA,names.arg = as.character(filt.down.genes.df$names),cex.names = .1,las=2)
  }
  #plot(gene.rotations.df,pch=19,cex=1.5,xlab='PC2',ylab='PC3',col=all.genes.col)
  #filt.gene.rotations.df=subset(gene.rotations.df,PC3>0)
  filt.gene.rotations.df=gene.rotations.df
  filt.gene.rotations.df=filt.gene.rotations.df[with(filt.gene.rotations.df,order(PC3,decreasing = T)),]
  #filt.gene.rotations.df=head(filt.gene.rotations.df,n = 170)
  gene.names=rownames(filt.gene.rotations.df)
  gene.load.col.str=ifelse(gene.names %in%  sex.markers.vec, 'red','lightblue')
  barplot(height  =  filt.gene.rotations.df$PC3,col = gene.load.col.str ,names.arg = gene.names,cex.names = .2,pch=19,border = NA,las=2)
  #samples.pca=prcomp(t(log.rpkm.df),retx=T,center=T)
  samples.pca.summary=summary(samples.pca)$importance
  #barplot(100*samples.pca.summary[2,],main='Components explained variance',ylab='%')
  #abline(h=100/length(colnames(samples.pca.summary)),b=1)
  first.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,first.pc]))*100,2)
  second.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,second.pc]))*100,2)
  out.results.list[['pca.obj']]= samples.pca
  samples.pca.scores=samples.pca$x
  samples=rownames(samples.pca.scores)
  pc.scores.meta.df=data.frame(samples.pca.scores,meta.data)
  x.points=as.numeric(as.character(pc.scores.meta.df[,first.pc]))
  y.points=as.numeric(as.character(pc.scores.meta.df[,second.pc]))
  if(log.pca.scores){
    pseudo.x=min(x.points[x.points!=0])/2
    pseudo.y=min(y.points[y.points!=0])/2
    x.points=log2(x.points+pseudo.x)
    y.points=log2(y.points+pseudo.y)
  }
  #abbrev.std.stages.col=c('R'='green','LR'='green4','ET'='blue','LT'='blue4','ES'='deeppink3','S'='red')
  abbrev.std.stages.col=c('R'='cyan','LR'='green4','ET'='blue4','T'='orange','ES'='purple','S'='firebrick')
  #pca.col.str=get.abbrev.std.stage.cols(in.stages.vec = as.character(pc.scores.meta.df$development.stage))
  pca.col.str=get.color.list.for.pheatmap(in.vec  = as.character(pc.scores.meta.df$markers.cluster.groups))
  mod.col.str=ifelse(as.character(pc.scores.meta.df$mRFP1.expr)=='Y','red','darkgreen')
  #mod.col.str=ifelse(as.character(pc.scores.meta.df$mRFP1.expr)=='Y','red',pca.col.str)
  #pca.col.str=mod.col.str
  pca.col.str=pca.col.str[as.character(pc.scores.meta.df$markers.cluster.groups)]
  sc.point.shape=ifelse(as.character(pc.scores.meta.df$mRFP1.expr)=='Y',17,19)
  #PCbiplot(PC = samples.pca,meta.df=pc.scores.meta.df,col.str=pca.col.str,gene.annotation.df = gene.annotation.df,x = first.pc,y = second.pc)
  first.pc.explained.var=paste('(Explained var: ',paste(first.pc.explained.var,'%)',sep=''),sep='')
  second.pc.explained.var=paste('(Explained var: ',paste(second.pc.explained.var,'%)',sep=''),sep='')
  first.pc.lab=paste(first.pc,first.pc.explained.var,'')
  second.pc.lab=paste(second.pc,second.pc.explained.var,'')
  x.lim=c(min(x.points,y.points),max(x.points,y.points))
  y.lim=c(min(x.points,y.points),max(x.points,y.points))
  temp.pc.scores.meta.df=pc.scores.meta.df
  plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,
       ylab=second.pc.lab,main=title.str,pch=sc.point.shape,cex=2.5)
  marker.expr.vec=as.vector(as.numeric(rpkm.df[marker,samples]))
  marker_gene_rpkm_=marker.expr.vec
  show(marker_gene_rpkm_)
  col.gradient.vec=colorRampPalette(c('green','red'))
  abline(h=0,b=1,col = "lightgray")
  abline(v=0,b=1,col = "lightgray")
  plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,
       ylab=second.pc.lab,main=title.str,pch=sc.point.shape,cex=2.5)
  sample.ggplot<-ggplot(data=pc.scores.meta.df,
                        aes(x=PC2,y=PC3,shape=mRFP1.expr))+geom_point(mapping=aes(colour=mRFP1.rpkm,size=3.5))+scale_colour_gradient2(
                          low = "royalblue4", mid = "white",
                                       high = "firebrick3", 
                          midpoint = ((max(marker_gene_rpkm_)+min(marker_gene_rpkm_))/2),
                                       space = "Lab", na.value = "midnightblue", 
                          guide = "colourbar",
                                       limits=c(min(marker_gene_rpkm_), 
                                                max(marker_gene_rpkm_)))
  print(sample.ggplot)
  plot.new()
  #legend('center',legend=rep('',2),fill=unique(pca.col.str),cex=1.5,box.lty = 0,border = NA)
  legend.str=sub.populations.grp.order.vec[sort(unique(names(pca.col.str)))]
  legend.vec=sort(sub.populations.grp.order.names.vec[unique(names(pca.col.str))])
  legend.str=paste('SP.',as.character(legend.vec),sep='')
  legend('center',legend=rep('',times = length(legend.str)),fill=pca.col.str[names(legend.vec)],cex=2.0,box.lty = 0,border = NA)
  plot.new()
  legend('center',legend=rep('',times=2),cex=3.0,box.lty = 0,border = NA,pch=c(17,19),col=pca.col.str['SP.2'])
  #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,cex=1.0,box.lty = 0,border = NA)
  #x.points=as.numeric(as.character(pc.scores.meta.df[,first.pc]))
  #y.points=as.numeric(as.character(pc.scores.meta.df[,second.pc]))
  total.no.pc=dim(samples.pca.scores)[2]
  pca.pair.mat=samples.pca.scores
  if(total.no.pc>=6){
    pca.pair.mat=samples.pca.scores[,1:5]
  }
  if(plot.pairs){
    pairs(pca.pair.mat,cex=1,pch=19,col=pca.col.str,main=title.str)
    #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill)
  }
  return(out.results.list)
}


PCbiplot <- function(PC, x="PC1", y="PC2",meta.df,col.str,gene.annotation.df) {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  meta.df=meta.df[rownames(data),]
  data=data.frame(data,meta.df)
  data$col.str=col.str
  plot(x=as.numeric(data[,x]),y=as.numeric(data[,y]),xlab='',ylab='',t='n')
  points(x =as.numeric(data[,x]) ,y = as.numeric(data[,y]),pch=19,col=data$col.str,cex=3.0)
  #abline(h=0,b=1,col = "lightyellow")
  #abline(v=0,b=1,col = "yellow")
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  var.labels=as.character(gene.annotation.df[datapc$varnames,'Symbol'])
  text(x =datapc$v1,y= datapc$v2,labels = var.labels,cex=.5,col='blue')
  arrows(x0=0, y0=0, x1 = datapc$v1, y1 = datapc$v2,length = .04,col = 'lightgray',angle=10,lwd=0.1)
}


#Uses PCA to find most variable gene
get.most.var.genes.pca=function(rpkm.df,meta.df,no.cells.detected.in=2,rpkm.cut.off=1){
  rpkm.df=filter.rpkm.less.than.cutoff(df = rpkm.df,rpkm.cutoff =rpkm.cut.off,no.samples = no.cells.detected.in )
  pca.obj=prcomp(t(rpkm.df))
  pca.scores.df=pca.obj$x
  temp.list=select.gene.with.sig.component.weights(pca.obj =pca.obj )
  genes.vec=unique(as.character(unlist(temp.list[['pc.loadings.list']])))
  rpkm.df=filter.none.expressed.samples(df = filter.none.expressed.genes(input.data = rpkm.df[genes.vec,]))
  meta.df=meta.df[colnames(rpkm.df),]
  plot.pca(in.rpkm.df =rpkm.df,meta.data = meta.df,title.str = 'Test',var.cut.off = 0,log.pca.scores = F)
}


highlight.sexual.markers.expression.in.pca.plot=function(in.rpkm.df,trans=F,meta.data,
                                                         title.str='Test pca plot',
                                                         var.cut.off=0,first.pc='PC1',
                                                         second.pc='PC2',log.pca.scores=F,
                                                         plot.pairs=T,log.rpkm=F,markers.vec=c('PF3D7_1317200'),markers.annotation.df){
  rpkm.df=filter.none.expressed.samples(in.rpkm.df)
  #rpkm.df=in.rpkm.df
  in.rpkm.df=rpkm.df
  pseudo.value=min(rpkm.df[rpkm.df!=0])/2
  samples=as.character(colnames(rpkm.df))
  meta.data=meta.data[samples,]
  meta.data=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.data)
  #meta.data$development.stage=factor(meta.data$development.stage,levels = factor(c('R','LR','ET','T','ES','S')))
  meta.data$no.of.detected.genes=as.numeric(as.character(apply(rpkm.df,2,nnzero)))
  rpkm.df=if(trans){t(rpkm.df)} else{rpkm.df}
  #var.cut.off=as.numeric(quantile(as.numeric(apply(rpkm.df,1,sd)),.5))
  rpkm.df=filter.non.variable.rows(rpkm.df,cut.off=var.cut.off)
  no.samples=dim(rpkm.df)[2]
  no.genes=dim(rpkm.df)[1]
  out.results.list=list()
  if(log.rpkm){
    log.rpkm.df=rpkm.df
    log.rpkm.df=log2(log.rpkm.df+pseudo.value)
    rpkm.df=log.rpkm.df
  }
  samples.pca=prcomp(t(rpkm.df),retx=T,center=T,scale.=T)
  samples.pca.summary=summary(samples.pca)$importance
  first.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,first.pc]))*100,2)
  second.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,second.pc]))*100,2)
  samples.pca.scores=samples.pca$x
  samples=rownames(samples.pca.scores)
  pc.scores.meta.df=data.frame(samples.pca.scores,meta.data)
  x.points=as.numeric(as.character(pc.scores.meta.df[,first.pc]))
  y.points=as.numeric(as.character(pc.scores.meta.df[,second.pc]))
  first.pc.explained.var=paste('(Explained var: ',paste(first.pc.explained.var,'%)',sep=''),sep='')
  second.pc.explained.var=paste('(Explained var: ',paste(second.pc.explained.var,'%)',sep=''),sep='')
  first.pc.lab=paste(first.pc,first.pc.explained.var,'')
  second.pc.lab=paste(second.pc,second.pc.explained.var,'')
  len.markers=length(markers.vec)
  for(m in 1:len.markers){
    marker=markers.vec[m]
    marker.expr.vec=as.vector(as.numeric(rpkm.df[marker,samples]))
    if(all(is.na(marker.expr.vec))){
      next
    }
    marker.name=marker
    #marker.name=markers.annotation.df[marker,'gene.name']
    title.str=paste(marker,'(',marker.name,')')
    marker_gene_rpkm_=marker.expr.vec
    pc.scores.meta.df$log2.rpkm=marker.expr.vec
    sample.ggplot<-ggplot(data=pc.scores.meta.df,
                          aes(x=PC2,y=PC3,shape=mRFP1.expr))+xlab(first.pc.lab)+ylab(second.pc.lab)
    sample.ggplot<-sample.ggplot+
      geom_point(mapping=aes(colour=log2.rpkm,size=6))+scale_colour_gradient2(
      low = "royalblue4", mid = "white",
      high = "firebrick3", 
      midpoint = ((max(marker_gene_rpkm_)+min(marker_gene_rpkm_))/2),
      space = "Lab", na.value = "midnightblue", 
      guide = "colourbar",
      limits=c(min(marker_gene_rpkm_), 
               max(marker_gene_rpkm_)))+ggtitle(title.str)+
      theme_bw()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "gray"),
            plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "cm"))+
      guides(size=FALSE)
    print(sample.ggplot)
  }
}

filter.rpkm.less.than.cutoff=function(df,rpkm.cutoff=0,no.samples=1){
  df=df[rownames(df)[rowSums(df>rpkm.cutoff)>=no.samples],]
  return(df)
}


filter.count.less.than.cutoff=function(df,count.cutoff=0,no.samples=1){
  df[df<count.cutoff]=0
  df$sample.detected.in.counts=as.numeric(as.character(apply(df,1,nnzero)))
  df=df[which(df$sample.detected.in.counts>=no.samples),]
  out.df=subset(df,select = -sample.detected.in.counts)
  return(out.df)
}


#Filters data.frame based on rpkm cut off and sample no.
reset.rpkm.less.than.cut.off=function(df,rpkm.cutoff){
  df[df<rpkm.cutoff]=0
  df=filter.none.expressed.samples(filter.none.expressed.genes(input.data))
  return(df)
}


filter.non.variable.rows=function(df,cut.off){
  df=data.frame(df)
  df[,'sd']=as.numeric(as.character(apply(df,1,sd)))
  out.df=df[which(df[,'sd']>=cut.off),]
  out.df=data.frame(subset(out.df,select=-sd))
  return(out.df)
}


plot.samples.correlations.heatmap = function(in.rpkm.df,title.str='Test plot',in.meta.df,
                                             filter.non.var.rows=F,in.method='spearman',
                                             tree.cut.off=3,annotation_legend = T,show_rownames = T,
                                             show_colnames = T,cellwidth=1,cellheight = 1,
                                             in.col.ramp=colorRampPalette(c("gray",'red'))(1000)){
  rpkm.df=in.rpkm.df
  meta.df=in.meta.df
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df$development.stage=factor(meta.df$development.stage,levels = factor(c('R','LR','ET','T','ES','S')))
  #temp.raw.rpkm.df=filter.none.expressed.genes(rpkm.df)
  temp.raw.rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df=rpkm.df,sample.count = 2)
  temp.rpkm.mat=as.matrix(temp.raw.rpkm.df)
  #pseudo.value=min(as.numeric(temp.rpkm.mat[temp.rpkm.mat!=0]))/2
  pseudo.value=1
  if(filter.non.var.rows){
    cut.off=as.numeric(quantile(as.numeric(apply(temp.rpkm.mat,1,sd)),.3))
    cut.off=1
    temp.rpkm.mat=as.matrix(filter.none.expressed.samples(filter.non.variable.rows(df=temp.rpkm.mat,cut.off=cut.off)))
    temp.raw.rpkm.df=filter.none.expressed.samples(filter.non.variable.rows(df = temp.raw.rpkm.df,cut.off = cut.off))
  }
  temp.rpkm.mat=cor(log2(temp.rpkm.mat+pseudo.value),method = in.method)
  #temp.rpkm.mat=cor(temp.rpkm.mat,method = in.method)
  meta.df=meta.df[colnames(temp.rpkm.mat),]
  #col.factor=as.character(meta.df[,'timepoint'])
  col.factor=as.character(meta.df[,'development.stage'])
  #col.factor=as.character(meta.df[,'parent.development.stage'])
  #batch.col.factor=as.character(meta.df[,'matching.bulks'])
  batch.col.factor=as.character(meta.df[,'batch'])
  #batch.col.factor=as.character(meta.df[,'top.sub.grp'])
  batch.col.list=get.col.factor(col.factor = batch.col.factor)
  #col.factor.list=get.col.factor(col.factor)
  col.factor.str=get.abbrev.std.stage.cols(in.stages.vec = col.factor)
  #col.factor.str=col.factor.list$col.str
  #col.factor.legend.str=col.factor.list$legend.str
  #col.factor.legend.col=col.factor.list$legend.cols
  #heatmap.2(x =temp.rpkm.mat,trace = 'none',cexRow = .4,cexCol = .4,ColSideColors = col.factor.str, main =title.str,dendrogram='col')
  #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,cex=.8)
  col.annot.col.df=subset(meta.df,select = c(development.stage))
  corr.dist.hclust=hclustfunc(dist(x =  temp.rpkm.mat,method = 'euclidean'))
  cut.tree.corr.dist.hclust=cutree(tree = corr.dist.hclust,k = tree.cut.off)
  #tree.annotation.df=data.frame(group=paste('SP:',as.character(as.numeric(cut.tree.corr.dist.hclust)),sep = ''))
  tree.annotation.df=data.frame(group=paste('SP:',
                                            as.character(as.numeric(cut.tree.corr.dist.hclust)),sep = ''))
  group.col.vec=get.color.rainbow.list.for.pheatmap(in.vec = as.character(tree.annotation.df$group))
  col.annot.col.list=list(group=group.col.vec,development.stage=abbrev.std.stages.col.vec)
  rownames(tree.annotation.df) = names(cut.tree.corr.dist.hclust)
  tree.annotation.df=data.frame(t(subset(t(tree.annotation.df),select = rownames(col.annot.col.df))))
  #col.annot.col.df=col.annot.col.df[rownames(tree.annotation.df),]
  #tree.annotation.df=tree.annotation.df[rownames(col.annot.col.df),]
  col.annot.col.df=data.frame(col.annot.col.df,tree.annotation.df)
  #col.annot.col.df=data.frame(col.annot.col.df)
  #pheatmap(mat = temp.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = 1,cellheight = .5,treeheight_row = 1,fontsize_row = 1,fontsize_col = 1,cutree_cols = tree.cut.off,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,annotation_legend = annotation_legend,show_rownames = show_rownames,show_colnames = show_colnames,border_color = NA)
  col.pelette=in.col.ramp
  pheatmap(mat = temp.rpkm.mat,main =title.str,
           color = col.pelette,
           cellwidth = cellwidth,cellheight = cellheight,fontsize_row = 1,
           fontsize_col = 1,cutree_cols = tree.cut.off,annotation_col = col.annot.col.df,
           annotation_colors = col.annot.col.list,legend = T,
           annotation_legend = annotation_legend,show_rownames = show_rownames,
           show_colnames = show_colnames,border_color = NA)
  pheatmap(mat = temp.rpkm.mat,main =title.str,color = col.pelette,
           cellwidth = cellwidth,cellheight = cellheight,fontsize_row = 1,
           fontsize_col = 1,annotation_col = col.annot.col.df,
           annotation_colors = col.annot.col.list,legend = F,
           annotation_legend = annotation_legend,show_rownames = show_rownames,
           show_colnames = show_colnames,border_color = NA)
  # pheatmap(mat = temp.rpkm.mat,main =title.str,color = col.pelette,
  #          cellwidth = cellwidth,cellheight = cellheight,fontsize_row = 1,
  #          fontsize_col = 1,annotation_col = col.annot.col.df,
  #          annotation_colors = col.annot.col.list,legend = F,
  #          annotation_legend = annotation_legend,show_rownames = show_rownames,
  #          show_colnames = show_colnames,border_color = NA)
  # 
  #heatmap.2(x =temp.rpkm.mat,trace = 'none',cexRow = .4,cexCol = .4,ColSideColors = col.factor.str, main =title.str,dendrogram='col',labCol = '',margins = c(10,10),col=bluered(10000),key.xlab = '',key.ylab = '',key.title = 'Spearman',RowSideColors = batch.col.list$col.str )
  #heatmap.2(x =temp.rpkm.mat,trace = 'none',cexRow = .4,cexCol = .4,ColSideColors = col.factor.str, main ='',dendrogram='col',labCol = '',margins = c(10,10),col=bluered(10000),key.xlab = '',key.ylab = '',key.title = 'Spearman' )
  #plot.new()
  #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,cex=.8)
  #plot.new()
  #legend('topright',legend=rep('',times = length(std.stages.col.abbrev.legend.str)),fill=std.stages.col.legend.fill,cex=1.5,box.lty = 0)
  #legend('left', legend =batch.col.list$legend.str, fill=batch.col.list$legend.col, cex=.5,box.lty = 0)
  #heatmap.2(x =as.matrix(1-temp.rpkm.mat),trace = 'none',cexRow = .4,cexCol = .4,ColSideColors = col.factor.str, main =title.str,dendrogram='col',labRow = '',margins = c(10,10) ,col=bluered(1000))
  #legend('topright', legend =col.factor.legend.str, fill=col.factor.legend.col, cex=1)
  #legend('left', legend =batch.col.list$legend.str, fill=batch.col.list$legend.col, cex=1)
  #plot.new()
  #legend('center', legend =rep('',length(col.factor.legend.str)), fill=col.factor.legend.col, cex=1.5,box.lty = 0)
  #heatmap.2(x =as.matrix(1-temp.rpkm.mat),trace = 'none',cexRow = .4,cexCol = .4,ColSideColors = col.factor.str, main ='',dendrogram='col',labCol='',labRow = '',margins = c(10,10) ,key.xlab = '',key.ylab = '',key.title = '',col=bluered(1000))
  #plot.new()
  #legend('center', legend =rep('',times=length(col.factor.legend.str)), fill=col.factor.legend.col, cex=1.5,box.lty = 0)
  logged.temp.raw.rpkm.mat=as.matrix(x = log10(temp.raw.rpkm.df+pseudo.value))
  #logged.temp.raw.rpkm.mat=as.matrix(x = temp.raw.rpkm.df)
  #logged.temp.raw.rpkm.mat=as.matrix(x = temp.raw.rpkm.df)
  #cor.dendrogram=as.dendrogram(hclust(d = as.dist(1-cor(logged.temp.raw.rpkm.mat,method = in.method))))
  #heatmap.2(x =logged.temp.raw.rpkm.mat,trace = 'none',cexRow = .4,cexCol = .4,ColSideColors = col.factor.str, main =title.str,labRow = '',labCol = '',margins = c(10,10),key.xlab = 'Z-score',key.ylab = '',key.title = '',col=bluered(10000),scale = 'row')
  #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,cex=.8)
  #pheatmap(mat =logged.temp.raw.rpkm.mat,main =title.str,color = bluered(10000),cellwidth = 1,cellheight = 1,treeheight_row = 1,fontsize_row = 2,fontsize_col = 2,cutree_cols = 4,annotation_col = col.annot.col.df)
  #pheatmap(mat =logged.temp.raw.rpkm.mat,main ='',color = bluered(10000),cellwidth = 1,cellheight = 1,treeheight_row = 1,fontsize_row = 1,fontsize_col = 1,cutree_cols = 4,annotation_col = col.annot.col.df)
  sup.pop.vec=as.character(data.frame(t(subset(t(tree.annotation.df), select = rownames(meta.df))))$group)
  meta.df$sub.group=sup.pop.vec
  out.list=list(rpkm=temp.raw.rpkm.df,meta=meta.df,cor.mat=temp.rpkm.mat)
  return(out.list)
  }


plot.samples.idc.expression.heatmap = function(in.rpkm.df,title.str='Test plot',in.meta.df,marker.df){
  col.pelette=colorRampPalette(c("darkgray","gray","white",
                                 "yellow","orange",'red'))(5)
  rpkm.df=in.rpkm.df
  meta.df=in.meta.df
  markers.cluster.groups.vec=as.character(meta.df$markers.cluster.groups)
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df$development.stage=factor(meta.df$development.stage,
                                   levels = factor(c('R','LR','ET','T','ES','S')))
  markers.cluster.groups.list=get.color.list.for.pheatmap(in.vec = markers.cluster.groups.vec)
  temp.raw.rpkm.df=filter.none.expressed.genes(rpkm.df)
  temp.rpkm.mat=as.matrix(temp.raw.rpkm.df)
  #pseudo.value=min(as.numeric(temp.rpkm.mat[temp.rpkm.mat!=0]))/2
  #temp.rpkm.mat=cor(log10(temp.rpkm.mat+pseudo.value),method = 'spearman')
  #temp.rpkm.mat=cor(temp.rpkm.mat,method = in.method)
  meta.df=meta.df[colnames(temp.rpkm.mat),]
  col.factor=as.character(meta.df[,'development.stage'])
  batch.col.factor=as.character(meta.df[,'batch'])
  #batch.col.factor=as.character(meta.df[,'top.sub.grp'])
  batch.col.list=get.col.factor(col.factor = batch.col.factor)
  #col.factor.list=get.col.factor(col.factor)
  col.factor.str=get.abbrev.std.stage.cols(in.stages.vec = col.factor)
  col.annot.col.df=subset(meta.df,select = development.stage)
  col.annot.col.df=subset(meta.df,select =markers.cluster.groups)
  #temp.annotation.cols=
  temp.annotation.cols=unique(as.character(col.annot.col.df[,1]))
  len.temp.annotation.cols=length(temp.annotation.cols)
  col.annot.col.list=list()
  for(m in 1:len.temp.annotation.cols){
    m
    #col.annot.col.list=append(col.annot.col.list,temp.annotation.cols[m]=markers.cluster.groups.list[temp.annotation.cols[m]],after = length(col.annot.col.list)) 
  }
  #col.annot.col.list=list(development.stage=abbrev.std.stages.col.vec)
  row.col.df=marker.df[rownames(temp.rpkm.mat),]
  row.col.df=subset(row.col.df,select=category)
  pheatmap(mat = log2(temp.rpkm.mat+1),main =title.str,
           color = col.pelette,
           cellwidth = 5,cellheight = 5,clustering_distance_cols = 'euclidean',
           clustering_distance_rows = 'correlation',treeheight_row = 1,
           fontsize_row = 2,fontsize_col = 2,cutree_cols = 2,
           annotation_col = col.annot.col.df, legend = T,annotation_legend = F,
           show_rownames = F,show_colnames = F,annotation_row = row.col.df)
  pheatmap(mat = log2(temp.rpkm.mat+1),main =title.str,
           color = col.pelette,
           cellwidth = 6,cellheight = 6,clustering_distance_cols = 'euclidean',
           clustering_distance_rows = 'correlation',treeheight_row = 1,
           fontsize_row = 2,fontsize_col = 2,cutree_cols = 2,
           annotation_col = col.annot.col.df,legend = T,annotation_legend = T,
           show_rownames = F,show_colnames = F,annotation_row = row.col.df)
  #heatmap.2(x =temp.rpkm.mat,trace = 'none',cexRow = .4,cexCol = .4,ColSideColors = col.factor.str, main =title.str,dendrogram='col',labCol = '',margins = c(10,10),col=bluered(10000),key.xlab = '',key.ylab = '',key.title = 'Spearman',RowSideColors = batch.col.list$col.str )
  #heatmap.2(x =temp.rpkm.mat,trace = 'none',cexRow = .4,cexCol = .4,ColSideColors = col.factor.str, main ='',dendrogram='col',labCol = '',margins = c(10,10),col=bluered(10000),key.xlab = '',key.ylab = '',key.title = 'Spearman' )
  plot.new()
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,cex=.8)
  plot.new()
  legend('topright',legend=rep('',times = length(std.stages.col.abbrev.legend.str)),fill=std.stages.col.legend.fill,cex=1.5,box.lty = 0)
  logged.temp.raw.rpkm.mat=as.matrix(x = log2(temp.raw.rpkm.df+1))
  pheatmap(mat =logged.temp.raw.rpkm.mat,main =title.str,
           color = col.pelette,
           cellwidth = 5,cellheight = 5,treeheight_row = 1,fontsize_row = 2,
           fontsize_col = 2,cutree_cols = 4,annotation_col = col.annot.col.df)
  pheatmap(mat =logged.temp.raw.rpkm.mat,main =title.str,
           color = col.pelette,
           cellwidth = 5,cellheight = 5,treeheight_row = 1,fontsize_row = 2,
           fontsize_col = 2,cutree_cols = 4,annotation_col = col.annot.col.df,
           annotation_legend = F)
  #pheatmap(mat =logged.temp.raw.rpkm.mat,main ='',color = bluered(10000),cellwidth = 1,cellheight = 1,treeheight_row = 1,fontsize_row = 1,fontsize_col = 1,cutree_cols = 4,annotation_col = col.annot.col.df)
  out.list=list(rpkm=temp.raw.rpkm.df,meta=meta.df,cor.mat=temp.rpkm.mat)
  return(out.list)
}


plot.samples.idc.expression.heatmap.per.sub.population = function(in.rpkm.df,title.str='Test plot',
                                                                  in.meta.df,idc.marker.df,
                                                                  all.plos.one.markers.df,
                                                                  in.col.ramp=colorRampPalette(c("gray","white",'red'))(100)){
  genes.col.vec=get.color.rainbow.list.for.pheatmap(as.character(all.plos.one.markers.df$category))
  #genes.col.vec=get.distinct.colour.list.for.pheatmap(in.vec = as.character(all.plos.one.markers.df$category))
  genes.col.vec=gene.cat.col.names
  col.pelette=in.col.ramp
  rpkm.df=in.rpkm.df
  rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  meta.df=in.meta.df
  markers.cluster.groups.vec=as.character(meta.df$markers.cluster.groups)
  expr.markers.vec=intersect(rownames(idc.marker.df),rownames(rpkm.df))
  test.df=all.plos.one.markers.df[expr.markers.vec,]
  expr.markers.vec=rownames(subset(test.df,category!='none'))
  temp.df=subset(test.df,category!='none')
  rpkm.df=rpkm.df[expr.markers.vec,]
  show(dim(rpkm.df))
  meta.df$sub.pop.order=sub.populations.grp.order.names.vec[markers.cluster.groups.vec]
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  markers.cluster.groups.list=get.color.list.for.pheatmap(in.vec = markers.cluster.groups.vec)
  #markers.cluster.groups.list=get.color.rainbow.list.for.pheatmap(in.vec = markers.cluster.groups.vec)
  ordered.sub.pop.col.list=unlist(lapply(names(markers.cluster.groups.list),function(list.item){
    list.item.col=markers.cluster.groups.list[list.item]
    return(list.item.col)
  }))
  names(ordered.sub.pop.col.list)=sub.populations.grp.order.vec[names(ordered.sub.pop.col.list)]
  temp.raw.rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  temp.rpkm.mat=as.matrix(temp.raw.rpkm.df)
  pseudo.value=min(as.numeric(temp.rpkm.mat[temp.rpkm.mat!=0]))/2
  log.rpkm.mat=log2(temp.rpkm.mat+1)
  meta.df=meta.df[colnames(log.rpkm.mat),]
  col.annot.col.df=subset(meta.df,select =c(markers.cluster.groups,sub.pop.order,
                                            development.stage,batch))
  row.annotation.col.df=all.plos.one.markers.df[rownames(log.rpkm.mat),]
  ordered.sub.pop.names=sort(unique(markers.cluster.groups.vec))
  #ordered.sub.pop.meta.df=meta.df[with(meta.df,order(markers.cluster.groups)), ]
  ordered.sub.pop.meta.df=meta.df[with(meta.df,order(sub.pop.order)),]
  ordered.sub.pop.samples.names=rownames(ordered.sub.pop.meta.df)
  row.annotation.col.df=subset(row.annotation.col.df,select = category)
  gene.categories.names=unique(row.annotation.col.df$category)
  sorted.plos.one.markers.df=all.plos.one.markers.df[rownames(log.rpkm.mat),]
  sorted.plos.one.markers.df=sorted.plos.one.markers.df[with(sorted.plos.one.markers.df,
                                                             order(category)),]
  sorted.log.rpkm.mat=log.rpkm.mat[,ordered.sub.pop.samples.names]
  sorted.log.rpkm.mat=sorted.log.rpkm.mat[rownames(sorted.plos.one.markers.df),]
  plasmodb.id=rownames(sorted.log.rpkm.mat)
  gene.names.symbols.str=as.character(idc.marker.df[plasmodb.id,'X.Gene.Name.or.Symbol.'])
  gene.names.symbols.str=ifelse(gene.names.symbols.str=='null',plasmodb.id,gene.names.symbols.str)
  gene.names.categories.str=as.character(all.plos.one.markers.df[plasmodb.id,'categories'])
  #rownames(sorted.log.rpkm.mat)=gene.names.symbols.str
  gene.category.col=genes.col.vec[unique(all.plos.one.markers.df[rownames(log.rpkm.mat),'category'])]
  print(gene.category.col)
  #col.list=list(category=gene.category.col,
                #markers.cluster.groups=markers.cluster.groups.list[unique(col.annot.col.df$markers.cluster.groups)],
                #sub.pop.order=ordered.sub.pop.col.list,development.stage=abbrev.std.stages.col.vec)
  col.list=list(category=gene.category.col,
                markers.cluster.groups=markers.cluster.groups.list[unique(col.annot.col.df$markers.cluster.groups)],
                development.stage=abbrev.std.stages.col.vec)
  final.col.annot.col.df=subset(col.annot.col.df,
                                select=c(markers.cluster.groups,development.stage))
  col.gap.indices=c()
  col.names=colnames(sorted.log.rpkm.mat)
  len.col.names=length(col.names)
  for(m in 1:len.col.names){
    if(m<len.col.names){
      first.sample=col.names[m]
      sec.sample=col.names[m+1]
      first.sample.sp=as.character(col.annot.col.df[first.sample,'markers.cluster.groups'])
      sec.sample.sp=as.character(col.annot.col.df[sec.sample,'markers.cluster.groups'])
      if(first.sample.sp!=sec.sample.sp){
        col.gap.indices=append(col.gap.indices,m,length(col.gap.indices))
      }
    }
  }
  row.indices=c()
  row.names=rownames(sorted.log.rpkm.mat)
  len.row.names=length(row.names)
  for(n in 1:len.row.names){
    if(n<len.row.names){
      first.gene=row.names[n]
      sec.gene=row.names[n+1]
      first.sample.sp=as.character(all.plos.one.markers.df[first.gene,'category'])
      sec.sample.sp=as.character(all.plos.one.markers.df[sec.gene,'category'])
      if(first.sample.sp!=sec.sample.sp){
        row.indices=append(row.indices,n,length(row.indices))
      }
    }
  }
  pheatmap(mat =sorted.log.rpkm.mat,main =title.str,color = col.pelette,cellwidth = 1,
           cellheight = 1,fontsize_row = 2,fontsize_col = 2,
           annotation_col = final.col.annot.col.df,annotation_row = row.annotation.col.df,
           annotation_legend = T,annotation_colors = col.list,cluster_rows = T,cluster_cols = F,
           border_color = NA,labels_row =  gene.names.categories.str,fontsize = 2)
  pheatmap(mat =sorted.log.rpkm.mat,main =title.str,color = col.pelette,cellwidth = 2,
           cellheight = 5,fontsize_row = 2,fontsize_col = 2,annotation_col = final.col.annot.col.df,
           annotation_row = row.annotation.col.df,annotation_legend = F,annotation_colors = col.list,
           cluster_rows = F,cluster_cols = F,border_color = NA,show_rownames = T,show_colnames = F,
           gaps_col = col.gap.indices,gaps_row = row.indices,labels_row = gene.names.symbols.str,fontsize = 2)
  pheatmap(mat =sorted.log.rpkm.mat,main =title.str,color = col.pelette,cellwidth = 3,
           cellheight = 3,fontsize_row = 2,fontsize_col = 2,annotation_col = final.col.annot.col.df,
           annotation_row = row.annotation.col.df,annotation_legend = F,annotation_colors = col.list,
           cluster_rows = F,cluster_cols = F,gaps_row = row.indices,border_color = NA,
           show_rownames = T,show_colnames = F,fontsize = 2)
}


plot.samples.gene.clusters.heatmap = function(in.rpkm.df,title.str='Test plot',in.meta.df,var.genes.cut.off=0,gene.cluster.cut.off=2){
  rpkm.df=in.rpkm.df
  meta.df=in.meta.df
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df$development.stage=factor(meta.df$development.stage,levels = factor(c('R','LR','ET','LT','ES','S')))
  temp.raw.rpkm.df=filter.none.expressed.genes(rpkm.df)
  temp.rpkm.mat=as.matrix(temp.raw.rpkm.df)
  temp.rpkm.mat=as.matrix(filter.none.expressed.samples(filter.non.variable.rows(df=temp.rpkm.mat,cut.off=var.genes.cut.off)))
  temp.raw.rpkm.df=filter.none.expressed.samples(filter.non.variable.rows(df = temp.raw.rpkm.df,cut.off = var.genes.cut.off))
  pseudo.value=min(as.numeric(temp.rpkm.mat[temp.rpkm.mat!=0]))/2
  logged.temp.raw.rpkm.mat=as.matrix(x = log10(temp.raw.rpkm.df+pseudo.value))
  meta.df=meta.df[colnames(temp.raw.rpkm.df),]
  #col.factor=as.character(meta.df[,'timepoint'])
  col.factor=as.character(meta.df[,'development.stage'])
  #col.factor=as.character(meta.df[,'parent.development.stage'])
  #batch.col.factor=as.character(meta.df[,'matching.bulks'])
  batch.col.factor=as.character(meta.df[,'batch'])
  #batch.col.factor=as.character(meta.df[,'top.sub.grp'])
  batch.col.list=get.col.factor(col.factor = batch.col.factor)
  #col.factor.list=get.col.factor(col.factor)
  col.factor.str=get.abbrev.std.stage.cols(in.stages.vec = col.factor)
  gene.euclidean.dist=dist(x = logged.temp.raw.rpkm.mat)
  gene.euclidean.dist.hclust=hclust(d = gene.euclidean.dist)
  #gene.euclidean.dist.hclust.dendrogram=as.dendrogram(gene.euclidean.dist.hclust)
  gene.clusters =cutree(gene.euclidean.dist.hclust,k = as.numeric(gene.cluster.cut.off))
  gene.groups.df=data.frame(samples=names(gene.clusters),group=as.numeric(gene.clusters),row.names = names(gene.clusters))
  gene.clusters.cols <- brewer.pal(as.numeric(gene.cluster.cut.off), "Set3")
  out.list=list()
  gene.groups.list=split(gene.groups.df,f=gene.groups.df$group)
  gene.grps.names=names(gene.groups.list)
  len.gene.grps.names=length(gene.grps.names)
  for(m in 1:len.gene.grps.names){
    grp.name=gene.grps.names[m]
    grp.genes=rownames(gene.groups.list[[grp.name]])
    temp.title.str=paste('Gene.grp: ',grp.name,sep='')
    temp.df=logged.temp.raw.rpkm.mat[grp.genes,]
    out.list[[grp.name]]=temp.df
    heatmap.2(x =temp.df,trace = 'none',cexRow = .4,cexCol = .4,ColSideColors = col.factor.str, main =temp.title.str,labCol = '',margins = c(10,10),key.xlab = 'Z-score',key.ylab = '',key.title = '',col=bluered(10000),scale = 'row')
  }
  #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,cex=.8)
  #out.list=list(rpkm=temp.raw.rpkm.df,meta=meta.df,cor.mat=temp.rpkm.mat)
  return(out.list)
}


create.parental.stage=function(meta.df){
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  specific.stages=as.character(meta.df$development.stage)
  parental.stages=gsub(gsub(gsub(gsub(gsub(gsub(x = specific.stages,pattern = 'T',replacement = 'T'),
                                           pattern = 'LR',replacement = 'R'),
                                      pattern = 'R',replacement = 'R'),
                                 pattern = 'ES',replacement = 'S'),
                            pattern = 'S',replacement = 'S'),
                       pattern = 'ET',replacement = 'T')
  meta.df$parent.development.stage=parental.stages
  return(meta.df)
}


#Plots the overlap btw sc and pop 
plot.sc.pop.diff.expressed.venn.diagram=function(pop.deseq.res.list,scde.res.list){
  pop.deseq.names=names(pop.deseq.res.list)
  sc.scde.names=names(scde.res.list)
  intersect.names=intersect(pop.deseq.names,sc.scde.names)
  all.names=c(pop.deseq.names,sc.scde.names)
  len.all.names=length(all.names)
  len.intersect.names=length(intersect.names)
  for(m in 1:len.all.names){
    #intersect.name=intersect.names[m]
    temp.name=all.names[m]
    temp.name.parts=strsplit(x = temp.name,split = '_vs_')
    rev.name=paste(temp.name.parts[[1]][2],temp.name.parts[[1]][1],sep='_vs_')
    if(temp.name %in% pop.deseq.names | rev.name %in% pop.deseq.names & temp.name %in% sc.scde.names | rev.name %in% sc.scde.names){
      temp.scde.res.list=scde.res.list[[temp.name]]
      if(is.null(temp.scde.res.list)){
        temp.scde.res.list=scde.res.list[[rev.name]]
      }
      if(is.null(temp.scde.res.list)){
        next
      }
      temp.scde.res.df=temp.scde.res.list
      temp.pop.deseq.res.df=pop.deseq.res.list[[temp.name]]
      if(is.null(temp.pop.deseq.res.df)){
        temp.pop.deseq.res.df=pop.deseq.res.list[[rev.name]]
      }
      if(is.null(temp.pop.deseq.res.df)){
        next
      }
      label.pop=paste('Pop',temp.name,sep='.')
      label.sc=paste('SC',temp.name,sep='.')
      pop.sig.genes=rownames(temp.pop.deseq.res.df)
      sc.sig.genes=rownames(temp.scde.res.df)
      temp.list=list(population=pop.sig.genes,sc=sc.sig.genes)
      intersect.genes.venn <- VennFromSets(setList=temp.list)
      #text( x =.5 ,y=.5, labels = title(main =temp.name))
      plot(intersect.genes.venn,doWeights=F,type='circles')
      #plot(intersect.genes.venn,doWeights=F,type='circles')
    }
    else{
      next
    }
  }
}


#Adds the p-adjusted value to the scde results
add.padjust.value.to.scde.res=function(scde.res.list,get.sig.only=T){
  scde.res.names=names(scde.res.list)
  out.list=list()
  len.scde.res.names=length(scde.res.names)
  for(m in 1:len.scde.res.names){
    scde.res.name=scde.res.names[m]
    temp.scde.res.df=scde.res.list[[scde.res.name]]$results
    if(is.null(temp.scde.res.df)){
      next
    }
    no.genes=dim(temp.scde.res.df)[1]
    z.values=as.numeric(as.character(temp.scde.res.df[,'Z']))
    adj.z.values=as.numeric(as.character(temp.scde.res.df[,'cZ']))
    p.values=2*pnorm(-abs(z.values))
    p.adj.values=2*pnorm(-abs(adj.z.values))
    temp.scde.res.df$pvalue=p.values
    temp.scde.res.df$padj=p.adj.values
    if(get.sig.only){
      #temp.scde.res.df=temp.scde.res.df[which(abs(temp.scde.res.df$mle)>=4),]
      temp.scde.res.df=temp.scde.res.df[which(temp.scde.res.df$padj<=.05),]
    }
    out.list[[scde.res.name]]=temp.scde.res.df
  }
  return(out.list)
}


plot.samples.maximum.cor =function(rpkm.df,meta.df,filter.non.var.rows=F,title.str){
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(rpkm.df))
  meta.df=meta.df[colnames(rpkm.df),]
  temp.raw.rpkm.df=rpkm.df
  temp.rpkm.mat=as.matrix(temp.raw.rpkm.df)
  if(filter.non.var.rows){
    temp.rpkm.mat=as.matrix(filter.non.variable.rows(df=temp.rpkm.mat,cut.off=1))
  }
  temp.rpkm.mat=cor(temp.rpkm.mat,method='spearman')
  temp.rpkm.mat[is.na(temp.rpkm.mat)]=0
  col.factor=as.character(meta.df[,'development.stage'])
  #col.factor=as.character(meta.df[,'parent.development.stage'])
  batch.col.factor=as.character(meta.df[,'batch'])
  col.factor.list=get.col.factor(col.factor)
  batch.col.list=get.col.factor(col.factor = batch.col.factor)
  col.factor.str=col.factor.list$col.str
  col.factor.legend.str=col.factor.list$legend.str
  col.factor.legend.col=col.factor.list$legend.cols
  heatmap.2(x =temp.rpkm.mat,trace = 'none',cexRow = .4,cexCol = .4,ColSideColors = col.factor.str,RowSideColors = batch.col.list$col.str ,
            main =title.str,dendrogram='col',labRow = '',col =  )
  legend('left', legend =col.factor.legend.str, fill=col.factor.legend.col, cex=1)
}



plot.samples.cor.at.diff.variation=function(rpkm.df,meta.df){
  #genes.var=as.numeric(apply(rpkm.df,1,sd))
  #var.quantiles=as.numeric(quantile(genes.var))[1:3]
  cut.offs=seq(0,to = 200,by = 10)
  #cut.offs=seq(0,100,50)
  #cut.offs=var.quantiles
  for(m in 1:length(cut.offs)){
    cut.off=cut.offs[m]
    temp.rpkm.df=filter.non.variable.rows(df = rpkm.df,cut.off = cut.off)
    temp.rpkm.df=filter.none.expressed.samples(df = temp.rpkm.df)
    temp.meta.df=meta.df[colnames(temp.rpkm.df),]
    temp.title.str=paste('No. genes cut off: ',toString(cut.off),sep = '')
    show(temp.title.str)
    plot.new()
    text(.5,.5 ,labels = temp.title.str,cex = 1.5)
    plot.samples.correlations.heatmap(in.rpkm.df =temp.rpkm.df,title.str ='SC cor.' ,in.meta.df =temp.meta.df,filter.non.var.rows = T)
  }
}


plot.samples.cor.at.diff.no.detected.genes=function(rpkm.df,meta.df){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df$no.of.detected.genes=as.numeric(as.character(apply(rpkm.df,2,nnzero)))
  cut.offs=seq(100,1000,100)
  for(m in 1:length(cut.offs)){
    cut.off=cut.offs[m]
    temp.meta.df=meta.df[which(meta.df$no.of.detected.genes>=cut.off),]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.genes(input.data =temp.rpkm.df)
    temp.title.str=paste('No. genes cut off: ',toString(cut.off),sep = '')
    cellwidth=ifelse(cut.off>300,8,2)
    cellheight = ifelse(cut.off>300,8,2)
    plot.samples.correlations.heatmap(in.rpkm.df =temp.rpkm.df,title.str =temp.title.str ,in.meta.df =temp.meta.df,filter.non.var.rows = T,
                                      annotation_legend = F,show_colnames = F,show_rownames = F,cellwidth=cellwidth, cellheight=cellheight)
  }
}


# Plot sample's pca at different gene detection counts
plot.samples.pca.at.diff.no.detected.genes=function(rpkm.df,meta.df){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df$no.of.detected.genes=as.numeric(as.character(apply(rpkm.df,2,nnzero)))
  cut.offs=seq(100,1000,100)
  for(m in 1:length(cut.offs)){
    cut.off=cut.offs[m]
    temp.meta.df=meta.df[which(meta.df$no.of.detected.genes>=cut.off),]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.genes(input.data =temp.rpkm.df)
    show(dim(temp.rpkm.df))
    temp.title.str=paste('No. genes cut off: ',toString(cut.off),sep = '')
    show(temp.title.str)
    plot.pca(in.rpkm.df =temp.rpkm.df,meta.data = temp.meta.df,title.str =temp.title.str,var.cut.off = 0,plot.pairs = T,log.rpkm = T )
  }
}


plot.samples.tsne.at.diff.no.detected.genes=function(rpkm.df,meta.df){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df$no.of.detected.genes=as.numeric(as.character(apply(rpkm.df,2,nnzero)))
  cut.offs=seq(100,1000,100)
  for(m in 1:length(cut.offs)){
    cut.off=cut.offs[m]
    temp.meta.df=meta.df[which(meta.df$no.of.detected.genes>=cut.off),]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.genes(input.data =temp.rpkm.df)
    temp.title.str=paste('No. genes cut off: ',toString(cut.off),sep = '')
    #plot.pca(in.rpkm.df =temp.rpkm.df,meta.data = temp.meta.df,title.str =temp.title.str,var.cut.off = 0,plot.pairs = T,log.rpkm = T )
    plot.samples.pca.tsne.clustering(rpkm.df = temp.rpkm.df,meta.df =temp.meta.df,title.str =  temp.title.str,log.rpkms = T)
  }
}


#Removes the 'rna- prefix and -1 suffix from the P. falp gene ids
edit.gene.ids=function(rpkm.df){
  editted.ids =as.character(unlist(lapply(rownames(rpkm.df),function(gene.row){
    gene.row=as.character(gene.row)
    editted.id=sub(pattern='^rna_',replacement='',x=gene.row)
    editted.id=sub(pattern='-1$',replacement='',x=editted.id)
    editted.id=sub(pattern=':rRNA$',replacement='',x=editted.id)
    gene.row.parts=strsplit(x = gene.row,split = '\\+')[[1]]
    len.gene.row.parts=length(gene.row.parts)
    if(len.gene.row.parts>1){
      temp.edit.parts=c()
      for(m in 1:len.gene.row.parts){
        temp.gene.row.parts=gene.row.parts[m]
        editted.id=sub(pattern='^rna_',replacement='',x=temp.gene.row.parts)
        editted.id=sub(pattern='-1$',replacement='',x=editted.id)
        editted.id=sub(pattern=':rRNA$',replacement='',x=editted.id)
        temp.edit.parts=append(temp.edit.parts,editted.id,length(temp.edit.parts))
      }
      editted.id=paste(temp.edit.parts,collapse = '+')
    }
    return(editted.id)
  })))
  rownames(rpkm.df)=editted.ids
  return(rpkm.df)
}


#Edit the data frame samples names
edit.samples.df.names=function(df,trans=F){
  df=if (trans){t(df)} else{df}
  sample.names=as.character(colnames(df))
  no.samples=length(sample.names)
  editted.names=c()
  for(m in 1:no.samples){
    sample=sample.names[m]
    editted.name=sub('^X',replacement='',x=sample,ignore.case=F)
    #editted.name=paste('sample_',editted.name,sep='')
    editted.name=sub('\\.',replacement='_N_',x=editted.name)
    editted.name=sub('-',replacement='_N_',x=editted.name)
    editted.name=sub('0R',replacement='OR',x=editted.name)
    editted.name=sub('10_h',replacement='10h',x=editted.name)
    editted.name=sub('15_h',replacement='15h',x=editted.name)
    editted.name=sub('16_h',replacement='16h',x=editted.name)
    editted.name=sub('20_h',replacement='20h',x=editted.name)
    editted.name=sub('24_h',replacement='24h',x=editted.name)
    editted.name=sub('25_h',replacement='25h',x=editted.name)
    editted.name=sub('30_h',replacement='30h',x=editted.name)
    editted.name=sub('36_h',replacement='36h',x=editted.name)
    editted.name=sub('40_h',replacement='40h',x=editted.name)
    editted.name=sub('[sample_]*',replacement='sample_',x=editted.name)
    editted.names=append(editted.names,editted.name,length(editted.names))
  }
  colnames(df)=editted.names
  out.df=data.frame(subset(df,select=editted.names))
  colnames(out.df)=editted.names
  return(out.df)
}


plot.scatter.per.development.stage=function(rpkm.df,meta.df){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df$no.of.detected.genes=as.numeric(as.character(apply(rpkm.df,2,nnzero)))
  meta.list=split(meta.df,f =meta.df$development.stage)
  development.stages = names(meta.list)
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  development.stages = intersect(ordered.development.stages,development.stages)
  #col.factor.list=get.col.factor(development.stages)
  col.factor.list=get.std.stage.cols(in.stages.vec = development.stages)
  no.development.stage=length(development.stages)
  for (i in 1:no.development.stage){
    development.stage= development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(temp.rpkm.df))
    temp.col=col.factor.list$col.str[i]
    plot.new()
    category=toupper(gsub(pattern = '\\.',replacement = ' ',x = development.stage))
    category=paste(toupper(substring(text = category,1,1)),tolower(substring(text = category,2)),sep='')
    text(x=.5,y=.5,category,cex = 2.0)
    plot.df.scatter(rpkm.df =temp.rpkm.df,title.str =development.stage,col=temp.col)
  }
}
plot.scatter.per.replicate.bulk.sample=function(rpkm.df,meta.df){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df$no.of.detected.genes=as.numeric(as.character(apply(rpkm.df,2,nnzero)))
  meta.list=split(meta.df,f =meta.df$development.stage)
  development.stages = names(meta.list)
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  development.stages = intersect(ordered.development.stages,development.stages)
  #col.factor.list=get.col.factor(development.stages)
  col.factor.list=get.std.stage.cols(in.stages.vec = development.stages)
  no.development.stage=length(development.stages)
  for (i in 1:no.development.stage){
    development.stage= development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(temp.rpkm.df))
    temp.col=col.factor.list$col.str[i]
    #plot.new()
    #category=toupper(gsub(pattern = '\\.',replacement = ' ',x = development.stage))
    #category=paste(toupper(substring(text = category,1,1)),tolower(substring(text = category,2)),sep='')
    #text(x=.5,y=.5,development.stage,cex = 2.0)
    samples=colnames(temp.rpkm.df)
    temp.samples=unique(gsub(pattern='_five',x  = gsub(pattern = 'inhouse_',replacement = '',x = samples),replacement = ''))
    #show(temp.samples)
    len.unique.samples=length(temp.samples)
    for(m in 1:len.unique.samples){
      temp.sample=temp.samples[m]
      match.pattern=gsub(pattern = '_Vy',x=gsub(pattern = '_Mt',replacement = '',x = gsub(pattern = '_five',replacement = '',x = temp.sample)),replacement = '')
      show(match.pattern)
      #show(samples)
      #show(grepl(pattern = match.pattern,x = temp.samples))
      temp.sample.pair=samples[grepl(pattern = match.pattern,x = samples)]
      if(length(temp.sample.pair)>1){
        pair.rpkm.df=filter.none.expressed.genes(input.data = subset(temp.rpkm.df,select = temp.sample.pair))
        psuedo.value=min(pair.rpkm.df[pair.rpkm.df>0])/2
        log.pair.rpkm.df=log2(pair.rpkm.df+psuedo.value)
        spearman.cor=paste('Spearman(', round(x = cor(x = log.pair.rpkm.df[,1],y = log.pair.rpkm.df[,2],method = 'spearman'),2),')',collapse = '')
        smoothScatter(x=log.pair.rpkm.df[,1],y=log.pair.rpkm.df[,2],xlab=colnames(log.pair.rpkm.df)[1],ylab=colnames(log.pair.rpkm.df)[2],pch=19,cex=0.8,main=spearman.cor)
        #plot(x=log.pair.rpkm.df[,1],y=log.pair.rpkm.df[,2],xlab=colnames(log.pair.rpkm.df)[1],ylab=colnames(log.pair.rpkm.df)[2],pch=19,cex=0.8,main=spearman.cor)
        pearson.cor=paste('Pearson(', round(x = cor(x = log.pair.rpkm.df[,1],y = log.pair.rpkm.df[,2],method = 'pearson'),2),')',collapse = '')
        smoothScatter(x=log.pair.rpkm.df[,1],y=log.pair.rpkm.df[,2],xlab=colnames(log.pair.rpkm.df)[1],ylab=colnames(log.pair.rpkm.df)[2],pch=19,cex=0.8,main=pearson.cor)
        #plot(x=log.pair.rpkm.df[,1],y=log.pair.rpkm.df[,2],xlab=colnames(log.pair.rpkm.df)[1],ylab=colnames(log.pair.rpkm.df)[2],pch=19,cex=0.8,main=pearson.cor)
      }
    }
    #plot.df.scatter(rpkm.df =temp.rpkm.df,title.str =development.stage,col=temp.col)
  }
}

plot.per.development.stage.gene.variation=function(rpkm.df,meta.df){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df( meta.df = meta.df)
  meta.list=split(meta.df,f =meta.df$development.stage)
  development.stages = names(meta.list)
  #development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  col.factor.list=get.col.factor(development.stages)
  no.development.stage=length(development.stages)
  for (i in 1:no.development.stage){
    development.stage= development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    samples=rownames(temp.meta.df)
    len.samples=length(samples)
    if(!len.samples>1){
      next
    }
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=as.matrix(filter.none.expressed.samples(filter.none.expressed.genes(temp.rpkm.df)))
    gene.ids=rownames(temp.rpkm.df)
    gene.col=ifelse(grepl(pattern = '^ERCC',x = gene.ids),'blue','black')
    gene.mean.rpkm=as.numeric(rowMeans(temp.rpkm.df))
    gene.var.rpkm=as.numeric(rowSds(temp.rpkm.df))
    #gene.coeff.var=gene.var.rpkm/gene.mean.rpkm^2
    gene.coeff.var=gene.var.rpkm/gene.mean.rpkm
    #plot(NULL,xaxt = 'n',yaxt = 'n',log='xy',xlab='Mean rpkm',ylab='Coefficient of variation',main=development.stage)
    title.str=gsub(pattern = '\\.',' ',x = development.stage)
    mean.expression.weights=1/gene.mean.rpkm
    weighted.mean.expression=mean.expression.weights*gene.mean.rpkm
    plot(x = log10(gene.mean.rpkm),y =gene.coeff.var,pch=19,cex=.7 ,main=title.str,xlab='log10(Mean normalized counts))',ylab='Coefficient of variation',col=gene.col)
    #abline(v = 1,b = 1,col='lightgray')
  }
}


plot.chrom.mapping.rate=function(samples.chrom.mapping.df,title.str='Test plot'){
  colnames.vec=colnames(samples.chrom.mapping.df)
  species.name.lab=colnames.vec[1]
  samples=colnames.vec[2:length(colnames.vec)]
  len.samples=length(samples)
  chrom.list=split(samples.chrom.mapping.df,f=samples.chrom.mapping.df[,species.name.lab])
  species=names(chrom.list)
  len.species=length(species)
  out.mat=matrix(nrow =len.samples,ncol = len.species )
  for(m in 1:len.samples){
    sample=samples[m]
    temp.read.counts=c()
    for(n in 1:len.species){
      sp.name=species[n]
      temp.sp.chrom.df=chrom.list[[sp.name]]
      read.counts=as.numeric(temp.sp.chrom.df[,sample])
      temp.read.counts=append(temp.read.counts,sum(read.counts),length(temp.read.counts))
    }
    out.mat[m,]=temp.read.counts
  }
  rownames(out.mat)=samples
  colnames(out.mat)=species
  out.prop.mat=out.mat/rowSums(out.mat)
  sp.col.list=get.col.factor(col.factor = species)
  barplot2(height = t(out.mat+1),beside = T,log = 'y',las=2,main = title.str,col = sp.col.list$col.str,cex.names = .4)
  legend('left',legend=sp.col.list$legend.str,fill=sp.col.list$legend.col,cex = 1.0)
  barplot2(height = t(out.prop.mat),beside = F,las=2,col = sp.col.list$col.str,cex.names = .4,main = title.str)
  legend('left',legend=sp.col.list$legend.str,fill=sp.col.list$legend.col,cex = 1.0)
}


plot.chrom.mapping.rate.subsets=function(samples.chrom.mapping.df,breaks.interval=10,title.str='Test'){
  df=samples.chrom.mapping.df
  colnames.vec=colnames(df)
  species.col.name=colnames.vec[1]
  samples=colnames.vec[2:length(colnames.vec)]
  no.samples=length(samples)
  if(no.samples==0){
    show('You have no samples')
  }
  else{
    if(no.samples>breaks.interval){
      breaks=seq(1,to =no.samples,by = breaks.interval)
      no.breaks=length(breaks)
      for(m in 1:no.breaks){
        start=breaks[m]
        end=start+(breaks.interval-1)
        if(end<no.samples){
          temp.df=df[,c(species.col.name,samples[seq(start:end)])]
          plot.chrom.mapping.rate(temp.df,title.str = title.str)
        }
        else{
          temp.df=df[,c(species.col.name,samples[seq(start:no.samples)])]
          plot.chrom.mapping.rate(temp.df,title.str = title.str)
        }
      }
    }
    else{
      plot.chrom.mapping.rate(samples.chrom.mapping.df = samples.chrom.mapping.df,title.str = title.str)
    }
  }
}


#Plots the mapping profile in a stage specific manner
plot.samples.mapping.per.development.stage=function(meta.df,chrom.mapping.profile,breaks.interval =10){
  meta.list=split(meta.df,f=meta.df$development.stage)
  stages=names(meta.list)
  ordered.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  stages=intersect(ordered.stages,stages)
  len.stages=length(stages)
  for(n in 1:len.stages){
    stage=stages[n]
    temp.meta.df=meta.list[[stage]]
    species.col.lab=colnames(chrom.mapping.profile)[1]
    intersect.samples=c(species.col.lab,intersect(colnames(chrom.mapping.profile),rownames(temp.meta.df)))
    temp.mapping.profile.df=chrom.mapping.profile[,intersect.samples]
    title.str=gsub('\\.',' ', x = convert.to.title.case(in.str = stage))
    plot.chrom.mapping.rate.subsets(samples.chrom.mapping.df =  temp.mapping.profile.df,breaks.interval = breaks.interval,title.str = title.str )
  }
}
#Plots all the pairwise development stage scatters
plot.all.pairwise.scatter.per.development.stage=function(rpkm.df,meta.df){
  meta.list=split(meta.df,f =meta.df$development.stage)
  development.stages = names(meta.list)
  no.development.stage=length(development.stages)
  for (i in 1:no.development.stage){
    development.stage= development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(temp.rpkm.df))
    temp.meta.df=temp.meta.df[colnames(temp.rpkm.df),]
    all.pairwise.df=create.all.pairwise.scatter.data.frame(temp.rpkm.df)
    plot.df.scatter(rpkm.df =all.pairwise.df,title.str = development.stage )
  }
}


create.all.pairwise.scatter.data.frame=function(in.df){
  samples=colnames(in.df)
  samples.combn.mat=combn(samples,m = 2)
  no.pairs=dim(samples.combn.mat)[2]
  gene.ids=c()
  first.row=c()
  sec.row=c()
  for(n in 1:no.pairs){
    first.sample=samples.combn.mat[,n][1]
    sec.sample=samples.combn.mat[,n][2]
    temp.df=data.frame(first.sample=in.df[,first.sample],sec.sample=in.df[,sec.sample])
    rownames(temp.df)=rownames(in.df)
    temp.df=filter.none.expressed.genes(temp.df)
    no.genes=dim(temp.df)[1]
    if(is.null(dim(temp.df))){
      next
    }
    if(no.genes<1){
      next
    }
    temp.gene.ids=paste(paste(first.sample,sec.sample,sep='_'),rownames(temp.df),sep='_')
    gene.ids=append(gene.ids,temp.gene.ids,after = length(gene.ids))
    first.row=append(first.row,as.numeric(as.character(temp.df$first.sample)),after = length(first.row))
    sec.row=append(sec.row,as.numeric(as.character(temp.df$sec.sample)),after = length(sec.row))
  }
  out.df=data.frame(first.sample=as.numeric(first.row),sec.sample=as.numeric(sec.row))
  rownames(out.df)=gene.ids
  return(out.df)
}


plot.df.scatter=function(rpkm.df,title.str,col='black',pearson=F){
  samples=colnames(rpkm.df)
  len.samples=length(samples)
  rpkm.cor.mat=cor(rpkm.df,method='spearman')
  cor.vec=c()
  if(len.samples>2){
    sample.pairs=combn(x = samples,m = 2)
    no.pairs=dim(sample.pairs)[2]
    tracker=c()
    for(i in 1:no.pairs){
      sample.pair=sample.pairs[,i]
      first.sample=sample.pair[1]
      sec.sample=sample.pair[2]
      x.points=as.numeric(as.character(rpkm.df[,first.sample]))
      y.points=as.numeric(as.character(rpkm.df[,sec.sample]))
      temp.df=data.frame(first.sample=x.points,sec.sample=y.points)
      temp.df=filter.none.expressed.genes(input.data = temp.df)
      x.points=as.numeric(as.character(temp.df[,1]))
      y.points=as.numeric(as.character(temp.df[,2]))
      no.genes=dim(temp.df)[1]
      sample.one.detected.genes=length(x.points[x.points>0])
      sample.two.detected.genes=length(y.points[y.points>0])
      sample.one.drop.out.rate=round(100*((no.genes-sample.one.detected.genes)/no.genes),digits = 2)
      sample.two.drop.out.rate=round(100*((no.genes-sample.two.detected.genes)/no.genes),digits = 2) 
      #spearman.cor=round(cor(x =  x.points,y = y.points,method = 'spearman'),digits = 2)
      #temp.cor.mat=data.frame(cor(temp.df,method='spearman'))
      #spearman.cor=round(as.numeric(temp.cor.mat[1,2]),digits = 2)
      spearman.cor=round(as.numeric(rpkm.cor.mat[first.sample,sec.sample]),digits = 2)
      #if(spearman.cor<.5){
        #next
      #}
      #spearman.cor=round(cor(x =  x.points,y = y.points,method = 'spearman'),digits = 2)
      #spearman.cor=round(cor(x = x.points,y=y.points,method='spearman'),2)
      cor.vec=append(x = cor.vec,values =spearman.cor,after = length(cor.vec))
      main.title= paste ('Correlation: ',toString(spearman.cor),sep = '')
      pearson.cor=round(cor(x =  x.points,y = y.points,method = 'pearson'),digits = 2)
      pearson.main.title= paste ('Pearson cor: ',toString(pearson.cor),sep = '')
      sample.one.drop.out.rate.str=paste('[Dropout-rate: ',toString(sample.one.drop.out.rate),sep='')
      sample.one.drop.out.rate.str=paste(sample.one.drop.out.rate.str,'%]',sep='')
      sample.two.drop.out.rate.str=paste('[Dropout-rate: ',toString(sample.two.drop.out.rate),sep='')
      sample.two.drop.out.rate.str=paste(sample.two.drop.out.rate.str,'%]',sep='')
      first.sample=paste(sub(pattern = '^sample_',replacement = '',first.sample),'Normalized read counts',sep=' ')
      xlab.str=paste(first.sample,sample.one.drop.out.rate.str,sep=' ')
      sec.sample=paste(sub(pattern = '^sample_', replacement = '',sec.sample),'Normalized read counts',sep=' ')
      ylab.str=paste(sec.sample,sample.two.drop.out.rate.str,sep = ' ')
      #plot(x =log2(x.points+1),y=log2(y.points+1),pch=19,cex=.5,col=col,xlab='',ylab='')
      #plot(x =log2(x.points+1),y=log2(y.points+1),pch=19,cex=.5,xlab='',ylab='')
      #title(xlab = xlab.str,cex=.3,ylab =ylab.str,main= main.title,cex.lab=.6)
      geneScatterplot( x =x.points,y = y.points,xlab = xlab.str,ylab =ylab.str,col = '#00207040',main.title= main.title  )
      geneScatterplot(x =x.points ,y=y.points,main= pearson.main.title,xlab=xlab.str,ylab=ylab.str,col = '#00207040')
    }
  }
#   else{
#     
#    
#     
#     if(len.samples==2){
#       
#       first.sample=samples[1]
#       
#       sec.sample=samples[2]
#       
#       x.points=as.numeric(as.character(rpkm.df[,first.sample]))
#       
#       y.points=as.numeric(as.character(rpkm.df[,sec.sample]))
#       
#       no.genes=length(x.points)
#       
#       sample.one.detected.genes=length(x.points[x.points>0])
#       
#       sample.two.detected.genes=length(y.points[y.points>0])
#       
#       sample.one.drop.out.rate=round(100*((no.genes-sample.one.detected.genes)/no.genes),digits = 2)
#       
#       sample.two.drop.out.rate=round(100*((no.genes-sample.two.detected.genes)/no.genes),digits = 2) 
#       
#       #spearman.cor=round(cor(x =  x.points,y = y.points,method = 'spearman'),digits = 2)
#       
#       spearman.cor=round(as.numeric(rpkm.cor.mat[first.sample,sec.sample]),digits = 2)
#       
#       if(spearman.cor<.5){
#         next
#       }
#       
#       cor.vec=append(x = cor.vec,values =spearman.cor,after = length(cor.vec))
#       
#       main.title= paste ('Correlation: ',toString(spearman.cor),sep = '')
#       
#       pearson.cor=round(cor(x =  log2(x.points+.0001),y = log2(y.points+0.0001),method = 'pearson'),digits = 2)
#       
#       pearson.main.title= paste ('Pearson cor: ',toString(pearson.cor),sep = '')
#       
#       xlab.str=paste(first.sample,toString(sample.one.drop.out.rate),': ')
#       
#       ylab.str=paste(sec.sample,toString(sample.two.drop.out.rate),sep = ': ')
#       
#       final.title.str=main.title
#       
#       plot(x =log2(x.points+0.001) ,y=log2(y.points+0.001),main= final.title.str,xlab=xlab.str,ylab=ylab.str,pch=19,cex=.3,col=col)
#       
#       #plot(x =log2(x.points+0.001) ,y=log2(y.points+0.001),main= pearson.main.title,xlab=first.sample,ylab=sec.sample,pch=19,cex=.5)
#      
#     }
#     
#     else{
#       
#       Show('You need atleast 2 samples in a data frame')
#       
#     }
#     
#     
#     
#   }
#   
}


plot.ercc.df.scatter=function(sc.meta.df,ercc.rpkm.df,spike.info.df,title.str,plot.per.stage=T){
  spike.info.list=split(spike.info.df,f=spike.info.df$Spike.in)
  categories=names(spike.info.list)
  for (m in 1:length(categories)){
    category=categories[m]
    #plot.new()
    #text(x = .5,y=.5, labels =category,cex = 2.0 )
    temp.spike.df=spike.info.list[[category]]
    rpkm.df=ercc.rpkm.df[,rownames(temp.spike.df)]
    show(dim(rpkm.df))
    if (dim(rpkm.df)[2]<2){
      next
    }
    samples=colnames(rpkm.df)
    sample.pairs=combn(x = samples,m = 2)
    no.pairs=dim(sample.pairs)[2]
    cor.vec=c()
    pearson.cor.vec=c()
    for(i in 1:no.pairs){
      sample.pair=sample.pairs[,i]
      first.sample=sample.pair[1]
      sec.sample=sample.pair[2]
      first.sample.stage=as.character(sc.meta.df[first.sample,'development.stage'])
      sec.sample.stage=as.character(sc.meta.df[sec.sample,'development.stage'])
      if(plot.per.stage){
        if(!first.sample.stage==sec.sample.stage){
          next
        }
      }
      if(first.sample==sec.sample){
        next
      }
      x.points=as.numeric(as.character(rpkm.df[,first.sample]))
      y.points=as.numeric(as.character(rpkm.df[,sec.sample]))
      no.genes=length(x.points)
      sample.one.detected.genes=length(x.points[x.points>0])
      sample.two.detected.genes=length(y.points[y.points>0])
      if(sample.one.detected.genes==0 | sample.two.detected.genes==0){
        next
      }
      sample.one.drop.out.rate=round(100*((no.genes-sample.one.detected.genes)/no.genes),digits = 2)
      sample.two.drop.out.rate=round(100*((no.genes-sample.two.detected.genes)/no.genes),digits = 2) 
      temp.ercc.count.df=data.frame(first=x.points,sec=y.points)
      temp.ercc.count.df=filter.none.expressed.genes(temp.ercc.count.df)
      spearman.cor=round(cor(temp.ercc.count.df,method = 'spearman')[1,2],digits = 2)
      #spearman.cor=round(cor(x =  x.points,y = y.points,method = 'spearman'),digits = 2)
      cor.vec=append(x = cor.vec,values =spearman.cor,after = length(cor.vec))
      main.title= paste ('Spearman cor: ',toString(spearman.cor),sep = '')
      #pearson.cor=round(cor(x =  x.points,y = y.points,method = 'pearson'),digits = 2)
      pearson.cor=round(cor(temp.ercc.count.df,method = 'pearson')[1,2],digits = 2)
      pearson.cor.vec=append(pearson.cor.vec,pearson.cor,after = length(pearson.cor.vec))
      pearson.main.title= paste ('Pearson cor: ',toString(pearson.cor),sep = '')
      xlab.str=paste(first.sample,toString(sample.one.drop.out.rate),sep=': ')
      ylab.str=paste(sec.sample,toString(sample.two.drop.out.rate),sep = ': ')
      #smoothScatter(x =log10(x.points+0.01) ,y=log10(y.points+0.01),main= main.title,xlab=xlab.str,ylab=ylab.str,pch=19,cex=1.0)
      #abline(0,1,col='lightgray')
      first.sample.lab=paste(first.sample,paste('(',first.sample.stage,')',collapse = ''),sep='')
      sec.sample.lab=paste(sec.sample,paste('(',sec.sample.stage,')',collapse = ''),sep='')
      #plot(x =log10(x.points+0.01) ,y=log10(y.points+0.01),main= pearson.main.title,xlab=first.sample.lab,ylab=sec.sample.lab,pch=19,cex=1.0)
      #abline(0,1,col='lightgray')
    }
    title.str=paste('ERCC (dilution: ' ,sub(pattern = 'ERCC','',category),sep='')
    title.str=paste(title.str,')',sep='')
    #hist(cor.vec,breaks = seq(-1,to = 1,by = .1),xlab='Spearman correlation',main=title.str)
    hist(cor.vec,breaks = seq(-1,to = 1,by = .1),xlab='',main='',ylab='')
    #hist(pearson.cor.vec,breaks = seq(-1,to = 1,by = .1),xlab='Pearson correlation',main=title.str)
  }
}


plot.ercc.barplot.df=function(in.rpkm.df,spike.info.df,title.str='Test'){
  spike.info.list=split(spike.info.df,f=spike.info.df$Spike.in)
  categories=names(spike.info.list)
  for (m in 1:length(categories)){
    category=categories[m]
    temp.spike.df=spike.info.list[[category]]
    rpkm.df=in.rpkm.df[,rownames(temp.spike.df)]
    if (dim(rpkm.df)[2]<2){
      next
    }
    plot.new()
    text(x = .5,y=.5, labels =category,cex = 2.0 )
    ercc.names=rownames(rpkm.df)
    samples=colnames(rpkm.df)
    for (n in 1:dim(rpkm.df)[1]){
      spike.in.name=rownames(rpkm.df[n,])[1]
      spike.in.rpkm=as.numeric(as.character(rpkm.df[n,]))
      spike.in.filt.rpkm=spike.in.rpkm[spike.in.rpkm>0]
      if(length( spike.in.filt.rpkm)==0){
        next
      }
      rpkm.row=log10(spike.in.rpkm+1)
      barplot(rpkm.row,main = spike.in.name,ylab = 'log10(rpkm+1)',names.arg =samples,las=2 ,cex.names = .3)
    }
  }
}


get.ercc.rpkm=function(ercc.counts.df,ercc.len.df){
  samples.size=as.numeric(apply(ercc.counts.df,2,sum))
  ercc.len=as.numeric(ercc.len.df[rownames(ercc.counts.df),'chrom.len'])
  ercc.rpkm.df=data.frame(t(t(ercc.counts.df*1000000000)*samples.size)/ercc.len)
  return(ercc.rpkm.df)
}


plot.samples.pfa.hg.chrom.mapping.rate=function(sc.chrom.mapping.df,sc.meta.df){
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df)
  sc.chrom.mapping.list=split(sc.chrom.mapping.df,f = sc.chrom.mapping.df$sample_species)
  sc.chrom.mapping.hg.df=sc.chrom.mapping.list[['hg']]
  sc.chrom.mapping.hg.dim=dim(sc.chrom.mapping.hg.df)
  sc.chrom.mapping.hg.df=sc.chrom.mapping.hg.df[,1:sc.chrom.mapping.hg.dim[2]-1]
  sc.chrom.mapping.hg.df=sc.chrom.mapping.hg.df[,rownames(sc.meta.df)]
  sc.chrom.mapping.hg.sum.df=as.numeric(as.character(apply(sc.chrom.mapping.hg.df,2,sum)))
  sc.chrom.mapping.pfa.df=sc.chrom.mapping.list[['pfa']]
  sc.chrom.mapping.pfa.dim=dim(sc.chrom.mapping.pfa.df)
  sc.chrom.mapping.pfa.df=sc.chrom.mapping.pfa.df[,1:sc.chrom.mapping.pfa.dim[2]-1]
  sc.chrom.mapping.pfa.df=sc.chrom.mapping.pfa.df[,rownames(sc.meta.df)]
  sc.chrom.mapping.pfa.sum.df=as.numeric(as.character(apply(sc.chrom.mapping.pfa.df,2,sum)))
  sc.pfa.hg.mapped.reads.sum.df=data.frame(hg=sc.chrom.mapping.hg.sum.df,pfa=sc.chrom.mapping.pfa.sum.df)
  rownames(sc.pfa.hg.mapped.reads.sum.df)=rownames(sc.meta.df)
  total.mapped.reads=sc.chrom.mapping.hg.sum.df+sc.chrom.mapping.pfa.sum.df
  pfa.mapped.prop=100*(sc.chrom.mapping.pfa.sum.df/total.mapped.reads)
  hg.mapped.prop=100*(sc.chrom.mapping.hg.sum.df/total.mapped.reads)
  #barplot(height = t(as.matrix(sc.pfa.hg.mapped.reads.sum.df)),las=2,log = 'y',col=c('red','green'),border = F,,xaxt='none')
  ylim.vec=c(0,max(sc.pfa.hg.mapped.reads.sum.df)+1000000)
  # barplot2(height = t(log10(as.matrix(sc.pfa.hg.mapped.reads.sum.df))),
  #         las=2,col=c('green','red'),border = F,xaxt='none', beside = F)
  # 
  # 
  sc.meta.df=sc.meta.df[rownames(sc.pfa.hg.mapped.reads.sum.df),]
  sc.meta.list=split(sc.meta.df,sc.meta.df$development.stage)
  stages=names(sc.meta.list)
  for (stage in stages){
    temp.meta.df=sc.meta.list[[stage]]
    temp.sc.pfa.hg.mapped.reads.sum.df=sc.pfa.hg.mapped.reads.sum.df[rownames(temp.meta.df),]
    plot.samples.chrom.mapping.proportions.barplot.for.subset.samples(meta.df = temp.meta.df,
                                                                      chrom.mapping.df = temp.sc.pfa.hg.mapped.reads.sum.df,
                                                                      title.str = stage,breaks.interval = 20)
    # barplot2(height = t(log10(as.matrix(temp.sc.pfa.hg.mapped.reads.sum.df))),
    #          las=2,col=c('green','red'),border = F,xaxt='none', beside = T,main=stage)
  }
  plot.samples.mapping.proportions.barplot.for.subset.samples(df =sc.meta.df,cut.off = 0,title.str = 'Mapping profile',
                                                              breaks.interval = length(rownames(sc.pfa.hg.mapped.reads.sum.df)) )
  temp.df=data.frame(pfa=pfa.mapped.prop,hg=hg.mapped.prop)
  barplot(height = t(as.matrix(temp.df)),las=2,col=c('red','green'),border = F,xaxt='none')
  samples.names=rownames(sc.meta.df)
  rownames(temp.df)=samples.names
  samples.names.sorted=sort(samples.names)
  temp.df=temp.df[samples.names.sorted,]
  sc.meta.df[,'sc.chrom.mapping.pfa.sum.df']=sc.chrom.mapping.pfa.sum.df
  temp.meta.df=sc.meta.df[samples.names.sorted,]
  temp.meta.df[,'pfa.genome.mapping.proportion']=as.numeric(as.character(temp.df$pfa))
  temp.meta.df[,'hg.genome.mapping.proportion']=as.numeric(as.character(temp.df$hg))
  temp.mat=as.matrix(temp.df)
  pfa.temp.mat=as.matrix(temp.df[which(temp.df$pfa>=50),])
  hg.temp.mat=as.matrix(temp.df[which(temp.df$hg>50),])
  stages.lab=as.character(temp.meta.df[rownames(pfa.temp.mat),'development.stage'])
  #temp.barplot=barplot(t(pfa.temp.mat),main='High in Pfa',col=c('red','green'),ylab='% of reads',xlab='',xaxt='none',beside=F,ylim = c(0,115))
  #temp.barplot=barplot(t(pfa.temp.mat),main='High in Pfa',col=c('red','green'),ylab='% of reads',xlab='',beside=F,xaxt='none')
  temp.barplot=barplot(t(temp.mat),main='',col=c('red','green'),ylab='',xlab='',beside=F,xaxt='none',border = NA)
  temp.barplot=barplot(t(pfa.temp.mat),main='',col=c('red','green'),ylab='',xlab='',beside=F,xaxt='none',border = NA)
  abline(b = 1,h = 50)
  #text(temp.barplot, par("usr")[3], labels =stages.lab, adj = c(1.1,1.1), xpd = T, cex=.5)
  plot.new()
  #legend('topright',legend=c('Pfa','hg'),fill=c('red','green'),cex=.35)
  legend('center',legend=c('',''),fill=c('red','green'),cex=2.0,box.lty = 0,border = NA)
  #temp.barplot.2=barplot(t(hg.temp.mat),main='High in hg19',col=c('red','green'),ylab='% of reads',xlab='',xaxt='none',beside=F,ylim = c(0,115))
  temp.barplot.2=barplot(t(hg.temp.mat),main='',col=c('red','green'),ylab='',xlab='',xaxt='none',beside=F,border = NA)
  abline(b = 1,h = 50)
  stages.2.lab=as.character(temp.meta.df[rownames(hg.temp.mat),'development.stage'])
  #text(temp.barplot.2, par("usr")[3], labels = stages.2.lab,  adj = c(1.1,1.1), xpd = T, cex=.5)
  #legend('topright',legend=c('Pfa','hg'),fill=c('red','green'),cex=.35)
  plot.new()
  legend('center',legend=c('',''),fill=c('red','green'),cex=1.5)
  return(temp.meta.df)
}


plot.contamint.read.proportions=function(contaminant.df,title.str,genome.grp){
  sample.names=rownames(contaminant.df)
  #test.df=data.frame(total.reads=as.numeric(contaminant.df$Number_of_input_reads_),virus=as.numeric(contaminant.df$Number_of_reads_mapped_to_multiple_loci_)+as.numeric(contaminant.df$Uniquely_mapped_reads_number_),row.names = sample.names)
  test.df=data.frame(total.reads=as.numeric(contaminant.df$Number_of_input_reads_),clean.reads=as.numeric(contaminant.df$Number_of_input_reads_)-as.numeric(as.numeric(contaminant.df$Number_of_reads_mapped_to_multiple_loci_)+as.numeric(contaminant.df$Uniquely_mapped_reads_number_)),virus=as.numeric(contaminant.df$Number_of_reads_mapped_to_multiple_loci_)+as.numeric(contaminant.df$Uniquely_mapped_reads_number_),row.names = sample.names)
  barplot.var=barplot(as.matrix(t(test.df[,2:3])),main=title.str,col=c('blue','red'),xaxt='none',ylab='Read counts')
  text(barplot.var, par("usr")[3], labels = sample.names, srt = 90, adj = c(1.1,1.1), xpd = T, cex=.3)
  legend('topright',legend=c('Input reads',"'contaminants'"),fill=c('blue','red'))
  prop.mat=100*data.frame(clean=test.df[,2]/test.df[,1],contaminants=test.df[,3]/test.df[,1])
  barplot.var.2=barplot(as.matrix(t(prop.mat)),main=title.str,col=c('blue','red'),xaxt='none',ylab='Proportion(%)',ylim=c(0,110))
  text(barplot.var.2, par("usr")[3], labels = sample.names, srt = 90, adj = c(1.1,1.1), xpd = T, cex=.3)
  legend('topright',legend=c('Clean reads',"'contaminants'"),fill=c('blue','red'),cex=.4)
  abline(h = 50,b = 1)
}


plot.development.stage.gene.counts=function(rpkm.df,meta.df){
  rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df =meta.df )
  meta.df$development.stage=factor(x = meta.df$development.stage, levels = as.factor(c('R','LR','ET','T','ES','S')))
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  #development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(development.stages,names(meta.list))
  col.fact.list=get.col.factor(col.factor = development.stages)
  detected.genes.count.list=list()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
    no.detected.genes=as.numeric(as.character(apply(temp.rpkm.df,2,nnzero)))
    detected.genes.count.list[[development.stage]]=no.detected.genes
  }
  col.str=col.fact.list$col.str
  legend.str=col.fact.list$legend.str
  legend.col=col.fact.list$legend.col
  title.str='Gene counts (normalized count > 0 in atleast 1 SC)'
  title.str=''
  #boxplot(x = detected.genes.count.list,las=2,cex=.3,names = rep('',times = length(development.stages)),ylab=title.str)
  #boxplot(x = detected.genes.count.list,las=2,cex=.3,col = col.str,names = rep('',times = length(development.stages)),ylab=title.str)
  #boxplot2(x = detected.genes.count.list,las=2,cex=.3,col =abbrev.std.stages.col.vec[development.stages],names = development.stages,ylab=title.str)
  boxplot(x = detected.genes.count.list,ylab='',names=rep('',times=length(development.stages)),pch=19,col =abbrev.std.stages.col.vec[development.stages],cex=1.5)
  #boxplot2(x = detected.genes.count.list)
  meta.df=meta.df[colnames(rpkm.df),]
  no.of.genes.vec=as.numeric(apply(rpkm.df,2,nnzero))
  lowest.cut.off=200
  highest.cut.off=as.numeric(quantile(x = no.of.genes.vec, .993))
  highest.cut.off=1750
  meta.df$detected.genes.counts=no.of.genes.vec
  intervals=seq(from = 0,to = max(no.of.genes.vec)+100,by = 100)
  #hist(x = no.of.genes.vec,main = 'No of genes',breaks =intervals,xlab = 'Gene counts',ylab='Frequency' )
  #abline(v = lowest.cut.off,b=1,col='lightgray')
  #abline(v=highest.cut.off,b=1,col='lightgray')
  #plot.new()
  #legend('center',legend=rep('',times =length(legend.str)),fill=legend.col,box.lty = 0,cex=1.5)
  out.list=list(rpkm=rpkm.df,meta=meta.df)
  return(out.list)
}


plot.development.stage.gene.counts.btw.sc.and.bulk=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df)
  #sc.rpkm.df=filter.rpkm.less.than.cutoff(df=sc.rpkm.df,rpkm.cutoff = 1,no.samples = 1)
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$development.stage)
  pop.rpkm.df=filter.rpkm.less.than.cutoff(df = pop.rpkm.df,rpkm.cutoff = 5.0,
                                           no.samples = 1)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = pop.meta.df)
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$development.stage)
  pop.sc.development.stages=intersect(names(sc.meta.list),names(pop.meta.list))
  #development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  intersect.development.stages=intersect(development.stages,
                                         pop.sc.development.stages)
  development.stages=intersect.development.stages
  #col.fact.list=get.col.factor(col.factor = development.stages)
  col.str=c()
  detected.genes.count.list=list()
  for (i in 1:length(development.stages)){
    temp.detected.genes.count.list=list()
    development.stage=development.stages[i]
    sc.temp.meta.df=sc.meta.list[[development.stage]]
    sc.temp.rpkm.df=sc.rpkm.df[,rownames(sc.temp.meta.df)]
    sc.temp.rpkm.df=filter.none.expressed.genes(input.data = sc.temp.rpkm.df)
    #sc.temp.rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df = sc.temp.rpkm.df,sample.count = 3)
    sc.temp.rpkm.df=filter.rpkm.less.than.cutoff(df = sc.temp.rpkm.df,
                                                 no.samples = 1,rpkm.cutoff = 1)
    pop.temp.meta.df=pop.meta.list[[development.stage]]
    pop.temp.rpkm.df=filter.none.expressed.genes(input.data = 
                                                   pop.rpkm.df[,rownames(pop.temp.meta.df)])
    sc.gene.counts.vec=as.numeric(apply(sc.temp.rpkm.df,2,nnzero))
    pop.gene.counts.vec=as.numeric(apply(pop.temp.rpkm.df,2,nnzero))
    sc.total.genes.count=nrow(sc.temp.rpkm.df)
    sc.mean.gene.counts=floor(mean(sc.gene.counts.vec))
    pop.total.genes.count=nrow(pop.temp.rpkm.df)
    sc.pop.intersect.genes.vec=intersect(rownames(sc.temp.rpkm.df),
                          rownames(pop.temp.rpkm.df))
    sc.bulk.intersect.count=length(sc.pop.intersect.genes.vec)
    sc.bulk.intersect.percent=round(length(sc.pop.intersect.genes.vec)
                                    /pop.total.genes.count,4)*100
    pop.mean.genes.count=floor(mean(pop.gene.counts.vec))
    c(development.stage,sc.total.genes.count,sc.mean.gene.counts,
            sc.bulk.intersect.percent)
    sc.str=paste('S:',development.stage,':',sc.total.genes.count,':',
                 ':',sc.bulk.intersect.count,':',sc.bulk.intersect.percent,collapse = '')
    pop.str=paste('P:',development.stage,':',pop.total.genes.count,
                  ':',pop.mean.genes.count,collapse = '')
    print(c(development.stage,'bulk',nrow(pop.temp.rpkm.df),'sc',nrow(sc.temp.rpkm.df),
            'intersect',sc.bulk.intersect.count,'percentage',sc.bulk.intersect.percent,'mean',sc.mean.gene.counts))
    temp.detected.genes.count.list[[sc.str]]=sc.gene.counts.vec
    temp.detected.genes.count.list[[pop.str]]=pop.gene.counts.vec
    #detected.genes.count.list[[development.stage]]=temp.detected.genes.count.list
    detected.genes.count.list[[sc.str]]=sc.gene.counts.vec
    detected.genes.count.list[[pop.str]]=pop.gene.counts.vec
    col.str=append(col.str,values = c('red','blue'),after = length(col.str))
  }
  #boxplot(x = detected.genes.count.list,las=2,cex=.3,ylab='Fraction of genes detected in bulk-RNAseq (%)',col = col.str)
  boxplot(x = detected.genes.count.list,las=2,cex=1.3,
          outcol=col.str,border=T,frame.plot=F,ylim=c(0,5000),cex.axis=0.3,
          ylab='Number of genes detected',col=col.str,pch=19)
  #plot.new()
  #legend('center',fill=c('red','blue'),legend=c('SCs','Bulk'),lty = 0,border = NA,box.lty = 0,cex=2.5)
}
plot.development.stage.gene.counts.per.count.intervals=function(counts.df,meta.df,title.str ='Test'){
  meta.df=meta.df[colnames(counts.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df =meta.df )
  meta.df$development.stage=factor(x = meta.df$development.stage, levels = as.factor(c('R','LR','ET','T','ES','S')))
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  col.fact.list=get.col.factor(col.factor = development.stages)
  detected.genes.count.list=list()
  #gene.count.starts=c(1,10,100,500)
  gene.count.starts=c(1,5,10,100)
  len.gene.count.starts=length(gene.count.starts)
  for(m in 1:len.gene.count.starts){
    start.count=as.numeric(gene.count.starts[m])
    start.interval=ifelse(m>1,start.count+1,start.count)
    end.interval=ifelse(m==len.gene.count.starts,Inf,gene.count.starts[m+1])
    #temp.count.df=filter.genes.with.counts.btw.intervals(count.df=counts.df,min.count=start.interval,max.count = end.interval)
    temp.count.df=counts.df
    temp.count.df[temp.count.df<start.interval|temp.count.df>end.interval]=0
    temp.detected.genes.count.list=list()
    for (i in 1:length(development.stages)){
      development.stage=development.stages[i]
      temp.meta.df=meta.list[[development.stage]]
      temp.counts.df=temp.count.df[,rownames(temp.meta.df)]
      no.detected.genes=as.numeric(as.character(apply(temp.counts.df,2,nnzero)))
      temp.detected.genes.count.list[[development.stage]]=no.detected.genes
    }
    end.interval.lab=ifelse(start.interval==gene.count.starts[len.gene.count.starts],'above',end.interval)
    interval.lab=paste(c(start.interval,end.interval.lab),collapse='-')
    show(interval.lab)
    detected.genes.count.list[[interval.lab]]=temp.detected.genes.count.list
    boxplot(temp.detected.genes.count.list,names = rep('',times = length(names(temp.detected.genes.count.list))),main=interval.lab,pch=19,cex=.8,col = abbrev.std.stages.col.vec[names(temp.detected.genes.count.list)])
    #boxplot(temp.detected.genes.count.list,main='',names =NA)
  }
  return(detected.genes.count.list)
}
plot.development.stage.gene.expression.per.count.intervals=function(rpkm.df,meta.df,title.str ='Test'){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df =meta.df )
  meta.df$development.stage=factor(x = meta.df$development.stage, levels = as.factor(c('R','LR','ET','T','ES','S')))
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  col.fact.list=get.col.factor(col.factor = development.stages)
  detected.genes.count.list=list()
  #gene.count.starts=c(1,10,100,500)
  gene.rpkm.starts=c(0.01,.1,1)
  len.gene.rpkm.starts=length(gene.rpkm.starts)
  for(m in 1:len.gene.rpkm.starts){
    start.rpkm=as.numeric(gene.rpkm.starts[m])
    start.interval=ifelse(m==1,0,start.rpkm)
    end.interval=ifelse(m==1,start.rpkm,gene.rpkm.starts[m+1])
    if(m==len.gene.rpkm.starts){
      end.interval=ifelse(m==len.gene.rpkm.starts,Inf,gene.rpkm.starts[m+1])
    }
    temp.rpkm.df=filter.genes.with.rpkm.btw.intervals(rpkm.df = rpkm.df,min.rpkm = start.interval,max.rpkm =  end.interval)
    temp.detected.genes.rpkm.list=list()
    for (i in 1:length(development.stages)){
      development.stage=development.stages[i]
      temp.meta.df=meta.list[[development.stage]]
      samples.intersect=intersect(colnames(temp.rpkm.df),rownames(temp.meta.df))
      if(length(samples.intersect)==0){
        next
      }
      temp.rpkm.df=subset(temp.rpkm.df,select = rownames(temp.meta.df))
      temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
      no.detected.genes=as.numeric(as.character(apply(temp.rpkm.df,2,nnzero)))
      temp.detected.genes.rpkm.list[[development.stage]]=no.detected.genes
    }
    end.interval.lab=ifelse(is.infinite(end.interval),'above',end.interval)
    interval.lab=paste(c(start.interval,end.interval.lab),collapse='-')
    detected.genes.count.list[[interval.lab]]=temp.detected.genes.rpkm.list
    boxplot2(temp.detected.genes.rpkm.list,main=interval.lab)
  }
  return(detected.genes.count.list)
}
filter.genes.with.counts.btw.intervals=function(count.df,min.count=1,max.count=Inf){
  min.count=as.numeric(min.count)
  max.count=as.numeric(max.count)
  count.df[count.df<min.count|count.df>max.count]=0
  out.count.df=filter.none.expressed.samples(filter.none.expressed.genes(count.df))
  return(out.count.df)
}
filter.genes.with.rpkm.btw.intervals=function(rpkm.df,min.rpkm=1,max.rpkm=Inf){
  min.rpkm=as.numeric(min.rpkm)
  max.rpkm=as.numeric(max.rpkm)
  rpkm.df[rpkm.df<min.rpkm|rpkm.df>max.rpkm]=0
  out.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(rpkm.df))
  return(out.rpkm.df)
}
plot.development.stage.gene.detection.rate=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$development.stage)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$development.stage)
  pop.sc.development.stages=intersect(names(sc.meta.list),names(pop.meta.list))
  development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  intersect.development.stages=intersect(development.stages,pop.sc.development.stages)
  development.stages=intersect.development.stages
  col.fact.list=get.col.factor(col.factor = development.stages)
  detected.genes.count.list=list()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    sc.temp.meta.df=sc.meta.list[[development.stage]]
    sc.temp.rpkm.df=sc.rpkm.df[,rownames(sc.temp.meta.df)]
    sc.temp.rpkm.df=filter.none.expressed.genes(input.data = sc.temp.rpkm.df)
    pop.temp.meta.df=pop.meta.list[[development.stage]]
    pop.temp.rpkm.df=pop.rpkm.df[,rownames(pop.meta.df)]
    sensitivity.vec=get.sensitivity.btw.two.df(first.df=sc.temp.rpkm.df,sec.df=pop.temp.rpkm.df)
    detected.genes.count.list[[development.stage]]=sensitivity.vec
  }
  col.str=col.fact.list$col.str
  legend.str=col.fact.list$legend.str
  legend.col=col.fact.list$legend.col
  boxplot(x = detected.genes.count.list,las=2,cex=.3,ylab='Fraction of genes detected in bulk-RNAseq (%)',col = col.str)
  #boxplot(x = detected.genes.count.list,las=2,cex=.3,ylab='Fraction of genes detected in bulk-RNAseq (%)',col = col.str,names = rep('',times = length(col.str)))
  boxplot(x = detected.genes.count.list,las=2,cex=.3,ylab='',names = rep('',times = length(col.str)))
  #legend('topright',legend=legend.str,fill=legend.col,cex=1.5,box.lty  = 0)
  return(detected.genes.count.list)
}
plot.gene.counts.hist=function(rpkm.df,breaks=100,title.str='Test'){
  no.genes.per.sample=as.numeric(apply(rpkm.df,2,nnzero))
  title.str=convert.to.title.case(title.str)
  hist(x = no.genes.per.sample,main='',xlab='',ylab='')
}
get.sensitivity.btw.two.df=function(first.df,sec.df){
  first.samples=colnames(first.df)
  sec.samples=colnames(sec.df)
  len.first.samples=length(first.samples)
  len.sec.samples=length(sec.samples)
  sensitivity.vec=c()
  for(m in 1:len.first.samples){
    for(n in 1:len.sec.samples){
      temp.first.sample=first.samples[m]
      temp.sec.sample=sec.samples[n]
      temp.first.genes=rownames(first.df[which(first.df[,temp.first.sample]>0),])
      temp.sec.genes=rownames(sec.df[which(sec.df[,temp.sec.sample]>0),])
      detected.genes=intersect(temp.first.genes,temp.sec.genes)
      sensitivity=as.numeric(100*(length(detected.genes)/length(temp.sec.genes)))
      sensitivity.vec=append(sensitivity.vec,sensitivity,after = length(sensitivity.vec))
    }
  }
  return(sensitivity.vec)
}
plot.development.stage.sample.distrubution=function(rpkm.df,meta.df,rpkm.cut.off=0){
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  detected.genes.count.list=list()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
    create.density.distribution.plot.df(df =temp.rpkm.df,meta.df =temp.meta.df,sample.dist = T,title.str = development.stage,log = T,filter.zeros = T)
  }
}
#Sample and gene distribution
create.density.distribution.plot.df=function(df,meta.df,log=F,sample.dist=F,filter.zeros=F,title.str){
  in.meta.df=data.frame(t(subset(t(meta.df),select=colnames(df))))
  x.lab ='RPKM'
  y.lab='Density'
  main.label = title.str
  if(log){
    df[df=0]=.00001
    df=log2(df)
    x.lab = 'Log2(Normalized expression)'
  }
  if(sample.dist){
    main.label =title.str
    df=t(df)
  }
  density.list = apply(df, 1, function(x){x=x[x!=-Inf];x[x!=Inf];if(length(x)>=2){density(x)}})
  density.list.length=length(density.list)
  if(filter.zeros){
    density.list = apply(df, 1, function(x){x=x[x!=-Inf];x[x!=Inf];if(length(x)>=2){density(x)}})
  }
  n.samples = length(density.list)
  xlim = range(unlist(lapply(density.list, '[[', 'x')))
  ylim = range(unlist(lapply(density.list, '[[', 'y')))
  plot(x = xlim, y = ylim, main=main.label, xlab = x.lab,ylab=y.lab,type='n')
  colour.factor.list=get.col.factor(col.factor=as.character(in.meta.df$development.stage))
  colour.list=colour.factor.list[['col.str']]
  for(j.sample in 1:n.samples){
    lines(density.list[[j.sample]],col=colour.list[j.sample])
  }
}
plot.development.stage.cor.distribution=function(rpkm.df,meta.df,method='spearman'){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  development.stages=intersect(ordered.development.stages,ordered.development.stages)
  development.stages.col.list=get.col.factor(col.factor = development.stages)
  development.stages.col.str=development.stages.col.list$col.str
  detected.genes.count.list=list()
  #plot.new()
  #legend('center',legend=rep('',times = length(development.stages.col.list$legend.str)),fill=development.stages.col.list$legend.col,cex=1.5,)
  cor.list=list()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    if(dim(temp.rpkm.df)[2]<=3){
      next
    }
    temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
    cor.vec=get.cor.list(df =temp.rpkm.df,in.method = method)
    cor.list[[development.stage]]=cor.vec
    #vioplot(x = cor.vec,names =c(development.stage),col = 'white')
    #title(ylab="Spearman cor.")
    #boxplot(x =cor.vec,names = c(development.stage),las=2 )
    #abline(h = mean(cor.vec),b = 1,col='red')
    #abline(h = median(cor.vec),b = 1,col='blue')
    #legend('topright',legend=c('mean cor.','median cor'),fill=c('red','blue'))
  }
  boxplot(cor.list,las=2,names=rep('',times=length(names(cor.list))),xlab='',ylab='',main='',pch=19)
  boxplot2(cor.list,las=2,xlab='',ylab='',main='')
}
plot.development.stage.cv.distribution=function(rpkm.df,meta.df){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(ordered.development.stages,ordered.development.stages)
  development.stages.col.list=get.col.factor(col.factor = development.stages)
  development.stages.col.str=development.stages.col.list$col.str
  detected.genes.count.list=list()
  #plot.new()
  #legend('center',legend=rep('',times = length(development.stages.col.list$legend.str)),fill=development.stages.col.list$legend.col,cex=1.5,)
  cv.list=list()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    if(dim(temp.rpkm.df)[2]<=3){
      next
    }
    temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
    cv=as.numeric(apply(temp.rpkm.df,1,sd))/as.numeric(apply(temp.rpkm.df,1,mean))
    cv.list[[development.stage]]=cv
    #vioplot(x = cor.vec,names =c(development.stage),col = 'white')
    #title(ylab="Spearman cor.")
    #boxplot(x =cor.vec,names = c(development.stage),las=2 )
    #abline(h = mean(cor.vec),b = 1,col='red')
    #abline(h = median(cor.vec),b = 1,col='blue')
    #legend('topright',legend=c('mean cor.','median cor'),fill=c('red','blue'))
  }
  boxplot(cv.list,las=2,names=rep('',times=length(names(cv.list))),xlab='',ylab='',main='')
  #boxplot(cv.list,las=2,xlab='',ylab='',main='')
}
get.cor.list=function(df,in.method='spearman'){
  cor.mat=cor(df,method = in.method)
  samples=colnames(df)
  len.samples=length(samples)
  cor.vec=c()
  samples.combn=combn(samples,2)
  no.pairs=dim(samples.combn)[2]
  if(len.samples==2){
    cor.vec=append(cor.vec,cor.mat[1,2],length(cor.vec))
  }
  else if(len.samples>2){
    for(n in 1:no.pairs){
      temp.pair=samples.combn[,n]
      first.sample=temp.pair[1]
      sec.sample=temp.pair[2]
      temp.cor=as.numeric(cor.mat[first.sample,sec.sample])
      cor.vec=append(cor.vec,temp.cor,length(cor.vec))
    }
  }
  else{
    cor.vec=append(cor.vec,NA,length(cor.vec))
  }
  return(cor.vec)
}
plot.gene.body.cov.vs.rpkm=function(gene.body.df,rpkm.df){
  samples=colnames(rpkm.df)
  gene.body.df=gene.body.df[,samples]
  intersect.genes=intersect(rownames(gene.body.df),rownames(rpkm.df))
  gene.body.df=gene.body.df[intersect.genes,]
  rpkm.df=rpkm.df[intersect.genes,]
  len.intersect.genes=length(intersect.genes)
  len.samples=length(samples)
  row.ids=c()
  rpkms=c()
  coverages=c()
  for(m in 1:len.intersect.genes){
    for(n in 1:len.samples){
      gene=intersect.genes[m]
      sample=samples[n]
      rpkm=rpkm.df[gene,sample]
      covarage=gene.body.df[gene,sample]
      if(covarage>=.1){
        next
      }
      row.id=paste(gene,sample,sep = ':')
      row.ids=append(row.ids,row.id,length(row.ids))
      rpkms=c(rpkms,rpkm,length(rpkms))
      coverages=c(coverages,coverage,length(coverages))
    }
  }
  out.df=data.frame(coverage=coverages,normalized.counts=rpkms)
  rownames(out.df)=row.ids
  return(out.df)
}
plot.development.stage.cor.heatmap=function(rpkm.df,meta.df,in.method='spearman',euclidian=F,plot.var.genes=T){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(ordered.development.stages,development.stages)
  detected.genes.count.list=list()
  out.list=list()
  var.genes.list=list()
  samples.to.keep=c()
  genes.to.keep=c()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    #if(!development.stage=='ring'){
      #next
    #}
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    if(is.null(dim(temp.rpkm.df))){
      next
    }
    no.samples.cut.off=round(.2*dim(temp.rpkm.df)[2])
    temp.rpkm.df=filter.rpkm.less.than.cutoff(df =temp.rpkm.df,rpkm.cutoff = 0,no.samples =  no.samples.cut.off)
    temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
    temp.meta.df$no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    temp.meta.df=temp.meta.df[which(temp.meta.df$no.of.detected.genes>=200),]
    no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    temp.rpkm.df=temp.rpkm.df[,rownames(temp.meta.df)]
    out.list[[development.stage]]=list(rpkm=temp.rpkm.df,meta=temp.meta.df)
    samples.to.keep=append(samples.to.keep,rownames(temp.meta.df),length(samples.to.keep))
    genes.to.keep=append(genes.to.keep,rownames(temp.rpkm.df),length(genes.to.keep))
    lab.rows = as.character(temp.meta.df[colnames(temp.rpkm.df),'no.of.detected.genes'])
    dist.mat=matrix()
    if(euclidian){
      euclidian.dist.mat=as.matrix(dist(x = t(log2(temp.rpkm.df+1))))
      dist.mat= euclidian.dist.mat
      show('euclidian')
    }
    else{
      if(in.method=='spearman'){
        dist.mat=cor(x = temp.rpkm.df,method = in.method)
        show('spearman')
      }
      else{
        log.temp.rpkm.df=log2(temp.rpkm.df+1)
        show('pearson')
        dist.mat=cor(x = log.temp.rpkm.df,method = in.method)
      }
    }
    title.str=sub('\\.',' ',x = development.stage)
    dist.mat=round(dist.mat,digits = 2)
    temp.rpkm.mat=as.matrix(log2(temp.rpkm.df+1))
    heatmap.2(x = dist.mat,notecex=.1,cellnote = dist.mat,trace = 'none',main =title.str,cexRow = .3,cexCol = .3,las=2,margins = c(11,11),labCol = lab.rows,key.xlab = '',key.ylab = '',key.title = '',col = greenred(7))
    #heatmap.2(x = temp.rpkm.mat,trace = 'none',main =title.str,cexRow = .3,cexCol = .3,las=2,margins = c(11,11),labCol = lab.rows,key.xlab = '',key.ylab = '',key.title = '',col = greenred(7))
    #heatmap.2(x = temp.rpkm.mat,trace = 'none',main =title.str,cexRow = .3,cexCol = .3,las=2,margins = c(11,11),labCol = lab.rows,key.xlab = '',key.ylab = '',key.title = '',col = greenred(7))
    mat.samples.count=dim(temp.rpkm.mat)[2]
    if(plot.var.genes){
      if(mat.samples.count>5){
        var.cut.off=as.numeric(quantile(as.numeric(apply(temp.rpkm.df,1,var)))[2])
        temp.rpkm.df=filter.non.variable.rows(df =temp.rpkm.df,cut.off = var.cut.off )
        temp.rpkm.mat=as.matrix(log2(temp.rpkm.df+1))
        title.str=paste(title.str, '(var. genes)',sep=' ')
        gene.cor.mat=cor(t(temp.rpkm.mat),method = 'pearson')
        top.var.genes =rownames(temp.rpkm.mat)
        var.genes.list[[development.stage]]=top.var.genes
        heatmap.2(x = temp.rpkm.mat,trace = 'none',main ='',cexRow = .3,cexCol = .8,las=2,margins = c(13,13),key.xlab = '',key.ylab = '',key.title = '',col = greenred(100))
      }
    }
  }
  out.df=rpkm.df[unique(genes.to.keep),samples.to.keep]
  out.list[['all.samples.rpkm.df']]=out.df
  out.list[['var.genes.list']]=var.genes.list
  return(out.list)
}
plot.development.stage.normalized.heatmap=function(rpkm.df,meta.df,in.method='spearman',euclidian=F){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(ordered.development.stages,development.stages)
  detected.genes.count.list=list()
  out.list=list()
  samples.to.keep=c()
  genes.to.keep=c()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    #     if(!development.stage=='ring'){
    #       next
    #     }
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    no.samples.cut.off=round(.2*dim(temp.rpkm.df)[2])
    temp.rpkm.df=filter.rpkm.less.than.cutoff(df =temp.rpkm.df,rpkm.cutoff = 0,no.samples =  no.samples.cut.off)
    temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
    temp.meta.df$no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    temp.meta.df=temp.meta.df[which(temp.meta.df$no.of.detected.genes>=200),]
    no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    temp.rpkm.df=temp.rpkm.df[,rownames(temp.meta.df)]
    temp.norm.rpkm.df=t(apply(temp.rpkm.df, 1, function(x)(x-min(x))/(max(x)-min(x))))
    out.list[[development.stage]]=list(rpkm.norm=temp.norm.rpkm.df,meta=temp.meta.df)
    samples.to.keep=append(samples.to.keep,rownames(temp.meta.df),length(samples.to.keep))
    genes.to.keep=append(genes.to.keep,rownames(temp.rpkm.df),length(genes.to.keep))
    lab.rows = as.character(temp.meta.df[colnames(temp.rpkm.df),'no.of.detected.genes'])
    dist.mat=matrix()
    if(euclidian){
      euclidian.dist.mat=as.matrix(dist(x = t(log2(temp.rpkm.df+1))))
      dist.mat= euclidian.dist.mat
      show('euclidian')
    }
    else{
      if(in.method=='spearman'){
        dist.mat=cor(x = temp.rpkm.df,method = in.method)
        show('spearman')
      }
      else{
        temp.rpkm.df=log2(temp.rpkm.df+1)
        show('pearson')
        dist.mat=cor(x = temp.rpkm.df,method = in.method)
      }
    }
    title.str=sub('\\.',' ',x = development.stage)
    dist.mat=round(dist.mat,digits = 2)
    temp.rpkm.mat=as.matrix(log2(temp.rpkm.df+1))
    heatmap.2(x = dist.mat,notecex=.1,cellnote = dist.mat,trace = 'none',main =title.str,cexRow = .3,cexCol = .3,las=2,margins = c(11,11),labCol = lab.rows,key.xlab = '',key.ylab = '',key.title = '',col = greenred(7))
    #heatmap.2(x = temp.rpkm.mat,trace = 'none',main =title.str,cexRow = .3,cexCol = .3,las=2,margins = c(11,11),labCol = lab.rows,key.xlab = '',key.ylab = '',key.title = '',col = greenred(7))
    gene.spearman.cor.dist=dist(1-cor(t(temp.rpkm.mat),method = 'spearman'))
    gene.spearman.cor.dist.dendrogram=as.dendrogram(hclust(gene.spearman.cor.dist))
    gene.spearman.cor.dist=dist(1-cor(t(temp.rpkm.mat),method = 'spearman'))
    gene.pearson.cor.dist.dendrogram=as.dendrogram(hclust(as.dist(1-cor(t(temp.rpkm.mat),method = 'pearson'))))
    heatmap.2(x = temp.norm.rpkm.df,trace = 'none',main =title.str,cexRow = .3,cexCol = .3,las=2,margins = c(11,11),labCol = lab.rows,key.xlab = '',key.ylab = '',key.title = '',col = greenred(7),Rowv =  gene.spearman.cor.dist.dendrogram)
    #heatmap.2(x = temp.rpkm.mat,trace = 'none',main =title.str,cexRow = .3,cexCol = .3,las=2,margins = c(11,11),labCol = lab.rows,key.xlab = '',key.ylab = '',key.title = '',col = greenred(7),Rowv = gene.pearson.cor.dist.dendrogram)
  }
  out.df=rpkm.df[unique(genes.to.keep),samples.to.keep]
  out.list[['all.samples.rpkm.df']]=out.df
  return(out.list)
}
get.development.stage.variable.genes.brenneck=function(counts.df,meta.df,no.samples.cut.off,mean.fit.genes.cut.off=.3,mean.exp.quantile.cut.off=.5,plot=F){
  meta.df=meta.df[colnames(counts.df),]
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(ordered.development.stages,development.stages)
  var.genes.list=list()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.counts.df=counts.df[,rownames(temp.meta.df)]
    if(is.null(dim(temp.counts.df))){
      next
    }
    if(dim(temp.counts.df)[2]<10){
      next
    }
    title.str=convert.to.title.case(gsub(pattern = '\\.',' ',development.stage))
    temp.counts.df=filter.none.expressed.genes(input.data = temp.counts.df)
    deseq.sf=estimateSizeFactorsForMatrix(counts = temp.counts.df)
    temp.rpkm.df=t(t(temp.counts.df)/deseq.sf)
    temp.meta.df$no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    temp.rpkm.df=temp.rpkm.df[,rownames(temp.meta.df)]
    #all.genes.pca.df=temp.rpkm.df[variable.genes,]
    all.genes.pca.obj=prcomp(x = t(temp.rpkm.df),retx = T,center = T,scale. = T)
    all.genes.pca.obj.score.df=data.frame(all.genes.pca.obj$x)
    all.genes.pca.ggplot.scatter=ggplot(data = all.genes.pca.obj.score.df,aes(x = PC1,y=PC2))+geom_point()+ggtitle(label = title.str)
    all.genes.pca.ggplot.scatter+xlab('Log(mean expression)')
    all.genes.pca.ggplot.scatter+ylab('Sq. coeff. variation')
    #print( all.genes.pca.ggplot.scatter)
    genes=rownames(temp.rpkm.df)
    temp.rpkm.mat=as.matrix(temp.rpkm.df)
    means.rpkm <- rowMeans(temp.rpkm.mat )
    vars.rpkm <- rowVars( temp.rpkm.mat )
    cv2.rpkm <- vars.rpkm/means.rpkm^2
    fraction.detected.in=as.numeric(apply(temp.rpkm.mat,1,nnzero))/dim(temp.rpkm.mat)[2]
    norm.cv2.rpkm =cv2.rpkm*fraction.detected.in
    col.points <- "#00207040"
    col.fact=ifelse(fraction.detected.in<.25,'Fraction < .75','Fraction >= .75')
    temp.df=data.frame(log.mean.rpkm=log10(means.rpkm),cv=cv2.rpkm, norm.cv= norm.cv2.rpkm,doubled.fraction=2*fraction.detected.in,fraction.category=col.fact)
    rownames(temp.df)=rownames(temp.rpkm.df)
    cv.quantiles.cut.off=as.numeric(quantile(x = norm.cv2.rpkm,.90))
    variable.genes=rownames(subset(temp.df,doubled.fraction>1.0 & norm.cv >=cv.quantiles.cut.off))
    shape.factor=ifelse(genes %in% variable.genes, 'high_cv','low_cv')
    temp.df$cv.category=shape.factor
    ggplot.scatter=ggplot(data = temp.df,aes(x = log.mean.rpkm,y=cv,colour=cv.category))+geom_point(aes(size = doubled.fraction,shape=cv.category))+ggtitle(label = title.str)
    ggplot.scatter+xlab('Log(mean expression)')
    ggplot.scatter+ylab('Sq. coeff. variation')
    #print(ggplot.scatter)
    pca.df=temp.rpkm.df[variable.genes,]
    pca.obj=prcomp(x = t(pca.df),retx = T,center = T,scale. = T)
    pca.obj.score.df=data.frame(pca.obj$x)
    pca.ggplot.scatter=ggplot(data = pca.obj.score.df,aes(x = PC1,y=PC2))+geom_point()+ggtitle(label = title.str)
    pca.ggplot.scatter+xlab('Log(mean expression)')
    pca.ggplot.scatter+ylab('Sq. coeff. variation')
    print( pca.ggplot.scatter)
    var.genes.list[[title.str]]=variable.genes
  }
  sig.genes.venn <- VennFromSets(setList= var.genes.list)
  out.venn=sig.genes.venn
  plot(sig.genes.venn,doWeights = F,doEuler = 'doEuler')
  categories.names=names(var.genes.list)
  for(m in 1:length(categories.names)){
    categories.name=categories.names[m]
  }
  return(var.genes.list)
}
get.variable.genes.brenneck.from.normalized.expr=function(rpkm.df,counts.df,
                                                          mean.fit.genes.cut.off=.3,
                                                          mean.exp.quantile.cut.off=.95,title.str='Test'){
  counts.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = counts.df))
  samples=colnames(counts.df)
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(rpkm.df[,samples]))
  deseq.sf=estimateSizeFactorsForMatrix(counts = counts.df)
  temp.counts.mat=as.matrix(counts.df)
  #temp.rpkm.df=t(t(counts.df)/deseq.sf)
  temp.rpkm.df=rpkm.df
  col.points <- "#00207040"
  temp.rpkm.mat=as.matrix(temp.rpkm.df)
  means.rpkm <- rowMeans(temp.rpkm.mat )
  vars.rpkm <- rowVars( temp.rpkm.mat )
  cv2.rpkm <- vars.rpkm / means.rpkm^2
  minMeanForFit <- unname(quantile(means.rpkm[ which( cv2.rpkm > mean.fit.genes.cut.off ) ],  mean.exp.quantile.cut.off) )
  useForFit <- means.rpkm >= minMeanForFit
  fit <- glmgam.fit(cbind( a0 = 1, a1tilde = 1/means.rpkm[useForFit] ),cv2.rpkm[useForFit] )
  xi <- mean(1 / deseq.sf )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"] - xi )
  col.points <- "#00207040"
  # Prepare the plot (scales, grid, labels, etc.)
  psia1theta <- mean( 1 / deseq.sf ) + a1 * mean( deseq.sf )
  minBiolDisp <- .5^2
  m <- ncol(temp.counts.mat)
  cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
  testDenom <- ( means.rpkm * psia1theta + means.rpkm^2 * cv2th ) / ( 1 + cv2th/m )
  p <- 1 - pchisq( vars.rpkm * (m-1) / testDenom, m-1 )
  padj <- p.adjust( p, "BH" )
  sig <- padj < .1
  sig[is.na(sig)] <- FALSE
  genes.test.res=data.frame(genes=names(sig),pval=as.character(sig))
  rownames(genes.test.res)=names(sig)
  variable.genes.vec=rownames(subset(genes.test.res,pval==T))
  # Prepare the plot (scales, grid, labels, etc.)
  plot( NULL, xaxt="n", yaxt="n",log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 8 ),xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)",main=title.str)
  axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",expression(10^4), expression(10^5) ) )
  axis( 2, 10^(-2:1), c( "0.01", "0.1", "1", "10" ), las=2 )
  abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
  # Add the data points
  points( means.rpkm, cv2.rpkm, pch=20, cex=1.0, col=ifelse(padj < .1, "#C0007090",col.points ))
  # Plot the fitted curve
  xg <- 10^seq( -2, 6, length.out=1000 )
  lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
  # Plot quantile lines around the fit
  df <- ncol(temp.counts.mat) - 1
  lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .975, df ) / df,
         col="#FF000080", lwd=2, lty="dashed" )
  lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .025, df ) / df,
         col="#FF000080", lwd=2, lty="dashed" )
  return(variable.genes.vec)
}
get.variable.genes.brenneck.from.counts.expr=function(counts.df,
                                                      mean.fit.genes.cut.off=.3,
                                                      mean.exp.quantile.cut.off=.5,
                                                      title.str='Test'){
  counts.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = counts.df))
  samples=colnames(counts.df)
  deseq.sf=estimateSizeFactorsForMatrix(counts = counts.df)
  temp.counts.mat=as.matrix(counts.df)
  temp.rpkm.df=t(t(counts.df)/deseq.sf)
  #temp.rpkm.df=rpkm.df
  col.points <- "#00207040"
  temp.rpkm.mat=as.matrix(temp.rpkm.df)
  means.rpkm <- rowMeans(temp.rpkm.mat )
  vars.rpkm <- rowVars( temp.rpkm.mat )
  cv2.rpkm <- vars.rpkm / means.rpkm^2
  minMeanForFit <- unname(quantile(means.rpkm[ which( cv2.rpkm > mean.fit.genes.cut.off ) ],  mean.exp.quantile.cut.off) )
  useForFit <- means.rpkm >= minMeanForFit
  fit <- glmgam.fit(cbind( a0 = 1, a1tilde = 1/means.rpkm[useForFit] ),cv2.rpkm[useForFit] )
  xi <- mean(1 / deseq.sf )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"] - xi )
  col.points <- "#00207040"
  # Prepare the plot (scales, grid, labels, etc.)
  psia1theta <- mean( 1 / deseq.sf ) + a1 * mean( deseq.sf )
  minBiolDisp <- .5^2
  m <- ncol(temp.counts.mat)
  cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
  testDenom <- ( means.rpkm * psia1theta + means.rpkm^2 * cv2th ) / ( 1 + cv2th/m )
  p <- 1 - pchisq( vars.rpkm * (m-1) / testDenom, m-1 )
  padj <- p.adjust( p, "BH" )
  sig <- padj < .1
  sig[is.na(sig)] <- FALSE
  genes.test.res=data.frame(genes=names(sig),pval=as.character(sig))
  rownames(genes.test.res)=names(sig)
  variable.genes.vec=rownames(subset(genes.test.res,pval==T))
  # Prepare the plot (scales, grid, labels, etc.)
  plot( NULL, xaxt="n", yaxt="n",log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 8 ),xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)",main=title.str)
  axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",expression(10^4), expression(10^5) ) )
  axis( 2, 10^(-2:1), c( "0.01", "0.1", "1", "10" ), las=2 )
  abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
  # Add the data points
  points( means.rpkm, cv2.rpkm, pch=20, cex=1.0, col=ifelse(padj < .1, "#C0007090",col.points ))
  # Plot the fitted curve
  xg <- 10^seq( -2, 6, length.out=1000 )
  lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
  # Plot quantile lines around the fit
  df <- ncol(temp.counts.mat) - 1
  lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .975, df ) / df,
         col="#FF000080", lwd=2, lty="dashed" )
  lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .025, df ) / df,
         col="#FF000080", lwd=2, lty="dashed" )
  return(variable.genes.vec)
}
get.top.cluster.params=function(counts.df,meta.df){
  meta.df=meta.df[colnames(counts.df),]
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(ordered.development.stages,development.stages)
  out.list=list()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.counts.df=counts.df[,rownames(temp.meta.df)]
    if(is.null(dim(temp.counts.df))){
      next
    }
    if(dim(temp.counts.df)[2]<10){
      next
    }
    title.str=convert.to.title.case(gsub(pattern = '\\.',' ',development.stage))
    col.points <- "#00207040"
    sample.counts=dim(temp.counts.df)[2]
    temp.counts.df=filter.none.expressed.genes(input.data = temp.counts.df)
    deseq.sf=estimateSizeFactorsForMatrix(counts = temp.counts.df)
    temp.rpkm.mat=t(t(temp.counts.df)/deseq.sf)
    temp.rpkm.df=data.frame(temp.rpkm.mat)
    temp.meta.df=temp.meta.df[colnames(temp.rpkm.df),]
    #temp.rpkm.df=filter.non.variable.rows(df = temp.rpkm.df,cut.off = 1)
    #temp.rpkm.df[temp.rpkm.df==0]=0.000001
    #temp.filt.rpkm.df=filter.non.variable.rows(df = temp.rpkm.df,cut.off = 1000)
    #temp.rpkm.df=temp.filt.rpkm.df
    #temp.log.rpkm.df=log10(temp.rpkm.df)
    temp.rpkm.pca=prcomp(x = t(temp.rpkm.df),retx = T,center = T,scale. = T)
    temp.rpkm.pca.scores=temp.rpkm.pca$x[,1:5]
    #temp.rpkm.pca.scores[temp.rpkm.pca.scores==0]=.000001
    #temp.rpkm.pca.log.scores=log10(temp.rpkm.pca.scores)
    temp.meta.df$no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    intern <- clValid(obj =temp.rpkm.pca.scores, nClust = 2:5, method = 'average',clMethods=c("hierarchical", "kmeans", "diana", "som", "model", "sota", "pam", "clara", "agnes"),validation="internal",metric = "euclidean")
    cluster.stability <- clValid(obj =temp.rpkm.pca.scores, nClust = 2:5, method = 'average',clMethods=c("hierarchical", "kmeans", "diana", "som", "model", "sota", "pam", "clara", "agnes"),validation ="stability",metric = "euclidean")
    show(optimalScores(object = intern))
    show(optimalScores(object = cluster.stability))
    plot(intern)
    plot(cluster.stability)
    #plot(intern, measure=c("APN","AD","ADM"),legend=T)
    internal.clusters=clusters(intern)
    plot(as.dendrogram(internal.clusters, use.modularity = T),cex=.4,main=title.str)
    rect.hclust(internal.clusters, k=4, border="lightblue")
    groups <- cutree(internal.clusters, k = 4)
    groups.df=data.frame(samples=names(groups),group=as.numeric(groups),row.names = names(groups))
    out.list[[development.stage]]=groups.df
    #plot(nClusters(intern),measures(intern,"Dunn")[,,1],type="n",axes=F, xlab="",ylab="")
    #legend("center", clusterMethods(intern), col=1:6, lty=1:6, pch=paste(1:6))
    #POSSIBLE MODIFICATIONS
    #temp.filt.rpkm.df=filter.non.variable.rows(df = temp.rpkm.df,cut.off = 1000)
    #pseudo.temp.rpkm.pca=prcomp(x = t(temp.rpkm.df),retx = T,center = T)
    #pca.summary=summary(pseudo.temp.rpkm.pca)
    #temp.rpkm.pca.scores=pseudo.temp.rpkm.pca$x
    #pseudo.temp.rpkm.pca.scores=temp.rpkm.pca.scores
    #pseudo.temp.rpkm.pca.scores[pseudo.temp.rpkm.pca.scores==0]=0.000001
    #log.pseudo.temp.rpkm.pca.scores=log10(pseudo.temp.rpkm.pca.scores)
  }
  return(out.list)
}
get.top.cluster.params.from.normalized.exprs=function(rpkm.df,meta.df,title.str){
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = rpkm.df))
  temp.rpkm.mat=as.matrix(rpkm.df)
  meta.df=meta.df[colnames(rpkm.df),]
  temp.rpkm.df=data.frame(temp.rpkm.mat)
  temp.meta.df=meta.df[colnames(temp.rpkm.df),]
  #temp.rpkm.df=filter.non.variable.rows(df = temp.rpkm.df,cut.off = 1)
  #temp.rpkm.df[temp.rpkm.df==0]=0.000001
  #temp.filt.rpkm.df=filter.non.variable.rows(df = temp.rpkm.df,cut.off = 1000)
  #temp.rpkm.df=temp.filt.rpkm.df
  #temp.log.rpkm.df=log10(temp.rpkm.df)
  temp.rpkm.pca=prcomp(x = t(temp.rpkm.df),retx = T,center = T,scale. = T)
  temp.rpkm.pca.scores=temp.rpkm.pca$x[,1:5]
  #temp.rpkm.pca.scores[temp.rpkm.pca.scores==0]=.000001
  #temp.rpkm.pca.log.scores=log10(temp.rpkm.pca.scores)
  temp.meta.df$no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
  intern <- clValid(obj =temp.rpkm.pca.scores, nClust = 2:5, method = 'average',clMethods=c("hierarchical", "kmeans", "diana", "som", "model", "sota", "pam", "clara", "agnes"),validation="internal",metric = "euclidean")
  cluster.stability <- clValid(obj =temp.rpkm.pca.scores, nClust = 2:5, method = 'average',clMethods=c("hierarchical", "kmeans", "diana", "som", "model", "sota", "pam", "clara", "agnes"),validation ="stability",metric = "euclidean")
  show(optimalScores(object = intern))
  show(optimalScores(object = cluster.stability))
  plot(intern)
  plot(cluster.stability)
  #plot(intern, measure=c("APN","AD","ADM"),legend=T)
  internal.clusters=clusters(intern)
  plot(as.dendrogram(internal.clusters, use.modularity = T),cex=.4,main=title.str)
  rect.hclust(internal.clusters, k=4, border="lightblue")
  groups <- cutree(internal.clusters, k = 4)
  groups.df=data.frame(samples=names(groups),group=as.numeric(groups),row.names = names(groups))
}
get.development.stage.heterogeneity=function(counts.df,meta.df){
  meta.df=meta.df[colnames(counts.df),]
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(ordered.development.stages,development.stages)
  coeff.var.list=list()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.counts.df=counts.df[,rownames(temp.meta.df)]
    if(is.null(dim(temp.counts.df))){
      next
    }
    if(dim(temp.counts.df)[2]<6){
      next
    }
    title.str=convert.to.title.case(gsub(pattern = '\\.',' ',development.stage))
    col.points <- "#00207040"
    sample.counts=dim(temp.counts.df)[2]
    temp.counts.df=filter.none.expressed.genes(input.data = temp.counts.df)
    deseq.sf=estimateSizeFactorsForMatrix(counts = temp.counts.df)
    temp.rpkm.df=t(t(temp.counts.df)/deseq.sf)
    temp.meta.df=temp.meta.df[colnames(temp.rpkm.df),]
    temp.meta.df$no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    #temp.rpkm.df=filter.non.variable.rows(temp.rpkm.df,100)
    temp.rpkm.df=temp.rpkm.df[,rownames(temp.meta.df)]
    gene.mean.rpkm=rowMeans(temp.rpkm.df)
    gene.var.rpkm=rowVars(temp.rpkm.df)
    fraction.detected.in=as.numeric(apply(temp.rpkm.df,1,nnzero))/dim(temp.rpkm.df)[2]
    gene.cv=gene.var.rpkm/(gene.mean.rpkm^2)
    temp.rpkm.df=data.frame(temp.rpkm.df)
    coeff.var.list[[title.str]]=gene.cv
    temp.gene.param.df=data.frame(cv=as.numeric(gene.cv),Detection=as.numeric(fraction.detected.in),log.mean=log10(gene.mean.rpkm))
    rownames(temp.gene.param.df)=rownames(temp.rpkm.df)
    cv.cut.off=as.numeric(quantile(gene.cv,.25))
    #show(cv.cut.off)
    #detection.cut.off=.4
    #variable.genes.df=subset(temp.gene.param.df,cv>=cv.cut.off & Detection>=detection.cut.off)
    variable.genes.df=subset(temp.gene.param.df,cv>=cv.cut.off)
    temp.rpkm.df=data.frame(temp.rpkm.df)
    #plot.all.scatter.in.one.axis(df =temp.rpkm.df,title.str = title.str )
    #temp.rpkm.df= filter.rpkm.less.than.cutoff(df = temp.rpkm.df,rpkm.cutoff = .1,no.samples = round(sample.counts*detection.cut.off))
    no.detected.in = as.numeric(as.character(apply(temp.rpkm.df,1,nnzero)))
    temp.rpkm.df$no.detected.in=no.detected.in
    temp.rpkm.df=subset(temp.rpkm.df,no.detected.in>=2)
    temp.rpkm.df=data.frame(subset(temp.rpkm.df,select=-no.detected.in))
    #temp.rpkm.df= filter.genes.not.expressed.in.all.samples(df =temp.rpkm.df)
    #temp.rpkm.df=temp.rpkm.df[intersect(rownames(variable.genes.df),rownames(temp.rpkm.df)),]
    cv.plot.point.col=ifelse(rownames(temp.gene.param.df) %in% rownames(temp.rpkm.df) ,'Pass','Fail')
    temp.gene.param.df$Category=cv.plot.point.col
    ggplot.scatter=ggplot(data = temp.gene.param.df,aes(x = log.mean,y=cv,colour=Category))+geom_point(aes(size = Detection))+ggtitle(label = title.str)+ylab('CV^2')+xlab('Log (mean expression)')
    #print(ggplot.scatter)
    #plot(x=x.points,y=y.points,main=title.str,ylab='Log10 (CV)',xlab='Log10 (average expression)',pch=19,col=cv.plot.point.col)
    if(dim(temp.rpkm.df)[1]<5){
      next
    }
    temp.rpkm.df[temp.rpkm.df==0]=.000000001
    plot.samples.heterogeneity(temp.rpkm.df,title.str)
    no.cluster.levels=3
    spearman.hclust=hclust(as.dist(1-cor(temp.rpkm.df,method='spearman')))
    for(m in 1:no.cluster.levels){
      cluster.cut.off=m+1
      #plot(spearman.hclust,cex=.4,main=paste(title.str,cluster.cut.off,sep=':')) 
      groups <- cutree(spearman.hclust, k = cluster.cut.off)
      groups.df=data.frame(samples=names(groups),group=as.numeric(groups),row.names = names(groups))
      #rect.hclust(spearman.hclust, k=cluster.cut.off, border="lightblue")
      samples.col.list=get.col.factor(col.factor = as.character(groups.df[colnames(temp.rpkm.df),]$group))
      temp.pca.obj=prcomp(log10(t(temp.rpkm.df)),retx = T,scale. = T,center = T)
      temp.pca.score=temp.pca.obj$x[,1:5]
      #euclidian.dist=dist(x = log10(temp.rpkm.df))
      samples.col=samples.col.list$col.str
      #pairs(x =temp.pca.score,cex = 1.0, col = samples.col,pch=19)
      samples.names=gsub(pattern = 'sample_',replacement = '',rownames(temp.pca.score))
      plot.progress=function(x){
        #plot(x,t='n',main='',xlab='',ylab='')
        plot(x,t='n',main=title.str,xlab='first dimension',ylab='second dimension')
        #points(x,pch=19,col=samples.col,cex=3.0)
        points(x,pch=19,col=samples.col,cex=3)
        #show(samples.names)
        #text(x = x,labels=samples.names,cex=.8)
        legend('topright',legend=samples.col.list$legend.str,fill=samples.col.list$legend.cols,cex=.5)
      }
      #tsne.out.mat=tsne(X =temp.pca.score,epoch_callback = plot.progress,perplexity = 5)
      #tsne.out.mat=tsne(X =euclidian.dist,epoch_callback = plot.progress,perplexity = 5,epoch = 200)
      #tsne.out.mat=tsne(X =log10(t(temp.rpkm.df)),epoch_callback = plot.progress,perplexity = 4)
    }
  }
  #boxplot(coeff.var.list,names = rep('',times=length(names(coeff.var.list))))
}
get.development.stage.heterogeneity.based.on.normalized.counts=function(rpkm.df,meta.df,show_rownames =F,annotation_legend=T, fontsize_row = 2 ,show_colnames =T){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df =meta.df )
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  development.stages=intersect(ordered.development.stages,development.stages)
  out.list=list()
  coeff.var.list=list()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    if(is.null(dim(temp.rpkm.df))){
      next
    }
    temp.rpkm.df=filter.none.expressed.genes(input.data =temp.rpkm.df )
    temp.rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df =temp.rpkm.df,sample.count = 2 )
    if(dim(temp.rpkm.df)[2]<6){
      next
    }
    title.str=convert.to.title.case(gsub(pattern = '\\.',' ',development.stage))
    col.points <- "#00207040"
    sample.counts=dim(temp.rpkm.df)[2]
    temp.meta.df=temp.meta.df[colnames(temp.rpkm.df),]
    temp.meta.df$no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    variance.cut.off=as.numeric(quantile(as.numeric(apply(temp.rpkm.df,1,sd)),probs = .3))
    show(variance.cut.off)
    #variance.cut.off=1
    temp.rpkm.df=filter.non.variable.rows(temp.rpkm.df,variance.cut.off)
    show(dim(temp.rpkm.df))
    temp.rpkm.df=temp.rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=data.frame(temp.rpkm.df)
    if(dim(temp.rpkm.df)[1]<5){
      next
    }
    pseudo.value=(min(as.numeric(temp.rpkm.df[temp.rpkm.df!=0])))/2
    show(min(temp.rpkm.df))
    #temp.rpkm.df=temp.rpkm.df+pseudo.value
    temp.out.df=plot.samples.heterogeneity(temp.rpkm.df = temp.rpkm.df,meta.df=temp.meta.df,title.str,show_rownames = show_rownames,add.other.categories = F,annotation_legend = annotation_legend, fontsize_row = fontsize_row,show_colnames = show_colnames,show_bach_var = T,sample.annotation.categories.vec =c('batch','spike','timepoint','Rosetting','in.stage.sub.grp') )
    temp.out.df[['rpkm']]=temp.rpkm.df
    out.list[[development.stage]]=temp.out.df
  }
  return(out.list)
}
get.development.stage.heterogeneity.based.on.normalized.counts.batch.corrected.expr=function(rpkm.df,meta.df,hetorogeneity.out.list,show_rownames =F,annotation_legend=T, fontsize_row = 2 ,show_colnames =T){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df =meta.df )
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  development.stages=intersect(ordered.development.stages,development.stages)
  out.list=list()
  coeff.var.list=list()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    if(is.null(dim(temp.rpkm.df))){
      next
    }
    temp.rpkm.df=filter.none.expressed.genes(input.data =temp.rpkm.df )
    #temp.rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df =temp.rpkm.df,sample.count = 2 )
    if(dim(temp.rpkm.df)[2]<6){
      next
    }
    title.str=convert.to.title.case(gsub(pattern = '\\.',' ',development.stage))
    col.points <- "#00207040"
    sample.counts=dim(temp.rpkm.df)[2]
    temp.meta.df=temp.meta.df[colnames(temp.rpkm.df),]
    temp.meta.df$no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    #variance.cut.off=as.numeric(quantile(as.numeric(apply(temp.rpkm.df,1,sd)),probs = .3))
    #show(variance.cut.off)
    #variance.cut.off=1
    #temp.rpkm.df=filter.non.variable.rows(temp.rpkm.df,variance.cut.off)
    show(dim(temp.rpkm.df))
    temp.rpkm.df=temp.rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=data.frame(temp.rpkm.df)
    genes.vec=rownames(hetorogeneity.out.list[[development.stage]][['rpkm']])
    if(dim(temp.rpkm.df)[1]<5){
      next
    }
    temp.rpkm.df=temp.rpkm.df[genes.vec,]
    pseudo.value=(min(as.numeric(temp.rpkm.df[temp.rpkm.df!=0])))/2
    #temp.rpkm.df=temp.rpkm.df+pseudo.value
    show(dim(temp.rpkm.df))
    temp.out.df=plot.samples.heterogeneity(temp.rpkm.df = temp.rpkm.df,meta.df=temp.meta.df,title.str,show_rownames = show_rownames,add.other.categories = F,annotation_legend = annotation_legend, fontsize_row = fontsize_row,show_colnames = show_colnames,show_bach_var = T,sample.annotation.categories.vec =c('batch','spike','timepoint','Rosetting','in.stage.sub.grp'),log.rpkm = T )
    temp.out.df[['rpkm']]=temp.rpkm.df
    out.list[[development.stage]]=temp.out.df
  }
  return(out.list)
}
plot.samples.heterogeneity.based.on.normalized.counts=function(rpkm.df,meta.df,log.rpkm=F,title.str='Test'){
  rpkm.df=filter.none.expressed.samples(filter.non.variable.rows(rpkm.df,cut.off = 1))
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  if(log.rpkm){
    pseudo.value=min(rpkm.df[rpkm.df!=0])/2
    #pseudo.value=1
    rpkm.df=log2(rpkm.df+pseudo.value)
  }
  title.str=convert.to.title.case(title.str)
  col.factor.vec=as.character(meta.df$development.stage)
  #samples.col=get.col.factor(col.factor = col.factor.vec)
  samples.col=get.abbrev.std.stage.cols(in.stages.vec = col.factor.vec)
  #heatmap.2(x = as.matrix(rpkm.df),main=paste(title.str,' (Euclidean)',sep=''),trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = 'Log2 (expression)',key.xlab = '',key.ylab='',ColSideColors = samples.col,labRow = F)
  col.annotation.df=subset(meta.df, select=c(development.stage,batch,spike,Rosetting))
  annotation.col.list=list(development.stage=abbrev.std.stages.col.vec)
  #pheatmap(mat = as.matrix(rpkm.df),main=paste(title.str,' (Euclidean)',sep=''),color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),clustering_distance_cols = 'euclidean',clustering_distance_rows = 'correlation',treeheight_row = 1,annotation_col = col.annotation.df,cellwidth = 2,cellheight = .07,annotation_legend = F,show_rownames = F,show_colnames = F)
  pheatmap(mat = as.matrix(rpkm.df),color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),clustering_distance_cols = 'euclidean',clustering_distance_rows = 'correlation',treeheight_row = 1,annotation_col = col.annotation.df,cellwidth = 4,cellheight = .07,annotation_legend = F,show_rownames = F,show_colnames = F,annotation_colors = annotation.col.list)
  #heatmap.2(x = cor(rpkm.df,method='spearman'),main=paste(title.str,' (Spearman)',sep=''),trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = 'Correlation',key.xlab = '',key.ylab='',ColSideColors = samples.col)
  pheatmap(mat = cor(rpkm.df,method='spearman'),color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),treeheight_row = 1,annotation_col = col.annotation.df,cellwidth = 3.5,cellheight = 3.5,annotation_legend = F,show_rownames = F,show_colnames = F,annotation_colors = annotation.col.list)
  #heatmap.2(x = cor(rpkm.df,method='pearson'),main=paste(title.str,' (Pearson)',sep=''),trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = 'Correlation',key.xlab = '',key.ylab='',ColSideColors = samples.col)
  pheatmap(mat = cor(rpkm.df,method='pearson'),color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),treeheight_row = 1,annotation_col = col.annotation.df,cellwidth = 3.5,cellheight = 3.5,annotation_legend = F,show_rownames = F,show_colnames = F,annotation_colors = annotation.col.list)
}
plot.pf.200.10k.samples.heterogeneity.based.on.normalized.counts=function(rpkm.df,meta.df,log.rpkm=F,title.str='Test'){
  rpkm.df=filter.none.expressed.samples(filter.non.variable.rows(rpkm.df,cut.off = 1))
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  if(log.rpkm){
    pseudo.value=min(rpkm.df[rpkm.df!=0])/2
    #pseudo.value=1
    rpkm.df=log2(rpkm.df+pseudo.value)
  }
  title.str=convert.to.title.case(title.str)
  col.factor.vec=as.character(meta.df$development.stage)
  #samples.col=get.col.factor(col.factor = col.factor.vec)
  samples.col=get.abbrev.std.stage.cols(in.stages.vec = col.factor.vec)
  #heatmap.2(x = as.matrix(rpkm.df),main=paste(title.str,' (Euclidean)',sep=''),trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = 'Log2 (expression)',key.xlab = '',key.ylab='',ColSideColors = samples.col,labRow = F)
  col.annotation.df=subset(meta.df, select=c(development.stage,batch,spike,Rosetting))
  annotation.col.list=list(development.stage=abbrev.std.stages.col.vec)
  #pheatmap(mat = as.matrix(rpkm.df),main=paste(title.str,' (Euclidean)',sep=''),color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),clustering_distance_cols = 'euclidean',clustering_distance_rows = 'correlation',treeheight_row = 1,annotation_col = col.annotation.df,cellwidth = 2,cellheight = .07,annotation_legend = F,show_rownames = F,show_colnames = F)
  pheatmap(mat = as.matrix(rpkm.df),color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),clustering_distance_cols = 'euclidean',clustering_distance_rows = 'correlation',treeheight_row = 1,annotation_col = col.annotation.df,cellwidth = 2,cellheight = .07,annotation_legend = F,show_rownames = F,show_colnames = F,annotation_colors = annotation.col.list)
  pheatmap(mat = cor(rpkm.df,method='spearman'),color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),treeheight_row = 1,annotation_col = col.annotation.df,cellwidth = 2,cellheight = 2,annotation_legend = F,show_rownames = F,show_colnames = F,annotation_colors = annotation.col.list)
  pheatmap(mat = cor(rpkm.df,method='pearson'),color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),treeheight_row = 1,annotation_col = col.annotation.df,cellwidth = 2,cellheight = 2,annotation_legend = F,show_rownames = F,show_colnames = F,annotation_colors = annotation.col.list)
  plot.new()
  legend('center',legend =std.stages.col.abbrev.legend.str, fill = std.stages.col.legend.fill,cex = 1.5,box.lty = 0)
}
get.sig.genes.from.scde=function(scde.res,return.genes=F,out.folder='~/rosetting/data_files/'){
  groups=names(scde.res)
  len.groups=length(groups)
  out.list=list()
  for(m in 1:len.groups){
    temp.group=groups[m]
    temp.res=scde.res[[temp.group]]$results
    #temp.res=temp.res[order(temp.res$adj.p.value),]
    sig.temp.res=subset(temp.res,adj.p.value<=.05)
    sig.temp.res=sig.temp.res[which(abs(sig.temp.res$Z)>=2),]
    sig.temp.res=sig.temp.res[order(sig.temp.res$Z),]
    no.sig.genes=dim(sig.temp.res)[1]
    if(no.sig.genes>=1){
      if(return.genes){
        #temp.sig.genes.vec=c(rownames(head(sig.temp.res,10)),rownames(tail(sig.temp.res,10)))
        temp.sig.genes.vec=rownames(sig.temp.res)
        out.list[[temp.group]]=temp.sig.genes.vec
      }
      else{
        out.list[[temp.group]]=sig.temp.res
        out.filename=paste(out.folder,paste(gsub(pattern = '\\.',replacement = '',x = temp.group),'.csv',sep=''),sep = '')
        write.table(x = sig.temp.res,sep='\t',quote = F,col.names = NA,file=out.filename)
      }
    }
  }
  return(out.list)
}
plot.scde.sig.genes.counts=function(scde.res.list){
  sig.genes.list=get.sig.genes.from.scde(scde.res =scde.res.list ,return.genes = T)
  names.vec=as.character(unlist(lapply(strsplit(names(sig.genes.list),split = '_vs_'),function(list.item){
    ordered.sp=sub.populations.grp.order.names.vec[list.item]
    ordered.sp=as.character(paste('SP.',ordered.sp,sep = ''))
    ordered.sp=sort(ordered.sp)
    ordered.sp=paste(ordered.sp,collapse  = '_vs_')
    return(ordered.sp)
  })))
  show(names.vec)
  names(sig.genes.list)=names.vec
  sig.genes.vec=as.numeric(lapply(sig.genes.list,length))
  temp.df=data.frame(sig.genes=sig.genes.vec,sp=names.vec)
  rownames(temp.df)=names.vec
  temp.df=temp.df[order(temp.df$sp),]
  barplot.lab.vec=gsub(pattern = 'SP.',replacement = '',x = names(sig.genes.list))
  barplot.lab.vec=names(sig.genes.list)
  show(barplot.lab.vec)
  #ylim.vec=c(0,max(sig.genes.vec)+10)
  barplot2(height = as.numeric(temp.df$sig.genes),names.arg =rownames(temp.df),xlab = '',ylab='',main = '',las=2,border = NA,,col = 'green')
  barplot2(height = as.numeric(temp.df$sig.genes),names.arg = rep('',length(rownames(temp.df))),xlab = '',ylab='',main = '',las=2,border = NA,col = 'blue')
  barplot2(height = as.numeric(temp.df$sig.genes),names.arg = rownames(temp.df),xlab = '',ylab='',main = '',las=2,cex.names = .2,border = NA,col='black')
}
#Plots SCDE significantly expressed genes
plot.scde.sig.diff.expr.genes.heatmap=function(rpkm.df,meta.df,scde.res.list,markers.df){
  col.pelette=colorRampPalette(c("darkgray","gray","white","yellow","orange",'red'))(5)
  psuedo.value=min(rpkm.df[rpkm.df>0])/2
  rpkm.df=filter.none.expressed.samples(df = filter.none.expressed.genes(input.data = rpkm.df))
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  sub.pop.col.list=get.color.list.for.pheatmap(in.vec = 
                                                 as.character(meta.df$markers.cluster.groups))
  sig.genes.list=get.sig.genes.from.scde(scde.res = scde.res.list,return.genes = T)
  groups.names=names(sig.genes.list)
  len.grp.names=length(groups.names)
  for(m in 1:len.grp.names){
    temp.grp=groups.names[m]
    temp.sig.genes.vec=intersect(rownames(markers.df),
                                 sig.genes.list[[temp.grp]])
    if(length(temp.sig.genes.vec)<=1){
      next
    }
    categories.vec=as.character(strsplit(x = temp.grp,split = '_vs_')[[1]])
    #categories.vec=paste('sub.cluster.',categories.vec,sep = '')
    temp.meta.df=subset(meta.df,sex.comit %in% categories.vec)
    #temp.meta.df=subset(meta.df,markers.cluster.groups %in% categories.vec)
    temp.rpkm.df=subset(rpkm.df,select=rownames(temp.meta.df))
    temp.sig.genes.vec=append(x = temp.sig.genes.vec,values = c("PF3D7_0416500", "PF3D7_0422300" ,"PF3D7_1343000" ,
                                                                "PF3D7_1008000"),
                              after=length(temp.sig.genes.vec))
    temp.sig.genes.vec=unique(temp.sig.genes.vec)
    temp.rpkm.mat=log2(as.matrix(temp.rpkm.df[temp.sig.genes.vec,])+psuedo.value)
    title.str=convert.to.title.case(in.str = gsub(pattern = 'grp.',replacement = '',x = temp.grp))
    rownames(temp.rpkm.mat)=gsub(pattern = '^gene_',replacement = '',x = rownames(temp.rpkm.mat))
    #col.annotation.df=subset(meta.df,select=c(PCA.class,post.induction.time,All.cell.types.pca.grp))
    col.annotation.df=subset(meta.df,select=c(markers.cluster.groups,development.stage,batch))
    #col.annotation.df=subset(meta.df,select=c(spike))
    col.annotation.df=col.annotation.df[order(col.annotation.df$sex.comit),]
    #temp.rpkm.mat=temp.rpkm.mat[,rownames(col.annotation.df)]
    col.list=list(in.stage.sub.grp=get.color.rainbow.list.for.pheatmap(in.vec = meta.df$in.stage.sub.grp),markers.cluster.groups=sub.pop.col.list,development.stage=abbrev.std.stages.col.vec)
    col.gap.indices=c()
    col.names=colnames(temp.rpkm.mat)
    len.col.names=length(col.names)
    for(m in 1:len.col.names){
      if(m<len.col.names){
        first.sample=col.names[m]
        sec.sample=col.names[m+1]
        first.sample.mrfp.expr=as.character(meta.df[first.sample,'markers.cluster.groups'])
        sec.sample.mrfp.expr=as.character(meta.df[sec.sample,'markers.cluster.groups'])
        if( first.sample.mrfp.expr!=sec.sample.mrfp.expr){
          col.gap.indices=append(col.gap.indices,m,length(col.gap.indices))
        }
      }
    }
    pheatmap(mat = temp.rpkm.mat,main=title.str,
             color=col.pelette,cellwidth = 2,cellheight = 5,fontsize_row = 1,
             fontsize_col = 2,border_color = NA,annotation_col = col.annotation.df,
             annotation_colors = col.list)
    pheatmap(mat = temp.rpkm.mat,main=title.str,show_colnames = F,legend = F,
             color=col.pelette,cellwidth = 3,cellheight = 5,fontsize_row = 3,
             fontsize_col = 2,border_color = NA,annotation_col = col.annotation.df,
             annotation_colors = col.list)
  }
}
get.sig.genes.from.scde.balanced.list=function(scde.balanced.res.list){
  scde.categories.list=names(scde.balanced.res.list)
  len.scde.categories.list=length(scde.categories.list)
  out.list=list()
  for(m in 1:len.scde.categories.list){
    scde.category=scde.categories.list[m]
    temp.scde.balanced.res.list=scde.balanced.res.list[[scde.category]]
    scde.category.names=names(temp.scde.balanced.res.list)
    len.scde.category.names=length(scde.category.names)
    temp.sig.genes.list=list()
    for(n in 1:len.scde.category.names){
      scde.category.name=scde.category.names[n]
      temp.scde.res.list=temp.scde.balanced.res.list[[scde.category.name]]
      scde.names=names(temp.scde.res.list)[1]
      if(is.null(scde.names)){
        next
      }
      temp.scde.res.df=temp.scde.res.list[[scde.names]]$results
      temp.sig.genes.list[[scde.category.name]]=rownames(subset(temp.scde.res.df,p.value <=.05))
      temp.sig.genes=get.sig.genes.from.scde(scde.res = temp.scde.res.list)
    }
    out.list[[scde.category]]=temp.sig.genes.list
  }
  return(out.list)
}
plot.development.stage.subtypes=function(counts.df,meta.df,scde.res.list){
  meta.df=meta.df[colnames(counts.df),]
  sub.groups=names(scde.res.list)
  len.sub.groups=length(sub.groups)
  for (i in 1:len.sub.groups){
    sub.group=sub.groups[i]
    temp.scde.res=scde.res.list[[sub.group]]
    diff.expressed.genes=c()
    temp.scde.genes.res.df=temp.scde.res$results
    temp.scde.sig.genes.df=subset(temp.scde.genes.res.df,adj.p.value<=.05)
    no.sig.genes=dim(temp.scde.sig.genes.df)[1]
    if(no.sig.genes>=10){
      diff.expressed.genes=append(diff.expressed.genes,rownames(temp.scde.sig.genes.df),length(diff.expressed.genes))
    }
    else{
      temp.scde.sig.genes.df=subset(temp.scde.genes.res.df,p.value<=.05)
      if(sub.group=='schizont_1'){
        cut.off=as.numeric(quantile(abs(as.numeric(temp.scde.sig.genes.df$Z)),.95))
        temp.scde.sig.genes.df=subset(temp.scde.genes.res.df,abs(Z)>=cut.off)
      }
      diff.expressed.genes=append(diff.expressed.genes,rownames(temp.scde.sig.genes.df),length(diff.expressed.genes))
    }
    temp.samples=rownames(temp.scde.res$o.ifm)
    temp.counts.df=filter.none.expressed.genes(counts.df[,temp.samples])
    deseq.sf=estimateSizeFactorsForMatrix(counts = temp.counts.df)
    temp.rpkm.df=t(t(temp.counts.df)/deseq.sf)
    temp.meta.df=meta.df[colnames(temp.rpkm.df),]
    temp.meta.df$no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    temp.rpkm.df[temp.rpkm.df==0]=.000000001
    gene.ids=intersect(diff.expressed.genes,rownames(temp.rpkm.df))
    gene.cols=ifelse(as.numeric(temp.scde.genes.res.df[gene.ids,'Z'])<0,'red','green')
    temp.rpkm.df=temp.rpkm.df[gene.ids,]
    title.str=gsub(pattern = '\\.',replacement = ' ',convert.to.title.case(in.str = sub.group))
    heatmap.2(x = as.matrix(log10(temp.rpkm.df)),main=paste(title.str,' (Euclidian)',sep=''),trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = 'Log (expression)',key.xlab = '',key.ylab='',RowSideColors = gene.cols)
    #heatmap.2(x = as.matrix(log10(temp.rpkm.df)),main='',trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = '',key.xlab = '',key.ylab='',labRow = '',labCol='',RowSideColors = gene.cols)
    }
  plot.new()
  legend('center',legend=rep('',times =2),fill=c('red','green'),box.lty = 0,cex=1.5)
}
plot.all.scatter.in.one.axis=function(df,title.str='Test'){
  samples=colnames(df)
  combn.mat=combn(samples,m = 2)
  dim.combn.coln=dim(combn.mat)[2]
  #plot(x=NULL,t = 'n',xaxt='n',yaxt = 'n',main=title.str,xlim=c(-5,log10(max(df))),ylim=c(-5,log10(max(df))),xlab='Log10(normalized counts)',ylab='Log10(normalized counts)')
  all.x.axis.points=c()
  all.y.axis.points=c()
  for(n in 1:dim.combn.coln){
    temp.samples=as.character(combn.mat[,n])
    temp.df=filter.none.expressed.genes(input.data[,temp.samples])
    x.points=as.numeric(temp.df[,temp.samples[1]])
    x.points=log10(ifelse(x.points==0,.00001,x.points))
    y.points=as.numeric(temp.df[,temp.samples[2]])
    y.points=log10(ifelse(y.points==0,.00001,y.points))
    all.x.axis.points=append(all.x.axis.points,x.points,length(all.x.axis.points))
    all.y.axis.points=append(all.y.axis.points,y.points,length(all.y.axis.points))
  }
  smoothScatter(x = all.x.axis.points,y=all.y.axis.points,pch = 19,cex = 1.5,main=title.str,xlab='Log10(normalized counts)',ylab='Log10(normalized counts)')
}
plot.samples.heterogeneity=function(temp.rpkm.df,meta.df,title.str,show_rownames =F,log.rpkm=T,sample.annotation.categories.vec=c('batch','spike','timepoint','Rosetting'),add.other.categories=T,annotation_legend=T,fontsize_row = 2,show_colnames =T,show_bach_var=T){
  #colnames(temp.rpkm.df) =gsub(pattern = 'sample_','',colnames(temp.rpkm.df))
  cor.dendrogram = as.dendrogram(hclust(as.dist(1-cor(temp.rpkm.df,method='spearman'))))
  #heatmap.2(x = as.matrix(log10(temp.rpkm.df)),main=paste(title.str,' (Euclidian)',sep=''),trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = 'Log (expression)',key.xlab = '',key.ylab='')
  #heatmap.2(x = cor(log10(temp.rpkm.df),method='spearman'),main=paste(title.str,' (Spearman)',sep=''),trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = 'Correlation',key.xlab = '',key.ylab='')
  out.list=list()
  grp.df=data.frame()
  if(log.rpkm){
    temp.rpkm.df=log2(temp.rpkm.df)
  }
  cor.mat=cor(temp.rpkm.df,method = 'spearman')
  #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), cluster_rows = F,cluster_cols = F,main='',cellwidth = 40, cellheight = 40,fontsize_row = 2,fontsize_col = 2)
  title.str=toupper(title.str)
  sample.annotation.df=subset(meta.df,select = sample.annotation.categories.vec)
  if(toupper(title.str)=='R'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 3)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    grp.df=ifelse(show_bach_var,data.frame(grp.df,sample.annotation.df),grp.df)
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cutree_cols = 3, cluster_rows = T,cluster_cols = T,treeheight_row = 1,main=title.str,cellwidth = 30, cellheight = 30,fontsize_row = fontsize_row,fontsize_col = 2,show_rownames = show_rownames,show_colnames =show_colnames )
  }
  if(toupper(title.str)=='LR'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 4)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    grp.batch.df=data.frame(grp.df,sample.annotation.df)
    grp.batch.df=grp.batch.df[rownames(grp.df),]
    #grp.batch.df =data.frame(grp.batch.df,grp.df)
    final.annotation.df=if(show_bach_var) grp.batch.df else grp.df
    final.col.list=list(in.stage.sub.grp=get.color.rainbow.list.for.pheatmap(in.vec =as.character(final.annotation.df$in.stage.sub.grp)),SP=get.color.rainbow.list.for.pheatmap(in.vec =as.character(final.annotation.df$SP)))
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000), cutree_cols  =4,cluster_rows = T,cluster_cols = T,main=title.str,treeheight_row = 1,cellwidth = 20, cellheight = 20,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = grp.df,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames )
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="PRGn")))(1000), cutree_cols  =4,cluster_rows = T,cluster_cols = T,main=title.str,treeheight_row = 1,cellwidth = 20, cellheight = 20,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = grp.df,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames )
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000), cutree_cols  =4,cluster_rows = T,cluster_cols = T,main=title.str,treeheight_row = 1,cellwidth = 20, cellheight = 20,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = final.annotation.df,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames,annotation_colors = final.col.list )
  }
  if(toupper(title.str)=='ET'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 3)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    grp.batch.df=data.frame(grp.df,sample.annotation.df)
    grp.batch.df=grp.batch.df[rownames(grp.df),]
    #grp.batch.df =data.frame(grp.batch.df,grp.df)
    final.annotation.df=if(show_bach_var) grp.batch.df else grp.df
    final.col.list=list(in.stage.sub.grp=get.color.rainbow.list.for.pheatmap(in.vec =as.character(final.annotation.df$in.stage.sub.grp)),SP=get.color.rainbow.list.for.pheatmap(in.vec =as.character(final.annotation.df$SP)))
    #grp.df=data.frame(grp.df,sample.annotation.df)
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000), cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 20,treeheight_row = 1, cellheight = 20,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = grp.df,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames )
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 4, name ="BrBG")))(4), cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 20,treeheight_row = 1, cellheight = 20,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = grp.df,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames )
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="BrBG")))(1000), cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 20,treeheight_row = 1, cellheight = 20,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = final.annotation.df,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames,annotation_colors = final.col.list )
  }
  if(toupper(title.str)=='T'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 5)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    grp.batch.df=data.frame(grp.df,sample.annotation.df)
    grp.batch.df=grp.batch.df[rownames(grp.df),]
    #grp.batch.df =data.frame(grp.batch.df,grp.df)
    final.annotation.df=if(show_bach_var) grp.batch.df else grp.df
    final.annotation.df=subset(final.annotation.df,select=c(batch,SP))
    final.col.list=list(in.stage.sub.grp=get.color.rainbow.list.for.pheatmap(in.vec =as.character(final.annotation.df$in.stage.sub.grp)),SP=get.color.rainbow.list.for.pheatmap(in.vec =as.character(final.annotation.df$SP)))
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cutree_cols  = 5, cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 7, cellheight = 7,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = grp.df,treeheight_row = 1,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames )
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 3, name ="BrBG")))(3),cutree_cols  = 5, cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 7, cellheight = 7,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = grp.df,treeheight_row = 1,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames )
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(1000),cutree_cols  = 5, cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 3, cellheight = 3,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = final.annotation.df,treeheight_row = 1,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames,annotation_colors = final.col.list )
  }
  if(toupper(title.str)=='ES'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 3)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    grp.batch.df=data.frame(grp.df,sample.annotation.df)
    grp.batch.df=grp.batch.df[rownames(grp.df),]
    #grp.batch.df =data.frame(grp.batch.df,grp.df)
    final.annotation.df=if(show_bach_var) grp.batch.df else grp.df
    final.col.list=list(in.stage.sub.grp=get.color.rainbow.list.for.pheatmap(in.vec =as.character(final.annotation.df$in.stage.sub.grp)),SP=get.color.rainbow.list.for.pheatmap(in.vec =as.character(final.annotation.df$SP)))
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cutree_cols  = 3, cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 30, cellheight = 30,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = grp.df,treeheight_row = 1,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames )
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="BrBG")))(10),cutree_cols  = 3, cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 30, cellheight = 30,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = grp.df,treeheight_row = 1,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames )
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cutree_cols  = 3, cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 30, cellheight = 30,fontsize_row = fontsize_row,fontsize_col = 2,annotation_col = final.annotation.df,treeheight_row = 1,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames,annotation_colors = final.col.list )
  }
  if(toupper(title.str)=='S'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 3)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    grp.batch.df=data.frame(grp.df,sample.annotation.df)
    grp.batch.df=grp.batch.df[rownames(grp.df),]
    #grp.batch.df =data.frame(grp.batch.df,grp.df)
    final.annotation.df=if(show_bach_var) grp.batch.df else grp.df
    final.col.list=list(in.stage.sub.grp=get.color.rainbow.list.for.pheatmap(in.vec =as.character(final.annotation.df$in.stage.sub.grp)),SP=get.color.rainbow.list.for.pheatmap(in.vec =as.character(final.annotation.df$SP)))
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000), cluster_rows = T,cluster_cols = T,main='',cellwidth = 10, cellheight = 10,fontsize_row = 2,fontsize_col = 2,annotation_legend = F,cutree_rows = 3,annotation_col  = grp.df)
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),treeheight_row = 1, cluster_rows = T,main=title.str,cellwidth = 7, cellheight = 7,cluster_cols = 3,fontsize_row = fontsize_row,fontsize_col = 2,cutree_cols  = 6,annotation_col = grp.df,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames )
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 5, name ="BrBG")))(3),treeheight_row = 1, cluster_rows = T,main=title.str,cellwidth = 7, cellheight = 7,cluster_cols = 3,fontsize_row = fontsize_row,fontsize_col = 2,cutree_cols  = 3,annotation_col = grp.df,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames )
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),treeheight_row = 1, cluster_rows = T,main=title.str,cellwidth = 7, cellheight = 7,cluster_cols = 3,fontsize_row = fontsize_row,fontsize_col = 2,cutree_cols  = 3,annotation_col = final.annotation.df,show_rownames = show_rownames,border_color = NA,annotation_legend = annotation_legend,show_colnames =show_colnames,annotation_colors = final.col.list )
  }
  out.list=list(rpkm.mat=cor.mat,grp.df=grp.df)
  #pheatmap.2(x = cor(log10(temp.rpkm.df),method='pearson'),main=paste(title.str,' (Pearson)',sep=''),trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = 'Correlation',key.xlab = '',key.ylab='')
  return(out.list)
}
plot.all.samples.heterogeneity=function(temp.rpkm.df,title.str){
  #colnames(temp.rpkm.df) =gsub(pattern = 'sample_','',colnames(temp.rpkm.df))
  cor.dendrogram = as.dendrogram(hclust(as.dist(1-cor(temp.rpkm.df,method='spearman'))))
  #heatmap.2(x = as.matrix(log10(temp.rpkm.df)),main=paste(title.str,' (Euclidian)',sep=''),trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = 'Log (expression)',key.xlab = '',key.ylab='')
  #heatmap.2(x = cor(log10(temp.rpkm.df),method='spearman'),main=paste(title.str,' (Spearman)',sep=''),trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = 'Correlation',key.xlab = '',key.ylab='')
  cor.mat=cor(log10(temp.rpkm.df),method = 'spearman')
  #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), cluster_rows = F,cluster_cols = F,main='',cellwidth = 40, cellheight = 40,fontsize_row = 2,fontsize_col = 2)
  title.str=toupper(title.str)
  if(toupper(title.str)=='R'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 3)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cutree_cols = 3, cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 30, cellheight = 30,fontsize_row = 2,fontsize_col = 2,show_colnames =show_colnames )
  }
  if(toupper(title.str)=='LR'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 3)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000), cutree_cols  =3,cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 20, cellheight = 20,fontsize_row = 2,fontsize_col = 2,annotation_col = grp.df,show_rownames = F,show_colnames = F)
  }
  if(toupper(title.str)=='ET'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 3)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000), cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 20, cellheight = 20,fontsize_row = 2,fontsize_col = 2,annotation_col = grp.df,show_rownames = F,show_colnames = F)
  }
  if(toupper(title.str)=='T'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 5)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cutree_cols  = 5, cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 7, cellheight = 7,fontsize_row = 2,fontsize_col = 2,annotation_col = grp.df,treeheight_row = 1,show_rownames = F,show_colnames = F)
  }
  if(toupper(title.str)=='ES'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 3)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cutree_cols  = 3, cluster_rows = T,cluster_cols = T,main=title.str,cellwidth = 30, cellheight = 30,fontsize_row = 2,fontsize_col = 2,annotation_col = grp.df,treeheight_row = 1,show_rownames = F,show_colnames = F)
  }
  if(toupper(title.str)=='S'){
    grp.hclust=hclustfunc(x = distfunc(x = cor.mat))
    grp.tree=cutree(grp.hclust,k = 3)
    grp.df=data.frame(row.names = names(grp.tree),SP=as.character(grp.tree))
    #pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000), cluster_rows = T,cluster_cols = T,main='',cellwidth = 10, cellheight = 10,fontsize_row = 2,fontsize_col = 2,annotation_legend = F,cutree_rows = 3,annotation_col  = grp.df)
    pheatmap(mat =cor.mat,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000), cluster_rows = T,main=title.str,cellwidth = 10, cellheight = 10,cluster_cols = 3,fontsize_row = 2,fontsize_col = 2,annotation_legend = T,cutree_cols  = 3,annotation_col = grp.df,treeheight_row = 1,show_rownames = F,show_colnames = F)
  }
  #pheatmap.2(x = cor(log10(temp.rpkm.df),method='pearson'),main=paste(title.str,' (Pearson)',sep=''),trace='none',margins = c(10,10),col=bluered(1000),cexRow = .3, cexCol = .3,key.title = 'Correlation',key.xlab = '',key.ylab='')
}
#Uses PCA loadings to pick the most variable genes
get.development.stage.variable.genes.pca=function(counts.df,meta.df,no.samples.cut.off,plot=F){
  meta.df=meta.df[colnames(counts.df),]
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(ordered.development.stages,development.stages)
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    title.str=convert.to.title.case(gsub(pattern = '\\.',' ',development.stage))
    temp.meta.df=meta.list[[development.stage]]
    temp.counts.df=counts.df[,rownames(temp.meta.df)]
    if(is.null(dim(temp.counts.df))){
      next
    }
    temp.counts.df=filter.none.expressed.genes(input.data = temp.counts.df)
    deseq.sf=estimateSizeFactorsForMatrix(counts = temp.counts.df)
    temp.rpkm.df=t(t(temp.counts.df)/deseq.sf)
    temp.meta.df$no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    temp.rpkm.df=temp.rpkm.df[,rownames(temp.meta.df)]
    col.points <- rep("#00207040",dim(temp.rpkm.df)[1])
    if(dim(temp.rpkm.df)[2]>=10){
      temp.meta.df=temp.meta.df[colnames(temp.rpkm.df),]
      temp.pca=prcomp(x = t(temp.rpkm.df),retx = T,scale. = T)
      samples.pca.summary=summary(temp.pca)$importance
      cut.off=100/dim(temp.pca$x)[1]
      #barplot2(100*samples.pca.summary[2,],ylab='%',las=2,main=title.str)
      #abline(h =cut.off,b = 1,col='lightgray' )
      sig.loading.genes.list=select.gene.with.sig.component.weights(pca.obj = temp.pca)
      sig.genes=unique(as.character(unlist(sig.loading.genes.list$pc.loadings.list$PC1)))
      temp.genes=rownames(temp.rpkm.df)
      col.points=ifelse(temp.genes %in% sig.genes,'blue',col.points)
      temp.rpkm.mat=as.matrix(temp.rpkm.df)
      means.rpkm <- rowMeans(temp.rpkm.mat )
      vars.rpkm <- rowVars( temp.rpkm.mat )
      cv2.rpkm <- vars.rpkm / means.rpkm^2
      # Prepare the plot (scales, grid, labels, etc.)
      #pairs(x = temp.pca$x[,1:5],cex = 1.0,pch=19)
      no.samples=dim(temp.rpkm.df)[2]
      fraction.detected.in=as.numeric(apply(temp.rpkm.df,1,nnzero))/no.samples
      temp.df=data.frame(log.mean=log10(means.rpkm),log.coeff.var=log10(cv2.rpkm),fraction.detected.in=fraction.detected.in,cex.size=2*fraction.detected.in)
      #plot(log10(means.rpkm), cv2.rpkm, pch=19, cex=.5, col=col.points,xlab = "Log(mean normalized read count)", ylab = "squared coefficient of variation (CV^2)",main=title.str)
      ggplot.scatter=ggplot(data = temp.df,aes(x = log.mean,y=log.coeff.var))+geom_point(aes(size = cex.size))+ggtitle(label = title.str)
      ggplot.scatter+xlab('Log(mean expression)')
      ggplot.scatter+ylab('Sq. coeff. variation')
      print(ggplot.scatter)
      create.density.distribution.plot.df(df =temp.rpkm.df,meta.df = temp.meta.df,sample.dist = T,title.str =  title.str,log = T)
      #plot.samples.pca.tsne.clustering(rpkm.df =temp.rpkm.df,meta.df = temp.meta.df,title.str = title.str, log.rpkms = T)
    }
  }
}
convert.to.title.case = function(in.str){
  out.str=paste(toupper(substring(text = in.str,1,1)),tolower(substring(text = in.str,2)),sep='')
  return(out.str)
}
select.gene.with.sig.component.weights=function(pca.obj){
  out.list=list()
  samples.pca=pca.obj
  pca.importance.df=data.frame(t(summary(samples.pca)$importance))
  important.pcs.cut.off =1/dim(pca.importance.df)[1]
  important.pcs.df=pca.importance.df[which(as.character(as.numeric(pca.importance.df$'Proportion.of.Variance'))>important.pcs.cut.off),]
  no.of.filt.pcs=dim(important.pcs.df)[1]
  pca.rotation=samples.pca$rotation
  pca.rotation.cut.off=1/nrow(pca.rotation)
  pca.rotation.squared=pca.rotation^2
  #pca.rotation.squared=abs(pca.rotation)
  genes.with.high.weight.test=pca.rotation.squared>pca.rotation.cut.off
  pc.names=colnames(pca.rotation.squared)
  #pc.high.weight.genes=list(names<-pc.names)
  pc.sig.gene.loadings.list=list()
  for(m in 1:no.of.filt.pcs){
    pc=pc.names[m]
    temp.df=data.frame(genes=rownames(pca.rotation.squared),row.names=rownames(pca.rotation.squared),weight.test=genes.with.high.weight.test[,pc])
    temp.df=temp.df[which(temp.df[,'weight.test']==T),]
    pc.sig.gene.loadings.list[[pc]]=rownames(temp.df)
  }
  out.list[['pc.loadings.list']]=pc.sig.gene.loadings.list
  out.list[['no.pc.cut.off']]=no.of.filt.pcs 
  return(out.list)
}
plot.development.stage.tsne.clustering=function(rpkm.df,meta.df,log.rpkm=F){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(ordered.development.stages,development.stages)
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    if(development.stage=="late.trophozoite"){
      next
    }
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    no.samples.cut.off=round(.2*dim(temp.rpkm.df)[2])
    #temp.rpkm.df=filter.rpkm.less.than.cutoff(df =temp.rpkm.df,rpkm.cutoff = 0,no.samples =  no.samples.cut.off)
    temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
    temp.meta.df$no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    temp.meta.df=temp.meta.df[which(temp.meta.df$no.of.detected.genes>=200),]
    temp.rpkm.df=temp.rpkm.df[,rownames(temp.meta.df)]
    no.of.detected.genes=as.numeric(apply(temp.rpkm.df,2,nnzero))
    cex.points.size.fact=no.of.detected.genes/sum(no.of.detected.genes)
    multiplication.fact=2/max(cex.points.size.fact)
    cex.points.size=cex.points.size.fact*multiplication.fact
    if(log.rpkm){
      temp.rpkm.df=log2(temp.rpkm.df+1)
    }
    tsne.out.mat <- tsne(X = t(temp.rpkm.df))
    rownames(tsne.out.mat)=rownames(temp.meta.df)
    col.factor.list=get.col.factor(col.factor = as.character(temp.meta.df$development.stage))
    title.str=paste(development.stage,'t-sne plot',sep=' ')
    plot(x = tsne.out.mat[,1],y = tsne.out.mat[,2],pch=19,main=title.str,xlab='First dimension',ylab='Second dimension',cex=cex.points.size,col='blue')
  }
}
plot.samples.tsne.clustering=function(rpkm.df,meta.df,title.str='Samples t-sne plot',log.rpkms=F,parental.stage=F){
  #var.cut.off=1
  #var.cut.off=as.numeric(quantile(as.numeric(apply(rpkm.df,1,sd),.50)))
  #rpkm.df=filter.none.expressed.samples(filter.non.variable.rows(df =rpkm.df, cut.off =var.cut.off ))
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  #meta.df$development.stage=factor(meta.df$development.stage,levels = factor(c('R','LR','ET','T','ES','S')))
  samples.names=colnames(rpkm.df)
  min.rpkm=min(rpkm.df[rpkm.df!=0])
  no.detected.genes=as.numeric(apply(rpkm.df,2,nnzero))
  points.sizes.fact=no.detected.genes/sum(no.detected.genes)
  points.sizes=4/max(points.sizes.fact)*(points.sizes.fact)
  development.stages=meta.df[samples.names,'development.stage']
  if(parental.stage){
    development.stages=meta.df[samples.names,'parent.development.stage']
  }
  #development.stages=meta.df[samples.names,'batch']
  #samples.col.list=get.col.factor(col.factor = development.stages)
  samples.col=get.abbrev.std.stage.cols(in.stages.vec =development.stages)
  plot.progress=function(x){
    plot(x,t='n',main=title.str,xlab='first dimension',ylab='second dimension')
    #plot(x,t='n',main='',xlab='',ylab='')
    points(x,pch=19,col=samples.col,cex=1.5)
    #points(x,pch=19,col=samples.col,cex=1.5)
    #legend('topright',legend=samples.col.list$legend.str,fill=samples.col.list$legend.cols,cex=.8)
    legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,cex=.8)
  }
  if(log.rpkms){
    pseudo.value=min(rpkm.df[rpkm.df!=0])/2
    rpkm.df=log2(rpkm.df+pseudo.value)
  }
  tsne.out.mat=tsne(X = t(rpkm.df),perplexity = 15,
                    epoch_callback = plot.progress)
  colnames(tsne.out.mat)=c('first','second')
  rownames(tsne.out.mat)=samples.names
  return(tsne.out.mat)
}
plot.samples.tsne.clustering=function(rpkm.df,meta.df,title.str='Samples t-sne plot',
                                      log.rpkms=F,pseudo.value=1,col.factor='grp',
                                      cex.points=1.5,in.epoch=200){
  samples=intersect(colnames(rpkm.df),rownames(meta.df))
  rpkm.df=subset(rpkm.df,select=samples)
  meta.df=meta.df[colnames(rpkm.df),]
  samples.names=colnames(rpkm.df)
  col.factor.vec=as.character(meta.df[,col.factor])
  col.fact.map=get.rainbow.color.vec(col.factor.vec)
  col.vec=col.fact.map[col.factor.vec]
  plot.progress=function(x){
    plot(x,t='n',main=title.str,xlab='first dimension',ylab='second dimension',
         frame.plot =F)
    points(x,pch=19,col=col.vec,cex=cex.points)
    legend('topright',legend=names(col.fact.map),fill=col.fact.map[names(col.fact.map)],
           cex=.8,border = NA,box.lty = 0)
  }
  if(log.rpkms){
    rpkm.df=log2(rpkm.df+pseudo.value)
  }
  tsne.out.mat=tsne(X = t(rpkm.df),epoch = in.epoch,
                    epoch_callback = plot.progress)
  colnames(tsne.out.mat)=c('first','second')
  rownames(tsne.out.mat)=samples.names
  return(tsne.out.mat)
}
plot.samples.api.transcriptional.signature=function(rpkm.df,meta.df,title.str='Samples t-sne plot',log.rpkms=F,parental.stage=F){
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(rpkm.df))
  meta.df=meta.df[colnames(rpkm.df),]
  samples.names=colnames(rpkm.df)
  tf.expressed.vec=rownames(rpkm.df)
  samples.tf.signature.vec=c()
  len.samples.names=length(samples.names)
  len.tf.expressed.vec=length(tf.expressed.vec)
  sig.vec=c()
  no.tf.detected=c()
  for(m in 1:len.samples.names){
    sample=samples.names[m]
    temp.tf.signature=c()
    for(n in 1:len.tf.expressed.vec){
      tf=tf.expressed.vec[n]
      tf.xpr=as.numeric(rpkm.df[tf,sample])
      temp.tf.class=ifelse(tf.xpr>=10,1,0)
      temp.tf.signature=append(temp.tf.signature,temp.tf.class,length(temp.tf.signature))
    }
    temp.tf.sig.str=paste(temp.tf.signature,collapse =',')
    sig.vec=append(sig.vec,temp.tf.sig.str,length(sig.vec))
    tf.detected=length(temp.tf.signature[temp.tf.signature!=0])
    no.tf.detected=append(no.tf.detected,tf.detected,length(no.tf.detected))
  }
  min.rpkm=min(rpkm.df[rpkm.df!=0])
  no.detected.genes=as.numeric(apply(rpkm.df,2,nnzero))
  points.sizes.fact=no.detected.genes/sum(no.detected.genes)
  points.sizes=4/max(points.sizes.fact)*(points.sizes.fact)
  development.stages=meta.df[samples.names,'development.stage']
  if(parental.stage){
    development.stages=meta.df[samples.names,'parent.development.stage']
    development.stages=meta.df[samples.names,'parent.development.stage']
  }
  development.stages=as.character(no.tf.detected)
  samples.col.list=get.col.factor(col.factor = development.stages)
  samples.col=samples.col.list$col.str
  plot.progress=function(x){
    #plot(x,t='n',main=title.str,xlab='first dimension',ylab='second dimension')
    plot(x,t='n',main='',xlab='',ylab='')
    points(x,pch=19,col=samples.col,cex=points.sizes)
    legend('topright',legend=samples.col.list$legend.str,fill=samples.col.list$legend.cols,cex=.8)
  }
  if(log.rpkms){
    pseudo.value=min(rpkm.df[rpkm.df!=0])/2
    rpkm.df=rpkm.df+pseudo.value
  }
  #tsne.out.mat=tsne(X = t(rpkm.df),epoch_callback = plot.progress)
  #colnames(tsne.out.mat)=c('first','second')
  #rownames(tsne.out.mat)=samples.names
  #return(tsne.out.mat)
}
plot.samples.api.transcriptional.signature.heatmap=function(rpkm.df,meta.df,title.str='Samples t-sne plot',log.rpkms=F,parental.stage=F){
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(rpkm.df))
  meta.df=meta.df[colnames(rpkm.df),]
  samples.names=colnames(rpkm.df)
  tf.expressed.vec=rownames(rpkm.df)
  samples.tf.signature.vec=c()
  len.samples.names=length(samples.names)
  len.tf.expressed.vec=length(tf.expressed.vec)
  sig.vec=c()
  no.tf.detected=c()
  for(m in 1:len.samples.names){
    sample=samples.names[m]
    temp.tf.signature=c()
    for(n in 1:len.tf.expressed.vec){
      tf=tf.expressed.vec[n]
      tf.xpr=as.numeric(rpkm.df[tf,sample])
      temp.tf.class=ifelse(tf.xpr>=10,1,0)
      temp.tf.signature=append(temp.tf.signature,temp.tf.class,length(temp.tf.signature))
    }
    temp.tf.sig.str=paste(temp.tf.signature,collapse =',')
    sig.vec=append(sig.vec,temp.tf.sig.str,length(sig.vec))
    tf.detected=length(temp.tf.signature[temp.tf.signature!=0])
    no.tf.detected=append(no.tf.detected,tf.detected,length(no.tf.detected))
  }
  if(log.rpkms){
    pseudo.value=min(rpkm.df[rpkm.df!=0])/2
    rpkm.df=rpkm.df+pseudo.value
  }
}
plot.samples.pca.tsne.clustering=function(rpkm.df,meta.df,title.str='Samples t-sne plot',
                                          log.rpkms =F,parental.stage=F,var.cut.off=0,
                                          sample.count.cut.off=1){
  temp.out.list=list()
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df$development.stage=factor(meta.df$development.stage,
                                   levels = factor(c('R','LR','ET','T','ES','S')))
  no.detected.genes=as.numeric(apply(rpkm.df,2,nnzero))
  rpkm.df=filter.none.expressed.samples(filter.genes.not.expressed.in.a.number.of.samples(
    df  = rpkm.df, sample.count = sample.count.cut.off))
  var.cut.off=as.numeric(quantile(as.numeric(apply(rpkm.df,1,sd)),.3))
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(
    df = filter.non.variable.rows(rpkm.df,cut.off = var.cut.off)))
  temp.out.list[['rpkm']]=rpkm.df
  #rpkm.df=filter.genes.with.zero.variance(df = rpkm.df)
  #rpkm.df=filter.non.variable.rows(df = rpkm.df,cut.off = var.cut.off)
  meta.df=meta.df[colnames(rpkm.df),]
  samples.names=colnames(rpkm.df)
  points.sizes.fact=no.detected.genes/sum(no.detected.genes)
  points.sizes=5/max(points.sizes.fact)*(points.sizes.fact)
  development.stages=as.character(meta.df[samples.names,
                                          'development.stage'])
  #development.stages=meta.df[samples.names,'protocol']
  if(parental.stage){
    development.stages=meta.df[samples.names,'parent.development.stage']
  }
  #ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  sample.cluster.grps=as.character(meta.df[samples.names,'markers.cluster.groups'])
  #sample.cluster.grps=as.character(meta.df[samples.names,'batch'])
  #samples.col.list=get.col.factor(col.factor = sample.cluster.grps)
  samples.col=get.abbrev.std.stage.cols(in.stages.vec =development.stages)
  #samples.col.list=get.col.factor(col.factor = development.stages)
  #samples.col.list=get.color.list.for.pheatmap(in.vec =sample.cluster.grps )
  #samples.col.list=get.color.brewer.list(in.vec=sample.cluster.grps)
  #samples.col.list=get.col.factor(col.factor =sample.cluster.grps)
  #samples.col=samples.col.list[sample.cluster.grps]
  #samples.col=samples.col.list$col.str
  sc.pch.shape.vec=ifelse(as.character(meta.df$mRFP1.expr)=='Y',17,19)
  plot.progress=function(x){
    #plot(x,t='n',main='',xlab='',ylab='')
    plot(x,t='n',main=title.str,xlab='first dimension',
         ylab='second dimension',frame.plot=F)
    #show(title.str)
    #points(x,pch=19,col=samples.col,cex=points.sizes)
    #points(x,pch=19,col=samples.col,cex=3)
    points(x,pch=sc.pch.shape.vec,col=samples.col,cex=1.5)
    text(x = x,labels=samples.names,cex=.1)
    #text(x = x,labels=sample.cluster.grps,cex=.3)
    #text(x = x,labels=development.stages,cex=.3)
    #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,cex=.8,box.lty = 0,border = NA)
    #legend('topright',legend=samples.col.list$legend.str,fill=samples.col.list$legend.cols,cex=.5)
    #legend('topright',legend=names(samples.col.list),fill=as.character(samples.col.list),cex=.5)
  }
  if(log.rpkms){
    pseudo.value=min(rpkm.df[rpkm.df!=0])/2
    #pseudo.value=1
    rpkm.df=rpkm.df+pseudo.value
    rpkm.df=log2(rpkm.df)
  }
  #plot.new()
  #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,cex=.8)
  #legend('center',legend=samples.col.list$legend.str,fill=samples.col.list$legend.cols,cex=.5)
  rpkm.pca=prcomp(t(rpkm.df),center = T,scale. = T)
  samples.pca.summary=summary(rpkm.pca)$importance
  explained.var.df=data.frame(explained.var=samples.pca.summary[2,])
  mean.explained.var=1/dim(explained.var.df)[1]
  explained.var.df=subset(explained.var.df,explained.var>=mean.explained.var)
  last.pc= dim(explained.var.df)[1]
  rpkm.pca.scores=rpkm.pca$x[,1:5]
  #rpkm.pca.scores=rpkm.pca$x
  #tsne.out.mat=tsne(X =rpkm.pca.scores,initial_dims = 2,epoch = 500,max_iter = 1500,
                    #epoch_callback = plot.progress)
  set.seed(seed = 1143)
  tsne.out.mat=tsne(X =t(rpkm.df),initial_dims = 5,epoch = 100,
                    max_iter = 1000,epoch_callback = plot.progress)
  colnames(tsne.out.mat)=c('first','second')
  rownames(tsne.out.mat)=samples.names
  #plot.new()
  #legend('center',legend=rep('',length(std.stages.col.abbrev.legend.str)),fill=std.stages.col.legend.fill,cex=1.5,box.lty = 0)
  #legend('center',legend=rep('',length(names(samples.col.list))),fill=as.character(samples.col.list),cex=1.5,box.lty=0)
  #plot.new()
  #legend('center',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,cex=1.5,box.lty = 0)
  #legend('center',legend=names(samples.col.list),fill=as.character(samples.col.list),cex=.5)
  return(temp.out.list)
}

plot.tsne.object=function(tsne.out.mat,col.factor.str){
  samples.col.list=get.col.factor(col.factor = col.factor.str)
  plot(tsne.out.mat[,1], tsne.out.mat[,2],xlab='first dimension',ylab='second dimension', main='Samples pca and t-sne plot',col=samples.col.list$col.str,pch=19)
  plot.new()
  legend('topright',legend=samples.col.list$legend.str,fill=samples.col.list$legend.cols,cex=1.5)
}
#Plots high correlatiing samples and returns their meta and rpkms as list
plot.high.cor.samples.per.development.stage=function(rpkm.df,meta.df){
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  detected.genes.count.list=list()
  samples=c()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
    temp.cor.mat=cor(x = temp.rpkm.df,method = 'spearman')
    temp.cor.mat[is.na(temp.cor.mat)]=0
    cor.vec=as.numeric(as.character(melt(temp.cor.mat)$value))
    cor.vec=unique(cor.vec[cor.vec!=1])
    cor.vec=sort(cor.vec)
    cut.off = tail(cor.vec,n = 3)[1]
    temp.cor.melt.df=melt(temp.cor.mat)
    temp.cor.melt.df=temp.cor.melt.df[which(temp.cor.melt.df$value>=cut.off & temp.cor.melt.df$value!=1),]
    samples=append(x = samples,unique(c(as.character(temp.cor.melt.df$X1),as.character(temp.cor.melt.df$X2))),after = length(samples))
  }
  out.rpkm.df=rpkm.df[,samples]
  out.rpkm.df=filter.none.expressed.genes(input.data = out.rpkm.df)
  out.rpkm.df=filter.none.expressed.samples(df = out.rpkm.df)
  out.meta.df=meta.df[colnames(out.rpkm.df),]
  plot.samples.correlations.heatmap(in.rpkm.df = out.rpkm.df,title.str = 'Top cor. SCs',in.meta.df =out.meta.df,filter.non.var.rows = T)
  out.list=list(rpkm=out.rpkm.df,meta=out.meta.df)
  return(out.list)
}
plot.development.stage.pca=function(rpkm.df,meta.df){
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(ordered.development.stages,development.stages)
  detected.genes.count.list=list()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
    title.str=gsub(pattern = '\\.',' ',x = development.stage)
    plot.pca(in.rpkm.df =temp.rpkm.df,meta.data = temp.meta.df,title.str = title.str,var.cut.off = 1)
  }
}
filter.samples.with.higher.cor.within.stage=function(rpkm.df,meta.df,plot=F,pearson.cor=F){
  meta.df=meta.df[colnames(rpkm.df),]
  samples.name=rownames(meta.df)
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  col.fact.list=get.col.factor(col.factor = development.stages)
  rpkm.cor.mat=cor(x = rpkm.df,method = 'spearman')
  if(pearson.cor){
    rpkm.pearson.cor.mat=cor(x = log10(rpkm.df+1))
    rpkm.cor.mat=rpkm.pearson.cor.mat
  }
  samples.to.keep=c()
  pass.counter=c()
  per.sample.cor.list=list()
  sig.adj.p.values=c()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    fail.counter=c()
    dev.stage.pass=c()
    temp.meta.df=meta.list[[development.stage]]
    temp.samples=rownames(temp.meta.df)
    temp.p.values.vec=c()
    no.samples=length(temp.samples)
    for(s in 1:no.samples){
      temp.sample=temp.samples[s]
      temp.sample.cor=list()
      other.temp.samples=setdiff(temp.samples,temp.sample)
      within.cor=as.numeric(rpkm.cor.mat[temp.sample,other.temp.samples])
      temp.sample.cor[['within.cor']]=within.cor
      without.samples=setdiff(samples.name,temp.samples)
      without.cor=as.numeric(rpkm.cor.mat[temp.sample,without.samples])
      temp.sample.cor[['without.cor']]=without.cor
      per.sample.cor.list[[temp.sample]]=temp.sample.cor
      temp.wilcox.test.p.value=format(wilcox.test(x=within.cor,y=without.cor,paired=F,alternative = 'greater')$p.value,digits=3)
      temp.p.values.vec=append(temp.p.values.vec,temp.wilcox.test.p.value,after = length(temp.p.values.vec))
    }
    adj.p.values=p.adjust(temp.p.values.vec,method = 'BH',n = length(temp.p.values.vec))
    sig.test.vec=ifelse(adj.p.values<=.05,'pass','fail')
    sig.adj.p.values=append(sig.adj.p.values,adj.p.values[sig.test.vec=='pass'],after = length(sig.adj.p.values))
    proportion.passed=100*(length(sig.test.vec[sig.test.vec=='pass'])/no.samples)
    pass.counter=append(pass.counter,proportion.passed,after=length(pass.counter))
    samples.to.keep=append(samples.to.keep,temp.samples[sig.test.vec=='pass'],after = length(samples.to.keep))
  }
  out.rpkm.df=rpkm.df[,samples.to.keep]
  out.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = out.rpkm.df))
  out.meta.df=meta.df[colnames(out.rpkm.df),]
  out.list=list(rpkm=out.rpkm.df,meta=out.meta.df)
  if(plot){
    barplot(pass.counter,ylab = 'SCs with higher cor. within (%)',col =col.fact.list$col.str )
    len.samples.to.keep=length(samples.to.keep)
    len.cor.list=length(per.sample.cor.list)
    plot.cor.list=list()
    for(m in 1:len.samples.to.keep){
      sample.to.keep=samples.to.keep[m]
      cor.list=per.sample.cor.list[[sample.to.keep]]
      within.cor=as.numeric(as.character(cor.list$within.cor))
      without.cor=as.numeric(as.character(cor.list$without.cor))
      wilcox.test.p.value=format(sig.adj.p.values[m],digits = 3)
      temp.boxplot=boxplot(list(within= within.cor,without=without.cor),col = c('blue','red'),names = c('within','without'),main=sample.to.keep)
      label.str=paste('Pvalue: ',wilcox.test.p.value,sep='')
      text(x=1.5,y=max(temp.boxplot$stats),labels=label.str,cex=.7)
    }
  }
  return(out.list)
}
#Gets correlation from non-filtered pairwise samples
get.cor.from.df=function(df){
  cor.df=data.frame(cor(df,method='spearman'))
  return(cor.df)
}
#Gets correlation from filtered pairwise samples
get.pairwise.filt.cor.from.df=function(df){
  samples=colnames(df)
  no.samples=length(samples)
  out.mat=matrix(nrow =no.samples,ncol = no.samples)
  for(m in 1:no.samples){
    for(n in 1:no.samples){
      first.sample=samples[m]
      sec.sample=samples[n]
      temp.df=df[,c(first.sample,sec.sample)]
      colnames(temp.df)=c(first.sample,sec.sample)
      temp.df=filter.none.expressed.genes(temp.df)
      if(dim(temp.df)[1]<5){
        out.mat[m,n]=NA
        next
      }
      samples.cor=cor(x = as.numeric(as.character(temp.df[,first.sample])), y = as.numeric(as.character(temp.df[,sec.sample])),method='spearman')
      if(is.na(samples.cor)){
        out.mat[m,n]=NA
        next
      }
      out.mat[m,n]=samples.cor
    }
  }
  out.df=data.frame(out.mat)
  rownames(out.df)=samples
  colnames(out.df)=samples
  return(out.df)
}
get.cor.per.development.stage=function(rpkm.df,meta.df,filter.pairs=F){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  len.development.stages=length(development.stages)
  out.list=list()
  for (i in 1:len.development.stages){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(temp.rpkm.df))
    development.stage.samples=colnames(temp.rpkm.df)
    temp.meta.df=temp.meta.df[development.stage.samples,]
    temp.cor.df=data.frame()
    if(filter.pairs){
      temp.cor.df=get.pairwise.filt.cor.from.df(df = temp.rpkm.df)
    }
    else{
      temp.cor.df=get.cor.from.df(df = temp.rpkm.df)
    }
    temp.list=list(cor=temp.cor.df,meta=temp.meta.df)
    out.list[[development.stage]]=temp.list
  }
  return(out.list)
}
correct.batch.effect.per.development.stage=function(rpkm.df,meta.df,log.rpkm=F,par.prior=T){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  #development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  len.development.stages=length(development.stages)
  out.list=list()
  for (i in 1:len.development.stages){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    if(dim(temp.meta.df)[1]<=1){
      next
    }
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.genes(temp.rpkm.df)
    development.stage.samples=colnames(temp.rpkm.df)
    temp.meta.df=temp.meta.df[development.stage.samples,]
    pheno.df=subset(temp.meta.df,select=c(batch,spike))
    if(log.rpkm){
      pseudo.val=min(temp.rpkm.df[temp.rpkm.df>0])/2
      log.rpkm.df=log2(temp.rpkm.df+pseudo.val)
      temp.rpkm.df=log.rpkm.df
    }
    batch = pheno.df$batch
    spike=pheno.df$spike
    batch.vec=paste(batch,spike,sep = '.')
    pheno.df$batch=batch.vec
    pheno.df=subset(pheno.df,select=batch)
    no.batches=length(unique(batch.vec))
    if(no.batches>=2){
      modcombat = model.matrix(~1, data=pheno.df)
      batch.corrected.log.rpkm.df = ComBat(dat=temp.rpkm.df, batch=batch.vec, mod=modcombat, par.prior=par.prior, prior.plot=F)
      temp.rpkm.df=batch.corrected.log.rpkm.df
    }
    batch.corrected.rpkm.df=data.frame()
    if(log.rpkm){
      batch.corrected.rpkm.df=data.frame(2^temp.rpkm.df)
    }
    out.list[[development.stage]]=data.frame(batch.corrected.rpkm.df)
  }
  out.df=create.data.frames.from.df.list(df.list = out.list,counts = F)
  meta.df=meta.df[colnames(out.df),]
  final.out.list=list(meta=meta.df,rpkm=out.df)
  return(final.out.list)
}
correct.batch.effect=function(rpkm.df,meta.df,log.rpkm=F,par.prior=T){
  rpkm.df=filter.none.expressed.genes(rpkm.df)
  meta.df=meta.df[colnames(rpkm.df),]
  if(log.rpkm){
    pseudo.val=min(rpkm.df[rpkm.df>0])/2
    log.rpkm.df=log2(rpkm.df+pseudo.val)
    rpkm.df=log.rpkm.df
  }
  meta.df=meta.df[colnames(rpkm.df),]
  pheno.df=subset(meta.df,select=c(batch,spike))
  batch = pheno.df$batch
  spike=pheno.df$spike
  batch.vec=paste(batch,spike,sep = '.')
  pheno.df$batch=batch.vec
  pheno.df=subset(pheno.df,select=batch)
  no.batches=length(unique(batch.vec))
  batch.corrected.rpkm.df=data.frame()
  if(no.batches>=2){
    modcombat = model.matrix(~1, data=pheno.df)
    temp.rpkm.df= ComBat(dat=rpkm.df, batch=batch.vec, mod=modcombat, par.prior=par.prior, prior.plot=F)
    batch.corrected.rpkm.df=data.frame(temp.rpkm.df)
  }
  if(log.rpkm){
    batch.corrected.rpkm.df=2^batch.corrected.rpkm.df
  }
  meta.df=meta.df[colnames(batch.corrected.rpkm.df),]
  final.out.list=list(meta=meta.df,rpkm=batch.corrected.rpkm.df)
  return(final.out.list)
}
temp.get.cor.from.df=function(df){
  samples=colnames(df)
  samples.combn.mat=combn(x = samples,m = 2)
  apply(samples.combn.mat,2,function(sample.pair){
    first.sample=sample.pair[1]
    sec.sample=sample.pair[2]
    first.sample.rpkm=as.numeric(as.character(df[,first.sample]))
    sec.sample.rpkm=as.numeric(as.character(df[,sec.sample]))
    temp.df=data.frame(first.sample=first.sample.rpkm,sec.sample=sec.sample.rpkm)
    temp.df=filter.none.expressed.genes(temp.df)
    temp.cor.mat=cor(temp.df,method='spearman')
    #cor.score=as.numeric(c(temp.cor.mat[first.sample,sec.sample]))
    #show(cor.score)
  })
}
get.top.cor.score=function(cor.df){
  samples=colnames(cor.df)
  len.samples=length(samples)
  samples.tracker=list()
  samples.pairs=c()
  cor.list=c()
  samples.combn=combn(samples,m = 2)
  no.pairs=dim(samples.combn)[2]
  for(m in 1:no.pairs){
    sample.pair=samples.combn[,m]
    sample.one=sample.pair[1]
    sample.two=sample.pair[2]
    cor.score=cor.df[sample.one,sample.two]
    out.pair=paste(sample.one,sample.two,sep = '_and_')
    samples.pairs=append(samples.pairs,out.pair,length(samples.pairs))
    cor.list=append( cor.list,cor.score,length( cor.list))
  }
  temp.df=data.frame(pairs=samples.pairs,cor=cor.list)
  rownames(temp.df)=samples.pairs
  max.cor=max(cor.list)
  top.pairs=as.character(temp.df[which(temp.df$cor==max.cor),'pairs'])
  top.pairs=strsplit(x = top.pairs,split = '_and_')[[1]]
  out.list=list(top.pairs=top.pairs,max.cor=max.cor)
  return(out.list)
}
get.top.pair.cor.samples.per.development.stage=function(rpkm.df,meta.df){
  meta.df=meta.df[colnames(rpkm.df),]
  samples.name=rownames(meta.df)
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(ordered.development.stages,names(meta.list))
  len.development.stages=length(development.stages)
  samples.to.keep=c()
  for (i in 1:len.development.stages){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    if(dim(temp.rpkm.df)[2]<3){
      next
    }
    development.stage.samples=colnames(temp.rpkm.df)
    development.stage.samples.len=length(development.stage.samples)
    temp.rpkm.df=filter.none.expressed.genes(temp.rpkm.df)
    temp.meta.df=temp.meta.df[development.stage.samples,]
    temp.cor.df=cor(temp.rpkm.df,method='spearman')
    #temp.cor.df=cor(temp.rpkm.df,method='pearson')
    top.cor.score.list=get.top.cor.score(cor.df = temp.cor.df)
    max.cor=top.cor.score.list$max.cor
    top.cor.pair=top.cor.score.list$top.pairs
    samples.to.keep=append(samples.to.keep,top.cor.pair,length(samples.to.keep))
  }
  samples.to.keep=unique(samples.to.keep)
  out.rpkm.df=rpkm.df[,samples.to.keep]
  out.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = out.rpkm.df))
  out.meta.df=meta.df[colnames(out.rpkm.df),]
  out.list=list(rpkm=out.rpkm.df,meta=out.meta.df)
  return(out.list)
}
get.top.cor.samples.per.development.stage=function(rpkm.df,meta.df,cor.cut.off=.2,in.method='spearman'){
  meta.df=meta.df[colnames(rpkm.df),]
  samples.name=rownames(meta.df)
  cor.cut.off=as.numeric(cor.cut.off)
  meta.list=split(meta.df,f=meta.df$development.stage)
  all.samples.cor=cor(x = rpkm.df,method = in.method)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=intersect(ordered.development.stages,development.stages)
  len.development.stages=length(development.stages)
  samples.to.keep=c()
  for (i in 1:len.development.stages){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
    development.stage.samples=colnames(temp.rpkm.df)
    development.stage.samples.len=length(development.stage.samples)
    temp.meta.df=temp.meta.df[development.stage.samples,]
    #temp.cor.df=cor(temp.rpkm.df,method = in.method)
    temp.cor.df=all.samples.cor[development.stage.samples,development.stage.samples]
    top.pair.list=get.top.cor.score(cor.df =temp.cor.df )
    top.pair.vec=top.pair.list[["top.pairs"]]
    samples.to.keep=append(samples.to.keep,top.pair.vec,length(samples.to.keep))
    non.top.samples=setdiff(development.stage.samples,top.pair.vec)
    non.development.stage.samples=setdiff(samples.name,development.stage.samples)
    for(m in 1:length(non.top.samples)){
      non.top.sample=non.top.samples[m]
      non.top.cor.with.top.pair=temp.cor.df[non.top.sample,top.pair.vec]
      non.top.cor.with.other.non.tops=all.samples.cor[non.top.sample,setdiff(non.top.samples,non.top.sample)]
      if(all(non.top.cor.with.top.pair>=cor.cut.off)){
        mean.cor=mean(non.top.cor.with.other.non.tops)
        median.cor=median(non.top.cor.with.other.non.tops)
        if(median.cor>=cor.cut.off){
          samples.to.keep=append(samples.to.keep,non.top.sample,length(samples.to.keep))
        }
      }
    }
  }
  samples.to.keep=unique(samples.to.keep)
  out.rpkm.df=rpkm.df[,samples.to.keep]
  out.rpkm.df=filter.none.expressed.genes(input.data = out.rpkm.df)
  out.meta.df=meta.df[colnames(out.rpkm.df),]
  out.list=list(rpkm=out.rpkm.df,meta=out.meta.df)
  return(out.list)
}
plot.top.pair.per.stage.scatter.stage=function(rpkm.df,meta.df){
  meta.df=meta.df[colnames(rpkm.df),]
  meta.list=split(meta.df,f =meta.df$development.stage)
  development.stages = names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages = intersect(ordered.development.stages,development.stages)
  col.factor.list=get.col.factor(development.stages)
  no.development.stage=length(development.stages)
  for (i in 1:no.development.stage){
    development.stage= development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    first.sample=colnames(temp.rpkm.df)[1]
    sec.sample=colnames(temp.rpkm.df)[2]
    temp.col=col.factor.list$col.str[i]
    #plot.new()
    category=toupper(gsub(pattern = '\\.',replacement = ' ',x = development.stage))
    category=paste(toupper(substring(text = category,1,1)),tolower(substring(text = category,2)),sep='')
    #text(x=.5,y=.5,category,cex = 2.0)
    x.points=as.numeric(temp.rpkm.df[,1])
    y.points=as.numeric(temp.rpkm.df[,2])
    temp.df=data.frame(first=x.points,second=y.points)
    temp.df=filter.none.expressed.genes(input.data = temp.df)
    x.points=as.numeric(as.character(temp.df[,1]))
    y.points=as.numeric(as.character(temp.df[,2]))
    no.genes=dim(temp.df)[1]
    sample.one.detected.genes=length(x.points[x.points>0])
    sample.two.detected.genes=length(y.points[y.points>0])
    sample.one.drop.out.rate=round(100*((no.genes-sample.one.detected.genes)/no.genes),digits = 2)
    sample.two.drop.out.rate=round(100*((no.genes-sample.two.detected.genes)/no.genes),digits = 2) 
    sample.one.drop.out.rate.str=paste('[Dropout-rate: ',toString(sample.one.drop.out.rate),sep='')
    sample.one.drop.out.rate.str=paste(sample.one.drop.out.rate.str,'%]',sep='')
    sample.two.drop.out.rate.str=paste('[Dropout-rate: ',toString(sample.two.drop.out.rate),sep='')
    sample.two.drop.out.rate.str=paste(sample.two.drop.out.rate.str,'%]',sep='')
    first.sample=paste(sub(pattern = '^sample_',replacement = '',first.sample),'Normalized read counts',sep=' ')
    xlab.str=paste(first.sample,sample.one.drop.out.rate.str,sep=' ')
    sec.sample=paste(sub(pattern = '^sample_', replacement = '',sec.sample),'Normalized read counts',sep=' ')
    ylab.str=paste(sec.sample,sample.two.drop.out.rate.str,sep = ' ')
    #plot(x =log2(x.points+1),y=lo
    #title.str=convert.to.title.case()
    geneScatterplot(x =x.points,y =y.points,xlab = xlab.str,ylab =ylab.str,col = '#00207040',main.title =category  )
    #geneScatterplot(x =x.points,y =y.points,xlab = '',ylab ='',col = '#00207040',main.title =''  )
  }
}
filter.samples.based.on.cor=function(rpkm.df,meta.df,plot=T,cor.cut.off=.2){
  #Filters samples based on their cor with vs without score
  temp.out.list=filter.samples.with.higher.cor.within.stage(rpkm.df =rpkm.df,meta.df =meta.df )
  filt.rpkm.df=temp.out.list$rpkm
  filt.meta.df=temp.out.list$meta
  #Filters samples with the highest cor within stage
  out.list=get.top.cor.samples.per.development.stage(rpkm.df =filt.rpkm.df,meta.df =filt.meta.df,cor.cut.off = cor.cut.off)
  return(out.list)
}
determine.no.genes.cut.off=function(rpkm.df,meta.df,plot=T){
  meta.df=meta.df[colnames(rpkm.df),]
  samples.name=rownames(meta.df)
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages=names(meta.list)
  rpkm.cor.mat=cor(x = rpkm.df,method = 'spearman')
  rpkm.cor.mat[is.na(rpkm.cor.mat)]=0
  detected.genes.count.list=list()
  out.samples=c()
  out.samples.dev.stage=c()
  out.samples.mean.cor.within=c()
  out.samples.mean.cor.without=c()
  out.samples.cor.within.pvalue=c()
  out.samples.cor.significant=c()
  out.samples.no.detected.genes=c()
  qc.pass.samples=c()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.samples=rownames(temp.meta.df)
    dev.stage.cor.within=c()
    dev.no.detected.genes=c()
    dev.stage.cor.pass=c()
    dev.stage.col=c()
    for(s in 1:length(temp.samples)){
      temp.sample=temp.samples[s]
      out.samples=append(out.samples,values = temp.sample,after = length(out.samples))
      out.samples.dev.stage=append(out.samples.dev.stage,values =development.stage,after = length(out.samples.dev.stage) )
      other.temp.samples=setdiff(temp.samples,temp.sample)
      within.cor=as.numeric(rpkm.cor.mat[temp.sample,other.temp.samples])
      without.samples=setdiff(samples.name,temp.samples)
      without.cor=as.numeric(rpkm.cor.mat[temp.sample,without.samples])
      out.samples.mean.cor.within=append(out.samples.mean.cor.within,values = mean(within.cor),after =length(out.samples.mean.cor.within) )
      dev.stage.cor.within=append(dev.stage.cor.within,values =mean(within.cor),after = length(dev.stage.cor.within))
      out.samples.mean.cor.without=append(out.samples.mean.cor.without,values = mean(without.cor),after = length(out.samples.mean.cor.without))
      wilcox.test.p.value=format(wilcox.test(x=within.cor,y=without.cor,paired=F,alternative = 'greater')$p.value,digits=3)
      out.samples.cor.within.pvalue=append(out.samples.cor.within.pvalue,values =wilcox.test.p.value,after = length(out.samples.cor.within.pvalue) )
      if(as.numeric(wilcox.test.p.value)<=0.05){
        out.samples.cor.significant=append(out.samples.cor.significant,values = 'yes',after = length(out.samples.cor.significant))
        dev.stage.col=append(dev.stage.col,values = 'blue',after=length(dev.stage.col))
        dev.stage.cor.pass=append(dev.stage.cor.pass,values = 'sig',after = length(dev.stage.cor.pass))
      }
      else{
        out.samples.cor.significant=append(out.samples.cor.significant,values = 'no',after = length(out.samples.cor.significant))
        dev.stage.col=append(dev.stage.col,values = 'red',after=length(dev.stage.col))
        dev.stage.cor.pass=append(dev.stage.cor.pass,values = 'not.sig',after = length(dev.stage.cor.pass))
      }
      no.of.detected.genes=as.numeric(as.character(temp.meta.df[temp.sample,'no.of.detected.genes']))
      dev.no.detected.genes=append(dev.no.detected.genes,values = no.of.detected.genes,after = length(dev.no.detected.genes))
      out.samples.no.detected.genes=append(out.samples.no.detected.genes,values =no.of.detected.genes,after = length(out.samples.no.detected.genes) )
    }
    if(plot){
      plot(x=dev.stage.cor.within,y=dev.no.detected.genes,main=development.stage,xlab='Cor. within',ylab='Gene counts',col=dev.stage.col,pch=19)
      abline(h = 200,b=1)
      legend('topright',legend=c('Higher cor. within','Lower cor. within'),fill=c('blue','red'),cex=.5)
    }
    dev.stage.df=data.frame(cor.within=dev.stage.cor.within,gene.count=dev.no.detected.genes,dev.stage.pass=dev.stage.cor.pass)
    rownames(dev.stage.df)=temp.samples
    dev.stage.df=dev.stage.df[which(dev.stage.df$dev.stage.pass=='sig'),]
    cor.cut.off=round(quantile(x = as.numeric(as.character(dev.stage.df$cor.within)))[['75%']])
    dev.stage.df=dev.stage.df[which(dev.stage.df$cor.within>=cor.cut.off),]
    no.genes.cut.off=round(quantile(x = as.numeric(as.character(dev.stage.df$gene.count)))[['75%']])
    dev.stage.df=dev.stage.df[which(dev.stage.df$gene.count>=no.genes.cut.off),]
    qc.pass.samples=append(qc.pass.samples,values = rownames(dev.stage.df),after = length(qc.pass.samples))
  }
  out.df=data.frame(development.stage=out.samples.dev.stage,mean.within.cor=out.samples.mean.cor.within,mean.without.cor=out.samples.mean.cor.without,significant=out.samples.cor.significant,detected.genes=out.samples.no.detected.genes)
  rownames(out.df)=out.samples
  col.factor.list=get.col.factor(col.factor = out.samples.dev.stage)
  col.str=col.factor.list$col.str
  qc.pass.rpkm.df=rpkm.df[,qc.pass.samples]
  qc.pass.meta.df=meta.df[qc.pass.samples,]
}
create.unique.pairwise.combn=function(input.vec,no.combn=2){
  #input.vec=sample(input.vec, size = length(input.vec))
  input.combn.mat=combn(x = input.vec,m = no.combn)
  no.combn=dim(input.combn.mat)[2]
  unique.list=c()
  col.to.keep=c()
  for(i in 1:no.combn){
    temp.pair=input.combn.mat[,i]
    if(length(intersect(temp.pair,unique.list))>0){
      next
    }
    else{
      col.to.keep=append(col.to.keep,values = i,after = length(col.to.keep))
      unique.list=append(unique.list,values =as.character(temp.pair),after = length(unique.list) )
    }
  }
  out.mat=input.combn.mat[,col.to.keep]
  return(out.mat)
}
create.aggregate.cells=function(sc.df,no.scs=2){
  samples=colnames(sc.df)
  no.genes=dim(sc.df)[1]
  aggregate.samples.list=chunk(in.vec =samples,chunk.sizes = no.scs )
  aggregate.samples.list.names=names(aggregate.samples.list)
  no.aggregates=length(aggregate.samples.list)
  out.mat=matrix(nrow = no.genes,ncol = no.aggregates)
  out.colnames=c()
  for(n in 1:no.aggregates){
    aggregate.samples.list.name=aggregate.samples.list.names[n]
    temp.aggregate.samples=c(unlist(aggregate.samples.list[aggregate.samples.list.name]))
    temp.aggregate.df=sc.df[,temp.aggregate.samples]
    temp.aggregate.rpkm=as.numeric(as.character(apply(temp.aggregate.df,1,mean)))
    out.colnames=append(out.colnames,values = paste(temp.aggregate.samples,collapse = ':'),after=length(out.colnames))
    out.mat[,n]=temp.aggregate.rpkm
  }
  out.df=data.frame(out.mat)
  colnames(out.df)=out.colnames
  rownames(out.df)=rownames(sc.df)
  return(out.df)
}
create.aggregate.per.stage=function(sc.rpkm.df,sc.meta.df,no.scs=2,parent.development.stage=F){
  sc.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(sc.rpkm.df))
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.list=split(sc.meta.df, f=sc.meta.df$development.stage)
  #development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  if(parent.development.stage){
    sc.meta.df=create.parental.stage(meta.df = sc.meta.df)
    sc.meta.list=split(sc.meta.df, f=sc.meta.df$parent.development.stage)
  }
  development.stages=names(sc.meta.list)
  aggregate.rpkm.list=list()
  col.factor=c()
  aggregate.col.names=c()
  stage.names=c()
  len.development.stage=length(development.stages)
  for (i in 1:len.development.stage){
    development.stage=development.stages[i]
    temp.meta.df=sc.meta.list[[development.stage]]
    temp.rpkm.df=sc.rpkm.df[,rownames(temp.meta.df)]
    temp.aggregate.rpkm.df=create.aggregate.cells(sc.df = temp.rpkm.df,no.scs = no.scs)
    aggregate.rpkm.list[[development.stage]]=temp.aggregate.rpkm.df
    temp.aggregate.names=paste(development.stage,paste('aggregate.one',seq(1,dim(temp.aggregate.rpkm.df)[2]),sep = '.'),sep = '.')
    stage.names=append(stage.names,rep(development.stage,times =dim(temp.aggregate.rpkm.df)[2]),length(stage.names))
    aggregate.col.names=append(aggregate.col.names,temp.aggregate.names,length(aggregate.col.names))
  }
  out.rpkm.df=as.data.frame(aggregate.rpkm.list)
  rownames(out.rpkm.df)=rownames(sc.rpkm.df)
  colnames(out.rpkm.df)=aggregate.col.names
  out.meta.df=data.frame(sample.name=aggregate.col.names,development.stage=stage.names)
  rownames(out.meta.df)= aggregate.col.names
  out.list=list(rpkm=out.rpkm.df,meta=out.meta.df)
  return(out.list)
}
plot.aggregate.per.development.stage.pca=function(rpkm.df,meta.df,no.scs=2,first.pc='PC1',second.pc='PC2',title.str='Test',log.samples.pca.scores=F,log.rpkm=F){
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = rpkm.df))
  #rpkm.df=filter.none.expressed.samples(df = rpkm.df)
  rpkm.df=filter.non.variable.rows(df = rpkm.df,100)
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.list=split(meta.df,f=meta.df$development.stage)
  #development.stages=names(meta.list)
  #development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  development.stages=intersect(names(meta.list),development.stages)
  aggregate.rpkm.list=list()
  col.factor=c()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    if(dim(temp.rpkm.df)[2]<2){
      next
    }
    temp.aggregate.rpkm.df=create.aggregate.cells(sc.df = temp.rpkm.df,no.scs = no.scs)
    aggregate.rpkm.list[[development.stage]]=temp.aggregate.rpkm.df
    col.factor=append(col.factor,values = rep(x = development.stage,times =dim(temp.aggregate.rpkm.df)[2]) ,after = length(col.factor))
  }
  aggregate.rpkm.df=as.data.frame(aggregate.rpkm.list)
  rownames(aggregate.rpkm.df)=rownames(rpkm.df)
  aggregate.rpkm.df=filter.non.variable.rows(df =aggregate.rpkm.df,cut.off = 1 )
  col.factor.list=get.col.factor(col.factor = col.factor)
  if(log.rpkm){
    pseudo.values=min(aggregate.rpkm.df[aggregate.rpkm.df!=0])/2
    aggregate.rpkm.df=aggregate.rpkm.df+pseudo.values
  }
  samples.pca=prcomp(t(aggregate.rpkm.df),retx=T,center=T,scale.=T)
  samples.pca.summary=summary(samples.pca)$importance
  #barplot(100*samples.pca.summary[2,],main='Components explained variance',ylab='%')
  #abline(h=100/length(colnames(samples.pca.summary)),b=1)
  first.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,first.pc]))*100,2)
  second.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,second.pc]))*100,2)
  samples.pca.scores=samples.pca$x
  x.points=as.numeric(as.character(samples.pca.scores[,first.pc]))
  y.points=as.numeric(as.character(samples.pca.scores[,second.pc]))
  if(log.samples.pca.scores){
    x.points=log2(x.points+1)
    y.points=log2(y.points+1)
  }
  std.stages.col.list
  #pca.col.list=get.col.factor(col.factor = col.factor)
  pca.col.str=abbrev.std.stages.col.vec[col.factor]
  first.pc.explained.var=paste('(Explained var: ',paste(first.pc.explained.var,'%)',sep=''),sep='')
  second.pc.explained.var=paste('(Explained var: ',paste(second.pc.explained.var,'%)',sep=''),sep='')
  first.pc.lab=paste(first.pc,first.pc.explained.var,'')
  second.pc.lab=paste(second.pc,second.pc.explained.var,'')
  x.lim=c(min(x.points,y.points),max(x.points,y.points))
  y.lim=c(min(x.points,y.points),max(x.points,y.points))
  plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,ylab=second.pc.lab,main=title.str,pch=19,cex=2)
  legend('center',legend = names(abbrev.std.stages.col.vec),fill=as.character(abbrev.std.stages.col.vec))
  #  
#   plot.new()
#   
#   plot.new()
#   
#   legend('center',legend = rep('',times=length(pca.col.list$legend.str)),fill=pca.col.list$legend.cols,box.lty = 0,cex=1.5)
#   
}
plot.aggregate.per.development.stage.pca.tsne=function(rpkm.df,meta.df,no.scs=2,first.pc='PC1',second.pc='PC2',title.str='Test',log.rpkm=F){
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = rpkm.df))
  meta.df=meta.df[colnames(rpkm.df),]
  aggregates.list=create.aggregate.per.stage(sc.rpkm.df =rpkm.df,sc.meta.df =meta.df,no.scs =  no.scs )
  plot.samples.pca.tsne.clustering(rpkm.df = aggregates.list$rpkm,meta.df = aggregates.list$meta,title.str = title.str,log.rpkms =log.rpkm )
}
plot.aggregate.per.development.stage.tsne=function(rpkm.df,meta.df,no.scs=2,first.pc='PC1',second.pc='PC2',title.str='Test',log.rpkm=F){
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = rpkm.df))
  meta.df=meta.df[colnames(rpkm.df),]
  aggregates.list=create.aggregate.per.stage(sc.rpkm.df =rpkm.df,sc.meta.df =meta.df,no.scs =  no.scs )
  plot.samples.pca.tsne.clustering(rpkm.df = aggregates.list$rpkm,meta.df = aggregates.list$meta,title.str = title.str,log.rpkms =log.rpkm )
}
plot.aggregate.per.development.stage.tsne.at.diff.aggregate=function(rpkm.df,meta.df,title.str='Average SCs ',log.rpkm=F){
  for(m in 3:3){
    temp.title.str=paste(title.str,paste(paste('(',m,sep=''),'SCs)',sep = ''),sep='')
    plot.aggregate.per.development.stage.tsne(rpkm.df  =rpkm.df, meta.df = meta.df,no.scs = m,title.str = temp.title.str,log.rpkm = log.rpkm )
  }
}
plot.aggregate.per.development.stage.pca.at.diff.aggregate=function(rpkm.df,meta.df,title.str){
  for(m in 2:5){
    temp.title.str=paste(title.str,paste(paste('(',m,sep=''),'SCs)',sep = ''),sep='')
    plot.aggregate.per.development.stage.pca(rpkm.df,meta.df,no.scs=m,first.pc='PC1',second.pc='PC2',title.str=temp.title.str,log.samples.pca.scores=F,log.rpkm = T)
  }
}
plot.aggregate.per.development.stage.cor.heatmap=function(rpkm.df,meta.df,no.scs=2,title.str='Test'){
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = rpkm.df))
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.list=split(meta.df,f=meta.df$development.stage)
  #development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=c("R" ,"LR", "ET", "LT" , "ES"  , "S")
  development.stages=intersect(development.stages,names(meta.list))
  aggregate.rpkm.list=list()
  col.factor=c()
  for (i in 1:length(development.stages)){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.rpkm.df=rpkm.df[,rownames(temp.meta.df)]
    temp.aggregate.rpkm.df=create.aggregate.cells(sc.df = temp.rpkm.df,no.scs = no.scs)
    aggregate.rpkm.list[[development.stage]]=temp.aggregate.rpkm.df
    col.factor=append(col.factor,values = rep(x = development.stage,times =dim(temp.aggregate.rpkm.df)[2]) ,after = length(col.factor))
  }
  aggregate.rpkm.df=as.data.frame(aggregate.rpkm.list)
  rownames(aggregate.rpkm.df)=rownames(rpkm.df)
  aggregate.rpkm.df=filter.non.variable.rows(df =aggregate.rpkm.df,cut.off = 1 )
  col.factor.list=get.col.factor(col.factor = col.factor)
  col.str=col.factor.list[['col.str']]
  legend.str=c(col.factor.list$legend.str)
  legend.col=c(col.factor.list$legend.col)
  aggregate.rpkm.cor.mat=cor(aggregate.rpkm.df,method = 'spearman')
  #aggregate.rpkm.cor.mat[is.na(aggregate.rpkm.cor.mat)]=0
  aggregate.rpkm.cor.dist.mat=dist(1-aggregate.rpkm.cor.mat)
  aggregate.rpkm.cor.dist.hclust=hclust(aggregate.rpkm.cor.dist.mat, method = "complete")
  #plot(aggregate.rpkm.cor.dist.hclust)
  aggregate.rpkm.cor.dist.dendrogram=as.dendrogram(aggregate.rpkm.cor.dist.hclust)
  #heatmap.2(x =aggregate.rpkm.cor.mat,margins = c(13,13),trace = 'none',ColSideColors = col.str,labRow = '',labCol = '',dendrogram = 'col',main = title.str,Colv =aggregate.rpkm.cor.dist.dendrogram,key.xlab = '',key.ylab = '',key.title = 'Cor. distribution')
  heatmap.2(x =aggregate.rpkm.cor.mat,margins = c(13,13),trace = 'none',ColSideColors = col.str,labRow = '',labCol = '',dendrogram = 'col',main = title.str,key.xlab = '',key.ylab = '',key.title = 'Cor. distribution',col=bluered(1000))
  heatmap.2(x =aggregate.rpkm.cor.mat,margins = c(13,13),trace = 'none',ColSideColors = col.str,labRow = '',labCol = '',dendrogram = 'col',main = '',key.xlab = '',key.ylab = '',key.title = '',col=bluered(1000))
  aggregate.rpkm.mat=as.matrix(log10(aggregate.rpkm.df+1))
  logged.aggregate.rpkm.mat=as.matrix(log10(aggregate.rpkm.df+1))
  show(dim(aggregate.rpkm.df))
  heatmap.2(x =logged.aggregate.rpkm.mat,margins = c(13,13),trace = 'none',ColSideColors = col.str,labRow = '',labCol = '',main = '',key.xlab = '',key.ylab = '',key.title = '',Colv = aggregate.rpkm.cor.dist.dendrogram,col=bluered(1000))
  plot.new()
  legend('topright',legend=legend.str,fill=legend.col)
  plot.new()
  legend('center',legend=rep('',length(legend.str)),fill=legend.col,cex=1.5,box.lty=0)
  return(aggregate.rpkm.df)
}
plot.aggregate.per.development.stage.cor.heatmap.at.diff.aggregate=function(rpkm.df,meta.df,title.str){
  aggregates.list=list()
  for(m in 2:5){
    temp.title.str=paste(title.str,paste(paste('(',m,sep=''),'SCs)',sep = ''),sep='')
    temp.aggregate.rpkm.df=plot.aggregate.per.development.stage.cor.heatmap(rpkm.df =rpkm.df,meta.df =meta.df,no.scs = m ,title.str =temp.title.str)
    temp.name=paste('average.scs.',m,sep = '')
    aggregates.list[[temp.name]]=temp.aggregate.rpkm.df
  }
  return(aggregates.list)
}
plot.aggregated.samples.heatmap=function(aggregates.list){
  average.sc.no=names(aggregates.list)
  len.average.sc.no=length(average.sc.no)
  for (i in 1:len.average.sc.no){
    temp.average.sc=average.sc.no[i]
    temp.aggregate.df=aggregates.list[[temp.average.sc]]
    sample.names=colnames(temp.aggregate.df)
    temp.stages=gsub(pattern = 'sample_.*',replacement = '',x = sample.names)
    temp.stages=gsub(pattern = '\\.$','',x = temp.stages)
    stages=temp.stages
    samples.col.list=get.col.factor(col.factor = stages)
    samples.col.str=samples.col.list$col.str
    temp.cor.mat=cor(temp.aggregate.df,method = 'spearman')
    title.str=temp.average.sc
    heatmap.2(x =temp.cor.mat,ColSideColors = samples.col.str ,main='',
              trace='none',dendrogram = 'col',key.ylab = '',key.xlab = '' ,
              key.title = '',margins = c(5,5),labCol = '',labRow = '')
    plot.new()
    legend('center',legend=samples.col.list$legend.str,
           fill=samples.col.list$legend.col,cex=1.5)
    plot.new()
    legend('center',legend=rep('',length(samples.col.list$legend.str)),
           fill=samples.col.list$legend.col,cex=1.5,box.lty = 0)
  }
}
plot.aggregated.samples.pca=function(aggregates.list){
  average.sc.no=names(aggregates.list)
  len.average.sc.no=length(average.sc.no)
  for (i in 1:len.average.sc.no){
    temp.average.sc=average.sc.no[i]
    temp.aggregate.df=aggregates.list[[temp.average.sc]]
    sample.names=colnames(temp.aggregate.df)
    temp.stages=gsub(pattern = 'sample_.*',replacement = '',x = sample.names)
    temp.stages=gsub(pattern = '\\.$','',x = temp.stages)
    stages=temp.stages
    samples.col.list=get.col.factor(col.factor = stages)
    samples.col.str=samples.col.list$col.str
    temp.aggregate.df=filter.non.variable.rows(df = temp.aggregate.df,cut.off = 10)
    #gene.var.vec=as.numeric(as.character(apply(temp.aggregate.df,1,var)))
    #plot(density(gene.var.vec),xlab='Gene variance',main='Gene variance distribution')
    pca.obj=prcomp(x = t(temp.aggregate.df),retx = T,center = T,scale. = T)
    pca.scores=pca.obj$x
    pairs(x = pca.scores[,1:5],col=samples.col.str,pch=19,cex=1.2)
    legend('topright',legend=samples.col.list$legend.str,fill=samples.col.list$legend.col,cex=1.5)
    plot.new()
    legend('center',legend=rep('',length(samples.col.list$legend.str)),fill=samples.col.list$legend.col,cex=1.5,box.lty = 0)
  }
}
plot.aggregate.samples.pca.at.varying.cut.off=function(rpkm.df,meta.df){
  for(m in 2:5){
    cut.off=m
    plot.aggregate.samples.pca(rpkm.df = rpkm.df,meta.df = meta.df,no.scs =cut.off )
  }
}
plot.aggregate.samples.pca=function(rpkm.df,meta.df,no.scs=2){
  aggregate.rpkm.df=create.aggregate.cells(sc.df = rpkm.df,no.scs = no.scs)
  aggregate.rpkm.df=filter.non.variable.rows(aggregate.rpkm.df,cut.off = 1)
  samples.pca=prcomp(t(aggregate.rpkm.df),retx=T,center=T,scale.=T)
  samples.pca.summary=summary(samples.pca)$importance
  barplot(100*samples.pca.summary[2,],main='Components explained variance',ylab='%')
  abline(h=100/length(colnames(samples.pca.summary)),b=1)
  first.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,1]))*100,2)
  second.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,2]))*100,2)
  samples.pca.scores=samples.pca$x
  x.points=as.numeric(as.character(samples.pca.scores[,1]))
  y.points=as.numeric(as.character(samples.pca.scores[,2]))
  plot(x=x.points,y=y.points,xlab='PC1',ylab='PC2',main='Test',pch=19,cex=.5)
  abline(h=0,b=1,col = "lightgray")
  abline(v=0,b=1,col = "lightgray")
}
plot.markers.pairwise.scatterplot=function(rpkm.df,meta.df,markers.df){
  markers.vec=intersect(rownames(rpkm.df),rownames(markers.df))
  len.markers.vec=length(markers.vec)
  meta.df=meta.df[colnames(rpkm.df),]
  samples.col.list=get.col.factor(col.factor = as.character(meta.df$development.stage))
  if(len.markers.vec>2){
    marker.combn=combn(x = markers.vec,2)
    len.marker.combn=dim(marker.combn)[2]
    for(n in 1:len.marker.combn){
      temp.pair.vec=marker.combn[,n]
      first.marker=temp.pair.vec[1]
      sec.marker=temp.pair.vec[2]
      categories.str=convert.to.title.case(in.str = paste(as.character(markers.df[first.marker,'category']),as.character(markers.df[sec.marker,'category']),sep = '_vs_'))
      first.gene.name=ifelse(markers.df[first.marker,'X.Gene.Name.or.Symbol.']=='null',first.marker,as.character(markers.df[first.marker,'X.Gene.Name.or.Symbol.']))
      sec.gene.name=ifelse(markers.df[sec.marker,'X.Gene.Name.or.Symbol.']=='null',sec.marker,as.character(markers.df[sec.marker,'X.Gene.Name.or.Symbol.']))
      x.points=log10(as.numeric(as.character(rpkm.df[first.marker,]))+1)
      y.points=log10(as.numeric(as.character(rpkm.df[sec.marker,]))+1)
      plot(x=x.points,y=y.points,xlab=first.gene.name,ylab=sec.gene.name,main=categories.str,pch=19,cex=1.5,col=samples.col.list$col.str)
      legend('topright',legend=samples.col.list$legend.str,fill=samples.col.list$legend.col,cex=1.0)
    }
  }
  else{
    first.marker=markers.vec[1]
    sec.marker=markers.vec[2]
    categories.str=convert.to.title.case(in.str = paste(as.character(markers.df[first.marker,'category']),as.character(markers.df[sec.marker,'category']),sep = '_vs_'))
    first.gene.name=ifelse(as.character(markers.df[first.marker,'X.Gene.Name.or.Symbol.'])=='null',first.marker,as.character(markers.df[first.marker,'X.Gene.Name.or.Symbol.']))
    sec.gene.name=ifelse(as.character(markers.df[sec.marker,'X.Gene.Name.or.Symbol.'])=='null',sec.marker,as.character(markers.df[sec.marker,'X.Gene.Name.or.Symbol.']))
    x.points=log10(as.numeric(as.character(rpkm.df[first.marker,]))+1)
    y.points=log10(as.numeric(as.character(rpkm.df[sec.marker,]))+1)
    plot(x=x.points,y=y.points,ylab=sec.gene.name,xlab=first.gene.name,main='',pch=19,cex=1.5,col=samples.col.list$col.str)
    legend('topright',legend=samples.col.list$legend.str,fill=samples.col.list$legend.col,cex=1.0)
  }
}
chunk <- function(in.vec,chunk.sizes){
  len.in.vec=length(in.vec)
  in.vec=sample(in.vec,size =len.in.vec,replace = F)
  chunk.index=seq(1,len.in.vec,by =chunk.sizes)
  len.chunk.index=length(chunk.index)
  out.list=list()
  for(m in 1:len.chunk.index){
    d=chunk.index[m]
    start.id=d
    end.in=d+chunk.sizes-1
    sub.vec=in.vec[start.id:end.in]
    sub.vec=sub.vec[!is.na(sub.vec)]
    if(length(sub.vec)<chunk.sizes){
      next
    }
    list.item.name=paste('chunk',d,sep ='_')
    out.list[[list.item.name]]=sub.vec
  }
  return(out.list)
} 
chunk.non.randomly <- function(in.vec,chunk.sizes){
  len.in.vec=length(in.vec)
  chunk.index=seq(1,len.in.vec,by =chunk.sizes)
  len.chunk.index=length(chunk.index)
  out.list=list()
  for(m in 1:len.chunk.index){
    d=chunk.index[m]
    start.id=d
    end.in=d+chunk.sizes-1
    sub.vec=in.vec[start.id:end.in]
    sub.vec=sub.vec[!is.na(sub.vec)]
    list.item.name=paste('chunk',d,sep ='_')
    out.list[[list.item.name]]=sub.vec
  }
  return(out.list)
} 
#Filters none expressed samples
#Filters none expressed samples
filter.none.expressed.samples=function(input.data){
  samples.count=apply(input.data,2,Matrix::nnzero) 
  filtered.samples=names(samples.count[samples.count>0])
  return(input.data[,filtered.samples])
}
#Filters genes with zero variance across columns(samples)
filter.genes.with.zero.variance=function(df){ 
  genes.sd=apply(df,1,sd)
  filt.genes=names(genes.sd[genes.sd>0])
  return(df[filt.genes,])
}
filter.none.expressed.genes=function(input.data){
  genes.detection.count=apply(input.data,1,Matrix::nnzero)
  filt.genes.detection.count=names(genes.detection.count[genes.detection.count>0])
  return(input.data[filt.genes.detection.count,])
}
#Abbriviates the stages of development in the meta file
change.development.stage.to.abbrev.stage.meta.df=function(meta.df){
  stages.vec=as.character(meta.df$development.stage)
  meta.df$development.stage=gsub(gsub(gsub(gsub(gsub(gsub(x = stages.vec,pattern = 'late.trophozoite',
                                                          replacement = 'T'),pattern = 'late.ring',
                                                     replacement = 'LR'),pattern = 'ring',replacement = 'R'),
                                           pattern = 'early.schizont',replacement = 'ES'),
                                      pattern = 'schizont',replacement = 'S'),pattern = 'early.trophozoite',
                                 replacement = 'ET')
  return(meta.df)
}
#Runs the SC samples and filters based on QC parameters
filter.sc.basic.qc=function(sc.rpkm.df,sc.meta.df,no.detected.genes=200,rpkm.cut.off=1,samples.cut.off=1,filt.by.rpkm=F,no.unique.reads=10000,proportion.per.stage=.2,filter.genes.proportions.per.stage=F,filter.genes.proportions.for.all=F,all.proportions=.2){
  in.rpkm.df= filter.none.expressed.samples(filter.none.expressed.genes(input.data = sc.rpkm.df))
  in.meta.df=sc.meta.df[colnames(in.rpkm.df),]
  in.meta.df=in.meta.df[which(as.numeric(as.character(in.meta.df$Uniquely_mapped_reads_number_))>=no.unique.reads),]
  in.rpkm.df=in.rpkm.df[,rownames(in.meta.df)]
  if(filt.by.rpkm){
    in.rpkm.df=filter.none.expressed.samples(filter.rpkm.less.than.cutoff(df = sc.rpkm.df,rpkm.cutoff =rpkm.cut.off,no.samples = 1 ))
  }
  #Filter based no. detected genes
  genes.count=as.numeric(as.character(apply(in.rpkm.df,2,nnzero)))
  tran.in.rpkm.df=data.frame(t(in.rpkm.df))
  tran.in.rpkm.df$no.of.detected.genes=genes.count
  tran.in.rpkm.df=tran.in.rpkm.df[which(tran.in.rpkm.df$no.of.detected.genes>=no.detected.genes),]
  in.rpkm.df=data.frame(t(subset(tran.in.rpkm.df,select=-no.of.detected.genes)))
  in.rpkm.df=filter.none.expressed.genes(input.data = in.rpkm.df)
  out.meta.df=in.meta.df[colnames(in.rpkm.df),]
  if(filter.genes.proportions.per.stage){
    meta.list=split(out.meta.df,f=out.meta.df$development.stage)
    temp.stages=names(meta.list)
    len.temp.stages=length(temp.stages)
    genes.to.keep=c()
    samples.to.keep=c()
    for(m in 1:len.temp.stages){
      temp.stage=temp.stages[m]
      temp.meta.df=meta.list[[temp.stage]]
      temp.samples=rownames(temp.meta.df)
      temp.rpkm.df= in.rpkm.df[,temp.samples]
      temp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = temp.rpkm.df))
      temp.samples=colnames(temp.rpkm.df)
      #detection.limit.cut.off=round(proportion.per.stage*length(temp.samples))
      detection.limit.cut.off=3
      temp.rpkm.df=filter.rpkm.less.than.cutoff(df =temp.rpkm.df,no.samples = detection.limit.cut.off )
      temp.genes.count=as.numeric(as.character(apply(temp.rpkm.df,2,nnzero)))
      temp.tran.in.rpkm.df=data.frame(t(temp.rpkm.df))
      temp.tran.in.rpkm.df$no.of.detected.genes=temp.genes.count
      temp.tran.in.rpkm.df=temp.tran.in.rpkm.df[which(temp.tran.in.rpkm.df$no.of.detected.genes>=no.detected.genes),]
      temp.rpkm.df=data.frame(t(subset(temp.tran.in.rpkm.df,select=-no.of.detected.genes)))
      temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
      temp.genes.to.keep=rownames(temp.rpkm.df)
      genes.to.keep=append(genes.to.keep,temp.genes.to.keep,after = length(genes.to.keep))
      temp.samples.to.keep=colnames(temp.rpkm.df)
      samples.to.keep=append(samples.to.keep,temp.samples.to.keep,length(samples.to.keep))
    }
    genes.to.keep=unique(genes.to.keep)
    in.rpkm.df=in.rpkm.df[genes.to.keep,samples.to.keep]  
    out.meta.df=out.meta.df[samples.to.keep,]
  }
  if(filter.genes.proportions.for.all){
    no.samples=dim(in.rpkm.df)[2]
    no.samples.cut.off=round(all.proportions*no.samples)
    no.samples.cut.off=3
    show(no.samples.cut.off)
    in.rpkm.df$no.samples.detected=as.numeric(apply(in.rpkm.df,1,nnzero))
    in.rpkm.df=subset(in.rpkm.df,no.samples.detected>=no.samples.cut.off)
    in.rpkm.df=subset(in.rpkm.df,select=-no.samples.detected)
    #in.rpkm.df=filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff =rpkm.cut.off,no.samples = no.samples.cut.off)
    #in.rpkm.df=filter.none.expressed.samples(in.rpkm.df)
    temp.genes.count=as.numeric(as.character(apply(in.rpkm.df,2,nnzero)))
    temp.tran.in.rpkm.df=data.frame(t(in.rpkm.df))
    temp.tran.in.rpkm.df$no.of.detected.genes=temp.genes.count
    temp.tran.in.rpkm.df=temp.tran.in.rpkm.df[which(temp.tran.in.rpkm.df$no.of.detected.genes>=no.detected.genes),]
    in.rpkm.df=data.frame(t(subset(temp.tran.in.rpkm.df,select=-no.of.detected.genes)))
    in.rpkm.df=filter.none.expressed.genes(input.data = in.rpkm.df)
    out.meta.df=out.meta.df[colnames(in.rpkm.df),]
  }
  out.list=list(rpkm=in.rpkm.df,meta=out.meta.df)
  return(out.list)
  #Possible filters:
  #cor.out.list=filter.samples.based.on.cor(rpkm.df =in.rpkm.df,meta.df = in.meta.df,plot = F )
}
filter.sc.counts.basic.qc=function(sc.count.df,sc.meta.df,no.detected.genes=200,count.cut.off=1,no.samples=1,no.unique.reads=10000,proportion.per.stage=.2,filter.genes.proportions.per.stage=F,filter.genes.proportions.for.all=F,all.proportions=.2){
  in.count.df= filter.none.expressed.samples(filter.none.expressed.genes(input.data = sc.count.df))
  in.meta.df=sc.meta.df[colnames(in.count.df),]
  in.meta.df$count.unique.reads=as.numeric(colSums(in.count.df))
  in.meta.df=in.meta.df[which(in.meta.df$count.unique.reads>=no.unique.reads),]
  in.count.df=in.count.df[,rownames(in.meta.df)]
  in.count.df=filter.none.expressed.samples(filter.count.less.than.cutoff(df = in.count.df,count.cutoff =count.cut.off,no.samples = no.samples ))
  #Filter based on no.detected genes
  genes.count=as.numeric(as.character(apply(in.count.df,2,nnzero)))
  tran.in.count.df=data.frame(t(in.count.df))
  tran.in.count.df$no.of.detected.genes=genes.count
  tran.in.count.df=tran.in.count.df[which(tran.in.count.df$no.of.detected.genes>=no.detected.genes),]
  in.count.df=data.frame(t(subset(tran.in.count.df,select=-no.of.detected.genes)))
  in.count.df=filter.none.expressed.genes(input.data = in.count.df)
  out.meta.df=in.meta.df[colnames(in.count.df),]
  if(filter.genes.proportions.per.stage){
    meta.list=split(out.meta.df,f=out.meta.df$development.stage)
    temp.stages=names(meta.list)
    len.temp.stages=length(temp.stages)
    genes.to.keep=c()
    samples.to.keep=c()
    for(m in 1:len.temp.stages){
      temp.stage=temp.stages[m]
      temp.meta.df=meta.list[[temp.stage]]
      temp.samples=rownames(temp.meta.df)
      temp.count.df= in.count.df[,temp.samples]
      temp.count.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = temp.count.df))
      temp.samples=colnames(temp.count.df)
      detection.limit.cut.off=round(proportion.per.stage*length(temp.samples))
      temp.count.df=filter.count.less.than.cutoff(df =temp.count.df,no.samples = detection.limit.cut.off )
      temp.genes.count=as.numeric(as.character(apply(temp.count.df,2,nnzero)))
      temp.tran.in.count.df=data.frame(t(temp.count.df))
      temp.tran.in.count.df$no.of.detected.genes=temp.genes.count
      temp.tran.in.count.df=temp.tran.in.count.df[which(temp.tran.in.count.df$no.of.detected.genes>=no.detected.genes),]
      temp.count.df=data.frame(t(subset(temp.tran.in.count.df,select=-no.of.detected.genes)))
      temp.count.df=filter.none.expressed.genes(input.data = temp.count.df)
      temp.genes.to.keep=rownames(temp.count.df)
      genes.to.keep=append(genes.to.keep,temp.genes.to.keep,after = length(genes.to.keep))
      temp.samples.to.keep=colnames(temp.count.df)
      samples.to.keep=append(samples.to.keep,temp.samples.to.keep,length(samples.to.keep))
    }
    genes.to.keep=unique(genes.to.keep)
    in.count.df=in.count.df[genes.to.keep,samples.to.keep]  
    out.meta.df=out.meta.df[samples.to.keep,]
  }
  if(filter.genes.proportions.for.all){
    no.samples=dim(in.count.df)[2]
    no.samples.cut.off=round(all.proportions*no.samples)
    in.count.df=filter.count.less.than.cutoff(df = in.count.df,count.cutoff =count.cut.off,no.samples = no.samples.cut.off)
    in.count.df=filter.none.expressed.samples(in.count.df)
    temp.genes.count=as.numeric(as.character(apply(in.count.df,2,nnzero)))
    temp.tran.in.count.df=data.frame(t(in.count.df))
    temp.tran.in.count.df$no.of.detected.genes=temp.genes.count
    temp.tran.in.count.df=temp.tran.in.count.df[which(temp.tran.in.count.df$no.of.detected.genes>=no.detected.genes),]
    in.count.df=data.frame(t(subset(temp.tran.in.count.df,select=-no.of.detected.genes)))
    in.count.df=filter.none.expressed.genes(input.data = in.count.df)
    out.meta.df=out.meta.df[colnames(in.count.df),]
  }
  out.list=list(count=in.count.df,meta=out.meta.df)
  return(out.list)
}
create.data.frames.from.df.list=function(df.list,counts=F){
  df.names=names(df.list)
  genes=unique(unlist((lapply(df.list,rownames))))
  len.df.names=length(df.names)
  out.list=list()
  samples.names=c()
  for(n in 1:len.df.names){
    df.list.name=df.names[n]
    temp.df=df.list[[df.list.name]]
    intersect.genes=intersect(rownames(temp.df),genes)
    temp.df=temp.df[genes,]
    temp.df[is.na(temp.df)]=ifelse(counts,0,0.00)
    out.list=append(out.list,temp.df,after=length(out.list))
    samples.names=append(samples.names,colnames(temp.df),after = length(samples.names))
  }
  out.df=data.frame(out.list)
  rownames(out.df)=genes
  colnames(out.df)=samples.names
  return(out.df)
}
plot.pairwise.detection.heatmap=function(df){
  samples=colnames(df)
  len.samples=length(samples)
  for(m in 1:len.samples){
    for(n in 1:len.samples){
      first.sample=samples[m]
      sec.sample=samples[n]
      temp.df=filter.genes.not.expressed.in.all.samples(df[,c(first.sample,sec.sample)])
    }
  }
}
filter.sc.based.on.cor=function(rpkm.df,meta.df,cor.cut.off=.2){
  cor.within.passed.list=filter.samples.with.higher.cor.within.stage(rpkm.df =rpkm.df,meta.df =meta.df )
  out.list=get.top.cor.samples.per.development.stage(rpkm.df =cor.within.passed.list$rpkm,meta.df = cor.within.passed.list$meta,cor.cut.off = cor.cut.off )
  return(out.list)
}
plot.markers.mean.expression.barplot.per.timepoint=function(rpkm.df,meta.df,markers.df,log.rpkm=T){
  rpkm.df=filter.none.expressed.genes( df = rpkm.df)
  psuedo.value=min(rpkm.df[rpkm.df>0])/2
  log.rpkm.df=log2(rpkm.df+psuedo.value)
  if(log.rpkm){
    rpkm.df=log.rpkm.df
  }
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df$development.stage=factor(meta.df$development.stage,level=c('R','LR','ET','T','ES','S'))
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes  )
  markers.df=markers.df[markers.vec,]
  in.rpkm.df=rpkm.df[markers.vec,]
  #in.rpkm.df[is.na(in.rpkm.df)]=0
  in.rpkm.df=filter.none.expressed.genes(in.rpkm.df)
  meta.list=split(meta.df,f = meta.df$development.stage)
  development.stages=names(meta.list)
  expressed.markers=rownames(in.rpkm.df)
  no.expressed.markers=length(expressed.markers)
  for(m in 1:no.expressed.markers){
    marker=expressed.markers[m]
    categories.name=as.character(markers.df[marker,'category'])
    marker.lab=as.character(markers.df[marker,'X.Gene.Name.or.Symbol.'])
    marker.lab=ifelse(marker.lab=='null',marker,marker.lab)
    values.vec=c()
    mean.vec=c()
    se.vec=c()
    col.vec=c()
    for (i in 1:length(development.stages)){
      development.stage=development.stages[i]
      temp.meta.df=meta.list[[development.stage]]
      temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
      peaking.stages=as.character(unlist(strsplit(x = as.character(markers.peaking.stages.list[categories.name]),split = ',')))
      if(marker=="PF3D7_0417200"){
        peaking.stages=c('T')
      }
      temp.col=ifelse(development.stage %in% peaking.stages , 'blue','gray')
      col.vec=append(col.vec,temp.col,after = length(col.vec))
      #temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
      if(is.null(dim(temp.rpkm.df))){
        next
      }
      rpkm=as.numeric(melt(temp.rpkm.df[marker,])$value)
      #expr.list[[development.stage]]=rpkm
      expr.vec=rpkm
      expr.stats=describe(expr.vec)
      values.vec=append(x = values.vec,values =i ,after = length(values.vec))
      mean.vec=append(x = mean.vec,values = expr.stats$mean,after = length(mean.vec))
      se.vec=append(x = se.vec,values =expr.stats$se ,after = length(se.vec))
    }
    marker.stats.df=data.frame(values=values.vec,mean=mean.vec,se=se.vec)
    rownames(marker.stats.df)=development.stages
    width.size.vec=c(rep(.1,times = dim(marker.stats.df)[1]))
    editted.error.bars(stats=marker.stats.df,bars=T,main=marker.lab,ylab = 'rpkm',xlab='',eyes = F,col.str =col.vec,width.size.vec = width.size.vec)
  }
}
plot.markers.mean.expression.barplot.per.parental.timepoint=function(rpkm.df,meta.df,markers.df,log.rpkm=T){
  rpkm.df=filter.none.expressed.genes( df = rpkm.df)
  psuedo.value=min(rpkm.df[rpkm.df>0])/2
  log.rpkm.df=log2(rpkm.df+psuedo.value)
  if(log.rpkm){
    rpkm.df=log.rpkm.df
  }
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df=create.parental.stage(meta.df = meta.df)
  meta.df$parent.development.stage=factor(meta.df$parent.development.stage,level=c('R','T','S'))
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes  )
  markers.df=markers.df[markers.vec,]
  in.rpkm.df=rpkm.df[markers.vec,]
  #in.rpkm.df[is.na(in.rpkm.df)]=0
  in.rpkm.df=filter.none.expressed.genes(in.rpkm.df)
  meta.list=split(meta.df,f = meta.df$parent.development.stage)
  development.stages=names(meta.list)
  expressed.markers=rownames(in.rpkm.df)
  no.expressed.markers=length(expressed.markers)
  for(m in 1:no.expressed.markers){
    marker=expressed.markers[m]
    categories.name=as.character(markers.df[marker,'category'])
    marker.lab=as.character(markers.df[marker,'X.Gene.Name.or.Symbol.'])
    marker.lab=ifelse(marker.lab=='null',marker,marker.lab)
    values.vec=c()
    mean.vec=c()
    se.vec=c()
    col.vec=c()
    for (i in 1:length(development.stages)){
      development.stage=development.stages[i]
      temp.meta.df=meta.list[[development.stage]]
      temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
      peaking.stages=as.character(unlist(strsplit(x = as.character(markers.peaking.stages.list[categories.name]),split = ',')))
      if(marker=="PF3D7_0417200"){
        peaking.stages=c('T')
      }
      temp.col=ifelse(development.stage %in% peaking.stages , 'blue','gray')
      col.vec=append(col.vec,temp.col,after = length(col.vec))
      #temp.rpkm.df=filter.none.expressed.genes(input.data = temp.rpkm.df)
      if(is.null(dim(temp.rpkm.df))){
        next
      }
      rpkm=as.numeric(melt(temp.rpkm.df[marker,])$value)
      #expr.list[[development.stage]]=rpkm
      expr.vec=rpkm
      expr.stats=describe(expr.vec)
      values.vec=append(x = values.vec,values =i ,after = length(values.vec))
      mean.vec=append(x = mean.vec,values = expr.stats$mean,after = length(mean.vec))
      se.vec=append(x = se.vec,values =expr.stats$se ,after = length(se.vec))
    }
    marker.stats.df=data.frame(values=values.vec,mean=mean.vec,se=se.vec)
    rownames(marker.stats.df)=development.stages
    editted.error.bars(stats=marker.stats.df,bars=T,main=marker.lab,ylab = 'rpkm',xlab='',eyes = F,col.str =col.vec)
  }
}
plot.markers.mean.expression.dotplot.per.timepoint=function(rpkm.df,meta.df,
                                                            markers.df,label.main=T){
  meta.df=meta.df[colnames(rpkm.df),]
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes  )
  markers.df=markers.df[markers.vec,]
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df[is.na(in.rpkm.df)]=0
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  meta.df=meta.df[colnames(in.rpkm.df),]
  markers.df=markers.df[rownames(in.rpkm.df),]
  markers.list=split(markers.df,f=markers.df$category)
  expressed.markers.vec=rownames(markers.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=split(meta.df,f = meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  intersect.development.stages = intersect(ordered.development.stages,development.stages)
  len.intersect.development.stages=length(intersect.development.stages)
  markers.mean.expression=c()
  log.marker.mean.expression.vec=c()
  markers.stages.vec=c()
  markers.names=c()
  fraction.detected.in.vec=c()
  samples.category.vec=c()
  for(m in 1:len.expressed.markers.vec){
    marker=expressed.markers.vec[m]
    marker.stage=as.character(markers.df[marker,'category'])
    marker.name=as.character(markers.df[marker,'X.Gene.Name.or.Symbol.'])
    marker.gene.name=ifelse(marker.name=='null',marker,marker.name)
    marker.name =ifelse(marker.name=='null',marker,paste(marker,marker.gene.name,sep=':'))
    for( n in 1:len.intersect.development.stages){
      intersect.development.stage=intersect.development.stages[n]
      temp.meta.df=meta.list[[intersect.development.stage]]
      stage.samples=rownames(temp.meta.df)
      stage.rpkm.vec=as.numeric(as.character(in.rpkm.df[marker,stage.samples]))
      marker.mean.expression=mean(stage.rpkm.vec)
      marker.median.expression=median(stage.rpkm.vec)
      log.marker.mean.expression=log2(marker.mean.expression+0.01)
      #log.marker.mean.expression=log.marker.mean.expression/length(stage.rpkm.vec)
      #log.marker.mean.expression=log10(marker.median.expression+1)
      log.marker.mean.expression.vec=append(log.marker.mean.expression.vec,log.marker.mean.expression,length(log.marker.mean.expression.vec))
      fraction.detected.in=length(stage.rpkm.vec[stage.rpkm.vec>0])/length(stage.samples)
      markers.mean.expression=append(markers.mean.expression,marker.mean.expression,after = length(marker.mean.expression))
      markers.stages.vec=append(markers.stages.vec,marker.stage,length(markers.stages.vec))
      fraction.detected.in.vec=append(fraction.detected.in.vec,fraction.detected.in,after = length(fraction.detected.in.vec))
      markers.names=append(markers.names,marker.name,length(markers.names))
      samples.category.vec=append(samples.category.vec,intersect.development.stage,length(samples.category.vec))
    }
  }
  markers.dotplot.df=data.frame(names=markers.names,stages=markers.stages.vec,
                                mean.rpkm=markers.mean.expression,
                                fraction.detected=fraction.detected.in.vec,
                                log.mean.expression=log.marker.mean.expression.vec,
                                samples.category=samples.category.vec)
  markers.dotplot.df.list=split(markers.dotplot.df,f=markers.dotplot.df$stages)
  markers.dotplot.list.names=names(markers.dotplot.df.list)
  len.markers.dotplot.list.names=length(markers.dotplot.list.names)
  for(l in 1:len.markers.dotplot.list.names){
    markers.dotplot.list.name=markers.dotplot.list.names[l]
    temp.markers.dotplot.df=markers.dotplot.df.list[[markers.dotplot.list.name]]
    temp.markers.dotplot.df$ordered.samples.categories <- 
      factor(temp.markers.dotplot.df$samples.category, 
             as.character(c('ring','late.ring','early.trophozoite',
                            'late.trophozoite','early.schizont','schizont')))
    no.markers=unique(as.character(temp.markers.dotplot.df[,'names']))
    len.no.markers=length(no.markers)
    if(len.no.markers<=10){
      point.size=as.numeric(as.character(temp.markers.dotplot.df$fraction.detected))
      if(label.main)
      {
        temp.ggplot=ggplot(temp.markers.dotplot.df,aes(x = ordered.samples.categories, y = names)) 
        temp.ggplot+ geom_point(aes(size=log.mean.expression,colour=fraction.detected)) 
        temp.ggplot+ ggtitle(label = markers.dotplot.list.name)+ylab(label = '')
        temp.ggplot+xlab(label = '') 
        temp.ggplot+theme(axis.text.x = element_text(angle = 90, hjust = 1))
      }
      else{
        temp.ggplot=ggplot(temp.markers.dotplot.df, 
                           aes(x = ordered.samples.categories, y = names)) + geom_point(aes(size=log.mean.expression,colour=fraction.detected)) + ggtitle(label = '')+ylab(label = '')+xlab(label = '') +theme(axis.text.x =element_blank())
      }
      print(temp.ggplot)
    }
    else{
      chunked.vec=chunk.non.randomly(in.vec = no.markers,chunk.sizes = 10)
      names.chunked.vec=names(chunked.vec)
      len.names.chunked.vec=length(names.chunked.vec)
      for(ch in 1:len.names.chunked.vec ){
        temp.names.chunked=names.chunked.vec[ch]
        temp.chunked.vec=chunked.vec[[temp.names.chunked]]
        markers.chunked.df=temp.markers.dotplot.df[which(temp.markers.dotplot.df[,'names'] %in% temp.chunked.vec),]
        point.size=as.numeric(as.character(markers.chunked.df$fraction.detected))
        temp.ggplot=ggplot(markers.chunked.df, aes(x = ordered.samples.categories, y = names)) + geom_point(aes(size=log.mean.expression,colour=fraction.detected)) + ggtitle(label = markers.dotplot.list.name)+ylab(label = '')+xlab(label = '') +theme(axis.text.x = element_text(angle = 90, hjust = 1))
        print(temp.ggplot)
      }
    }
  }
  return(markers.dotplot.df)
}
plot.variant.genes.mean.expression.dotplot.per.timepoint=function(rpkm.df,meta.df,
                                                            markers.df,label.main=T){
  meta.df=meta.df[colnames(rpkm.df),]
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),
                        y =all.genes  )
  markers.df=markers.df[markers.vec,]
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df[is.na(in.rpkm.df)]=0
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  meta.df=meta.df[colnames(in.rpkm.df),]
  markers.df=markers.df[rownames(in.rpkm.df),]
  markers.list=split(markers.df,f=markers.df$category)
  expressed.markers.vec=rownames(markers.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.list=split(meta.df,f = meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("R" ,"LR","ET","T" ,"ES", "S")
  intersect.development.stages = intersect(ordered.development.stages,
                                           development.stages)
  len.intersect.development.stages=length(intersect.development.stages)
  markers.mean.expression=c()
  log.marker.mean.expression.vec=c()
  markers.stages.vec=c()
  markers.names=c()
  fraction.detected.in.vec=c()
  samples.category.vec=c()
  for(m in 1:len.expressed.markers.vec){
    marker=expressed.markers.vec[m]
    marker.stage=as.character(markers.df[marker,'category'])
    marker.name=as.character(markers.df[marker,'X.Gene.Name.or.Symbol.'])
    marker.gene.name=ifelse(marker.name=='null',marker,marker.name)
    marker.name =ifelse(marker.name=='null',marker,paste(marker,marker.gene.name,sep=':'))
    marker.name =marker
    for( n in 1:len.intersect.development.stages){
      intersect.development.stage=intersect.development.stages[n]
      temp.meta.df=meta.list[[intersect.development.stage]]
      stage.samples=rownames(temp.meta.df)
      stage.rpkm.vec=as.numeric(as.character(in.rpkm.df[marker,stage.samples]))
      marker.mean.expression=mean(stage.rpkm.vec)
      marker.median.expression=median(stage.rpkm.vec)
      log.marker.mean.expression=log2(marker.mean.expression+0.01)
      #log.marker.mean.expression=log.marker.mean.expression/length(stage.rpkm.vec)
      #log.marker.mean.expression=log10(marker.median.expression+1)
      log.marker.mean.expression.vec=append(log.marker.mean.expression.vec,
                                            log.marker.mean.expression,
                                            length(log.marker.mean.expression.vec))
      fraction.detected.in=length(stage.rpkm.vec[stage.rpkm.vec>0])/length(stage.samples)
      markers.mean.expression=append(markers.mean.expression,marker.mean.expression,
                                     after = length(marker.mean.expression))
      markers.stages.vec=append(x = markers.stages.vec,values = marker.stage,after =
                                  length(markers.stages.vec))
      fraction.detected.in.vec=append(fraction.detected.in.vec,
                                      fraction.detected.in,after = length(fraction.detected.in.vec))
      markers.names=append(markers.names,marker.name,length(markers.names))
      samples.category.vec=append(samples.category.vec,
                                  intersect.development.stage,
                                  length(samples.category.vec))
    }
  }
  markers.dotplot.df=data.frame(names=markers.names,stages=markers.stages.vec,
                                mean.rpkm=markers.mean.expression,
                                fraction.detected=fraction.detected.in.vec,
                                log.mean.expression=log.marker.mean.expression.vec,
                                samples.category=samples.category.vec)
  show(head(markers.dotplot.df))
  markers.dotplot.df.list=split(markers.dotplot.df,f=markers.dotplot.df$stages)
  markers.dotplot.list.names=names(markers.dotplot.df.list)
  len.markers.dotplot.list.names=length(markers.dotplot.list.names)
  for(l in 1:len.markers.dotplot.list.names){
    markers.dotplot.list.name=markers.dotplot.list.names[l]
    temp.markers.dotplot.df=markers.dotplot.df.list[[markers.dotplot.list.name]]
    temp.markers.dotplot.df$ordered.samples.categories <- 
      factor(temp.markers.dotplot.df$samples.category, 
             as.character(c('R','LR','ET',
                            'T','ES','S')))
    no.markers=unique(as.character(temp.markers.dotplot.df[,'names']))
    len.no.markers=length(no.markers)
    if(len.no.markers<=10){
      point.size=as.numeric(as.character(temp.markers.dotplot.df$fraction.detected))
      if(label.main)
      {
        temp.ggplot=ggplot(temp.markers.dotplot.df,aes(x = ordered.samples.categories, y = names)) 
        temp.ggplot+ geom_point(aes(size=log.mean.expression,colour=fraction.detected)) 
        temp.ggplot+ ggtitle(label = markers.dotplot.list.name)+ylab(label = '')
        temp.ggplot+xlab(label = '') 
        temp.ggplot+theme(axis.text.x = element_text(angle = 90, hjust = 0))
      }
      else{
        temp.ggplot=ggplot(temp.markers.dotplot.df, 
                           aes(x = ordered.samples.categories, y = names)) 
        temp.ggplot+ geom_point(aes(size=log.mean.expression,colour=fraction.detected)) 
        temp.ggplot+ ggtitle(label = '')
        temp.ggplot+ylab(label = '')
        temp.ggplot+xlab(label = '') 
        temp.ggplot+theme(axis.text.x =element_blank(angle = 0, hjust = 1))
      }
      print(temp.ggplot)
    }
    else{
      chunked.vec=chunk.non.randomly(in.vec = no.markers,chunk.sizes = 15)
      names.chunked.vec=names(chunked.vec)
      len.names.chunked.vec=length(names.chunked.vec)
      for(ch in 1:len.names.chunked.vec ){
        temp.names.chunked=names.chunked.vec[ch]
        temp.chunked.vec=chunked.vec[[temp.names.chunked]]
        markers.chunked.df=
          temp.markers.dotplot.df[which(temp.markers.dotplot.df[,'names'] %in% temp.chunked.vec),]
        point.size=as.numeric(as.character(markers.chunked.df$fraction.detected))
        temp.ggplot=ggplot(markers.chunked.df, aes(x = ordered.samples.categories, y = names)) + geom_point(aes(size=log.mean.expression,colour=fraction.detected)) + ggtitle(label = markers.dotplot.list.name)+ylab(label = '')+xlab(label = '') +theme(axis.text.x = element_text(angle = 90, hjust = 1))
        print(temp.ggplot)
      }
    }
  }
  return(markers.dotplot.df)
}
plot.variant.genes.expression.dotplot.per.sample=function(rpkm.df,meta.df,
                                                                  markers.df,label.main=T){
  meta.df=meta.df[colnames(rpkm.df),]
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),
                        y =all.genes  )
  markers.df=markers.df[markers.vec,]
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df[is.na(in.rpkm.df)]=0
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  meta.df=meta.df[colnames(in.rpkm.df),]
  markers.df=markers.df[rownames(in.rpkm.df),]
  markers.list=split(markers.df,f=markers.df$category)
  expressed.markers.vec=rownames(markers.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.list=split(meta.df,f = meta.df$development.stage)
  development.stages=names(meta.list)
  #ordered.development.stages=c("R" ,"LR","ET","T" ,"ES", "S")
  ordered.development.stages=c("R" ,"LR","ET")
  intersect.development.stages = intersect(ordered.development.stages,
                                           development.stages)
  meta.df=order.by.target.development.stage.meta.df(meta.df = meta.df,
                                                    order.vec = intersect.development.stages)
  show(dim(meta.df))
  in.rpkm.df=subset(in.rpkm.df,select=rownames(meta.df))
  var.expr.vec.list=list()
  var.expr.list=apply(in.rpkm.df,2,function(col.vec){
    temp.expr.vec=sort(as.numeric(col.vec),decreasing = T)[1:5]
    genes=unlist(lapply(temp.expr.vec,function(x){
      temp.indices.vec=which(col.vec==x)
      indices.vec=ifelse(length(temp.indices.vec)>1,
                              temp.indices.vec,temp.indices.vec[1])
      return(col.vec[indices.vec])
    }))
    return(names(genes))
  })
  all.sc.vec=colnames(in.rpkm.df)
  len.all.sc.vec=length(all.sc.vec)
  fraction.vec=numeric()
  rpkm.vec=numeric()
  top.var.vec=as.character()
  for(m in 1:len.all.sc.vec){
    sc.name=all.sc.vec[m]
    gene.ids=var.expr.list[,sc.name]
    temp.expr.vec=in.rpkm.df[gene.ids,sc.name]
    names(temp.expr.vec)=gene.ids
    barplot.col=abbrev.std.stages.col.vec[meta.df[sc.name,'development.stage']]
    title.str=paste('SC_',meta.df[sc.name,'development.stage'],'_',m,sep='')
    barplot2(temp.expr.vec+0.01,las=2,ylab = 'RPKM',col=barplot.col,
             horiz = F,space=c(0,.5),
             cex.names = .5,border = NA,main = title.str)
    sample.var.vec=in.rpkm.df[expressed.markers.vec,sc.name]
    names(sample.var.vec)=expressed.markers.vec
    top.var.rpkm=max(sample.var.vec)
    other.var.rpkm=sample.var.vec[which(sample.var.vec!=top.var.rpkm)]
    top.var.name=names(sample.var.vec[which(sample.var.vec==top.var.rpkm)])
    top.var.vec=append(top.var.vec,top.var.name,after = length(top.var.vec))
    other.var.sum=ifelse(sum(other.var.rpkm)>0,sum(other.var.rpkm),1)
    fraction.rpkm=top.var.rpkm/other.var.sum
    fraction.vec=append(fraction.vec,values = fraction.rpkm,after = length(fraction.vec))
    rpkm.vec=append(rpkm.vec,max(sample.var.vec),after = length(rpkm.vec))
  }
  col.vec=abbrev.std.stages.col.vec[meta.df[all.sc.vec,'development.stage']]
  xlim.vec=c(-3,14)   
  plot(x=log2(fraction.vec),y=log2(rpkm.vec),ylab='log2(RPKM)',
       pch=19,col=col.vec,cex=1.5,
       frame.plot=T,xlab='Log-fold expression ratio',main='')
  plot(x=log2(fraction.vec),y=log2(rpkm.vec),ylab='log2(RPKM)',
       pch=19,cex=1.5,xlim=xlim.vec,
       frame.plot=T,xlab='Log-fold expression ratio',main='')
  show(length(fraction.vec))
  text(x =log2(fraction.vec) ,y=log2(rpkm.vec),labels = top.var.vec,cex = .3)
  out.df=data.frame(row.names = all.sc.vec,rpkm=rpkm.vec,dominant.var=top.var.vec,
                    fraction.of.total.var=fraction.vec)
  return(out.df)
}
plot.markers.clustering=function(rpkm.df,meta.df,markers.df){
  meta.df=meta.df[colnames(rpkm.df),]
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes)
  markers.df=markers.df[markers.vec,]
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.genes(in.rpkm.df)
  pseudo.rpkm=min(in.rpkm.df[in.rpkm.df!=0])/2
  show(pseudo.rpkm)
  in.rpkm.df=in.rpkm.df+pseudo.rpkm
  expr.markers.vec=rownames(in.rpkm.df)
  markers.df=markers.df[expr.markers.vec,]
  markers.list=split(markers.df,f=markers.df$category)
  marker.categories=names(markers.list)
  len.marker.categories=length(marker.categories)
  for(m in 1:len.marker.categories){
    marker.category=marker.categories[m]
    temp.markers.df=markers.list[[marker.category]]
    temp.rpkm.df=in.rpkm.df[rownames(temp.markers.df),]
    if(dim(temp.rpkm.df)[1]<2){
      next
    }
    markers.col.fact.vec=markers.df[expr.markers.vec,'category']
    markers.col.fact.list=get.col.factor(col.factor =markers.col.fact.vec )
    samples.col.fact.vec=as.character(meta.df[colnames(temp.rpkm.df),'development.stage'])
    samples.col.fact.list=get.col.factor(col.factor =samples.col.fact.vec )
    title.str=convert.to.title.case(in.str = marker.category)
    heatmap.2(x = as.matrix(log10(temp.rpkm.df)),ColSideColors = samples.col.fact.list$col.str,trace='none',margins=c(10,10),main=title.str,col=bluered(1000),labCol = '',cexRow = .5)
    legend('topright',legend=samples.col.fact.list$legend.str,fill=samples.col.fact.list$legend.col)
    annotation.meta.df=meta.df[colnames(temp.rpkm.df),]
    pheatmap(mat = as.matrix(log10(temp.rpkm.df)),main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),fontsize_row = 2,fontsize_col = 2,cellwidth = 2,cellheight = 2,annotation_col = annotation.meta.df)
  }
}
get.mean.development.stage.markers.score.per.cell=function(rpkm.df,meta.df,markers.df,title.str='Test',log.rpkm=T){
  rpkm.df=filter.none.expressed.samples(filter.genes.not.expressed.in.a.number.of.samples(rpkm.df,sample.count = 2))
  all.genes=rownames(rpkm.df)
  expr.markers.vec=intersect(rownames(rpkm.df),rownames(markers.df))
  markers.df=markers.df[expr.markers.vec,]
  markers.list=split(markers.df,f=markers.df$category)
  marker.categories=names(markers.list)
  len.marker.categories=length(marker.categories)
  no.samples=dim(rpkm.df)[2]
  out.mat=matrix(nrow = len.marker.categories,ncol = no.samples)
  row.names.vec=c()
  for(m in 1:len.marker.categories){
    marker.category=marker.categories[m]
    temp.markers.df=markers.list[[marker.category]]
    temp.markers.vec=rownames(temp.markers.df)
    marker.category.score.vec=c()
    if(log.rpkm){
      marker.category.score.vec=get.mean.gene.signature.score.per.cell(rpkm.df = rpkm.df, genes.sig.vec =temp.markers.vec,log.rpkm = T,quant.norm = T )
    }
    else{
      marker.category.score.vec=get.mean.gene.signature.score.per.cell(rpkm.df = rpkm.df, genes.sig.vec =temp.markers.vec,log.rpkm = F )
    }
    out.mat[m,]=as.numeric(marker.category.score.vec)
    marker.category=gsub(pattern = ' ',replacement = '_',x = marker.category)
    row.names.vec=append(row.names.vec,marker.category,after = length(row.names.vec))
  }
  out.df=data.frame(out.mat)
  rownames(out.df)=row.names.vec
  colnames(out.df)=colnames(rpkm.df)
  out.df=data.frame(t(subset(t(out.df),select=-none)))
  meta.df=meta.df[colnames(out.df),]
  col.side.col.str=get.std.stage.cols(as.character(meta.df$development.stage))
  out.mat=as.matrix(out.df)
  #out.mat=(out.mat-rowMeans(out.mat))/rowSds(out.mat)
  #samples.dendrogram=as.dendrogram(hclust(dist(x = t(out.mat)),method = 'complete'))
  heatmap.2(x = out.mat,margins = c(10,10),ColSideColors = col.side.col.str,main=title.str,trace='none',cexRow = .5,cexCol = .3,col=bluered(10000),scale='row')
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
  pheatmap(mat = out.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),cellwidth = 3,cellheight = 5,fontsize_row = 2,fontsize_col = 2)
  heatmap.2(x = cor(out.mat,method = 'spearman'),margins = c(10,10),ColSideColors = col.side.col.str,main=title.str,trace='none',cexRow = .3,cexCol = .3,col=bluered(10000),key.title = 'Correlation')
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
  pheatmap(mat = cor(out.mat,method = 'spearman'),main=title.str,trace='none',color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),cellwidth = 3,cellheight = 3,fontsize_row = 2,fontsize_col = 2)
  heatmap.2(x = as.matrix(dist(t(out.mat))),margins = c(10,10),ColSideColors = col.side.col.str,main=title.str,trace='none',cexRow = .3,cexCol = .3,col=bluered(10000),key.title = 'Euclidean dist.')
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
  pheatmap(mat = as.matrix(dist(t(out.mat))),main=title.str,trace='none',color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),cellwidth = 3,cellheight = 3,fontsize_row = 2,fontsize_col = 2)
  return(out.df)
}
get.development.stage.markers.score.per.category.per.cell=function(rpkm.df,meta.df,markers.df,title.str='Test',log.rpkm = F){
  rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(rpkm.df,sample.count = 2)
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(rpkm.df))
  all.genes=rownames(rpkm.df)
  expr.markers=intersect(all.genes,rownames(markers.df))
  markers.df=markers.df[expr.markers,]
  markers.list=split(markers.df,f=markers.df$category)
  marker.categories=names(markers.list)
  len.marker.categories=length(marker.categories)
  samples=colnames(rpkm.df)
  no.samples=length(samples)
  for(n in 1:len.marker.categories){
    temp.marker.category=marker.categories[n]
    temp.marker.df=markers.list[[temp.marker.category]]
    marker.category.score.df=data.frame()
    genes.sig.vec=rownames(temp.marker.df)
    if(length(genes.sig.vec)<2){
      next
    }
    if(log.rpkm){
      marker.category.score.df=get.gene.signature.score.per.cell(rpkm.df = rpkm.df, genes.sig.vec =genes.sig.vec,log.rpkm = T )
    }
    else{
      marker.category.score.df=get.gene.signature.score.per.cell(rpkm.df = rpkm.df, genes.sig.vec =genes.sig.vec,log.rpkm = F )
    }
    #marker.category.score.df=filter.non.variable.rows(df = marker.category.score.df,cut.off = .5)
    gene.cols.list=get.col.factor(col.factor = as.character(markers.df[rownames(marker.category.score.df),'category']))
    #temp.title.str=paste(title.str,'(',temp.marker.category,')',collapse='')
    temp.title.str=convert.to.title.case(temp.marker.category)
    out.mat=as.matrix(marker.category.score.df)
    samples=colnames(marker.category.score.df)
    stages.vec=as.character(meta.df[samples,]$development.stage)
    col.side.col.str=get.std.stage.cols(in.stages.vec=stages.vec)
    #heatmap.2(x = out.mat,margins = c(13,13),ColSideColors = col.side.col.str,main=temp.title.str,trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),key.title = 'Euclidean dist.',scale='row')
    heatmap.2(x = out.mat,margins = c(13,13),ColSideColors = col.side.col.str,main=temp.title.str,trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),key.title = 'Euclidean dist.')
    legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
    #heatmap.2(x = cor(out.mat,method = 'spearman'),margins = c(10,10),ColSideColors = col.side.col.str,main=title.str,trace='none',cexRow = .3,cexCol = .3,col=bluered(10000),key.title = 'Correlation')
    #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .5)
    #heatmap.2(x = as.matrix(dist(t(out.mat),method = 'euclidian')),margins = c(10,10),ColSideColors = col.side.col.str,main=temp.title.str,trace='none',cexRow = .3,cexCol = .3,col=bluered(10000),key.title = 'Euclidian distance')
    #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
    pca.obj=prcomp(t(out.mat))
    #plot(x=as.numeric(pca.obj$x[,1]),y=as.numeric(pca.obj$x[,2]),pch=19,col = col.side.col.str,cex=1.5,main=temp.title.str,xlab='PC1',ylab='PC2')
    #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
    #pairs(x=pca.obj$x[,1:5],pch=19,col = col.side.col.str,cex=.5)
  }
  #out.mat=matrix(nrow = len.marker.categories,ncol = no.samples)
  #plot.samples.tsne.clustering(rpkm.df =out.mat,meta.df =meta.df, title.str =title.str)
  return(out.mat)
}
get.development.stage.markers.score.per.cell=function(rpkm.df,meta.df,markers.df,title.str='Test',log.rpkm = F){
  rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(rpkm.df,sample.count = 2)
  rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(rpkm.df))
  all.genes=rownames(rpkm.df)
  samples=colnames(rpkm.df)
  no.samples=length(samples)
  #out.mat=matrix(nrow = len.marker.categories,ncol = no.samples)
  expr.markers.vec=intersect(all.genes,rownames(markers.df))
  markers.df=markers.df[expr.markers.vec,]
  all.markers=rownames(markers.df)
  len.markers=length(all.markers)
  marker.category.score.df=data.frame()
  if(log.rpkm){
    marker.category.score.df=get.gene.signature.score.per.cell(rpkm.df = rpkm.df, genes.sig.vec =all.markers,log.rpkm = T )
  }
  else{
    marker.category.score.df=get.gene.signature.score.per.cell(rpkm.df = rpkm.df, genes.sig.vec =all.markers,log.rpkm = F )
  }
  #marker.category.score.df=filter.non.variable.rows(df = marker.category.score.df,cut.off = .5)
  gene.cols.list=get.col.factor(col.factor = as.character(markers.df[rownames(marker.category.score.df),'category']))
  title.str=convert.to.title.case(title.str)
  out.mat=as.matrix(marker.category.score.df)
  samples=colnames(marker.category.score.df)
  stages.vec=as.character(meta.df[samples,]$development.stage)
  col.side.col.str=get.std.stage.cols(in.stages.vec=stages.vec)
  heatmap.2(x = out.mat,margins = c(13,13),ColSideColors = col.side.col.str,main=title.str,trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),RowSideColors = gene.cols.list$col.str,key.title = 'Euclidean dist.',scale='row')
  annotation.meta.df=subset(meta.df[colnames(out.mat),],select = development.stage )
  pheatmap(mat = out.mat,main=title.str,col=bluered(10000),cellwidth = 1,cellheight = 1,cutree_rows = 6,fontsize_row = 2,fontsize_col = 2)
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
  plot.new()
  legend('left',legend=gene.cols.list$legend.str,fill=gene.cols.list$legend.col,box.lty = 0,cex = .8)
  heatmap.2(x = cor(out.mat,method = 'spearman'),margins = c(10,10),ColSideColors = col.side.col.str,main=title.str,trace='none',cexRow = .3,cexCol = .3,col=bluered(10000),key.title = 'Correlation')
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
  heatmap.2(x = as.matrix(dist(t(out.mat),method = 'euclidian')),margins = c(10,10),ColSideColors = col.side.col.str,main=title.str,trace='none',cexRow = .3,cexCol = .3,col=bluered(10000),key.title = 'Euclidian distance')
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
  pca.obj=prcomp(t(out.mat))
  plot(x=as.numeric(pca.obj$x[,1]),y=as.numeric(pca.obj$x[,2]),pch=19,col = col.side.col.str,cex=1.5,main=title.str,xlab='PC1',ylab='PC2')
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
  pairs(x=pca.obj$x[,1:5],pch=19,col = col.side.col.str,cex=.5)
  #plot.samples.tsne.clustering(rpkm.df =out.mat,meta.df =meta.df, title.str =title.str)
  return(out.mat)
}
get.development.stage.clustering.markers=function(rpkm.df,meta.df,markers.df,title.str='Test',log.rpkm = F,in.sc.cut.off=6,gene.cluster.cut.off=5,get.quantile.score=T){
  rpkm.df=filter.none.expressed.samples(filter.genes.not.expressed.in.a.number.of.samples(rpkm.df,sample.count = 2))
  in.rpkm.df= rpkm.df
  #rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(rpkm.df))
  all.genes=rownames(rpkm.df)
  samples=colnames(rpkm.df)
  no.samples=length(samples)
  #out.mat=matrix(nrow = len.marker.categories,ncol = no.samples)
  expr.markers.vec=intersect(all.genes,rownames(markers.df))
  markers.df=markers.df[expr.markers.vec,]
  all.markers=rownames(markers.df)
  len.markers=length(all.markers)
  marker.category.score.df=data.frame()
  if(get.quantile.score){
    if(log.rpkm){
      marker.category.score.df=get.gene.signature.score.per.cell(rpkm.df = rpkm.df, genes.sig.vec =all.markers,log.rpkm = T )
    }
    else{
      marker.category.score.df=get.gene.signature.score.per.cell(rpkm.df = rpkm.df, genes.sig.vec =all.markers,log.rpkm = F )
    }
  }
  else{
    show('Un -scored for markers')
    marker.category.score.df=rpkm.df[all.markers,]
    if(log.rpkm){
      pseudo.value=as.numeric(min(marker.category.score.df[marker.category.score.df!=0]))/2
      marker.category.score.df=log2(marker.category.score.df+pseudo.value)
    }
    else{
      marker.category.score.df=marker.category.score.df
    }
  }
  genes.to.use=c()
  #marker.category.score.df=filter.non.variable.rows(df = marker.category.score.df,cut.off = .5)
  gene.cols.list=get.col.factor(col.factor = as.character(markers.df[rownames(marker.category.score.df),'category']))
  title.str=convert.to.title.case(title.str)
  out.mat=as.matrix(marker.category.score.df)
  samples=colnames(marker.category.score.df)
  stages.vec=as.character(meta.df[samples,]$development.stage)
  col.side.col.str=get.std.stage.cols(in.stages.vec=stages.vec)
  gene.dist=dist(out.mat)
  gene.dist.mat=as.matrix(gene.dist)
  intern <- clValid(obj =gene.dist.mat, nClust = 2:6, method = 'complete',clMethods=c("hierarchical", "kmeans", "diana", "som", "model", "sota", "pam", "clara", "agnes"),validation="internal",metric = "euclidean")
  #internal.clusters=clusters(intern)
  internal.clusters=hclust(d = gene.dist)
  gene.groups <- cutree(internal.clusters, k = 3)
  gene.groups.df=data.frame(genes=names(gene.groups),group=paste('Gene grp: ',as.character(as.numeric(gene.groups)),sep=''),row.names = names(gene.groups ))
  gene.grp.col.side.col.list=get.col.factor(col.factor = as.character(gene.groups.df$group))
  #all.markers.mat=log10(as.matrix(rpkm.df[all.markers,]))
  #heatmap.2(all.markers.mat,margins = c(10,10),main='All markers expression',trace='none',cexRow = .3,cexCol = .3,col=bluered(10000),key.title = '')
  #heatmap.2(gene.dist.mat,margins = c(10,10),ColSideColors = gene.grp.col.side.col.list$col.str,RowSideColors = gene.cols.list$col.str,main='Genes euclidean clusters',trace='none',cexRow = .3,cexCol = .3,col=bluered(10000),key.title = 'Genes euclidian dist.')
  #plot.new()
  #legend('left',legend=gene.grp.col.side.col.list$legend.str,fill=gene.grp.col.side.col.list$legend.col,box.lty = 0,cex = .8)
  gene.groups.list=split(gene.groups.df,f=gene.groups.df$group)
  gene.groups.names=names( gene.groups.list)
  len.gene.groups.names=length(gene.groups.names)
  for(n in 1:len.gene.groups.names){
    gene.groups.name=gene.groups.names[n]
    temp.title.str=convert.to.title.case(gene.groups.name)
    temp.gene.groups.df=gene.groups.list[[gene.groups.name]]
    temp.genes.vec=rownames(temp.gene.groups.df)
    temp.out.mat=out.mat[temp.genes.vec,]
    temp.gene.cols.list=get.col.factor(col.factor = as.character(markers.df[temp.genes.vec,'category']))
    heatmap.2(x = temp.out.mat,margins = c(13,13),ColSideColors = col.side.col.str,main=temp.title.str,trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),key.title = 'Euclidean dist.',scale='row',RowSideColors = temp.gene.cols.list$col.str,dendrogram = 'col')
    #plot.new()
    legend('left',legend=temp.gene.cols.list$legend.str,fill=temp.gene.cols.list$legend.col,box.lty = 0,cex = .8)
    if(gene.groups.name!='Gene grp: 1'){
      genes.to.use=append(genes.to.use,temp.genes.vec,length(genes.to.use))
    }
    else{
      temp.gene.grp.out.mat=out.mat[temp.genes.vec,]
      gene.clust.grps=hclustfunc(x = distfunc(x = temp.gene.grp.out.mat))
      gene.grps.tree=cutree(tree = gene.clust.grps,k = 30)
      gene.col.list=get.col.factor(col.factor = as.numeric(gene.grps.tree))
      gene.grps.df=data.frame(name=names(gene.grps.tree),group=as.numeric(gene.grps.tree))
      rownames(gene.grps.df)=names(gene.grps.tree)
      #heatmap.2(x = temp.gene.grp.out.mat,main='Sub-group one genes',trace='none',col=bluered(10000),key.title = 'Euclidean dist.',scale='row',cexRow = .3,cexCol = .3,dendrogram = 'col',RowSideColors = gene.col.list$col.str)
      #legend('topright',legend =gene.col.list$legend.str,fill=gene.col.list$legend.col )
      row.annotation.col.df=subset(gene.grps.df,select = group)
      grp.one.genes.vec=rownames(subset(gene.grps.df,group!=1))
      genes.to.use=append(genes.to.use,grp.one.genes.vec,length(genes.to.use))
      #temp.gene.grp.out.mat=temp.gene.grp.out.mat[rownames(subset(gene.grps.df,group==1)),]
      temp.gene.grp.out.mat=temp.gene.grp.out.mat[grp.one.genes.vec,]
      #heatmap.2(x = temp.gene.grp.out.mat,main='Sub-group one selected genes',trace='none',col=bluered(10000),key.title = 'Euclidean dist.',scale='row',cexRow = .3,cexCol = .3, dendrogram = 'col')
      #legend('topright',legend =gene.col.list$legend.str,fill=gene.col.list$legend.col )
      #pheatmap(mat = temp.gene.grp.out.mat,main='Gene group one',col=bluered(10000),key.title = 'Euclidean dist.',treeheight_row = 1,fontsize_row = 1,fontsize_col = 1,cutree_rows = 5,cellwidth = 2,cellheight = 2,annotation_row =row.annotation.col.df,kmeans_k = 5 )
    }
  }
  final.out.mat=as.matrix(marker.category.score.df[genes.to.use,])
  sc.final.out.dist=dist(t(final.out.mat))
  sc.final.out.dist.mat=as.matrix(sc.final.out.dist)
  sc.final.samples.hclust=hclustfunc(x = sc.final.out.dist)
  sc.final.genes.hclust=hclustfunc(x = distfunc(x =final.out.mat))
  sc.final.groups <- cutree(sc.final.samples.hclust, k = in.sc.cut.off)
  sc.final.groups.df=data.frame(samples=names(sc.final.groups),sp=paste('SC grp: ',as.character(as.numeric(sc.final.groups)),sep=''),row.names = names(sc.final.groups))
  #sc.final.grp.legend.str=sort(unique(as.character(sc.final.groups.df$sp)))
  #sc.final.col.list=brewer.pal(n = as.numeric(in.sc.cut.off),'Set3')
  sc.final.col.list=get.color.brewer.list(in.vec = as.character(sc.final.groups.df$sp))
  #sc.final.genes.groups=cutree(tree = sc.final.genes.hclust,as.numeric(gene.cluster.cut.off))
  gene.categories.vec=as.character(markers.df[rownames(final.out.mat),'category'])
  gene.cols.list=get.col.factor(col.factor = gene.categories.vec)
  title.str=convert.to.title.case(title.str)
  samples=colnames(final.out.mat)
  stages.vec=as.character(meta.df[samples,]$development.stage)
  col.side.col.str=get.std.stage.cols(in.stages.vec=stages.vec)
  heatmap.2(x = final.out.mat,margins = c(13,13),ColSideColors = col.side.col.str,trace='none',main='Final gene list',cexRow = .1,cexCol = .3,col=bluered(10000),dendrogram = 'col',RowSideColors = gene.cols.list$col.str,scale='row', labRow = gene.categories.vec)
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
  legend('left',legend =gene.cols.list$legend.str,fill=gene.cols.list$legend.col ,cex=.8,box.lty = 0)
  heatmap.2(x = final.out.mat,margins = c(13,13),ColSideColors = sc.final.col.list$col.str,trace='none',main='Final gene list',cexRow = .1,cexCol = .3,col=bluered(10000),dendrogram = 'col',scale='row', RowSideColors = gene.cols.list$col.str,labCol = paste('Sub.grp.',as.character(as.numeric(sc.final.groups)),sep=''))
  legend('topright',legend = sc.final.col.list$legend.str,fill = sc.final.col.list$legend.col,box.lty = 0,cex = 1.0)
  annotation.df=subset(sc.final.groups.df,select = colnames(final.out.mat))
  pheatmap(mat=final.out.mat,cellwidth = 3,cellheight = 2,main='Final gene list',fontsize_row  = 2,fontsize_col = 2,treeheight_row = 1,cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',clustering_distance_rows ='euclidean')
  pheatmap(mat=sc.final.out.dist.mat,cellwidth = 3,cellheight = 2,main='Final gene list',fontsize_row  = 2,fontsize_col = 2,treeheight_row = 1,cutree_cols = in.sc.cut.off,cutree_rows = in.sc.cut.off ,clustering_distance_rows = 'correlation',clustering_distance_cols = 'correlation')
  meta.df=meta.df[colnames(final.out.mat),]
  sc.sub.pop.vec=gsub(pattern = ':',replacement = '',as.character(sc.final.groups.df[rownames(meta.df),'sp']))
  sc.sub.pop.vec=gsub(pattern = ' ',replacement = '.',x = sc.sub.pop.vec)
  meta.df$markers.cluster.groups=sc.sub.pop.vec
  in.markers.rpkm.df=subset(in.rpkm.df,select = rownames(meta.df))
  in.markers.rpkm.df=in.markers.rpkm.df[genes.to.use,]
  pseudo.value=min(in.markers.rpkm.df[in.markers.rpkm.df>0])/2
  in.markers.log.rpkm.df=log2(in.markers.rpkm.df+pseudo.value)
  pca.obj=prcomp(t(in.markers.log.rpkm.df))
  pca.score=pca.obj$x
  sample.col.side.col.list=get.color.brewer.list(in.vec = as.character(as.numeric(sc.final.groups)))
  plot(x=as.numeric(pca.obj$x[,1]),y=as.numeric(pca.obj$x[,2]),col=sample.col.side.col.list$col.str,pch=19,cex=2.5,main='PCA plots',xlab='PC1',ylab='PC2')
  legend('topright',legend=sample.col.side.col.list$legend.str,fill=sample.col.side.col.list$legend.col,box.lty = 0,cex = .8)
  pairs(x=pca.obj$x[,1:5],pch=19,col = sample.col.side.col.list$col.str,cex=.5)
  pseudo.value=min(in.rpkm.df[in.rpkm.df>0])/2
  in.log.rpkm.df=log2(in.rpkm.df+pseudo.value)
  plot.samples.pca.tsne.clustering(rpkm.df = in.log.rpkm.df,meta.df = meta.df,title.str = 'All genes')
  plot.samples.pca.tsne.clustering(rpkm.df =in.log.rpkm.df[genes.to.use,],meta.df =meta.df, title.str ='Markers only')
  out.mat=final.out.mat
  return(out.mat)
}
get.development.stage.markers.sample.clusters=function(rpkm.df,meta.df,markers.df,
                                                       title.str='Test',log.rpkm = F,
                                                       in.sc.cut.off=8,gene.cluster.cut.off=5,
                                                       pop.rpkm.df,pop.meta.df,
                                                       get.quantile.score=T,
                                                       markers.vec,
                                                       in.col.ramp=colorRampPalette(c("gray","red"))(1000)){
  gene.cols.list=get.col.factor(col.factor = 
                                  as.character(markers.df$category))
  raw.rpkm.df=rpkm.df
  rpkm.df=filter.none.expressed.samples(filter.genes.not.expressed.in.a.number.of.samples(rpkm.df,sample.count = 2))
  pseudo.val=min(rpkm.df[rpkm.df!=0])/2
  pseudo.val=1
  rpkm.df=rpkm.df+pseudo.val
  in.rpkm.df= rpkm.df
  pop.rpkm.df=filter.none.expressed.samples(
    filter.genes.not.expressed.in.a.number.of.samples(pop.rpkm.df,sample.count = 2))
  #rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(rpkm.df))
  all.genes=rownames(rpkm.df)
  samples=colnames(rpkm.df)
  no.samples=length(samples)
  #out.mat=matrix(nrow = len.marker.categories,ncol = no.samples)
  expr.markers.vec=intersect(all.genes,rownames(markers.df))
  markers.df=markers.df[expr.markers.vec,]
  all.markers=rownames(markers.df)
  len.markers=length(all.markers)
  marker.category.score.df=data.frame()
  if(get.quantile.score){
    if(log.rpkm){
      marker.category.score.df=get.gene.signature.score.per.cell(
        rpkm.df = rpkm.df, genes.sig.vec =all.markers,log.rpkm = T )
    }
    else{
      marker.category.score.df=get.gene.signature.score.per.cell(
        rpkm.df = rpkm.df, genes.sig.vec =all.markers,log.rpkm = F )
    }
  }
  else{
    show('Un -scored for markers')
    marker.category.score.df=rpkm.df[all.markers,]
    if(log.rpkm){
      #pseudo.value=as.numeric(min(marker.category.score.df[marker.category.score.df!=0]))/2
      marker.category.score.df=log10(marker.category.score.df)
    }
    else{
      marker.category.score.df=marker.category.score.df
    }
  }
  genes.to.use=intersect(markers.vec,expr.markers.vec)
  final.out.mat=as.matrix(marker.category.score.df[genes.to.use,])
  sc.final.out.dist=dist(t(final.out.mat),method = 'euclidean')
  sc.final.out.dist.mat=as.matrix(sc.final.out.dist)
  sc.final.samples.hclust=hclustfunc(x = sc.final.out.dist)
  sc.final.genes.hclust=hclustfunc(x = distfunc(x =final.out.mat))
  sc.final.groups <- cutree(sc.final.samples.hclust, k = in.sc.cut.off)
  sc.final.groups.df=data.frame(samples=names(sc.final.groups),sp=paste('SC grp: ',as.character(as.numeric(sc.final.groups)),sep=''),row.names = names(sc.final.groups))
  temp.sc.meta.df=meta.df[rownames(sc.final.groups.df),]
  temp.sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = temp.sc.meta.df)
  temp.sc.meta.df=subset(temp.sc.meta.df,select=c(development.stage,kmeans.8.clusters))
  sc.final.groups.df$development.stage=as.character(temp.sc.meta.df$development.stage)
  sc.final.groups.df$kmeans.8.clusters=as.character(temp.sc.meta.df$kmeans.8.clusters)
  sc.final.groups.df=sc.final.groups.df[order(sc.final.groups.df$kmeans.8.clusters),]
  #sc.final.grp.legend.str=sort(unique(as.character(sc.final.groups.df$sp)))
  #sc.final.col.list=brewer.pal(n = as.numeric(in.sc.cut.off),'Set3')
  sc.final.col.list=get.color.brewer.list(in.vec = as.character(sc.final.groups.df$sp))
  #sc.final.stage.col.list=get.abbrev.std.stage.cols
  #sc.final.genes.groups=cutree(tree = sc.final.genes.hclust,as.numeric(gene.cluster.cut.off))
  gene.categories.vec=as.character(markers.df[rownames(final.out.mat),'category'])
  title.str=convert.to.title.case(title.str)
  samples=colnames(final.out.mat)
  stages.vec=as.character(meta.df[samples,]$development.stage)
  col.side.col.str=get.std.stage.cols(in.stages.vec=stages.vec)
  #heatmap.2(x = final.out.mat,margins = c(13,13),ColSideColors = col.side.col.str,trace='none',main='Final gene list',cexRow = .1,cexCol = .3,col=bluered(10000),dendrogram = 'col',RowSideColors = gene.cols.list$col.str,scale='row', labRow = gene.categories.vec)
  #legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
  #legend('left',legend =gene.cols.list$legend.str,fill=gene.cols.list$legend.col ,cex=.8,box.lty = 0)
  #heatmap.2(x = final.out.mat,margins = c(13,13),ColSideColors = sc.final.col.list$col.str,trace='none',main='Final gene list',cexRow = .1,cexCol = .3,col=bluered(10000),dendrogram = 'col',RowSideColors = gene.cols.list$col.str,scale='row', labRow = gene.categories.vec)
  #legend('topright',legend = sc.final.col.list$legend.str,fill = sc.final.col.list$legend.col,box.lty = 0,cex = .5)
  #legend('left',legend =gene.cols.list$legend.str,fill=gene.cols.list$legend.col ,cex=.3,box.lty = 0)
  #heatmap.2(x = final.out.mat,margins = c(13,13),ColSideColors = sc.final.col.list$col.str,trace='none',main='Final gene list',cexRow = .1,cexCol = .3,col=bluered(10000),dendrogram = 'col',scale='row', RowSideColors = gene.cols.list$col.str,labCol = paste('Sub.grp.',as.character(as.numeric(sc.final.groups)),sep=''))
  #legend('left',legend =gene.cols.list$legend.str,fill=gene.cols.list$legend.col ,cex=.4,box.lty = 0)
  #legend('topright',legend = sc.final.col.list$legend.str,fill = sc.final.col.list$legend.col,box.lty = 0,cex = .5)
  batch.df=subset(meta.df,select=c(batch,spike,timepoint))
  final.out.annotation.df=subset(sc.final.groups.df,
                                 select = c(kmeans.8.clusters,sp,development.stage))
  border.gaps=as.character(final.out.annotation.df$kmeans.8.clusters)
  len.border.gaps=length(border.gaps)
  border.gaps.indices=c()
  for(m in 1:len.border.gaps){
    if(m!=len.border.gaps & border.gaps[m]!=border.gaps[m+1]){
      border.gaps.indices=append(border.gaps.indices,values = m,
                                 after = length(border.gaps.indices))
    }
  }
  final.out.annotation.df=subset(sc.final.groups.df,
                                 select = c(sp,development.stage))
  batch.df=batch.df[rownames(final.out.annotation.df),]
  #final.out.annotation.df=data.frame(final.out.annotation.df,batch.df)
  #final.out.annotation.df=data.frame(final.out.annotation.df,batch.df)
  sub.grps.vec=as.character(sc.final.groups.df$sp)
  final.out.annotation.col.list=list(
    sp=get.color.list.for.pheatmap(in.vec = sub.grps.vec),
    development.stage=abbrev.std.stages.col.vec)
  #ordered.final.out.annotation.df=final.out.annotation.df
  #ordered.final.out.annotation.df$sp=factor(final.out.annotation.df$sp, levels =sub.grps.vec )
  #ordered.final.out.mat=final.out.mat[,rownames(ordered.final.out.annotation.df)]
  row.annotation.df=markers.df[rownames(final.out.mat),]
  row.annotation.df=subset(row.annotation.df,select = category)
  #pheatmap(mat=final.out.mat,cellwidth = 2,cellheight = 2,main='',fontsize_row  = 1,fontsize_col = 1,treeheight_row = 1,cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',clustering_distance_rows ='euclidean',annotation_col = final.out.annotation.df,annotation_colors = final.out.annotation.col.list,annotation_row = row.annotation.df,annotation_legend = F,color = bluered(10000),labels_row =as.character(row.annotation.df[rownames(final.out.mat),]),labels_row = F,labels_col = F)
  col.pelette=colorRampPalette(c("darkgray","gray","white",
                                 "yellow","orange",'red'))(5)
  col.pelette=in.col.ramp
  # pheatmap(mat=final.out.mat[,rownames(final.out.annotation.df)],cellwidth = 2,
  #          cellheight = 2,main='',
  #          fontsize_row  = 1,fontsize_col = 1,treeheight_row = 10,gaps_col = border.gaps.indices,cluster_cols = F,
  #          cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',
  #          clustering_distance_rows ='euclidean',
  #          annotation_col = final.out.annotation.df,
  #          annotation_colors = final.out.annotation.col.list,
  #          annotation_row = row.annotation.df,annotation_legend = F,
  #          color = col.pelette,show_rownames = F,show_colnames = F)
  # 
  pheatmap(mat=final.out.mat,cellwidth = 2,cellheight = 2,main='',
           fontsize_row  = 1,fontsize_col = 1,treeheight_row = 10,
           cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',
           clustering_distance_rows ='euclidean',
           annotation_col = final.out.annotation.df,
           annotation_colors = final.out.annotation.col.list,
           annotation_row = row.annotation.df,annotation_legend = F,
           color = col.pelette,show_rownames = F,show_colnames = F)
  pheatmap(mat=final.out.mat,cellwidth = 2,cellheight = 2,main='',
           fontsize_row  = 1,fontsize_col = 1,treeheight_row = 1,
           cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',
           clustering_distance_rows ='correlation',
           annotation_col = final.out.annotation.df,
           annotation_colors = final.out.annotation.col.list,
           annotation_row = row.annotation.df,annotation_legend = F,
           color = col.pelette,show_rownames = F,show_colnames = T)
  #pheatmap(mat=final.out.mat,cellwidth = 2,cellheight = 2,main='',fontsize_row  = 1,fontsize_col = 1,treeheight_row = 1,cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',clustering_distance_rows ='euclidean',annotation_col = final.out.annotation.df,annotation_colors = final.out.annotation.col.list,annotation_row = row.annotation.df,annotation_legend = T,color = bluered(10000),labels_row =as.character(row.annotation.df[rownames(final.out.mat),]) ,labels_row = F,labels_col = F)
  # pheatmap(mat=final.out.mat,cellwidth = 2,cellheight = 2,main='',
  #          fontsize_row  = 1,fontsize_col = 1,treeheight_row = 1,
  #          cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',
  #          clustering_distance_rows ='euclidean',
  #          annotation_col = final.out.annotation.df,
  #          annotation_colors = final.out.annotation.col.list,
  #          annotation_row = row.annotation.df,annotation_legend = T,
  #          color = bluered(10000),show_rownames = F,show_colnames = F)
  # 
  #pheatmap(mat=final.out.mat,cellwidth = 2,cellheight = 2,main='Final gene list',fontsize_row  = 2,fontsize_col = 2,treeheight_row = 1,cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',clustering_distance_rows ='euclidean',annotation_col = final.out.annotation.df,annotation_row = row.annotation.df,color = bluered(10000),labels_row = F,labels_col = F)
  # pheatmap(mat=final.out.mat,cellwidth = 2,cellheight = 2,main='Final gene list',
  #          fontsize_row  = 2,fontsize_col = 2,treeheight_row = 1,
  #          cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',
  #          clustering_distance_rows ='euclidean',
  #          annotation_col = final.out.annotation.df,
  #          annotation_row = row.annotation.df,color = bluered(10000),
  #          show_rownames = F,show_colnames = F)
  # 
  #pheatmap(mat=sc.final.out.dist.mat,cellwidth = 3,cellheight = 2,main='Final gene list',fontsize_row  = 2,fontsize_col = 2,treeheight_row = 1,cutree_cols = in.sc.cut.off,cutree_rows = in.sc.cut.off ,clustering_distance_rows = 'correlation',clustering_distance_cols = 'correlation',labels_row = F,labels_col = F)
  # pheatmap(mat=sc.final.out.dist.mat,cellwidth = 3,cellheight = 2,
  #          main='Final gene list',fontsize_row  = 2,fontsize_col = 2,
  #          treeheight_row = 1,cutree_cols = in.sc.cut.off,
  #          cutree_rows = in.sc.cut.off ,clustering_distance_rows = 'correlation',
  #          clustering_distance_cols = 'correlation',show_rownames = F,
  #          show_colnames = F)
  # 
  meta.df=meta.df[colnames(final.out.mat),]
  sc.sub.pop.vec=gsub(pattern = ':',replacement = '',as.character(sc.final.groups.df[rownames(meta.df),'sp']))
  sc.sub.pop.vec=gsub(pattern = ' ',replacement = '.',x = sc.sub.pop.vec)
  meta.df$markers.cluster.groups=sc.sub.pop.vec
  in.markers.rpkm.df=subset(in.rpkm.df,select = rownames(meta.df))
  in.markers.rpkm.df=in.markers.rpkm.df[genes.to.use,]
  pseudo.value=min(in.markers.rpkm.df[in.markers.rpkm.df>0])/2
  #in.markers.log.rpkm.df=log2(in.markers.rpkm.df+pseudo.value)
  in.markers.log.rpkm.df=log2(in.markers.rpkm.df)
  in.markers.log.rpkm.mat=as.matrix(in.markers.log.rpkm.df)
  # pheatmap(mat=in.markers.log.rpkm.mat,cellwidth = 3,cellheight = 2,
  #          main='Final gene list',fontsize_row  = 1,fontsize_col = 1,
  #          treeheight_row = 1,cutree_cols = in.sc.cut.off,
  #          clustering_distance_cols = 'euclidean',
  #          clustering_distance_rows ='euclidean',
  #          annotation_col = final.out.annotation.df,color = bluered(10000))
  # 
  #plot.samples.pca.tsne.clustering(rpkm.df =in.markers.rpkm.df,log.rpkms = F,
                                   #meta.df = meta.df )
  pca.obj=prcomp(t(in.markers.log.rpkm.df))
  pca.score=pca.obj$x
  #get.color.list.for.pheatmap(in.vec = sub.grps.vec))
  pca.col.str=as.character(meta.df[rownames(pca.score),]$markers.cluster.groups)
  sample.col.side.col.list=get.color.list.for.pheatmap(in.vec = pca.col.str)
  col.vec=as.character(sample.col.side.col.list[pca.col.str])
  # #get.color.list.for.pheatmap(in.vec = )
  # 
  # plot(x=as.numeric(pca.obj$x[,1]),y=as.numeric(pca.obj$x[,2]),col=col.vec,pch=19,cex=3.0,main='PCA plots',xlab='PC1',ylab='PC2')
  # 
  # text(x = pca.obj$x[,1],y= pca.obj$x[,2],labels = rownames(pca.obj$x),cex = .1)
  # 
  # legend('topright',legend=names(sample.col.side.col.list),fill=as.character(sample.col.side.col.list),box.lty = 0,cex = .8)
  # 
  # plot(x=as.numeric(pca.obj$x[,1]),y=as.numeric(pca.obj$x[,2]),col=col.vec,pch=19,cex=3.0,main='',xlab='',ylab='')
  # 
  # stage.col.vec=get.std.stage.cols(in.stages.vec = as.character(meta.df[rownames(pca.score),]$development.stage))
  # 
  # plot(x=as.numeric(pca.obj$x[,1]),y=as.numeric(pca.obj$x[,2]),col=stage.col.vec,pch=19,cex=3.0,main='',xlab='',ylab='')
  # 
  # legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .8)
  # 
  # plot(x=as.numeric(pca.obj$x[,1]),y=as.numeric(pca.obj$x[,2]),col=stage.col.vec,pch=19,cex=3.0,main='',xlab='',ylab='')
  # 
  # pairs(x=pca.obj$x[,4:7],pch=19,col = col.vec,cex=1.5)
  schizont.stage.sub.group.samples.df=filter.genes.not.expressed.in.a.number.of.samples(df = filter.none.expressed.genes( df = raw.rpkm.df[,schizont.stage.sub.group.samples.vec]),sample.count = 2)
  schizont.var.cut.off=as.numeric(quantile(x = as.numeric(apply(schizont.stage.sub.group.samples.df,1,sd)),probs = .3))
  schizont.stage.sub.group.samples.df=filter.non.variable.rows(df =schizont.stage.sub.group.samples.df,cut.off =schizont.var.cut.off)
  sub.grps.vec=as.character(sc.final.groups.df$sp)
  final.out.annotation.col.list=list(sp=get.color.list.for.pheatmap(in.vec = sub.grps.vec))
#   schizont.stage.sub.group.samples.cor.mat=cor(schizont.stage.sub.group.samples.df,method = 'spearman')
#   
#   pheatmap(mat=schizont.stage.sub.group.samples.cor.mat,cellwidth = 10,cellheight = 10,main='',fontsize_row  = 1,fontsize_col = 1,cutree_cols = 3,annotation_col = final.out.annotation.df,annotation_colors = final.out.annotation.col.list,annotation_legend = F,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
#   
#   pheatmap(mat=schizont.stage.sub.group.samples.cor.mat,cellwidth = 10,cellheight = 10,main='',fontsize_row  = 1,fontsize_col = 1,cutree_cols = 3,annotation_col = final.out.annotation.df,annotation_colors = final.out.annotation.col.list,annotation_legend = T,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
#  
#   schizont.stage.sub.group.samples.df=filter.genes.not.expressed.in.a.number.of.samples(df = filter.none.expressed.genes( df = rpkm.df[,schizont.stage.sub.group.samples.vec]),sample.count = 2)
#   
#   schizont.var.cut.off=as.numeric(quantile(x = as.numeric(apply(schizont.stage.sub.group.samples.df,1,sd)),probs = .3))
#   
#   schizont.stage.sub.group.samples.df=filter.non.variable.rows(df =schizont.stage.sub.group.samples.df,cut.off =schizont.var.cut.off)
#   
#   sub.grps.vec=as.character(sc.final.groups.df$sp)
#   
#   final.out.annotation.col.list=list(sp=get.color.list.for.pheatmap(in.vec = sub.grps.vec))
#   
#   schizont.stage.sub.group.samples.cor.mat=cor(schizont.stage.sub.group.samples.df,method = 'spearman')
#   
#   pheatmap(mat=schizont.stage.sub.group.samples.cor.mat,cellwidth = 10,cellheight = 10,main='',fontsize_row  = 1,fontsize_col = 1,cutree_cols = 3,annotation_col = final.out.annotation.df,annotation_colors = final.out.annotation.col.list,annotation_legend = F,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
#   
#   pheatmap(mat=schizont.stage.sub.group.samples.cor.mat,cellwidth = 10,cellheight = 10,main='',fontsize_row  = 1,fontsize_col = 1,cutree_cols = 3,annotation_col = final.out.annotation.df,annotation_colors = final.out.annotation.col.list,annotation_legend = T,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
#   
#   
#   late.ring.sample.vec=rownames(subset(meta.df,development.stage=='late.ring'))
#   
#   lr.stage.sub.group.samples.df=filter.genes.not.expressed.in.a.number.of.samples(df = filter.none.expressed.genes( df = raw.rpkm.df[,late.ring.sample.vec]),sample.count = 2)
#   
#   lr.var.cut.off=as.numeric(quantile(x = as.numeric(apply(lr.stage.sub.group.samples.df,1,sd)),probs = .3))
#   
#   lr.stage.sub.group.samples.df=filter.non.variable.rows(df =lr.stage.sub.group.samples.df,cut.off =lr.var.cut.off)
#   
#   sub.grps.vec=as.character(sc.final.groups.df$sp)
#   
#   final.out.annotation.col.list=list(sp=get.color.list.for.pheatmap(in.vec = sub.grps.vec))
#   
#   lr.stage.sub.group.samples.cor.mat=cor(lr.stage.sub.group.samples.df,method = 'spearman')
#   
#   pheatmap(mat=lr.stage.sub.group.samples.cor.mat,cellwidth = 10,cellheight = 10,main='',fontsize_row  = 1,fontsize_col = 1,cutree_cols = 3,annotation_col = final.out.annotation.df,annotation_colors = final.out.annotation.col.list,annotation_legend = F,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
#   
#   pheatmap(mat=lr.stage.sub.group.samples.cor.mat,cellwidth = 10,cellheight = 10,main='LR',fontsize_row  = 1,fontsize_col = 1,cutree_cols = 3,annotation_col = final.out.annotation.df,annotation_colors = final.out.annotation.col.list,annotation_legend = T,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
#   
  #get.development.stage.heterogeneity.based.on.normalized.counts(rpkm.df = schizont.stage.sub.group.samples.df,meta.df = meta.df)
  #pseudo.value=min(in.rpkm.df[in.rpkm.df>0])/2
  #in.log.rpkm.df=log2(in.rpkm.df+pseudo.value)
  in.log.rpkm.df=log2(in.rpkm.df)
  all.genes.pca.obj=prcomp(t(in.log.rpkm.df))
  all.genes.pca.obj.score=all.genes.pca.obj$x
  #get.color.list.for.pheatmap(in.vec = sub.grps.vec))
  #all.genes.pca.col.str=as.character(meta.df[rownames(all.genes.pca.obj.score),]$markers.cluster.groups)
  #all.genes.sample.col.side.col.list=get.color.list.for.pheatmap(in.vec = all.genes.pca.col.str)
  #all.genes.col.vec=as.character(all.genes.sample.col.side.col.list[all.genes.pca.col.str])
  #get.color.list.for.pheatmap(in.vec = )
  #plot(x=as.numeric(all.genes.pca.obj.score$x[,1]),y=as.numeric(all.genes.pca.obj.score$x[,2]),col=all.genes.col.vec,pch=19,cex=3.0,main='All genes PCA plots',xlab='PC1',ylab='PC2')
  #legend('topright',legend=names(all.genes.sample.col.side.col.list),fill=as.character(all.genes.sample.col.side.col.list),box.lty = 0,cex = .8)
  #plot(x=as.numeric(all.genes.pca.obj.score$x[,1]),y=as.numeric(all.genes.pca.obj.score$x[,2]),col=all.genes.col.vec,pch=19,cex=3.0,main='',xlab='',ylab='')
  #plot.new()
  #legend('center',legend=rep('',times = length(names(all.genes.sample.col.side.col.list))),fill=as.character(all.genes.sample.col.side.col.list),box.lty = 0,cex = 1.5)
  #pairs(x=all.genes.pca.obj.score$x[,1:5],pch=19,col = col.vec,cex=.5)
  #plot.samples.pca.tsne.clustering(rpkm.df = in.log.rpkm.df,meta.df = meta.df,title.str = 'All genes')
  #plot.samples.pca.tsne.clustering(rpkm.df =in.log.rpkm.df[genes.to.use,],meta.df =meta.df, title.str ='Markers only')
  #Compare to bulk-classes
  #bulk.intervals=seq(3,10)
  bulk.intervals=seq(9,9)
  len.bulk.intervals=length(bulk.intervals)
  for(d in 1:len.bulk.intervals){
    #pop.pseudo.val=as.numeric(min(pop.rpkm.df[pop.rpkm.df>0]))/2
    #log.pop.rpkm.df=log2(pop.rpkm.df+pop.pseudo.val)
    log.pop.rpkm.df=log2(pop.rpkm.df)
    pop.meta.df=pop.meta.df[colnames(log.pop.rpkm.df),]
    sub.grp.cut.off=bulk.intervals[d]
    log.pop.cor.mat=cor(log.pop.rpkm.df,method = 'spearman')
    log.pop.cor.dist.hclust=hclustfunc(x = distfunc(x = log.pop.cor.mat))
    pop.grp.list=cutree(tree =log.pop.cor.dist.hclust,sub.grp.cut.off )
    pop.grp.df=data.frame(samples=names(pop.grp.list),group=as.character(as.numeric(pop.grp.list)),row.names =names(pop.grp.list) )
    bulk.stages.vec=pop.meta.df[rownames(pop.grp.df),]
    bulk.stages.vec=as.character(bulk.stages.vec$development.stage)
    pop.grp.df$stage=bulk.stages.vec
    title.str=convert.to.title.case(in.str = paste('Bulk correlation (all genes, interval=',sub.grp.cut.off,')',sep=''))
    # pheatmap(mat=log.pop.cor.mat,cellwidth = 5,cellheight = 5,main= title.str,fontsize_row  = 2,fontsize_col = 2,treeheight_row = 1,cutree_cols = sub.grp.cut.off ,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
    # 
    temp.meta.list=split(meta.df,f = meta.df$markers.cluster.groups)
    sc.mean.rpkm.list=lapply(temp.meta.list,function(list.item){
      list.item=data.frame(list.item)
      temp.df=as.numeric(apply(subset(in.rpkm.df,select = rownames(list.item)),1,mean))
      return(temp.df)
    })
    sc.mean.rpkm.df=data.frame(sc.mean.rpkm.list)
    colnames(sc.mean.rpkm.df)=names(temp.meta.list)
    rownames(sc.mean.rpkm.df)=rownames(in.rpkm.df)
    intersect.genes=intersect(rownames(sc.mean.rpkm.df),rownames(pop.rpkm.df))
    sc.mean.rpkm.df=sc.mean.rpkm.df[intersect.genes,]
    first.samples=colnames(sc.mean.rpkm.df)
    pop.rpkm.df=pop.rpkm.df[intersect.genes,]
    second.samples=colnames(pop.rpkm.df)
    combn.df=data.frame(sc.mean.rpkm.df,pop.rpkm.df)
    colnames(combn.df)=c(first.samples,second.samples)
    cor.mat=data.frame(cor(combn.df,method = 'spearman'))
    # plot.cor.mat=as.matrix(cor.mat[first.samples,second.samples])
    # 
    bulk.annotation.df=data.frame(t(subset(t(pop.grp.df),select = second.samples)))
    bulk.annotation.df=subset(pop.grp.df,select = c(group,stage))
    # pheatmap(mat=plot.cor.mat,cellwidth = 6,cellheight = 10,main= 'Euclidean distance',fontsize_row  = 2,fontsize_col = 2,treeheight_row = 1 ,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean',annotation_col = bulk.annotation.df)
    # 
    # pheatmap(mat=plot.cor.mat,cellwidth = 6,cellheight = 10,main= 'Correlation distance',fontsize_row  = 2,fontsize_col = 2,treeheight_row = 1 ,clustering_distance_rows = 'correlation',clustering_distance_cols = 'correlation',annotation_col = bulk.annotation.df)
    # 
    #in.markers.rpkm.df=subset(in.rpkm.df,select = rownames(meta.df))
    #in.markers.rpkm.df=in.markers.rpkm.df[genes.to.use,]
    log.markers.pop.rpkm.df=log.pop.rpkm.df[genes.to.use,]
    pop.meta.df=pop.meta.df[colnames(log.markers.pop.rpkm.df),]
    log.markers.pop.cor.mat=cor(log.markers.pop.rpkm.df,method = 'spearman')
    log.markers.pop.cor.mat.dist.hclust=hclustfunc(x = distfunc(x = log.markers.pop.cor.mat))
    markers.pop.grp.list=cutree(tree =log.markers.pop.cor.mat.dist.hclust,sub.grp.cut.off )
    markers.pop.grp.df=data.frame(samples=names(markers.pop.grp.list),group=as.character(as.numeric(markers.pop.grp.list)),row.names =names(markers.pop.grp.list) )
    markers.bulk.stages.vec=pop.meta.df[rownames(markers.pop.grp.df),]
    markers.bulk.stages.vec=as.character(markers.bulk.stages.vec$development.stage)
    markers.pop.grp.df$stage=markers.bulk.stages.vec
    title.str=convert.to.title.case(in.str = paste('Bulk correlation (Marker genes, interval=',sub.grp.cut.off,')',sep=''))
    # pheatmap(mat=log.markers.pop.cor.mat,cellwidth = 5,cellheight = 5,main= title.str,fontsize_row  = 2,fontsize_col = 2,treeheight_row = 1,cutree_cols = sub.grp.cut.off ,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean')
    # 
    temp.meta.list=split(meta.df,f = meta.df$markers.cluster.groups)
    markers.sc.mean.rpkm.list=lapply(temp.meta.list,function(list.item){
      temp.df=subset(in.rpkm.df,select = rownames(list.item))
      temp.df=temp.df[genes.to.use,]
      temp.df=as.numeric(apply(temp.df,1,mean))
      return(temp.df)
    })
    marker.sc.mean.rpkm.df=data.frame(markers.sc.mean.rpkm.list)
    colnames(marker.sc.mean.rpkm.df)=names(temp.meta.list)
    rownames(marker.sc.mean.rpkm.df)=genes.to.use
    markers.intersect.genes=intersect(rownames(marker.sc.mean.rpkm.df),rownames(log.markers.pop.rpkm.df))
    marker.sc.mean.rpkm.df=marker.sc.mean.rpkm.df[markers.intersect.genes,]
    first.samples=colnames(marker.sc.mean.rpkm.df)
    pop.rpkm.df=pop.rpkm.df[markers.intersect.genes,]
    second.samples=colnames(pop.rpkm.df)
    combn.df=data.frame(marker.sc.mean.rpkm.df,pop.rpkm.df)
    colnames(combn.df)=c(first.samples,second.samples)
    cor.mat=data.frame(cor(combn.df,method = 'spearman'))
    plot.cor.mat=as.matrix(cor.mat[first.samples,second.samples])
    bulk.annotation.df=data.frame(t(subset(t(markers.pop.grp.df),select = second.samples)))
    bulk.annotation.df=subset(bulk.annotation.df,select = c(group,stage))
    # pheatmap(mat=plot.cor.mat,cellwidth = 6,cellheight = 10,main= 'Euclidean distance',fontsize_row  = 2,fontsize_col = 2,treeheight_row = 1 ,clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean',annotation_col = bulk.annotation.df)
    # 
    # pheatmap(mat=plot.cor.mat,cellwidth = 6,cellheight = 10,main= 'Correlation distance',fontsize_row  = 2,fontsize_col = 2,treeheight_row = 1 ,clustering_distance_rows = 'correlation',clustering_distance_cols = 'correlation',annotation_col = bulk.annotation.df)
  }
  out.list=list(markers.score=final.out.mat,meta.df=meta.df,rpkm=in.rpkm.df)
  return(out.list)
}
get.development.stage.markers.sample.clusters.based.on.limited.plots=function(rpkm.df,meta.df,markers.df,
                                                                              title.str='Test',log.rpkm = F,
                                                                              in.sc.cut.off=6,
                                                                              gene.cluster.cut.off=5,
                                                                              pop.rpkm.df,pop.meta.df,
                                                                              get.quantile.score=T,
                                                                              markers.vec,
                                                                              in.col.ramp=colorRampPalette(c("darkgreen",'red2'))(1000)){
  sub.pop.batch.col.list=get.color.rainbow.list.for.pheatmap(in.vec =as.character(meta.df$batch ))
  sub.pop.spike.col.list=get.color.rainbow.list.for.pheatmap(in.vec =as.character(meta.df$spike ))
  sub.pop.timepoint.col.list=get.color.rainbow.list.for.pheatmap(in.vec =as.character(meta.df$timepoint))
  sub.pop.rosetting.col.list=get.color.rainbow.list.for.pheatmap(in.vec =as.character(meta.df$Rosetting))
  #gene.categories.col.list=get.color.rainbow.list.for.pheatmap(in.vec =as.character(markers.df$category))
  #gene.categories.col.list=get.distinct.colour.list.for.pheatmap(in.vec = as.character(markers.df$category))
  #print(gene.categories.col.list)
  gene.categories.col.list=gene.cat.col.names
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  raw.rpkm.df=rpkm.df
  rpkm.df=filter.none.expressed.samples(filter.genes.not.expressed.in.a.number.of.samples(rpkm.df,sample.count = 2))
  pseudo.val=min(rpkm.df[rpkm.df!=0])/2
  rpkm.df=rpkm.df+pseudo.val
  in.rpkm.df= rpkm.df
  pop.rpkm.df=filter.none.expressed.samples(filter.genes.not.expressed.in.a.number.of.samples(pop.rpkm.df,sample.count = 2))
  #rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(rpkm.df))
  all.genes=rownames(rpkm.df)
  samples=colnames(rpkm.df)
  no.samples=length(samples)
  #out.mat=matrix(nrow = len.marker.categories,ncol = no.samples)
  expr.markers.vec=intersect(all.genes,rownames(markers.df))
  markers.df=markers.df[expr.markers.vec,]
  all.markers=rownames(markers.df)
  len.markers=length(all.markers)
  marker.category.score.df=data.frame()
  if(get.quantile.score){
    if(log.rpkm){
      marker.category.score.df=get.gene.signature.score.per.cell(rpkm.df = rpkm.df, 
                                                                 genes.sig.vec =all.markers,
                                                                 log.rpkm = T )
    }
    else{
      marker.category.score.df=get.gene.signature.score.per.cell(rpkm.df = rpkm.df, 
                                                                 genes.sig.vec =all.markers,
                                                                 log.rpkm = F )
    }
  }
  # else{
  #   
  #   show('Un -scored for markers')
  #   
  #   marker.category.score.df=rpkm.df[all.markers,]
  #   
  #   if(log.rpkm){
  #     
  #     #pseudo.value=as.numeric(min(marker.category.score.df[marker.category.score.df!=0]))/2
  #     
  #     marker.category.score.df=log2(marker.category.score.df)
  #     
  #   }
  #   
  #   else{
  #     
  #     marker.category.score.df=marker.category.score.df
  #   }
  #   
  # }
  genes.to.use=intersect(markers.vec,expr.markers.vec)
  final.out.mat=as.matrix(marker.category.score.df[genes.to.use,])
  sc.final.out.dist=dist(t(final.out.mat),method = 'euclidean')
  sc.final.out.dist.mat=as.matrix(sc.final.out.dist)
  sc.final.samples.hclust=hclustfunc(x = sc.final.out.dist)
  sc.final.genes.hclust=hclustfunc(x = distfunc(x =final.out.mat))
  sc.final.groups <- cutree(sc.final.samples.hclust, k = in.sc.cut.off)
  sc.final.groups.df=data.frame(samples=names(sc.final.groups),
                                sp=paste('SC grp: ',as.character(as.numeric(sc.final.groups)),sep=''),
                                row.names = names(sc.final.groups))
  sc.final.col.list=get.color.brewer.list(in.vec = as.character(sc.final.groups.df$sp))
  gene.categories.vec=as.character(markers.df[rownames(final.out.mat),'category'])
  gene.cols.list=get.col.factor(col.factor = gene.categories.vec)
  temp.gene.categories.col.list=gene.categories.col.list[unique(gene.categories.vec)]
  title.str=convert.to.title.case(title.str)
  samples=colnames(final.out.mat)
  stages.vec=as.character(meta.df[samples,]$development.stage)
  col.side.col.str=get.std.stage.cols(in.stages.vec=stages.vec)
  batch.df=subset(meta.df,select=c(batch,spike,timepoint,Rosetting,development.stage))
  #batch.df=subset(meta.df,select=c(sampling.time,development.stage))
  final.out.annotation.df=subset(sc.final.groups.df,select = sp)
  batch.df=batch.df[rownames(final.out.annotation.df),]
  final.out.annotation.df=data.frame(final.out.annotation.df,batch.df)
  sub.grps.vec=as.character(sc.final.groups.df$sp)
  batch.col.list=sub.pop.batch.col.list[unique(batch.df$batch)]
  spike.col.list=sub.pop.spike.col.list[unique(batch.df$spike)]
  tm.col.list=sub.pop.timepoint.col.list[unique(batch.df$timepoint)]
  rosetting.col.list=sub.pop.rosetting.col.list[unique(batch.df$Rosetting)]
  # final.out.annotation.col.list=list(sp=get.color.list.for.pheatmap(in.vec = sub.grps.vec),
  #                                    category=temp.gene.categories.col.list,
  #                                    batch=batch.col.list,spike=spike.col.list,
  #                                    timepoint=tm.col.list,
  #                                    Rosetting=rosetting.col.list)
  #abbrev.std.stages.col.vec
  final.out.annotation.col.list=list(sp=get.color.list.for.pheatmap(in.vec = sub.grps.vec),
                                     category=temp.gene.categories.col.list,
                                     batch=batch.col.list,spike=spike.col.list,timepoint=tm.col.list,
                                     Rosetting=rosetting.col.list,development.stage=abbrev.std.stages.col.vec)
  row.annotation.df=markers.df[rownames(final.out.mat),]
  row.annotation.df=subset(row.annotation.df,select = c(category,category.abbrev))
  sp.only.final.out.annotation.df=subset(final.out.annotation.df,select = c(sp,development.stage))
  pheatmap(mat=final.out.mat,cellwidth = 1,cellheight = 1.5,main='',fontsize_row  = 1,fontsize_col = 1,
           treeheight_row = 1,cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',
           clustering_distance_rows ='euclidean',annotation_col = sp.only.final.out.annotation.df,
           annotation_colors = final.out.annotation.col.list,annotation_row = row.annotation.df,
           annotation_legend = F,color = in.col.ramp,show_rownames = T,show_colnames = F,fontsize=3)
  pheatmap(mat=final.out.mat,cellwidth = 1,cellheight = 1.5,main='',fontsize_row  = 1,fontsize_col = 1,
           treeheight_row = 1,cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',
           clustering_distance_rows ='euclidean',annotation_col = final.out.annotation.df,
           annotation_colors = final.out.annotation.col.list,annotation_row = row.annotation.df,
           annotation_legend = F,color = in.col.ramp,show_rownames = T,show_colnames = F,fontsize=3)
  temp.df=subset(row.annotation.df,select = category)
  rownames(final.out.mat)=as.character(temp.df[rownames(final.out.mat),'category'])
  pheatmap(mat=final.out.mat,cellwidth = 1,cellheight = 1.5,main='',fontsize_row  = 1,fontsize_col = 1,
           treeheight_row = 1,cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',
           clustering_distance_rows ='euclidean',annotation_col = sp.only.final.out.annotation.df,
           annotation_colors = final.out.annotation.col.list,annotation_row = temp.df,
           annotation_legend = T,color = in.col.ramp,show_rownames = T,show_colnames = T,fontsize=3)
  meta.df=meta.df[colnames(final.out.mat),]
  sc.sub.pop.vec=gsub(pattern = ':',replacement = '',as.character(sc.final.groups.df[rownames(meta.df),'sp']))
  sc.sub.pop.vec=gsub(pattern = ' ',replacement = '.',x = sc.sub.pop.vec)
  meta.df$markers.cluster.groups=sc.sub.pop.vec
  in.markers.rpkm.df=subset(in.rpkm.df,select = rownames(meta.df))
  in.markers.rpkm.df=in.markers.rpkm.df[genes.to.use,]
  pseudo.value=min(in.markers.rpkm.df[in.markers.rpkm.df>0])/2
  #in.markers.log.rpkm.df=log2(in.markers.rpkm.df+pseudo.value)
  in.markers.log.rpkm.df=log2(in.markers.rpkm.df)
  in.markers.log.rpkm.mat=as.matrix(in.markers.log.rpkm.df)
  pca.obj=prcomp(t(in.markers.log.rpkm.df))
  pca.score=pca.obj$x
  samples.pca.summary=summary(pca.obj)$importance
  first.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,'PC1']))*100,2)
  second.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,'PC2']))*100,2)
  first.pc.explained.var=paste('(Explained var: ',paste(first.pc.explained.var,'%)',sep=''),sep='')
  second.pc.explained.var=paste('(Explained var: ',paste(second.pc.explained.var,'%)',sep=''),sep='')
  first.pc.lab=paste('PC1',first.pc.explained.var,'')
  second.pc.lab=paste('PC2',second.pc.explained.var,'')
  #get.color.list.for.pheatmap(in.vec = sub.grps.vec))
  pca.col.str=as.character(meta.df[rownames(pca.score),]$markers.cluster.groups)
  #sample.col.side.col.list=get.color.list.for.pheatmap(in.vec = pca.col.str)
  sample.col.side.col.list=get.color.rainbow.list.for.pheatmap(in.vec = pca.col.str)
  col.vec=as.character(sample.col.side.col.list[pca.col.str])
  #get.color.list.for.pheatmap(in.vec = )
  #plot(x=as.numeric(pca.obj$x[,1]),y=as.numeric(pca.obj$x[,2]),col=col.vec,pch=19,cex=3.0,main='',xlab='',ylab='')
  #plot(x=as.numeric(pca.obj$x[,1]),y=as.numeric(pca.obj$x[,2]),col=col.vec,pch=19,cex=3.0,main='',xlab=first.pc.lab,ylab=second.pc.lab)
  #stage.col.vec=get.std.stage.cols(in.stages.vec = as.character(meta.df[rownames(pca.score),]$development.stage))
  #schizont.stage.sub.group.samples.df=filter.genes.not.expressed.in.a.number.of.samples(df = filter.none.expressed.genes(input.data  = raw.rpkm.df[,schizont.stage.sub.group.samples.vec]),sample.count = 2)
  #schizont.var.cut.off=as.numeric(quantile(x = as.numeric(apply(schizont.stage.sub.group.samples.df,1,sd)),probs = .3))
  #schizont.stage.sub.group.samples.df=filter.non.variable.rows(df =schizont.stage.sub.group.samples.df,cut.off =schizont.var.cut.off)
  #sub.grps.vec=as.character(sc.final.groups.df$sp)
  #final.out.annotation.col.list=list(sp=get.color.list.for.pheatmap(in.vec = sub.grps.vec))
  #in.log.rpkm.df=log2(in.rpkm.df)
  #all.genes.pca.obj=prcomp(t(in.log.rpkm.df))
  #all.genes.pca.obj.score=all.genes.pca.obj$x
  out.list=list(markers.score=final.out.mat,meta.df=meta.df,rpkm=in.rpkm.df)
  return(out.list)
}
plot.smilr.clusters=
  function(simlr.res,meta.df,kmeans.centers=2,title.str='SIMILR 2D visualization'){
    meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
    samples=intersect(rownames(simlr.res$ydata),rownames(meta.df))
    #col.str=as.character(meta.df[samples,'markers.cluster.groups'])
    col.str=as.character(meta.df[samples,'development.stage'])
    #smilr.col.vec=get.color.list.for.pheatmap(in.vec  = col.str)
    #smilr.col.str=smilr.col.vec[col.str]
    smilr.col.str=abbrev.std.stages.col.vec[col.str]
    plot(simlr.res$ydata,cex=1.0,frame.plot=T,
         col=smilr.col.str,
         xlab="SIMLR component 1", ylab="SIMLR component 2",pch=19,
         main=title.str)
    #Plot the K-means results
    k.means.title.str='Kmeans clusters'
    dist.mat=as.matrix(dist(simlr.res$ydata))
    kmeans.cl=kmeans(dist.mat,centers = kmeans.centers)
    kmeans.col.vec=kmeans.cl$cluster[rownames(simlr.res$ydata)]
    kmean.clusters=kmeans.cl$cluster[rownames(simlr.res$ydata)]
    names(kmean.clusters)=rownames(simlr.res$ydata)
    #plot(simlr.res$ydata, col = kmeans.col.vec,pch=20,main=k.means.title.str,frame.plot=F,cex=1.0,
         #xlab="SIMLR component 1", ylab="SIMLR component 2")
    kmeans.col.str=length(unique(simlr.res$y$cluster))
    #points(kmeans.cl$centers, col = 1:kmeans.col.str, pch = 8, cex = 1)
    return(kmean.clusters)
  }
run.simlr.clustering.at.different.cluster.cut.off=
  function(expression.df,meta.df){
    clusters.cut.offs=seq(6,6)
    out.list=list()
    for (m in clusters.cut.offs ){
        show(paste('Running: ',m, ' clusters, ',sep = ))
        temp.title.str=paste('SIMLR(','Clusters:',m,')',sep='')
        out.list[[paste('clusters_',m,sep='')]]=
          run.simlr.clustering(rpkm.df=expression.df,meta.df=meta.df,norm = T,
                               number.cluster = m)
        show(paste('Finished: ',m, ' clusters, ',sep = ))
      }
    return(out.list)
  }
run.simlr.clustering=function(rpkm.df,meta.df,norm=T,plot=T,
                              number.cluster=2){
  pseudo.value=min(rpkm.df[rpkm.df>0])/2
  log.rpkm.df=log2(rpkm.df+pseudo.value)
  simlr.norm.res=SIMLR(X = log.rpkm.df,c = number.cluster,no.dim = 3,
                       normalize = norm)
  rownames(simlr.norm.res$ydata)=colnames(rpkm.df)
  simlr.unnorm.res=SIMLR(X = log.rpkm.df,c = number.cluster,normalize = F)
  rownames(simlr.unnorm.res$ydata)=colnames(rpkm.df)
  if(plot){
    plot.smilr.clusters(simlr.res = simlr.norm.res,meta.df = meta.df,
                        title.str = 'Log-expression norm SIMLR plot')
    plot.smilr.clusters(simlr.res = simlr.unnorm.res,meta.df = meta.df,
                        title.str = 'Log-expression un-norm SIMLR plot')
  }
  return(simlr.norm.res)
}
get.development.stage.markers.sample.clusters.at.different.cluster.cut.off=
  function(rpkm.df,meta.df,markers.df,title.str='Test',log.rpkm = F,
           in.sc.cut.off=6,gene.cluster.cut.off=5,get.quantile.score = T,
           pop.rpkm.df,pop.meta.df,markers.vec =markers.vec,
           in.col.ramp=colorRampPalette(c("darkgreen","red2"))(1000)){
    intervals=seq(8,8)
    #intervals=seq(3,8)
    out.list=list()
    len.intervals=length(intervals)
    for(m in 1:len.intervals){
      temp.interval=as.numeric(intervals[m])
      plot.new()
      text(x = .5,y=.5,labels =temp.interval,cex = 2.5 )
      temp.out.mat=matrix()
      if(log.rpkm){
        if(get.quantile.score){
          # temp.out.list=get.development.stage.markers.sample.clusters(
          #   rpkm.df=rpkm.df,meta.df=meta.df,markers.df=markers.df,title.str=title.str,
          #   log.rpkm = T,in.sc.cut.off=temp.interval,gene.cluster.cut.off=5,
          #   get.quantile.score = T,pop.rpkm.df=pop.rpkm.df,pop.meta.df=pop.meta.df,
          #   markers.vec =markers.vec,in.col.ramp = in.col.ramp)
          # 
          temp.out.list=get.development.stage.markers.sample.clusters(
            rpkm.df=rpkm.df,meta.df=meta.df,markers.df=markers.df,title.str=title.str,
            log.rpkm = T,in.sc.cut.off=temp.interval,get.quantile.score = T,pop.rpkm.df=pop.rpkm.df,
            pop.meta.df=pop.meta.df,markers.vec =markers.vec,in.col.ramp = in.col.ramp)
          temp.out.mat=temp.out.list
          #temp.out.mat=temp.out.list$markers.score
        }
        else{
          temp.out.list=get.development.stage.markers.sample.clusters(rpkm.df=rpkm.df,meta.df=meta.df,markers.df=markers.df,
                                                                      title.str=title.str,log.rpkm = T,
                                                                      in.sc.cut.off=temp.interval,gene.cluster.cut.off=5,
                                                                      get.quantile.score = F,pop.rpkm.df=pop.rpkm.df,
                                                                      pop.meta.df=pop.meta.df,markers.vec =markers.vec,
                                                                      in.col.ramp = in.col.ramp)
          temp.out.mat=temp.out.list
          #temp.out.mat=temp.out.list$markers.score
        }
      }
      else{
        if(get.quantile.score){
          temp.out.list=get.development.stage.markers.sample.clusters(rpkm.df=rpkm.df,meta.df=meta.df,markers.df=markers.df,
                                                                      title.str=title.str,log.rpkm = F,
                                                                      in.sc.cut.off=temp.interval,gene.cluster.cut.off=5,
                                                                      get.quantile.score = T,pop.rpkm.df=pop.rpkm.df,
                                                                      pop.meta.df=pop.meta.df,markers.vec =markers.vec,
                                                                      in.col.ramp = in.col.ramp)
          temp.out.mat=temp.out.list
          #temp.out.mat=temp.out.list$markers.score
        }
        else{
          temp.out.list=get.development.stage.markers.sample.clusters(rpkm.df=rpkm.df,meta.df=meta.df,markers.df=markers.df,
                                                                      title.str=title.str,log.rpkm = F,
                                                                      in.sc.cut.off=temp.interval,gene.cluster.cut.off=5,
                                                                      get.quantile.score = F,pop.rpkm.df=pop.rpkm.df,
                                                                      pop.meta.df=pop.meta.df,markers.vec =markers.vec,
                                                                      in.col.ramp = in.col.ramp)
          temp.out.mat=temp.out.list
          #temp.out.mat=temp.out.list$markers.score
        }
      }
      out.list[[paste('interval.',temp.interval,sep='')]]=temp.out.mat
    }
    return(out.list)
  }
# get.development.stage.markers.sample.clusters.at.different.cluster.cut.off=
#   function(rpkm.df,meta.df,markers.df,title.str='Test',log.rpkm = F,
#            in.sc.cut.off=6,gene.cluster.cut.off=5,get.quantile.score = T,
#            pop.rpkm.df,pop.meta.df,markers.vec =markers.vec ){
#   
#   intervals=seq(8,8)
#   
#   #intervals=seq(3,8)
#   
#   out.list=list()
#   
#   len.intervals=length(intervals)
#   
#   for(m in 1:len.intervals){
#     
#     temp.interval=as.numeric(intervals[m])
#     
#     plot.new()
#     
#     text(x = .5,y=.5,labels =temp.interval,cex = 2.5 )
#     
#     temp.out.mat=matrix()
#     
#     if(log.rpkm){
#       
#       if(get.quantile.score){
#         
#         temp.out.list=get.development.stage.markers.sample.clusters(rpkm.df=rpkm.df,meta.df=meta.df,
#                                                                     markers.df=markers.df,title.str=title.str,
#                                                                     log.rpkm = T,in.sc.cut.off=temp.interval,
#                                                                     gene.cluster.cut.off=5,get.quantile.score = T,
#                                                                     pop.rpkm.df=pop.rpkm.df,pop.meta.df=pop.meta.df,
#                                                                     markers.vec =markers.vec)
#         
#         temp.out.mat=temp.out.list
#         
#         #temp.out.mat=temp.out.list$markers.score
#       }
#       
#       else{
#         
#         temp.out.list=get.development.stage.markers.sample.clusters(rpkm.df=rpkm.df,meta.df=meta.df,markers.df=markers.df,title.str=title.str,log.rpkm = T,in.sc.cut.off=temp.interval,gene.cluster.cut.off=5,get.quantile.score = F,pop.rpkm.df=pop.rpkm.df,pop.meta.df=pop.meta.df,markers.vec =markers.vec)
#         
#         temp.out.mat=temp.out.list
#         
#         #temp.out.mat=temp.out.list$markers.score
#         
#       }
#       
#       
#     }
#     
#     else{
#       
#       if(get.quantile.score){
#         
#         temp.out.list=get.development.stage.markers.sample.clusters(rpkm.df=rpkm.df,meta.df=meta.df,markers.df=markers.df,title.str=title.str,log.rpkm = F,in.sc.cut.off=temp.interval,gene.cluster.cut.off=5,get.quantile.score = T,pop.rpkm.df=pop.rpkm.df,pop.meta.df=pop.meta.df,markers.vec =markers.vec)
#         
#         temp.out.mat=temp.out.list
#         
#         #temp.out.mat=temp.out.list$markers.score
#         
#       }
#       
#       else{
#         
#         temp.out.list=get.development.stage.markers.sample.clusters(rpkm.df=rpkm.df,meta.df=meta.df,markers.df=markers.df,title.str=title.str,log.rpkm = F,in.sc.cut.off=temp.interval,gene.cluster.cut.off=5,get.quantile.score = F,pop.rpkm.df=pop.rpkm.df,pop.meta.df=pop.meta.df,markers.vec =markers.vec)
#         
#         temp.out.mat=temp.out.list
#         
#         #temp.out.mat=temp.out.list$markers.score
#         
#       }
#     
#     }
#     
#     out.list[[paste('interval.',temp.interval,sep='')]]=temp.out.mat
#     
#   
#   }
#   
#   return(out.list)
# 
# }
get.development.stage.markers.bulk.samples.clusters.at.different.cluster.cut.off=function(pop.rpkm.df,pop.meta.df,markers.df,title.str='Test',log.rpkm = F,gene.cluster.cut.off=5,get.quantile.score = T,markers.vec){
  intervals=seq(6,8)
  #intervals=seq(3,8)
  out.list=list()
  len.intervals=length(intervals)
  for(m in 1:len.intervals){
    temp.interval=as.numeric(intervals[m])
    plot.new()
    text(x = .5,y=.5,labels =temp.interval,cex = 2.5 )
    temp.out.mat=matrix()
    if(log.rpkm){
      if(get.quantile.score){
        #temp.out.list=get.development.stage.markers.sample.clusters(rpkm.df=rpkm.df,meta.df=meta.df,markers.df=markers.df,title.str=title.str,log.rpkm = T,in.sc.cut.off=temp.interval,gene.cluster.cut.off=5,get.quantile.score = T,pop.rpkm.df=pop.rpkm.df,pop.meta.df=pop.meta.df,markers.vec =markers.vec)
        temp.out.mat=temp.out.list
      }
      else{
        #temp.out.list=get.development.stage.markers.sample.clusters(rpkm.df=rpkm.df,meta.df=meta.df,markers.df=markers.df,title.str=title.str,log.rpkm = T,in.sc.cut.off=temp.interval,gene.cluster.cut.off=5,get.quantile.score = F,pop.rpkm.df=pop.rpkm.df,pop.meta.df=pop.meta.df,markers.vec =markers.vec)
        temp.out.list=get.development.stage.markers.bulk.samples.clusters(pop.rpkm.df=pop.rpkm.df,pop.meta.df=pop.meta.df,markers.df=markers.df,title.str=title.str,log.rpkm = T,col.cut.off=temp.interval,gene.cluster.cut.off=5,get.quantile.score = F,markers.vec =markers.vec)
        #temp.out.mat=temp.out.list$markers.score
      }
    }
  }
  #return(out.list)
}
get.development.stage.markers.bulk.samples.clusters=function(pop.rpkm.df,pop.meta.df,markers.df=markers.df,title.str='Test',log.rpkm = T,col.cut.off=temp.interval,gene.cluster.cut.off=5,get.quantile.score = F,markers.vec =markers.vec){
  pop.rpkm.df=filter.none.expressed.genes(input.data =filter.none.expressed.samples(df = pop.rpkm.df) )
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pseudo.val=min(pop.rpkm.df[pop.rpkm.df>0])/2
  expr.markers.vec=intersect(rownames(markers.df),rownames(pop.rpkm.df))
  expr.markers.vec=intersect(expr.markers.vec,markers.vec)
  if(log.rpkm){
    pop.rpkm.df=log2(pop.rpkm.df+pseudo.val)
  }
  pop.rpkm.mat=as.matrix(pop.rpkm.df[expr.markers.vec,])
  row.annotation.df=subset(x = markers.df[expr.markers.vec,],select = category)
  col.annotation.df=subset(pop.meta.df,select=c(development.stage,sub.group))
  pheatmap(mat = pop.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2,fontsize_col = 2,cutree_cols = col.cut.off,cluster_rows = F,cellheight = 2,cellwidth = 5,annotation_legend = F,border_color = NA,annotation_row = row.annotation.df,annotation_col = col.annotation.df)
}
#Detect at different parameters
detect.sc.pop.cluster.relation.based.development.stage.markers.at.varying.cut.off=function(pop.marker.category.score.df,pop.meta.df,pop.gene.cluster.cut.off=5,sc.marker.category.score.df,sc.gene.cluster.cut.off = 3,sc.meta.df,markers.df,title.str='Test'){
  pop.cluster.cut.off=seq(6,6)
  #pop.cluster.cut.off=seq(6,6)
  sc.cluster.cut.off=seq(7,7)
  #sc.cluster.cut.off=seq(6,6)
  out.list=list()
  len.pop.cluster.cut.off=length(pop.cluster.cut.off)
  len.sc.cluster.cut.off=length(sc.cluster.cut.off)
  for(n in 1:len.pop.cluster.cut.off){
    for(m in 1:len.sc.cluster.cut.off){
      temp.sc.cluster.cut.off=sc.cluster.cut.off[m]
      temp.pop.cluster.cut.off=pop.cluster.cut.off[n]
      plot.new()
      text(x = .5, y=.3,labels = paste(c('Pop cluster cut off: ',temp.pop.cluster.cut.off),collapse = ''))
      text(x = .5, y=.8,labels = paste(c('SC cluster cut off: ',temp.sc.cluster.cut.off),collapse = ''))
      temp.out=detect.sc.pop.cluster.relation.based.development.stage.markers(pop.marker.category.score.df=pop.marker.category.score.df,pop.meta.df=pop.meta.df,pop.clusters.cut.off=temp.pop.cluster.cut.off,pop.gene.cluster.cut.off=5,sc.marker.category.score.df=sc.marker.category.score.df,sc.clusters.cut.off=temp.sc.cluster.cut.off,sc.gene.cluster.cut.off = 3,sc.meta.df=sc.meta.df,markers.df=markers.df,title.str=title.str)
      #detect.sc.pop.cluster.relation.based.development.stage.markers(pop.marker.category.score.df = all.pop.blood.grp.A.atleast.500.genes.100k.unique.including.batch.eight.reads.artefacts.filt.tn5.tpm.markers.score.mat,pop.clusters.cut.off = 6,pop.gene.cluster.cut.off = 9,sc.marker.category.score.df =all.sc.blood.group.A.atleast.500.genes.10k.unique.reads.artefacts.filt.tpm.markers.score.mat,sc.clusters.cut.off = 6,sc.gene.cluster.cut.off = 8,sc.meta.df = all.samples.artefacts.filt.meta.df,markers.df = all.plos.one.markers.df )
      lab.str=paste('sc.cut.off:',temp.sc.cluster.cut.off,',pop.cut.off:',temp.pop.cluster.cut.off,sep = '')
      out.list[[lab.str]]=temp.out
    }
  }
  return(out.list)
}
detect.sc.pop.cluster.relation.based.development.stage.markers=function(pop.marker.category.score.df,pop.meta.df,pop.clusters.cut.off=3,pop.gene.cluster.cut.off=5,sc.marker.category.score.df,sc.clusters.cut.off=3,sc.gene.cluster.cut.off = 3,sc.meta.df,markers.df,title.str='Test',sc.title.str='SC clusters'){
  out.list=list()
  pop.marker.category.score.df=data.frame(pop.marker.category.score.df)
  #pop.marker.category.score.df=filter.non.variable.rows(df = marker.category.score.df,cut.off = .5)
  pop.gene.cols.list=get.col.factor(col.factor = as.character(markers.df[rownames(pop.marker.category.score.df),'category']))
  title.str=convert.to.title.case(title.str)
  pop.out.mat=as.matrix(pop.marker.category.score.df)
  #pop.cor.mat=cor(pop.out.mat,method = 'spearman')
  pop.samples.dist=dist(t(pop.out.mat))
  pop.samples.dist.mat=as.matrix(pop.samples.dist)
  gene.clusters<- hclustfunc(distfunc(pop.out.mat))
  sample.clusters <- hclustfunc(distfunc(t(pop.out.mat)))
  gene.groups <- cutree(gene.clusters, as.numeric(pop.gene.cluster.cut.off))
  sample.groups <- cutree(sample.clusters, as.numeric(pop.clusters.cut.off))
  pop.groups.df=data.frame(genes=names(sample.groups),group=paste('Pop grp: ',as.character(as.numeric(sample.groups)),sep=''),row.names = names(sample.groups))
  out.list[['pop.grps']]=pop.groups.df
  samples.grp.lab=unique(as.character(pop.groups.df$group))
  pop.gene.cols <- brewer.pal(as.numeric(pop.gene.cluster.cut.off), "Pastel1")
  samples.cols <-rainbow(n = as.numeric(pop.clusters.cut.off))
  #heatmap.2(x = pop.out.mat,margins = c(13,13),main='Pop clusters',trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),key.title = 'Euclidean dist.',scale='row',RowSideColors=pop.gene.cols.list$col.str,ColSideColors=samples.cols[sample.groups],dendrogram = 'col')
  #plot.new()
  #legend('left',legend=pop.gene.cols.list$legend.str,fill=pop.gene.cols.list$legend.col,box.lty = 0,cex = .3)
  pop.meta.df=pop.meta.df[colnames(pop.out.mat),]
  pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = pop.meta.df)
  pop.timepoint.col=get.abbrev.std.stage.cols(in.stages.vec =as.character(pop.meta.df$development.stage))
  heatmap.2(x = pop.out.mat,margins = c(13,13),main='Pop. clusters',trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),ColSideColors = pop.timepoint.col,key.title = 'Euclidean dist.',scale='row',dendrogram = 'col',RowSideColors = pop.gene.cols.list$col.str)
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .4)
  #legend('left',legend=pop.gene.cols.list$legend.str,fill=pop.gene.cols.list$legend.col,box.lty = 0,cex = .3)
  pop.classified.grps.cols.list=brewer.pal(n = as.numeric(pop.clusters.cut.off),'Set3')
  heatmap.2(x = pop.out.mat,margins = c(13,13),main='Pop. re-grouped',trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),ColSideColors = pop.classified.grps.cols.list[sample.groups],key.title = 'Euclidean dist.',scale='row',dendrogram = 'col',RowSideColors = pop.gene.cols.list$col.str)
  legend('topright',legend=samples.grp.lab,fill=pop.classified.grps.cols.list[unique(sample.groups)],box.lty = 0,cex = .4)
  sc.out.mat=as.matrix(sc.marker.category.score.df)
  gene.dist=dist(sc.out.mat)
  gene.dist.mat=as.matrix(gene.dist)
  gene.dist.hclust=hclustfunc(x = gene.dist)
  sc.genes.intern <- clValid(obj =gene.dist.mat, nClust = 2:3, method = 'complete',clMethods=c("hierarchical", "kmeans", "diana", "som", "model", "sota", "pam", "clara", "agnes"),validation="internal",metric = "euclidean")
  #sc.genes.intern <- hclust(gene.dist)
  sc.genes.internal.clusters=clusters(sc.genes.intern)
  #gene.groups <- cutree(sc.genes.intern, k = 3)
  gene.groups <- cutree(sc.genes.internal.clusters, k = 3)
  gene.groups.df=data.frame(genes=names(gene.groups),group=paste('Gene grp: ',as.character(as.numeric(gene.groups)),sep=''),row.names = names(gene.groups ))
  gene.grp.col.side.col.list=get.col.factor(col.factor = as.character(gene.groups.df$group))
  #heatmap.2(gene.dist.mat,margins = c(10,10),RowSideColors = gene.grp.col.side.col.list$col.str,main='Genes euclidean clusters',trace='none',cexRow = .3,cexCol = .3,col=bluered(10000),key.title = 'Genes euclidian dist.',dendrogram = 'col')
  gene.groups.df=data.frame(genes=names(gene.groups),group=paste('Gene grp: ',as.character(as.numeric(gene.groups)),sep=''),row.names = names(gene.groups))
  gene.groups.list=split(gene.groups.df,f=gene.groups.df$group)
  gene.groups.names=names(gene.groups.list)
  len.gene.groups.names=length(gene.groups.names)
  genes.to.use=c()
  for(n in 1:len.gene.groups.names){
    gene.groups.name=gene.groups.names[n]
    temp.title.str=convert.to.title.case(gene.groups.name)
    temp.gene.groups.df=gene.groups.list[[gene.groups.name]]
    temp.genes.vec=rownames(temp.gene.groups.df)
    temp.out.mat=sc.out.mat[temp.genes.vec,]
    if(!gene.groups.name=='Gene grp: 1'){
      genes.to.use=append(genes.to.use,temp.genes.vec,length(genes.to.use))
    }
    #genes.to.use=append(genes.to.use,temp.genes.vec,length(genes.to.use))
    #else{
    temp.gene.cl=hclustfunc(distfunc(x = temp.out.mat))
    temp.gene.cl.grps=cutree(tree = temp.gene.cl,k = 3)
    temp.gene.cl.cols=brewer.pal(n = 3,name = 'Pastel1')
    grp.one.gene.groups.df=data.frame(genes=names(temp.gene.cl.grps),group=paste('Gene grp: ',as.character(as.numeric(temp.gene.cl.grps)),sep=''),row.names = names(temp.gene.cl.grps))
    if(gene.groups.name=='Gene grp: 1'){
      temp.grp.one.genes.vec=rownames(subset(grp.one.gene.groups.df,group!='Gene grp: 1'))
      genes.to.use=append(genes.to.use,temp.grp.one.genes.vec,length(genes.to.use))
    }
    grp.one.gene.groups.list=split(grp.one.gene.groups.df,f=grp.one.gene.groups.df$group)
    gene.groups.names=names(gene.groups.list)
    #heatmap.2(x = temp.out.mat,margins = c(13,13),main=temp.title.str,trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),key.title = 'Euclidean dist.',scale='row',RowSideColors = temp.gene.cl.cols[temp.gene.cl.grps],dendrogram = 'col')
    heatmap.2(x = temp.out.mat,margins = c(13,13),main=temp.title.str,trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),key.title = 'Euclidean dist.',RowSideColors = temp.gene.cl.cols[temp.gene.cl.grps],dendrogram = 'col')
    legend('left',legend=sort(unique(gene.groups.names)),fill=temp.gene.cl.cols, box.lty = 0,cex=.3)
    #}
    #temp.gene.cols.list=get.col.factor(col.factor = as.character(markers.df[temp.genes.vec,'category']))
    #heatmap.2(x = temp.out.mat,margins = c(13,13),main=temp.title.str,trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),key.title = 'Euclidean dist.',scale='row',RowSideColors = temp.gene.cols.list$col.str,dendrogram = 'col')
    #plot.new()
    #legend('left',legend=temp.gene.cols.list$legend.str,fill=temp.gene.cols.list$legend.col,box.lty = 0,cex = .3)
  }
  #genes.to.use=genes.to.use[1:10]
  genes.to.use=unique(genes.to.use)
  sc.final.out.mat=sc.out.mat[genes.to.use,]
  sc.final.markers.df=markers.df[genes.to.use,]
  sc.final.out.dist=dist(t(sc.final.out.mat))
  sc.final.out.dist.mat=as.matrix(sc.final.out.dist)
  sc.final.samples.hclust=hclustfunc(x = sc.final.out.dist)
  sc.final.genes.hclust=hclustfunc(x = distfunc(x =sc.final.out.mat))
  sc.final.groups <- cutree(sc.final.samples.hclust, k = as.numeric(sc.clusters.cut.off))
  sc.final.groups.df=data.frame(genes=names(sc.final.groups),sp=paste('SC grp: ',as.character(as.numeric(sc.final.groups)),sep=''),row.names = names(sc.final.groups))
  sc.final.grp.legend.str=sort(unique(as.character(sc.final.groups.df$sp)))
  sc.final.col.list=brewer.pal(n = as.numeric(sc.clusters.cut.off),'Set3')
  sc.final.genes.groups=cutree(tree = sc.final.genes.hclust,as.numeric(sc.gene.cluster.cut.off))
  sc.final.markers.gene.col.list=get.col.factor(col.factor = as.character(sc.final.markers.df$category))
  sc.gene.cols=brewer.pal(n =as.numeric(sc.gene.cluster.cut.off),'Pastel1' )
  #heatmap.2(x = sc.final.out.mat,margins = c(13,13),main='SC clusters',trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),ColSideColors = sc.final.col.list[sc.final.groups],key.title = 'Euclidean dist.',scale='row',dendrogram = 'col',RowSideColors = sc.final.markers.gene.col.list$col.str)
  gene.lab.str=as.character(sc.final.markers.df[rownames(sc.final.out.mat),]$category)
  heatmap.2(x = sc.final.out.mat,margins = c(13,13),main=sc.title.str,trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),ColSideColors = sc.final.col.list[sc.final.groups],key.title = 'Euclidean dist.',dendrogram = 'col',RowSideColors = sc.final.markers.gene.col.list$col.str)
  legend('bottomleft',legend=sc.final.markers.gene.col.list$legend.str,fill=sc.final.markers.gene.col.list$legend.col,box.lty = 0,cex = .5)
  legend('topright',legend=sc.final.grp.legend.str,fill=sc.final.col.list[unique(sc.final.groups)],box.lty = 0,cex = .5)
  temp.sc.annotation.df=subset(sc.final.groups.df,select = sp)
  row.sc.final.markers.df=subset(sc.final.markers.df,select=category)
  ann_colors=list(
    sp=factor(sc.final.col.list[sc.final.groups]),
    category=factor(sc.final.markers.gene.col.list$col.str)
    )
  pheatmap(mat  = sc.final.out.mat,main=sc.title.str,color=bluered(10000),cutree_cols = as.numeric(sc.clusters.cut.off),fontsize_row = 2,fontsize_col = 2,cellwidth = 1,cellheight = 1,clustering_distance_rows = 'correlation',clustering_distance_cols = 'euclidean',annotation_col = temp.sc.annotation.df,annotation_row = row.sc.final.markers.df,annotation_legend = T, treeheight_row = 1)
  col.sc.final.groups.df=subset(sc.final.groups.df,select=sp)
  test.col.list=get.col.factor(as.character(sc.final.groups.df$sp))
  #pheatmap(mat  = sc.final.out.mat,main='',color=bluered(10000),scale = 'row',cutree_cols = as.numeric(sc.clusters.cut.off),fontsize_row = 2,fontsize_col = 2,cellwidth = 3,cellheight = 3,annotation_col =  col.sc.final.groups.df)
  pheatmap(mat  = sc.final.out.mat,main='',color=bluered(10000),cutree_cols = as.numeric(sc.clusters.cut.off),fontsize_row = 2,fontsize_col = 2,cellwidth = 3,cellheight = 3,clustering_distance_rows = 'correlation',clustering_distance_cols = 'euclidean',annotation_colors  = test.col.list)
  plot.new()
  legend('bottomleft',legend=sc.final.markers.gene.col.list$legend.str,fill=sc.final.markers.gene.col.list$legend.col,box.lty = 0,cex = 1.5)
  legend('topright',legend=sc.final.grp.legend.str,fill=sc.final.col.list[unique(sc.final.groups)],box.lty = 0,cex = 1.5)
  final.sc.meta.df=sc.meta.df[colnames(sc.final.out.mat),]
  final.sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = final.sc.meta.df)
  timepoint.col=get.abbrev.std.stage.cols(in.stages.vec =as.character(final.sc.meta.df$development.stage))
  heatmap.2(x = sc.final.out.mat,margins = c(13,13),main='SC clusters',trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),ColSideColors = timepoint.col,key.title = 'Euclidean dist.',scale='row',dendrogram = 'col',RowSideColors = sc.final.markers.gene.col.list$col.str)
  legend('topright',legend=std.stages.col.abbrev.legend.str,fill=std.stages.col.legend.fill,box.lty = 0,cex = .4)
  batch.col.list=get.col.factor(as.character(final.sc.meta.df$batch))
  ###Validate the clusters
  #Randmize groups 10 times
  sc.grps.vec=as.character(sc.final.groups.df$group)
  all.re.ordered.cluster.euclidean.distance.vec=c()
  randomize.param=100
#   
#   for(s in 1:randomize.param){
#     
#     sc.final.groups.df$re.ordered.grp=sample(sc.grps.vec,size = length(sc.grps.vec))
#     
#     sc.final.re.ordered.groups.list=split(sc.final.groups.df,f=sc.final.groups.df$re.ordered.grp)
#     
#     sc.final.re.ordered.groups.names=names(sc.final.re.ordered.groups.list)
#     
#     len.sc.final.re.ordered.groups.names=length(sc.final.re.ordered.groups.names)
#     
#     re.ordered.cluster.euclidean.distance.vec=c()
#     
#     for(l in 1:len.sc.final.re.ordered.groups.names){
#       
#       temp.sc.final.re.ordered.groups.name=sc.final.re.ordered.groups.names[l]
#       
#       temp.sc.final.re.ordered.samples.vec=rownames(sc.final.re.ordered.groups.list[[temp.sc.final.re.ordered.groups.name]])
#       
#       temp.sc.final.re.ordered.samples.mat=combn(temp.sc.final.re.ordered.samples.vec,m = 2)
#       
#       for(f in 1:dim(temp.sc.final.re.ordered.samples.mat)[2]){
#         
#         sc.re.ordered.final.pair=as.character(temp.sc.final.re.ordered.samples.mat[,f])
#         
#         re.ordered.cluster.euclidean.distance.vec=append(re.ordered.cluster.euclidean.distance.vec,as.numeric(sc.final.out.dist.mat[sc.re.ordered.final.pair[1],sc.re.ordered.final.pair[2]]),after = length(re.ordered.cluster.euclidean.distance.vec))
#         
#       }
#       
#     }
#     
#     all.re.ordered.cluster.euclidean.distance.vec=append(all.re.ordered.cluster.euclidean.distance.vec,sum(re.ordered.cluster.euclidean.distance.vec),length(all.re.ordered.cluster.euclidean.distance.vec))
#     
#   }
#   
#   sc.final.groups.list=split(sc.final.groups.df,f=sc.final.groups.df$group)
#   
#   sc.final.groups.names=names(sc.final.groups.list)
# 
#   len.sc.final.groups.list=length(sc.final.groups.names)
#   
#   cluster.euclidean.distance.vec=c()
#   
#   for(m in 1:len.sc.final.groups.list){
#     
#     temp.sc.final.group.name=sc.final.groups.names[m]
#     
#     temp.sc.final.samples.vec=rownames(sc.final.groups.list[[temp.sc.final.group.name]])
#     
#     temp.sc.final.samples.mat=combn(temp.sc.final.samples.vec,m = 2)
#     
#     for(d in 1:dim(temp.sc.final.samples.mat)[2]){
#       
#       sc.final.pair=as.character(temp.sc.final.samples.mat[,d])
#       
#       cluster.euclidean.distance.vec=append(cluster.euclidean.distance.vec,as.numeric(sc.final.out.dist.mat[sc.final.pair[1],sc.final.pair[2]]),after = length(cluster.euclidean.distance.vec))
#       
#     }
#     
#   }
#   
#   #test.value=ifelse(all.re.ordered.cluster.euclidean.distance.vec>cluster.euclidean.distance.vec,0,1)
#   
#   test.value=c()
#   
#   len.all.re.ordered.cluster.euclidean.distance.vec=length(all.re.ordered.cluster.euclidean.distance.vec)
#   
#   for(a in 1:len.all.re.ordered.cluster.euclidean.distance.vec){
#     
#     temp.dist=as.numeric(all.re.ordered.cluster.euclidean.distance.vec[a])
#     
#     test.value=append(test.value,ifelse(temp.dist>cluster.euclidean.distance.vec,0,1),length(test.value))
#     
#   }
#  
#   p.value=sum(test.value)/randomize.param
#  
#   lab.p.value=format(p.value,digits = 3)
# 
#   title.str=paste('Pvalue: ',lab.p.value,sep='')
#   
#   #barplot2(height = c(log10(sum(cluster.euclidean.distance.vec)),log10(mean(all.re.ordered.cluster.euclidean.distance.vec))),names.arg = c('Clustered','Random'),las=2)
#   
#   euclidean.dist.list=list(Clusters=log10(cluster.euclidean.distance.vec),Random=log10(all.re.ordered.cluster.euclidean.distance.vec))
#   
#   euclidean.dist.boxplot=boxplot2(x = euclidean.dist.list,las=2)
#   
#   text(x=1.5,y=max(euclidean.dist.boxplot$stats)/2,labels=c(title.str),cex=1.0)
#   
#   euclidean.dist.boxplot=boxplot(x = euclidean.dist.list,las=2,names=rep('',times = 2))
#   
  #text(euclidean.dist.boxplot)
  #heatmap.2(x = sc.final.out.mat,margins = c(13,13),main='SC clusters',trace='none',cexRow = .1,cexCol = .3,col=bluered(10000),key.title = 'Euclidean dist.',scale='row',dendrogram = 'col',RowSideColors =sc.final.markers.gene.col.list$col.str ,ColSideColors = batch.col.list$col.str)
  #legend('topright',legend=batch.col.list$legend.str,fill=batch.col.list$legend.col,box.lty = 0,cex = .4)
  #plot.the.sub.cluster.distances(first.df=sc.final.out.mat,first.marker.sub.grps.df=sc.final.groups.df,sec.df=pop.out.mat,sec.marker.sub.grps.df=pop.groups.df)
  out.list[['sc.grps']]=sc.final.groups.df
  return(out.list)
}
detect.sc.pop.cluster.relation.based.development.stage.markers.at.diff.cluster.cut.off=function(pop.marker.category.score.df,pop.meta.df,pop.gene.cluster.cut.off=5,sc.marker.category.score.df,sc.meta.df,markers.df,title.str='Test',sc.title.str='SC clusters'){
  out.list=list()
  pop.cut.off=seq(6,6)
  sc.cut.off=seq(7,7)
  len.pop.cut.off=length(pop.cut.off)
  len.sc.cut.off=length(sc.cut.off)
  for(m in 1:len.pop.cut.off){
    temp.pop.cut.off=as.numeric(pop.cut.off[m])
    for(n in 1:len.sc.cut.off){
      temp.sc.cut.off=as.numeric(sc.cut.off[n])
      sc.title.str=convert.to.title.case(in.str = sc.title.str)
      cut.off.str=paste(c('(n=',temp.sc.cut.off,')'),collapse='')
      temp.sc.title.str=paste(sc.title.str,cut.off.str,sep='')
      temp.out.df=detect.sc.pop.cluster.relation.based.development.stage.markers(pop.marker.category.score.df=pop.marker.category.score.df,pop.meta.df=pop.meta.df,pop.clusters.cut.off =temp.pop.cut.off,pop.gene.cluster.cut.off=5,sc.marker.category.score.df=sc.marker.category.score.df,sc.clusters.cut.off=temp.sc.cut.off,sc.meta.df=sc.meta.df,markers.df=markers.df,title.str='Test',sc.title.str= temp.sc.title.str)
      list.label=paste('cluster.cut.off:',temp.sc.cut.off,sep='')
      out.list[[list.label]]=temp.out.df
    }
  }
  return(out.list)
}
plot.the.sub.cluster.distances=function(first.df,first.marker.sub.grps.df,sec.df,sec.marker.sub.grps.df){
  first.marker.sub.grps.list=split(first.marker.sub.grps.df,f=first.marker.sub.grps.df$group)
  sec.marker.sub.grps.list=split(sec.marker.sub.grps.df,f=sec.marker.sub.grps.df$group)
  first.marker.sub.grps.names=names(first.marker.sub.grps.list)
  len.first.marker.sub.grps.names=length(first.marker.sub.grps.names)
  sec.marker.sub.grps.names=names(sec.marker.sub.grps.list)
  len.sec.marker.sub.grps.names=length(sec.marker.sub.grps.names)
  dist.out.mat=matrix(nrow = len.first.marker.sub.grps.names,ncol =len.sec.marker.sub.grps.names )
  out.list=list()
  for(n in 1:len.first.marker.sub.grps.names){
    first.marker.sub.grps.name=first.marker.sub.grps.names[n]
    temp.first.marker.sub.grps.df=first.marker.sub.grps.list[[first.marker.sub.grps.name]]
    temp.first.df=subset(first.df,select = rownames(temp.first.marker.sub.grps.df))
    temp.dist.list=list()
    for(m in 1:len.sec.marker.sub.grps.names){
      sec.marker.sub.grps.name=sec.marker.sub.grps.names[m]
      temp.sec.marker.sub.grps.df=sec.marker.sub.grps.list[[sec.marker.sub.grps.name]]
      temp.sec.df=subset(sec.df,select = rownames(temp.sec.marker.sub.grps.df))
      intersect.genes.vec=intersect(rownames(temp.first.df),rownames(temp.sec.df))
      temp.filt.first.df=temp.first.df[intersect.genes.vec,]
      temp.filt.sec.df=temp.sec.df[intersect.genes.vec,]
      combn.df=data.frame(temp.filt.first.df,temp.filt.sec.df)
      combn.cor.dist.mat=cor(combn.df,method = 'spearman')
      temp.cor.dist.mat=combn.cor.dist.mat[colnames(temp.filt.first.df),colnames(temp.filt.sec.df)]
      if(dim(temp.cor.dist.mat)[1]<1 ||dim(temp.cor.dist.mat)[2]<1 ){
        dist.out.mat[n,m]=0
        next
      }
      else{
        cor.dist.vec=as.numeric(melt(temp.cor.dist.mat)$value)
        temp.dist.list[[sec.marker.sub.grps.name]]=cor.dist.vec
        dist.out.mat[n,m]=median(cor.dist.vec)
      }
    }
    if(length(temp.dist.list)<1){
      next
    }
    out.list[[first.marker.sub.grps.name]]=temp.dist.list
  }
  return(out.list)
}
calculate.proportion.samples.clustered.together=function(first.meta.df,sec.meta.df){
  intersect.samples=intersect(rownames(first.meta.df),rownames(sec.meta.df))
  sample.combn=combn(x = intersect.samples,m = 2)
  no.pairs=dim(sample.combn)[2]
  pass.test.counter=c()
  for(n in 1:no.pairs){
    temp.pair=sample.combn[,n]
    first.sample=temp.pair[1]
    sec.sample=temp.pair[2]
    batch.uncorrected.test=as.character(first.meta.df[first.sample,'markers.cluster.groups'])==as.character(first.meta.df[sec.sample,'markers.cluster.groups'])
    batch.corrected.test=as.character(first.meta.df[first.sample,'markers.cluster.groups'])==as.character(sec.meta.df[sec.sample,'markers.cluster.groups'])
    pass.test.counter=append(pass.test.counter,batch.uncorrected.test==batch.corrected.test,after = length(pass.test.counter))
  }
}
hclustfunc = function(x,method="complete"){
  hclust.out=hclust(x, method=method)
  return(hclust.out)
} 
distfunc <- function(x){
  dist.out=dist(x, method="euclidean")
  return(dist.out)
} 
get.std.stage.cols=function(in.stages.vec){
  col.vec=c()
  len.in.stages.vec=length(in.stages.vec)
  for(n in 1:len.in.stages.vec){
    temp.stage=as.character(in.stages.vec[n])
    if(is.null(std.stages.col.list[[temp.stage]])){
      show(temp.stage)
    }
    col.vec=append(col.vec,std.stages.col.list[[temp.stage]],after = length(col.vec))
  }
  return(col.vec)
}
get.abbrev.std.stage.cols=function(in.stages.vec){
  col.vec=c()
  len.in.stages.vec=length(in.stages.vec)
  for(n in 1:len.in.stages.vec){
    temp.stage=as.character(in.stages.vec[n])
    if(is.null(abbrev.std.stages.col.list[[temp.stage]])){
      show(temp.stage)
    }
    col.vec=append(col.vec,abbrev.std.stages.col.list[[temp.stage]],after = length(col.vec))
  }
  return(col.vec)
}
add.stage.col.to.meta.df=function(meta.df,abbrev=F){
  stages.vec=as.character(meta.df$development.stage)
  if(abbrev){
    meta.df$stage.col=get.abbrev.std.stage.cols(in.stages.vec = stages.vec)
  }
  else{
    meta.df$stage.col=get.std.stage.cols(in.stages.vec = stages.vec)
  }
  return(meta.df)
}
add.abbrev.stage.col.to.meta.df=function(meta.df){
  stages.vec=as.character(meta.df$development.stage)
  meta.df$stage.col=get.std.stage.cols(in.stages.vec = stages.vec)
  return(meta.df)
}
change.development.stage.to.abbrev.stage.meta.df=function(meta.df){
  stages.vec=as.character(meta.df$development.stage)
  meta.df$development.stage=gsub(gsub(gsub(gsub(gsub(gsub(x = stages.vec,pattern = 'late.trophozoite',
                                                          replacement = 'T'),pattern = 'late.ring',
                                                     replacement = 'LR'),pattern = 'ring',replacement = 'R'),
                                           pattern = 'early.schizont',replacement = 'ES'),
                                      pattern = 'schizont',replacement = 'S'),pattern = 'early.trophozoite',
                                 replacement = 'ET')
  return(meta.df)
}
get.distance.between.bulk.and.sc=function(pop.mat,pop.rpkm.df,pop.meta.df,sc.mat,sc.rpkm.df,sc.meta.df,title.str='Test',markers.vec,pop.clusters.cut.off=2,sc.clusters.cut.off=2){
  markers.vec=intersect(intersect(rownames(pop.rpkm.df),rownames(sc.rpkm.df)),markers.vec)
  pop.dist=dist(t(pop.mat))
  sc.dist=dist(t(sc.mat))
  pop.hclust=hclust(d =pop.dist )
  #ColSideColors = pop.col.side.col.list$col.str
  #heatmap.2(x = as.matrix(pop.dist),margins = c(10,10),main='Bulk euclidean dist.',trace='none',cexRow = .3,cexCol = .3,col=bluered(10000),key.title = 'Euclidian distance')
  intern <- clValid(obj =as.matrix(pop.dist), nClust = 2:6, method = 'complete',clMethods=c("hierarchical", "kmeans", "diana", "som", "model", "sota", "pam", "clara", "agnes"),validation="internal",metric = "euclidean")
  #plot(intern)
  internal.clusters=clusters(intern)
  #plot(as.dendrogram(internal.clusters, use.modularity = F),cex=.2,main=title.str)
  #rect.hclust(internal.clusters, k=6, border="lightblue")
  pop.groups <- cutree(internal.clusters, k = pop.clusters.cut.off)
  pop.groups.df=data.frame(samples=names(pop.groups),group=as.numeric(pop.groups),row.names = names(pop.groups))
  #heatmap.2(x = as.matrix(sc.dist),margins = c(10,10),main='SC euclidean dist.',trace='none',cexRow = .3,cexCol = .3,col=bluered(10000),key.title = 'Euclidian distance')
  sc.intern <- clValid(obj =as.matrix(sc.dist), nClust = 2:7, method = 'complete',clMethods=c("hierarchical", "kmeans", "diana", "som", "model", "sota", "pam", "clara", "agnes"),validation="internal",metric = "euclidean")
  #plot(intern)
  sc.internal.clusters=clusters(sc.intern)
  #plot(as.dendrogram(sc.internal.clusters, use.modularity = T),cex=.2,main=title.str)
  #rect.hclust(sc.internal.clusters, k=12, border="lightblue")
  sc.groups <- cutree(sc.internal.clusters, k = sc.clusters.cut.off)
  sc.groups.df=data.frame(samples=names(sc.groups),group=as.numeric(sc.groups),row.names = names(sc.groups))
  #pop.meta.df,sc.mat,sc.meta.df
  filt.pop.meta.df=data.frame(t(subset(t(pop.meta.df),select=rownames(pop.groups.df))))
  filt.pop.meta.df$clusters=as.character(pop.groups.df$group)
  filt.pop.meta.list=split(filt.pop.meta.df,f=filt.pop.meta.df$clusters)
  filt.pop.meta.names=names(filt.pop.meta.list)
  len.filt.pop.meta.names=length(filt.pop.meta.names)
  filt.sc.meta.df=data.frame(t(subset(t(sc.meta.df),select=rownames(sc.groups.df))))
  filt.sc.meta.df$clusters=as.character(sc.groups.df$group)
  filt.sc.meta.list=split(filt.sc.meta.df,f=filt.sc.meta.df$clusters)
  filt.sc.meta.names=names(filt.sc.meta.list)
  len.filt.sc.meta.names=length(filt.sc.meta.names)
  cor.dist.mat=matrix(nrow =len.filt.sc.meta.names,ncol = len.filt.pop.meta.names )
  euclidian.dist.mat=matrix(nrow =len.filt.sc.meta.names,ncol = len.filt.pop.meta.names )
  sc.sub.grp.names=c()
  pop.sub.grp.names=c()
  for(n in 1:len.filt.sc.meta.names){
    filt.sc.meta.name=filt.sc.meta.names[n]
    sc.temp.name=paste('SC.cluster.',filt.sc.meta.name,collapse = '')
    sc.sub.grp.names=append(sc.sub.grp.names,sc.temp.name,after = length(sc.sub.grp.names))
    temp.filt.sc.meta.df=filt.sc.meta.list[[filt.sc.meta.name]]
    #temp.sc.df=subset(sc.mat,select=rownames(temp.filt.sc.meta.df))
    temp.sc.df=filter.none.expressed.genes(subset(sc.rpkm.df,select=rownames(temp.filt.sc.meta.df)))
    #temp.sc.df= temp.sc.df[markers.vec,]
    temp.average.sc.df=data.frame(sc=as.numeric(apply(temp.sc.df,1,mean)),row.names = rownames(temp.sc.df))
    cor.dist.list=list()
    euclidian.dist.list=list()
    mean.cor.dist.vec=c()
    mean.euclidian.dist.vec=c()
    for(m in 1:len.filt.pop.meta.names){
      filt.pop.meta.name=filt.pop.meta.names[m]
      if(n==1){
        pop.temp.name=paste('Pop.cluster.',filt.pop.meta.name,collapse = '')
        pop.sub.grp.names=append(pop.sub.grp.names,pop.temp.name,after = length(pop.sub.grp.names))
      }
      temp.filt.pop.meta.df=filt.pop.meta.list[[filt.pop.meta.name]]
      #temp.pop.df=subset(pop.mat,select=rownames(temp.filt.pop.meta.df))
      temp.pop.df=filter.none.expressed.genes(subset(pop.rpkm.df,select=rownames(temp.filt.pop.meta.df)))
      #temp.pop.df=temp.pop.df[markers.vec,]
      temp.cor.dist=get.cor.between.two.df.cor(first.df = temp.average.sc.df,second.df = temp.pop.df)
      cor.dist.list[[filt.pop.meta.name]]=temp.cor.dist
      temp.euclidian.dist=get.euclidian.dist.between.two.df.cor(first.df = temp.average.sc.df,second.df = temp.pop.df)
      euclidian.dist.list[[filt.pop.meta.name]]=temp.euclidian.dist
      mean.cor.dist.vec=append(mean.cor.dist.vec,median(temp.cor.dist),after = length(mean.cor.dist.vec))
      mean.euclidian.dist.vec=append(mean.euclidian.dist.vec,median(temp.euclidian.dist),after = length(mean.euclidian.dist.vec))
    }
    cor.dist.mat[n,]=as.numeric(mean.cor.dist.vec)
    euclidian.dist.mat[n,]=as.numeric(mean.euclidian.dist.vec)
  }
  rownames(cor.dist.mat)=sc.sub.grp.names
  colnames(cor.dist.mat)=pop.sub.grp.names
  parameters.str=paste('SC(n=',sc.clusters.cut.off,'),','Pop(n=',pop.clusters.cut.off,')',collapse = '')
  title.str=parameters.str
  heatmap.2(x = as.matrix(cor.dist.mat),margins = c(10,10),main=title.str,trace='none',cexRow = .3,cexCol = .3,col=bluered(1000),key.title = 'Cor. distance')
  rownames(euclidian.dist.mat)=sc.sub.grp.names
  colnames(euclidian.dist.mat)=pop.sub.grp.names
  heatmap.2(x = as.matrix(euclidian.dist.mat),margins = c(10,10),main=title.str,trace='none',cexRow = .3,cexCol = .3,col=bluered(1000),key.title = 'Euclidean dist.')
}
get.distance.between.bulk.and.sc.at.different.cluster.cut.off=function(pop.dist.mat,pop.rpkm.df,pop.meta.df,sc.dist.mat,sc.rpkm.df,sc.meta.df,title.str='Test',markers.vec){
  sc.intervals=seq(3,19)
  #sc.intervals=seq(3,4)
  pop.intervals=seq(3,9)
  #pop.intervals=seq(3,4)
  len.sc.intervals=length(sc.intervals)
  len.pop.intervals=length(pop.intervals)
  pop.mat=as.matrix(pop.dist.mat)
  sc.mat=as.matrix(sc.dist.mat)
  for(n in 1:len.pop.intervals){
    for(m in 1:len.sc.intervals){
      sc.cut.off=sc.intervals[m]
      pop.cut.off=pop.intervals[n]
      get.distance.between.bulk.and.sc(pop.mat=pop.mat,pop.rpkm.df=pop.rpkm.df,pop.meta.df=pop.meta.df,sc.mat=sc.mat,sc.rpkm.df=sc.rpkm.df,sc.meta.df=sc.meta.df,title.str=title.str,markers.vec=markers.vec,pop.clusters.cut.off=pop.cut.off,sc.clusters.cut.off=sc.cut.off)
    }
  }
}
get.mean.gene.signature.score.per.cell=function(rpkm.df,genes.sig.vec,quant.norm=F,log.rpkm=F){
  if(log.rpkm){
    pseudo.rpkm=min(rpkm.df[rpkm.df!=0])/2
    rpkm.df=log10(rpkm.df+pseudo.rpkm)
  }
  quantile.norm.rpkm.df=quantile.normalisation(rpkm.df)
  mean.norm.rpkm.df=rpkm.df
  if(quant.norm){
    mean.norm.rpkm.df=quantile.norm.rpkm.df
  }
  #mean.norm.rpkm.df=get.mean.normalize.exprs(quantile.normalisation(df = rpkm.df))
  other.genes.vec=setdiff(rownames(rpkm.df),genes.sig.vec)
  samples=colnames(rpkm.df)
  len.samples=length(samples)
  out.scores.vec=c()
  for(n in 1:len.samples){
    sample=samples[n]
    other.genes.mean.norm.exprs=as.numeric(mean.norm.rpkm.df[other.genes.vec,sample])
    sig.genes.vec=as.numeric(mean.norm.rpkm.df[genes.sig.vec,sample])
    if(all(is.na(sig.genes.vec))){
      next
    }
    sig.genes.vec=as.numeric(ifelse(is.na(sig.genes.vec),0,sig.genes.vec))
    #sig.score=as.numeric(mean(sig.genes.vec)-mean(other.genes.mean.norm.exprs))
    sig.score=as.numeric(mean(sig.genes.vec))
    out.scores.vec=append(out.scores.vec,sig.score,after = length(out.scores.vec))
  }
  return(out.scores.vec)
}
get.gene.signature.score.per.cell=function(rpkm.df,genes.sig.vec,log.rpkm=F){
  if(log.rpkm){
    #pseudo.rpkm=min(rpkm.df[rpkm.df!=0])/2
    #rpkm.df=log10(rpkm.df+pseudo.rpkm)
    rpkm.df=log2(rpkm.df)
  }
  #mean.norm.rpkm.df=get.mean.normalize.exprs(rpkm.df)
  mean.norm.rpkm.df=quantile.normalisation(df = rpkm.df)
  genes.sig.vec=intersect(rownames(rpkm.df),genes.sig.vec)
  other.genes.vec=setdiff(rownames(rpkm.df),genes.sig.vec)
  samples=colnames(rpkm.df)
  len.samples=length(samples)
  len.genes.sig.vec=length(genes.sig.vec)
  out.mat=matrix(nrow = len.genes.sig.vec,ncol =len.samples )
  out.scores.vec=c()
  for(n in 1:len.samples){
    sample=samples[n]
    other.genes.mean.norm.exprs=mean(as.numeric(mean.norm.rpkm.df[other.genes.vec,sample]))
    len.genes.sig.vec=length(genes.sig.vec)
    for(m in 1:len.genes.sig.vec){
      sig.gene=genes.sig.vec[m]
      sig.gene.norm.exprs=as.numeric(mean.norm.rpkm.df[sig.gene,sample])
      if(is.na(sig.gene.norm.exprs)){
        show(c(sig.gene,sig.gene.norm.exprs))
      }
      sig.gene.norm.exprs=ifelse(is.na(sig.gene.norm.exprs),
                                 other.genes.mean.norm.exprs,sig.gene.norm.exprs)
      #sig.gene.score=sig.gene.norm.exprs
      sig.gene.score=sig.gene.norm.exprs-other.genes.mean.norm.exprs
      out.mat[m,n]=sig.gene.score
    }
  }
  out.df=data.frame(out.mat)
  colnames(out.df)=samples
  rownames(out.df)=genes.sig.vec
  return(out.df)
}
quantile.normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  df_final=data.frame(df_final)
  return(df_final)
}
plot.dominant.marker.per.timepoint=function(rpkm.df,meta.df,markers.df){
  rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  #rpkm.df=get.mean.normalize.exprs(rpkm.df =  rpkm.df)
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df$development.stage=factor(meta.df$development.stage,levels  = as.factor(c('R','LR','ET','T','ES','S')))
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes  )
  show(length(markers.vec))
  markers.df=markers.df[markers.vec,]
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df[is.na(in.rpkm.df)]=0
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  meta.df=meta.df[colnames(in.rpkm.df),]
  markers.df=markers.df[rownames(in.rpkm.df),]
  markers.list=split(markers.df,f=markers.df$category)
  expressed.markers.vec=rownames(markers.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=split(meta.df,f = meta.df$development.stage)
  development.stages=names(meta.list)
  #ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  #intersect.development.stages = intersect(ordered.development.stages,development.stages)
  intersect.development.stages= development.stages
  len.intersect.development.stages=length(intersect.development.stages)
  markers.mean.expression=c()
  log.marker.mean.expression.vec=c()
  markers.stages.vec=c()
  markers.names=c()
  fraction.detected.in.vec=c()
  samples.category.vec=c()
  markers.categories=names(markers.list)
  #markers.categories=c('Glycolytic pathway', 'Cytoplasmic Translation machinery','DNA replication','Merozoite Invasion')
  len.markers.categories=length(markers.categories)
  #categories.to.include=c('Glycolytic pathway', 'Cytoplasmic Translation machinery','DNA replication','Merozoite Invasion')
  #len.markers.categories=length(len.markers.categorie)
  for(n in 1:len.markers.categories){
    markers.category=markers.categories[n]
    #categories.to.skip=c('Deoxynucleotide synthesis','Mitochondrial','Ribonucleotide synthesis','TCA cycle','none')
    #categories.to.include=c('Glycolytic pathway', 'Cytoplasmic Translation machinery','DNA replication','Merozoite Invasion')
    categories.to.skip=c('none','Ribonucleotide synthesis')
    #categories.to.skip=c('Deoxynucleotide synthesis','Mitochondrial','Ribonucleotide synthesis','none')
    #categories.to.skip=c('none')
    if(markers.category %in% categories.to.skip){
      next
    }
    marker.category.peaking.stages=as.character(unlist(strsplit(markers.peaking.stages.list[[markers.category]],split  = ',')))
    temp.marker.df=markers.list[[markers.category]]
    temp.expressed.markers.vec=rownames(temp.marker.df)
    len.temp.expressed.markers.vec=length(temp.expressed.markers.vec)
    no.markers=length(temp.expressed.markers.vec)
    winning.category.vec=c()
    for(m in 1:len.temp.expressed.markers.vec){
      marker=temp.expressed.markers.vec[m]
      temp.markers.mean.expression=c()
      marker.stage=as.character(markers.df[marker,'category'])
      marker.name=as.character(markers.df[marker,'X.Gene.Name.or.Symbol.'])
      for( n in 1:len.intersect.development.stages){
        intersect.development.stage=intersect.development.stages[n]
        temp.meta.df=meta.list[[intersect.development.stage]]
        stage.samples=rownames(temp.meta.df)
        stage.rpkm.vec=as.numeric(as.character(in.rpkm.df[marker,stage.samples]))
        fraction.detected.in=length(stage.rpkm.vec[stage.rpkm.vec>0])/length(stage.samples)
        marker.mean.expression=mean(stage.rpkm.vec)*fraction.detected.in
        #marker.mean.expression=mean(stage.rpkm.vec)
        marker.median.expression=median(stage.rpkm.vec)
        temp.markers.mean.expression=append(temp.markers.mean.expression,marker.mean.expression,after = length(temp.markers.mean.expression))
      }
      top.mean.expr=max(temp.markers.mean.expression)
      winning.stage=intersect.development.stages[which(temp.markers.mean.expression==top.mean.expr)]
      winning.category.vec=append(winning.category.vec,winning.stage,length(winning.category.vec))
    }
    winning.df=data.frame(table(winning.category.vec))
    rownames(winning.df)=as.character(winning.df[,'winning.category.vec'])
    tran.winning.df=t(winning.df)
    non.winning.stages=setdiff(c('R','LR','ET','T','ES','S'),rownames(winning.df))
    len.non.winning.stages=length(non.winning.stages)
    freq.vec=c()
    for(l in 1:len.intersect.development.stages){
      winning.intersect.development.stages.score=intersect.development.stages[l]
      temp.vec=as.numeric(winning.df[winning.intersect.development.stages.score,][2])
      temp.vec[is.na(temp.vec)]=0
      freq.vec=append(freq.vec,temp.vec,length(freq.vec))
    }
    final.winning.df=data.frame(stage=intersect.development.stages,Freq=as.numeric(freq.vec))
    rownames(final.winning.df)=intersect.development.stages
    #non.detected.stages=setdiff(ordered.development.stages,rownames(winning.df))
    non.detected.stages=setdiff(intersect.development.stages,rownames(winning.df))
    len.non.detected.stages=length(non.detected.stages)
    winning.df[,'Freq']=as.numeric(winning.df[,'Freq']/no.markers)
    final.winning.df[,'Freq']=as.numeric(final.winning.df[,'Freq']/no.markers)
    #winning.df$winning.category.vec=factor(winning.df$winning.category.vec,as.factor(c('R','LR','ET','LT','ES','S')))
    #temp.markers.dotplot.df$ordered.samples.categories <- factor(temp.markers.dotplot.df$samples.category, as.character(c('ring','late.ring','early.trophozoite','late.trophozoite','early.schizont','schizont')))
    #winning.df$winning.category.vec.temp=factor(winning.df$winning.category.vec,as.factor(intersect(ordered.development.stages,as.character(winning.df$winning.category.vec))))
    intersect.names=intersect(c('R','LR','ET','LT','ES','S'),rownames(winning.df))
    winning.df=winning.df[intersect.names,]
    title.str=paste(markers.category,paste('(n=',no.markers,')',sep =''),sep='')
    names.label=rownames(final.winning.df)
    col.vec=ifelse(names.label %in% marker.category.peaking.stages,'lightblue','gray')
    #barplot(height = winning.df[,'Freq'],las=2,main =title.str ,names.arg =names.label,cex.names  = .8)
    barplot(height = final.winning.df[,'Freq'],main =title.str,cex.names  = .8,names.arg =names.label,col =col.vec,ylim=c(0,1),las=1,border = NA)
    #barplot(height = final.winning.df[,'Freq'],main =title.str,cex.names  = 1.5,names.arg =names.label,col =col.vec)
    #barplot(height = final.winning.df[,'Freq'],las=2,main ='',cex.names  = 1.5,names.arg =rep('',times=length(names.label)),col =col.vec ,ylim=c(0,1))
  }
  plot.new()
  legend('center',legend=rep('',times = 2),fill =c('lightblue','gray') ,box.lty = 0,cex = 2.5)
  return(markers.df)
}
plot.dominant.marker.per.sub.population=function(rpkm.df,meta.df,markers.df){
  sub.pop.col.list=get.color.list.for.pheatmap(in.vec = as.character(meta.df$markers.cluster.groups))
  rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  pseudo.value=min(rpkm.df[rpkm.df!=0])/2
  rpkm.df=rpkm.df+pseudo.value
  #rpkm.df=get.mean.normalize.exprs(rpkm.df =  rpkm.df)
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  #meta.df$markers.cluster.groups=factor(meta.df$markers.cluster.groups,levels  = as.factor(sub.populations.grp.names))
  meta.df$markers.cluster.groups=factor(meta.df$markers.cluster.groups,levels  = as.factor(names(sub.populations.grp.order.names.vec)))
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes  )
  markers.df=markers.df[markers.vec,]
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df[is.na(in.rpkm.df)]=0
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  meta.df=meta.df[colnames(in.rpkm.df),]
  markers.df=markers.df[rownames(in.rpkm.df),]
  markers.list=split(markers.df,f=markers.df$category)
  expressed.markers.vec=rownames(markers.df)
  meta.list=split(meta.df,f =meta.df$markers.cluster.groups)
  development.stages=names(meta.list)
  intersect.development.stages= development.stages
  len.intersect.development.stages=length(intersect.development.stages)
  markers.mean.expression=c()
  log.marker.mean.expression.vec=c()
  markers.stages.vec=c()
  markers.names=c()
  fraction.detected.in.vec=c()
  samples.category.vec=c()
  markers.categories=names(markers.list)
  #markers.categories=c('Glycolytic pathway', 'Cytoplasmic Translation machinery','DNA replication','Merozoite Invasion')
  #categories.to.include=c('Glycolytic pathway', 'Cytoplasmic Translation machinery','DNA replication','Merozoite Invasion')
  #categories.to.include=c('Actin myosin motors','Early ring transcripts' ,'Merozoite Invasion' ,'Proteasome')
  #markers.categories=intersect(markers.categories,categories.to.include)
  len.markers.categories=length(markers.categories)
  #len.markers.categories=length(len.markers.categorie)
  for(n in 1:len.markers.categories){
    markers.category=markers.categories[n]
    #categories.to.skip=c('Deoxynucleotide synthesis','Mitochondrial','Ribonucleotide synthesis','TCA cycle','none')
    #categories.to.include=c('Glycolytic pathway', 'Cytoplasmic Translation machinery','DNA replication','Merozoite Invasion')
    categories.to.skip=c('none','Ribonucleotide synthesis')
    #categories.to.skip=c('Deoxynucleotide synthesis','Mitochondrial','Ribonucleotide synthesis','none')
    if(markers.category %in% categories.to.skip){
      next
    }
    marker.category.peaking.stages=as.character(unlist(strsplit(markers.peaking.stages.list[[markers.category]],split  = ',')))
    temp.marker.df=markers.list[[markers.category]]
    temp.expressed.markers.vec=rownames(temp.marker.df)
    len.temp.expressed.markers.vec=length(temp.expressed.markers.vec)
    no.markers=length(temp.expressed.markers.vec)
    winning.category.vec=c()
    for(m in 1:len.temp.expressed.markers.vec){
      marker=temp.expressed.markers.vec[m]
      temp.markers.mean.expression=c()
      marker.stage=as.character(markers.df[marker,'category'])
      marker.name=as.character(markers.df[marker,'X.Gene.Name.or.Symbol.'])
      for( n in 1:len.intersect.development.stages){
        intersect.development.stage=intersect.development.stages[n]
        temp.meta.df=meta.list[[intersect.development.stage]]
        if(dim(temp.meta.df)[1]<2){
          next
        }
        stage.samples=rownames(temp.meta.df)
        stage.rpkm.vec=as.numeric(as.character(in.rpkm.df[marker,stage.samples]))
        fraction.detected.in=length(stage.rpkm.vec[stage.rpkm.vec>0])/length(stage.samples)
        #marker.mean.expression=mean(stage.rpkm.vec)*fraction.detected.in
        marker.mean.expression=mean(stage.rpkm.vec)
        #marker.median.expression=median(stage.rpkm.vec)
        temp.markers.mean.expression=append(temp.markers.mean.expression,marker.mean.expression,after = length(temp.markers.mean.expression))
      }
      top.mean.expr=max(temp.markers.mean.expression)
      winning.stage=intersect.development.stages[which(temp.markers.mean.expression==top.mean.expr)]
      winning.category.vec=append(winning.category.vec,winning.stage,length(winning.category.vec))
    }
    winning.df=data.frame(table(winning.category.vec))
    rownames(winning.df)=as.character(winning.df[,'winning.category.vec'])
    tran.winning.df=t(winning.df)
    non.winning.stages=setdiff(sub.populations.grp.names,rownames(winning.df))
    len.non.winning.stages=length(non.winning.stages)
    freq.vec=c()
    for(l in 1:len.intersect.development.stages){
      winning.intersect.development.stages.score=intersect.development.stages[l]
      temp.vec=as.numeric(winning.df[winning.intersect.development.stages.score,][2])
      temp.vec[is.na(temp.vec)]=0
      freq.vec=append(freq.vec,temp.vec,length(freq.vec))
    }
    final.winning.df=data.frame(stage=intersect.development.stages,Freq=as.numeric(freq.vec))
    rownames(final.winning.df)=intersect.development.stages
    #non.detected.stages=setdiff(ordered.development.stages,rownames(winning.df))
    non.detected.stages=setdiff(intersect.development.stages,rownames(winning.df))
    len.non.detected.stages=length(non.detected.stages)
    winning.df[,'Freq']=as.numeric(winning.df[,'Freq']/no.markers)
    final.winning.df[,'Freq']=as.numeric(final.winning.df[,'Freq']/no.markers)
    intersect.names=intersect(sub.populations.grp.names,rownames(winning.df))
    winning.df=winning.df[intersect.names,]
    title.str=paste(markers.category,paste('(n=',no.markers,')',sep =''),sep='')
    names.label=rownames(final.winning.df)
    bar.col.vec=sub.pop.col.list[names.label]
    col.vec=bar.col.vec
    #barplot(height = winning.df[,'Freq'],las=2,main =title.str ,names.arg =names.label,cex.names  = .8)
    barplot(height = final.winning.df[,'Freq'],main =title.str,cex.names  = .3,names.arg =rep('',times = length(names.label)),col =col.vec,ylim=c(0,1),las=1)
    #barplot(height = final.winning.df[,'Freq'],main ='',cex.names  = .3,names.arg =rep('',times = length(names.label)),col =col.vec,ylim=c(0,1),las=1)
  }
}
plot.marker.score.per.timepoint=function(rpkm.df,meta.df,markers.df){
  rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df = rpkm.df,sample.count = 3)
  #rpkm.df=quantile.normalisation(df = rpkm.df)
  rpkm.df=get.mean.normalize.exprs(rpkm.df = rpkm.df)
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df$development.stage=factor(meta.df$development.stage,levels  = as.factor(c('R','LR','ET','T','ES','S')))
  all.genes=rownames(rpkm.df)
  all.markers.vec=rownames(markers.df)
  markers.vec=intersect(x=all.markers.vec,y =all.genes  )
  non.marker.genes.vec=setdiff(x = all.genes,all.markers.vec)
  markers.df=markers.df[markers.vec,]
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df[is.na(in.rpkm.df)]=0
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  meta.df=meta.df[colnames(in.rpkm.df),]
  markers.df=markers.df[rownames(in.rpkm.df),]
  markers.list=split(markers.df,f=markers.df$category)
  expressed.markers.vec=rownames(markers.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=split(meta.df,f = meta.df$development.stage)
  development.stages=names(meta.list)
  #ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  #intersect.development.stages = intersect(ordered.development.stages,development.stages)
  intersect.development.stages= development.stages
  len.intersect.development.stages=length(intersect.development.stages)
  markers.mean.expression=c()
  log.marker.mean.expression.vec=c()
  markers.stages.vec=c()
  markers.names=c()
  fraction.detected.in.vec=c()
  samples.category.vec=c()
  markers.categories=names(markers.list)
  len.markers.categories=length(markers.categories)
  for(n in 1:len.markers.categories){
    all.samples.markers.score.list=list()
    markers.category=markers.categories[n]
    #categories.to.skip=c('Deoxynucleotide synthesis','Mitochondrial','Ribonucleotide synthesis','TCA cycle','none')
    categories.to.skip=c('none','Ribonucleotide synthesis')
    #categories.to.skip=c('Deoxynucleotide synthesis','Mitochondrial','Ribonucleotide synthesis','none')
    if(markers.category %in% categories.to.skip){
      next
    }
    temp.marker.df=markers.list[[markers.category]]
    temp.expressed.markers.vec=rownames(temp.marker.df)
    len.temp.expressed.markers.vec=length(temp.expressed.markers.vec)
    no.markers=length(temp.expressed.markers.vec)
    for( f in 1:len.intersect.development.stages){
        intersect.development.stage=intersect.development.stages[f]
        temp.meta.df=meta.list[[intersect.development.stage]]
        stage.samples=rownames(temp.meta.df)
        stage.marker.scores.vec=c()
        len.stage.samples=length(stage.samples)
        for(s in 1:len.stage.samples){
          stage.sample=stage.samples[s]
          sig.markers.expr=rpkm.df[temp.expressed.markers.vec,stage.sample]
          temp.markers.expr.mean=mean(sig.markers.expr)
          other.genes.expr.mean=mean(rpkm.df[non.marker.genes.vec,stage.sample])
          marker.category.score=temp.markers.expr.mean-other.genes.expr.mean
          random.genes.score.vec=c()
          for(h in 1:100){
            random.genes=sample(x = non.marker.genes.vec,size = len.temp.expressed.markers.vec)
            random.genes.expr.score=rpkm.df[random.genes,stage.sample]-other.genes.expr.mean
            #random.genes.score=random.genes.expr.mean-other.genes.expr.mean
            random.genes.score.vec=append(random.genes.score.vec,random.genes.expr.score,length(random.genes.score.vec))
          }
          var.test.res=var.test(x = sig.markers.expr,random.genes.score.vec)
          if(var.test.res$p.value <=0.05){
            stage.marker.scores.vec=append(stage.marker.scores.vec,marker.category.score,length(stage.marker.scores.vec))
          }
        }
        all.samples.markers.score.list[[intersect.development.stage]]=stage.marker.scores.vec
        #log.stage.marker.scores.vec=ifelse(stage.marker.scores.vec>0,log2(stage.marker.scores.vec),-log2(abs(stage.marker.scores.vec)))
        #all.samples.markers.score.list[[intersect.development.stage]]=log.stage.marker.scores.vec
    }
    boxplot2(all.samples.markers.score.list,main=markers.category)
  }
}
plot.two.corresponding.samples.scatter=function(first.df,sec.df){
  first.df=filter.none.expressed.genes(input.data = first.df)
  sec.df=filter.none.expressed.genes(input.data = sec.df)
  pseudo.value=min(min(first.df[first.df>0]),min(sec.df[sec.df>0]))/2
  samples=intersect(colnames(first.df),colnames(sec.df))
  len.samples=length(samples)
  genes=intersect(rownames(first.df),rownames(sec.df))
  first.df=log2(first.df[genes,]+pseudo.value)
  sec.df=log2(sec.df[genes,]+pseudo.value)
  for(n in 1:len.samples){
    sample=samples[n]
    x.points=as.numeric(as.character(first.df[,sample]))
    y.points=as.numeric(as.character(sec.df[,sample]))
    cor.score=round(cor(x.points,y = y.points,method = 'pearson'),2)
    title.str=paste(sample,cor.score,sep=': ')
    #geneScatterplot(x =x.points,y =  y.points,xlab = paste(sample,'rpkm',sep='_'),ylab = paste(sample,'deseq2',sep='_'),main.title = title.str,col = '#00207040')
    plot(x=x.points,y=y.points,pch=19,cex=.5,main=title.str,xlab=paste(sample,'rpkm',sep='_'),ylab = paste(sample,'rpm',sep='_'))
  }
}
plot.markers.expression.density=function(rpkm.df,markers.df,in.title.str){
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes  )
  markers.df=markers.df[markers.vec,]
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  create.markers.density.distribution.plot.df(df = in.rpkm.df,markers.df =markers.df,title.str = in.title.str,log = T)
}
#Sample and gene distribution
create.markers.density.distribution.plot.df=function(df,markers.df,log=F,sample.dist=F,filter.zeros=F,title.str){
  x.lab ='RPKM'
  y.lab='Density'
  main.label = title.str
  if(log){
    df=log2(df+1)
    x.lab = 'log2(RPKM+1)'
  }
  if(sample.dist){
    main.label =title.str
    df=t(df)
  }
  density.list = apply(df, 1, function(x){x=x[x!=-Inf];x[x!=Inf];if(length(x)>=2){density(x)}})
  density.list.length=length(density.list)
  if(filter.zeros){
    density.list = apply(df, 1, function(x){x=x[x!=-Inf];x[x!=Inf];if(length(x)>=2){density(x)}})
  }
  n.samples = length(density.list)
  xlim = range(unlist(lapply(density.list, '[[', 'x')))
  ylim = range(unlist(lapply(density.list, '[[', 'y')))
  plot(x = xlim, y = ylim, main=main.label, xlab = x.lab,ylab=y.lab,type='n')
  colour.factor.list=get.col.factor(col.factor=as.character(markers.df$category))
  colour.list=colour.factor.list[['col.str']]
  for(j.sample in 1:n.samples){
    lines(density.list[[j.sample]],col=colour.list[j.sample])
  }
}
encode.fraction.detected.in.category=function(fraction.detected.vec,interval=.1){
  intervals=seq(0,1,by = interval)
  category.vec=c()
  len.fraction.detected.vec=length(fraction.detected.vec)
  for(m in 1:len.fraction.detected.vec){
    fraction.detected=fraction.detected.vec[m]
  }
}
plot.surface.markers.mean.expression.dotplot.per.timepoint=function(rpkm.df,meta.df,markers.df){
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes  )
  markers.df=markers.df[markers.vec,]
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.genes(in.rpkm.df)
  meta.df=meta.df[colnames(in.rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df =meta.df )
  markers.df=markers.df[rownames(in.rpkm.df),]
  markers.list=split(markers.df,f=markers.df$category)
  expressed.markers.vec=rownames(markers.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=split(meta.df,f = meta.df$development.stage)
  development.stages=names(meta.list)
  #ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  intersect.development.stages = intersect(ordered.development.stages,development.stages)
  len.intersect.development.stages=length(intersect.development.stages)
  markers.mean.expression=c()
  log.marker.mean.expression.vec=c()
  markers.stages.vec=c()
  markers.names=c()
  fraction.detected.in.vec=c()
  samples.category.vec=c()
  for(m in 1:len.expressed.markers.vec){
    marker=expressed.markers.vec[m]
    marker.stage=as.character(markers.df[marker,'category'])
    marker.name=as.character(markers.df[marker,'X.Gene.Name.or.Symbol.'])
    marker.name=marker
    #marker.name=ifelse(marker.name=='null',marker,marker.name)
    for( n in 1:len.intersect.development.stages){
      intersect.development.stage=intersect.development.stages[n]
      #show(intersect.development.stage)
      temp.meta.df=meta.list[[intersect.development.stage]]
      stage.samples=rownames(temp.meta.df)
      stage.rpkm.vec=as.numeric(as.character(in.rpkm.df[marker,stage.samples]))
      marker.mean.expression=mean(stage.rpkm.vec)
      marker.median.expression=median(stage.rpkm.vec)
      log.marker.mean.expression=log10(marker.mean.expression+1)
      #log.marker.mean.expression=log10(marker.median.expression+1)
      log.marker.mean.expression.vec=append(log.marker.mean.expression.vec,log.marker.mean.expression,length(log.marker.mean.expression.vec))
      fraction.detected.in=length(stage.rpkm.vec[stage.rpkm.vec>0])/length(stage.samples)
      markers.mean.expression=append(markers.mean.expression,marker.mean.expression,after = length(marker.mean.expression))
      markers.stages.vec=append(markers.stages.vec,marker.stage,length(markers.stages.vec))
      fraction.detected.in.vec=append(fraction.detected.in.vec,fraction.detected.in,after = length(fraction.detected.in.vec))
      markers.names=append(markers.names,marker.name,length(markers.names))
      samples.category.vec=append(samples.category.vec,intersect.development.stage,length(samples.category.vec))
    }
  }
  markers.dotplot.df=data.frame(names=markers.names,stages=markers.stages.vec,mean.rpkm=markers.mean.expression,fraction.detected=fraction.detected.in.vec,log.mean.expression=log.marker.mean.expression.vec,samples.category=samples.category.vec)
  markers.dotplot.df.list=split(markers.dotplot.df,f=markers.dotplot.df$stages)
  markers.dotplot.list.names=names(markers.dotplot.df.list)
  len.markers.dotplot.list.names=length(markers.dotplot.list.names)
  for(l in 1:len.markers.dotplot.list.names){
    markers.dotplot.list.name=markers.dotplot.list.names[l]
    temp.markers.dotplot.df=markers.dotplot.df.list[[markers.dotplot.list.name]]
    temp.markers.dotplot.df$ordered.samples.categories <- factor(temp.markers.dotplot.df$samples.category, as.character(c('R','LR','ET','T','ES','S')))
    #temp.markers.dotplot.df=subset(temp.markers.dotplot.df,fraction.detected>=.3)
    no.markers=unique(as.character(temp.markers.dotplot.df[,'names']))
    len.no.markers=length(no.markers)
    if(len.no.markers<=10){
      point.size=as.numeric(as.character(temp.markers.dotplot.df$fraction.detected))
      title.str=convert.to.title.case(in.str = markers.dotplot.list.name)
      temp.ggplot=ggplot(temp.markers.dotplot.df, aes(x = ordered.samples.categories, y = names)) + geom_point(aes(size=log.mean.expression,colour=fraction.detected)) + ggtitle(label = title.str)+ylab(label = '')+xlab(label = '') +theme(axis.text.x = element_text(angle = 90, hjust = 1))
      #temp.ggplot=ggplot(temp.markers.dotplot.df, aes(x = ordered.samples.categories, y = names)) + geom_point(aes(size=log.mean.expression,colour=fraction.detected)) + ggtitle(label = '')+ylab(label = '')+xlab(label = '')
      print(temp.ggplot)
    }
    else{
      chunked.vec=chunk.non.randomly(in.vec = no.markers,chunk.sizes = 10)
      names.chunked.vec=names(chunked.vec)
      len.names.chunked.vec=length(names.chunked.vec)
      for(ch in 1:len.names.chunked.vec ){
        temp.names.chunked=names.chunked.vec[ch]
        temp.chunked.vec=chunked.vec[[temp.names.chunked]]
        markers.chunked.df=temp.markers.dotplot.df[which(temp.markers.dotplot.df[,'names'] %in% temp.chunked.vec),]
        point.size=as.numeric(as.character(markers.chunked.df$fraction.detected))
        title.str=convert.to.title.case(in.str = markers.dotplot.list.name)
        temp.ggplot=ggplot(markers.chunked.df, aes(x = ordered.samples.categories, y = names)) + geom_point(aes(size=log.mean.expression,colour=fraction.detected)) + ggtitle(label = title.str)+ylab(label = '')+xlab(label = '') +theme(axis.text.x = element_text(angle = 90, hjust = 1))
        #temp.ggplot=ggplot(markers.chunked.df, aes(x = ordered.samples.categories, y = names)) + geom_point(aes(size=log.mean.expression,colour=fraction.detected)) + ggtitle(label = '')+ylab(label = '')+xlab(label = '') +theme(axis.text.x = element_text(angle = 90, hjust = 1))
        print(temp.ggplot)
      }
    }
  }
  return(in.rpkm.df)
  #return(markers.dotplot.df)
}
plot.surface.markers.expression.heatmap.per.sample=function(rpkm.df,meta.df,markers.df){
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df[colnames(in.rpkm.df),])
  markers.df=markers.df[rownames(in.rpkm.df),]
  markers.col.factor=get.col.factor(col.factor = as.character(markers.df$category))
  markers.df$col=markers.col.factor$col.str
  markers.list=split(markers.df,f=markers.df$category)
  expressed.markers.vec=rownames(markers.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=split(meta.df,f = meta.df$development.stage)
  development.stages=names(meta.list)
  #ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  intersect.development.stages = intersect(ordered.development.stages,development.stages)
  len.intersect.development.stages=length(intersect.development.stages)
  for(m in 1:len.intersect.development.stages){
    intersect.development.stage=intersect.development.stages[m]
    temp.meta.df=meta.list[[intersect.development.stage]]
    if(intersect.development.stage=='R'){
      next
    }
    #temp.rpkm.df=filter.none.expressed.samples(filter.rpkm.less.than.cutoff(df = in.rpkm.df[,rownames(temp.meta.df)],rpkm.cutoff =1,no.samples = 1))
    temp.rpkm.df=filter.none.expressed.samples(df = subset(in.rpkm.df, select = rownames(temp.meta.df)))
    show(dim(temp.rpkm.df))
    ordered.markers.df=markers.df[rownames(temp.rpkm.df),]
    ordered.markers.df=ordered.markers.df[order(ordered.markers.df$category),]
    temp.rpkm.df=temp.rpkm.df[rownames(ordered.markers.df),]
    gene.col=as.character(markers.df[rownames(temp.rpkm.df),'col'])
    col.str=gene.col
    temp.rpkm.mat=as.matrix(log10(temp.rpkm.df+1))
    title.str=convert.to.title.case(in.str = gsub(pattern = '\\.',replacement = ' ',x = intersect.development.stage))
    #title.str=''
    #heatmap.2(x = temp.rpkm.mat,cellnote =round(temp.rpkm.mat,3),notecex=.1 ,trace='none',margins = c(13,13),main=title.str,dendrogram = 'col',RowSideColors = col.str,cexRow = .5,cexCol = .4,Rowv = NA,key.title = '',key.xlab = '',key.ylab = '',col=greenred(100))
    #heatmap.2(x = temp.rpkm.mat,cellnote =round(temp.rpkm.mat,3),notecex=.1 ,trace='none',margins = c(10,10),main=intersect.development.stage,dendrogram = 'col',cexRow = .5,cexCol = .4,key.title = 'log10(norm.counts+1)',key.xlab = '',key.ylab = '')
    row.annotation.df=subset(ordered.markers.df,select = category)
    #row.annotation.col.map=get.color.list.for.pheatmap(in.vec = as.character(row.annotation.df$category))
    row.annotation.col.map=surface.markers.col.list
    row.annotation.col.list=list(category=row.annotation.col.map)
    if(intersect.development.stage=='R'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,col=bluered(100),scale = 'column',fontsize_row = 2,fontsize_col = 2,cluster_rows = F,cluster_cols  = F,annotation_row = row.annotation.df,cellheight = 20,cellwidth = 20,annotation_colors = row.annotation.col.list,annotation_legend = F)
      next
      #pheatmap(mat = temp.rpkm.mat,main=title.str,col=bluered(100),scale = 'column',fontsize_row = 2,fontsize_col = 2,cluster_rows = F,cluster_cols  = F,annotation_row = row.annotation.df,cellheight = 20,cellwidth = 20,annotation_colors = row.annotation.col.list,annotation_legend = F)
    }
    if(intersect.development.stage=='LR'){
      pheatmap(mat = temp.rpkm.mat,main=title.str,col=bluered(100),scale = 'column',fontsize_row = 2,fontsize_col = 2,cluster_rows = F,cluster_cols = F,annotation_colors = row.annotation.col.list,annotation_row = row.annotation.df,cellheight = 2,cellwidth = 20,annotation_legend = F)
    }
    if(intersect.development.stage=='ET'){
      #next
      pheatmap(mat = temp.rpkm.mat,main=title.str,col=bluered(100),scale = 'column',fontsize_row = 2,fontsize_col = 2,cluster_rows = F,cluster_cols = F,annotation_colors = row.annotation.col.list,annotation_row = row.annotation.df,cellheight = 2,cellwidth = 20,annotation_legend = F)
    }
    if(intersect.development.stage=='T'){
      pheatmap(mat = temp.rpkm.mat,main=title.str,col=bluered(100),scale = 'column',fontsize_row = 1,fontsize_col = 1,cluster_rows = F,cluster_cols = F,annotation_colors = row.annotation.col.list,annotation_row = row.annotation.df,cellheight = 3,cellwidth = 5,annotation_legend = F)
    }
    if(intersect.development.stage=='ES'){
      pheatmap(mat = temp.rpkm.mat,main=title.str,col=bluered(100),scale = 'column',fontsize_row = 1,fontsize_col = 1,cluster_rows = F,cluster_cols = F,annotation_colors = row.annotation.col.list,annotation_row = row.annotation.df,cellheight = 20,cellwidth = 30,annotation_legend = F)
    }
    if(intersect.development.stage=='S'){
      pheatmap(mat = temp.rpkm.mat,main=title.str,col=bluered(100),scale = 'column',fontsize_row = 1,fontsize_col = 1,cluster_rows = F,cluster_cols = F,annotation_colors = row.annotation.col.list,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 10,annotation_legend = F)
    }
    #plot.new()
    legend.str=unique(as.character(markers.df[rownames(temp.rpkm.df),'category']))
    legend.col=unique(as.character(markers.df[rownames(temp.rpkm.df),'col']))
    #legend('topright',legend=legend.str,fill=legend.col,cex=1.0,box.lty = 0)
    plot.new()
    legend('center',legend=names(surface.markers.col.list),fill=surface.markers.col.list[names(surface.markers.col.list)],cex=1.5,box.lty = 0)
    show(names(row.annotation.col.list))
    #legend('topright',legend=legend.str,fill=legend.col,cex=1.0,box.lty = 0)
    #legend('center',legend=rep('',times = length(legend.str)),fill=legend.col,cex=1.5,box.lty = 0)
  }
}
plot.surface.markers.mean.expression.heatmap.per.sub.population=function(rpkm.df,meta.df,markers.df){
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df[colnames(in.rpkm.df),])
  markers.df=markers.df[rownames(in.rpkm.df),]
  markers.col.factor=get.col.factor(col.factor = as.character(markers.df$category))
  markers.df$col=markers.col.factor$col.str
  markers.list=split(markers.df,f=markers.df$category)
  expressed.markers.vec=rownames(markers.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=split(meta.df,f = meta.df$markers.cluster.groups)
  development.stages=names(meta.list)
  ordered.sub.populations=c('SP.1', 'SP.2' ,'SP.3' ,'SP.4' ,'SP.5', 'SP.6', 'SP.7', 'SP.8')
  intersect.development.stages = intersect(ordered.sub.populations,development.stages)
  len.intersect.development.stages=length(intersect.development.stages)
  out.mat=matrix(nrow = dim(in.rpkm.df)[1],ncol = len.intersect.development.stages)
  for(m in 1:len.intersect.development.stages){
    intersect.development.stage=intersect.development.stages[m]
    temp.meta.df=meta.list[[intersect.development.stage]]
    temp.rpkm.df=subset(in.rpkm.df, select = rownames(temp.meta.df))
    out.mat[,m]=as.numeric(apply(temp.rpkm.df,1,mean))
  }
  rownames(out.mat)=rownames(in.rpkm.df)
  colnames(out.mat)=intersect.development.stages
  ordered.markers.df=markers.df[rownames(in.rpkm.df),]
  ordered.markers.df=ordered.markers.df[order(ordered.markers.df$category),]
  out.mat=out.mat[rownames(ordered.markers.df),]
  pseudo.val=min(out.mat[out.mat>0])/2
  out.mat=log2(out.mat+pseudo.val)
  gene.col=as.character(markers.df[rownames(out.mat),'col'])
  col.str=gene.col
  title.str=convert.to.title.case(in.str = gsub(pattern = '\\.',replacement = ' ',x = 'Sub-populations'))
  row.annotation.df=subset(ordered.markers.df,select = category)
  #row.annotation.col.map=get.color.list.for.pheatmap(in.vec = as.character(row.annotation.df$category))
  row.annotation.col.map=surface.markers.col.list
  col.annotation.df=data.frame(sp=colnames(out.mat))
  rownames(col.annotation.df)=colnames(out.mat)
  col.col.list=get.color.list.for.pheatmap(in.vec = colnames(out.mat))
  annotation.col.list=list(sp=col.col.list,category=row.annotation.col.map)
  pheatmap(mat = out.mat,main=title.str,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_col = 5,cluster_rows = F,cluster_cols = T,annotation_colors = annotation.col.list,annotation_row = row.annotation.df,annotation_col =  col.annotation.df,cellheight = 2.5,cellwidth = 35,annotation_legend = F,show_rownames = F)
  plot.new()
  legend('top',legend=names(surface.markers.col.list),fill=surface.markers.col.list[names(surface.markers.col.list)],cex=1.5,box.lty = 0)
  legend('bottom',legend=rep('',times = length(names(surface.markers.col.list))),fill=surface.markers.col.list[names(surface.markers.col.list)],cex=2.5,box.lty = 0,border = NA)
}
plot.markers.expression.heatmap.per.stage.per.sample=function(rpkm.df,meta.df,markers.df){
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  meta.df=meta.df[colnames(in.rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=split(meta.df,f = meta.df$development.stage)
  development.stages=names(meta.list)
  #ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  intersect.development.stages = intersect(ordered.development.stages,development.stages)
  len.intersect.development.stages=length(intersect.development.stages)
  for(m in 1:len.intersect.development.stages){
    intersect.development.stage=intersect.development.stages[m]
    temp.meta.df=meta.list[[intersect.development.stage]]
    temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
    if(is.null(dim(temp.rpkm.df))){
      next
    }
    if(dim(temp.rpkm.df)[1]<=0){
      next
    }
    temp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = in.rpkm.df[,rownames(temp.meta.df)]))
    show(dim(temp.rpkm.df))
    if(dim(temp.rpkm.df)[1]<=0){
      next
    }
    pseudo.value=min(temp.rpkm.df[temp.rpkm.df!=0])/2
    temp.rpkm.mat=as.matrix(log10(temp.rpkm.df+pseudo.value))
    #title.str=convert.to.title.case(gsub(pattern = '\\.',replacement = ' ',x = intersect.development.stage))
    title.str=intersect.development.stage
    title.str=''
    #temp.rpkm.mat=head(temp.rpkm.mat,100)
    temp.markers.df=markers.df[rownames(temp.rpkm.mat),]
    temp.markers.df= temp.markers.df[order(temp.markers.df$category),]
    marker.category=sort(as.character(markers.df[rownames(temp.rpkm.mat),'category']))
    marker.col.factor=get.col.factor(col.factor = marker.category)
    row.annotation.df=subset(temp.markers.df,select=category)
    if(intersect.development.stage=='R'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,col=bluered(100),scale = 'column',fontsize_row = 2,fontsize_col = 2,cluster_rows = F,cluster_cols  = F,annotation_row = row.annotation.df,cellheight = 20,cellwidth = 20,annotation_legend = F)
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2,fontsize_col = 2,cluster_rows = F,cluster_cols  = T,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2.5,fontsize_col = 2.5,cluster_rows = F,cluster_cols  = T,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
    }
    if(intersect.development.stage=='LR'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2,fontsize_col = 2,cluster_rows = F,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row =2.5,fontsize_col =2.5,cluster_rows = F,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
    }
    if(intersect.development.stage=='ET'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2,fontsize_col = 2,cluster_rows = F,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2.5,fontsize_col = 2.5,cluster_rows = F,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
    }
    if(intersect.development.stage=='T'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 1,fontsize_col = 1,cluster_rows = F,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 4,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2.5,fontsize_col = 2.5,cluster_rows = F,cellheight = 5,cellwidth = 4,annotation_legend = F,border_color = NA)
    }
    if(intersect.development.stage=='ES'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 1,fontsize_col = 1,cluster_rows = F,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2.5,fontsize_col = 2.5,cluster_rows = F,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
    }
    if(intersect.development.stage=='S'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 1,fontsize_col = 1,cluster_rows = F,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 5,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2.5,fontsize_col = 2.5,cluster_rows = F,cellheight = 5,cellwidth = 5,annotation_legend = F,border_color = NA)
    }
    #heatmap.2(x = temp.rpkm.mat,cellnote =round(temp.rpkm.mat,3),notecex=.1 ,trace='none',margins = c(13,13),main=title.str,dendrogram = 'col',cexRow = .3,cexCol = .3,Rowv = rownames(temp.markers.df),key.title = 'Log10 (Expression)',key.xlab = '',key.ylab = '',col=bluered(1000),RowSideColors = marker.col.factor$col.str)
    #legend('left',legend=marker.col.factor$legend.str,fill=marker.col.factor$legend.col,cex=.3)
  }
}
plot.markers.expression.heatmap.per.stage.per.bulk.sample=function(rpkm.df,meta.df,markers.df){
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=rownames(markers.df),y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  meta.df=meta.df[colnames(in.rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=split(meta.df,f = meta.df$development.stage)
  development.stages=names(meta.list)
  #ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  intersect.development.stages = intersect(ordered.development.stages,development.stages)
  len.intersect.development.stages=length(intersect.development.stages)
  for(m in 1:len.intersect.development.stages){
    intersect.development.stage=intersect.development.stages[m]
    temp.meta.df=meta.list[[intersect.development.stage]]
    temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
    if(is.null(dim(temp.rpkm.df))){
      next
    }
    if(dim(temp.rpkm.df)[1]<=0){
      next
    }
    temp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = in.rpkm.df[,rownames(temp.meta.df)]))
    show(dim(temp.rpkm.df))
    if(dim(temp.rpkm.df)[1]<=0){
      next
    }
    pseudo.value=min(temp.rpkm.df[temp.rpkm.df!=0])/2
    temp.rpkm.mat=as.matrix(log10(temp.rpkm.df+pseudo.value))
    #title.str=convert.to.title.case(gsub(pattern = '\\.',replacement = ' ',x = intersect.development.stage))
    title.str=intersect.development.stage
    title.str=''
    #temp.rpkm.mat=head(temp.rpkm.mat,100)
    temp.markers.df=markers.df[rownames(temp.rpkm.mat),]
    temp.markers.df= temp.markers.df[order(temp.markers.df$category),]
    marker.category=sort(as.character(markers.df[rownames(temp.rpkm.mat),'category']))
    marker.col.factor=get.col.factor(col.factor = marker.category)
    row.annotation.df=subset(temp.markers.df,select=category)
    if(intersect.development.stage=='R'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,col=bluered(100),scale = 'column',fontsize_row = 2,fontsize_col = 2,cluster_rows = F,cluster_cols  = F,annotation_row = row.annotation.df,cellheight = 20,cellwidth = 20,annotation_legend = F)
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2,fontsize_col = 2,cluster_rows = F,cluster_cols  = T,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cluster_rows = F,cluster_cols  = T,cellheight = 1,cellwidth = 15,show_rownames = F,annotation_row = row.annotation.df,annotation_legend = F,border_color = NA)
    }
    if(intersect.development.stage=='LR'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2,fontsize_col = 2,cluster_rows = F,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row =2.5,fontsize_col =2.5,cluster_rows = F,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
    }
    if(intersect.development.stage=='ET'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2,fontsize_col = 2,cluster_rows = F,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2.5,fontsize_col = 2.5,cluster_rows = F,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
    }
    if(intersect.development.stage=='T'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 1,fontsize_col = 1,cluster_rows = F,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 4,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2.5,fontsize_col = 2.5,cluster_rows = F,cellheight = 5,cellwidth = 4,annotation_legend = F,border_color = NA)
    }
    if(intersect.development.stage=='ES'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 1,fontsize_col = 1,cluster_rows = F,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2.5,fontsize_col = 2.5,cluster_rows = F,cellheight = 5,cellwidth = 10,annotation_legend = F,border_color = NA)
    }
    if(intersect.development.stage=='S'){
      #pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 1,fontsize_col = 1,cluster_rows = F,annotation_row = row.annotation.df,cellheight = 5,cellwidth = 5,annotation_legend = F,border_color = NA)
      pheatmap(mat = temp.rpkm.mat,main=title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 2.5,fontsize_col = 2.5,cluster_rows = F,cellheight = 5,cellwidth = 5,annotation_legend = F,border_color = NA)
    }
    #heatmap.2(x = temp.rpkm.mat,cellnote =round(temp.rpkm.mat,3),notecex=.1 ,trace='none',margins = c(13,13),main=title.str,dendrogram = 'col',cexRow = .3,cexCol = .3,Rowv = rownames(temp.markers.df),key.title = 'Log10 (Expression)',key.xlab = '',key.ylab = '',col=bluered(1000),RowSideColors = marker.col.factor$col.str)
    #legend('left',legend=marker.col.factor$legend.str,fill=marker.col.factor$legend.col,cex=.3)
  }
}
plot.markers.expression.heatmap.per.sample=function(rpkm.df,meta.df,markers.df,title.str='Test',cellwidth = 5,cellheight = 8,fontsize_row = 2,fontsize_col = 2,show_rownames = T,show_colnames = T){
  markers.sub.pop.col.list=get.color.list.for.pheatmap(in.vec = meta.df$markers.cluster.groups)
  all.genes=rownames(rpkm.df)
  markers.vec=rownames(markers.df)
  markers.vec=intersect(x=markers.vec,y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  pseudo.value=min(in.rpkm.df[in.rpkm.df!=0])/2
  meta.df=meta.df[colnames(in.rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df=order.by.target.development.stage.meta.df(meta.df = meta.df,order.vec = std.stages.col.abbrev.legend.str)
  #meta.df=order.by.marker.sub.population.meta.df(meta.df = meta.df,order.vec = sub.populations.grp.order.names.vec)
  in.rpkm.df=subset(in.rpkm.df,select = rownames(meta.df))
  development.stages.col=get.col.factor(as.character(meta.df[,'development.stage']))
  #development.stages.col=get.col.factor(as.character(meta.df[,'SP']))
  #development.stages.col=get.col.factor(as.character(meta.df[,'markers.cluster.groups']))
  in.rpkm.mat=as.matrix(log2(in.rpkm.df+pseudo.value))
  col.annot.col.df=subset(meta.df,select = development.stage)
  #col.annot.col.df=subset(meta.df,select = markers.cluster.groups)
  #col.annot.col.list=list(development.stage=abbrev.std.stages.col.vec,markers.cluster.groups=markers.sub.pop.col.list)
  col.annot.col.list=list(development.stage=abbrev.std.stages.col.vec)
  pheatmap(mat = in.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = 4,cellheight = 4,fontsize_row = 2,fontsize_col = 2,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,annotation_legend = F,show_rownames = T,show_colnames = T)
  #pheatmap(mat = in.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = cellwidth,cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,annotation_legend = F,cluster_rows = F,cluster_cols = F,show_rownames = show_rownames,show_colnames = show_colnames,border_color = NA)
  pheatmap(mat = in.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = cellwidth,cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,annotation_legend = F,cluster_rows = F,cluster_cols = F,show_rownames = show_rownames,show_colnames = show_colnames,border_color = NA)
  sub.markers.df=markers.df[rownames(in.rpkm.mat),]
  #rownames(in.rpkm.mat)=as.character(sub.markers.df$Product_description)
  rownames(in.rpkm.mat)=as.character(sub.markers.df$gene.name)
  pheatmap(mat = in.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = cellwidth,cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,legend = T,annotation_legend = F,cluster_rows = F,cluster_cols = F,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,show_rownames = show_rownames,show_colnames = show_colnames,border_color = NA)
  return(expressed.markers.vec)
}
plot.sig.diff.sexual.markers.per.sample=function(rpkm.df,meta.df,
                                                 markers.vec,title.str='Test',cellwidth = 5,
                                                 cellheight = 8,fontsize_row = 2,
                                                 fontsize_col = 2,show_rownames = T,
                                                 show_colnames = T,
                                                 in.col.ramp=colorRampPalette(c("gray",'red'))(100)){
  markers.sub.pop.col.list=get.color.list.for.pheatmap(in.vec = meta.df$markers.cluster.groups)
  #col.pelette=colorRampPalette(c("darkgray","gray","white","yellow","orange",'red'))(5)
  col.pelette=in.col.ramp
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(x=markers.vec,y =all.genes  )
  out.list=list()
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  out.list[['rpkm']]=in.rpkm.df
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  pseudo.value=min(in.rpkm.df[in.rpkm.df!=0])/2
  meta.df=meta.df[colnames(in.rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df=meta.df[order(meta.df$mRFP1.expr),]
  #meta.df=meta.df[order(meta.df$putative.commited),]
  out.list[['meta']]=meta.df
  # meta.df=order.by.target.development.stage.meta.df(meta.df = meta.df,
  #                                                   order.vec = std.stages.col.abbrev.legend.str)
  # 
  #meta.df=order.by.marker.sub.population.meta.df(meta.df = meta.df,order.vec = sub.populations.grp.order.names.vec)
  in.rpkm.df=subset(in.rpkm.df,select = rownames(meta.df))
  #development.stages.col=get.col.factor(as.character(meta.df[,'development.stage']))
  #development.stages.col=get.col.factor(as.character(meta.df[,'SP']))
  #development.stages.col=get.col.factor(as.character(meta.df[,'markers.cluster.groups']))
  in.rpkm.mat=as.matrix(log2(in.rpkm.df+pseudo.value))
  col.annot.col.df=subset(meta.df,select = c(mRFP1.expr,
                                             putative.commited,
                                             development.stage))
  #col.annot.col.df=subset(meta.df,select = markers.cluster.groups)
  #col.annot.col.list=list(development.stage=abbrev.std.stages.col.vec,markers.cluster.groups=markers.sub.pop.col.list)
  #col.annot.col.df=subset(meta.df,select = mRFP1.expr)
  in.rpkm.mat=as.matrix(log2(in.rpkm.df+1))
  col.gap.indices=c()
  col.names=colnames(in.rpkm.mat)
  len.col.names=length(col.names)
  for(m in 1:len.col.names){
    if(m<len.col.names){
      first.sample=col.names[m]
      sec.sample=col.names[m+1]
      first.sample.mrfp.expr=as.character(col.annot.col.df[first.sample,'mRFP1.expr'])
      sec.sample.mrfp.expr=as.character(col.annot.col.df[sec.sample,'mRFP1.expr'])
      if(first.sample.mrfp.expr!=sec.sample.mrfp.expr){
        col.gap.indices=append(col.gap.indices,m,length(col.gap.indices))
      }
    }
  }
  col.annot.col.list=list(development.stage=abbrev.std.stages.col.vec)
  # pheatmap(mat = in.rpkm.mat,main =title.str,color = col.pelette,cellwidth = 4,cellheight = 4,
  #          fontsize_row = 2,fontsize_col = 2,annotation_col = col.annot.col.df,
  #          annotation_colors = col.annot.col.list,legend = T,annotation_legend = F,
  #          show_rownames = T,show_colnames = T,gaps_col = col.gap.indices)
  #pheatmap(mat = in.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = cellwidth,cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,annotation_legend = F,cluster_rows = F,cluster_cols = F,show_rownames = show_rownames,show_colnames = show_colnames,border_color = NA)
  # pheatmap(mat = in.rpkm.mat,main =title.str,color = col.pelette,cellwidth = cellwidth,
  #          cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,
  #          annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,
  #          annotation_legend = F,cluster_rows = F,cluster_cols = F,show_rownames = show_rownames,
  #          show_colnames = show_colnames,border_color = NA,gaps_col = col.gap.indices)
  #sub.markers.df=markers.df[rownames(in.rpkm.mat),]
  #rownames(in.rpkm.mat)=as.character(sub.markers.df$Product_description)
  #rownames(in.rpkm.mat)=as.character(sub.markers.df$gene.name)
  pheatmap(mat = in.rpkm.mat,main =title.str,color = col.pelette,cellwidth = cellwidth,
           cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,
           legend = T,annotation_legend = F,cluster_rows = F,cluster_cols = F,
           annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,
           show_rownames = show_rownames,show_colnames = show_colnames,border_color = NA,
           gaps_col = col.gap.indices)
  in.rpkm.mat=as.matrix(log2(in.rpkm.df+1))
  #col.annot.list=split(meta.df,f = meta.df$mRFP1.expr)
  col.annot.list=split(meta.df,f = meta.df$putative.commited)
  col.annot.names=names(col.annot.list)
  expression.list=list()
  for (gene in rownames(in.rpkm.mat)){
    temp.list=list()
    for (sb in col.annot.names){
      #name=ifelse(sb=='Y','Gametocytes','Asexual')
      name=ifelse(sb=='putative.committed','Gametocytes','Asexual')
      temp.list[[name]]=as.numeric(in.rpkm.mat[gene,rownames(col.annot.list[[sb]])])
      expression.list[[paste(gene,sb,sep = '')]]=as.numeric(in.rpkm.mat[gene,rownames(col.annot.list[[sb]])])
    }
    boxplot(temp.list,pch=19,col=c('green','purple'),frame.plot=F,
            main=gene,outline=T,
            border=c('black','black'),width=c(.2,.2))
  } 
  return(out.list)
}
plot.markers.expression.heatmap.per.sample.per.sub.population=function(rpkm.df,meta.df,markers.df,title.str='Test',cellwidth = 5,cellheight = 8,fontsize_row = 2,fontsize_col = 2,show_rownames = T,show_colnames = T){
  markers.sub.pop.col.list=get.color.list.for.pheatmap(in.vec = 
                                                         as.character(meta.df$markers.cluster.groups))
  col.pelette=colorRampPalette(c("darkgray","gray","white","yellow","orange",'red'))(3)
  all.genes=rownames(rpkm.df)
  markers.vec=rownames(markers.df)
  markers.vec=intersect(x=markers.vec,y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = in.rpkm.df))
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  pseudo.value=min(in.rpkm.df[in.rpkm.df!=0])/2
  meta.df=meta.df[colnames(in.rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  #meta.df=order.by.target.development.stage.meta.df(meta.df = meta.df,order.vec = std.stages.col.abbrev.legend.str)
  meta.df=order.by.marker.sub.population.meta.df(meta.df = meta.df,
                                                 order.vec = sub.populations.grp.order.names.vec)
  temp.mfrp.expr.vec=rownames(subset(meta.df,mRFP1.expr=='Y'))
  samples.order.vec=c(temp.mfrp.expr.vec,setdiff(rownames(meta.df),temp.mfrp.expr.vec))
  #in.rpkm.df=subset(in.rpkm.df,select = rownames(meta.df))
  in.rpkm.df=subset(in.rpkm.df,select = samples.order.vec)
  #development.stages.col=get.col.factor(as.character(meta.df[,'development.stage']))
  #development.stages.col=get.col.factor(as.character(meta.df[,'SP']))
  #development.stages.col=get.col.factor(as.character(meta.df[,'markers.cluster.groups']))
  #col.annot.col.df=subset(meta.df,select = c(markers.cluster.groups,development.stage,mRFP1.expr))
  col.annot.col.df=subset(meta.df,select = mRFP1.expr)
  in.rpkm.mat=as.matrix(log2(in.rpkm.df+pseudo.value))
  col.gap.indices=c()
  col.names=colnames(in.rpkm.mat)
  len.col.names=length(col.names)
  for(m in 1:len.col.names){
    if(m<len.col.names){
      first.sample=col.names[m]
      sec.sample=col.names[m+1]
      first.sample.sp=as.character(col.annot.col.df[first.sample,'markers.cluster.groups'])
      sec.sample.sp=as.character(col.annot.col.df[sec.sample,'markers.cluster.groups'])
      first.sample.mrfp.expr=as.character(col.annot.col.df[first.sample,'mRFP1.expr'])
      sec.sample.mrfp.expr=as.character(col.annot.col.df[sec.sample,'mRFP1.expr'])
      if(first.sample.mrfp.expr!=sec.sample.mrfp.expr){
        col.gap.indices=append(col.gap.indices,m,length(col.gap.indices))
      }
    }
  }
  col.annot.col.list=list(development.stage=abbrev.std.stages.col.vec,
                          markers.cluster.groups=markers.sub.pop.col.list)
  #col.annot.col.list=list(development.stage=abbrev.std.stages.col.vec)
  #pheatmap(mat = in.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = 4,cellheight = 4,fontsize_row = 2,fontsize_col = 2,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,annotation_legend = T,show_rownames = T,show_colnames = T,gaps_col = col.gap.indices)
  #pheatmap(mat = in.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = cellwidth,cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,annotation_legend = F,cluster_rows = F,cluster_cols = F,show_rownames = show_rownames,show_colnames = show_colnames,border_color = NA,gaps_col = col.gap.indices)
  pheatmap(mat = in.rpkm.mat,main =title.str,color = 
             col.pelette,cellwidth = cellwidth,
           cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,
           annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,
           annotation_legend = T,cluster_rows = F,cluster_cols = F,show_rownames = show_rownames,
           show_colnames = show_colnames,border_color = NA,gaps_col = col.gap.indices)
  sub.markers.df=markers.df[rownames(in.rpkm.mat),]
  #rownames(in.rpkm.mat)=as.character(sub.markers.df$Product_description)
  rownames(in.rpkm.mat)=as.character(sub.markers.df$gene.name)
  pheatmap(mat = in.rpkm.mat,main =title.str,color = 
             col.pelette,cellwidth = cellwidth,
           cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,
           legend = T,annotation_legend = F,cluster_rows = F,cluster_cols = F,
           annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,
           show_rownames = show_rownames,show_colnames = show_colnames,border_color = NA,
           gaps_col = col.gap.indices)
  return(expressed.markers.vec)
}
plot.markers.expression.barplot.per.sub.population.per.category=function(rpkm.df,meta.df,markers.df,
                                                                         title.str='Test'){
  markers.list=split(markers.df,f=markers.df$category)
  categories.names=names(markers.list)
  len.categories.names=length(categories.names)
  for(m in 1:len.categories.names){
    temp.categories.name=categories.names[m]
    # ignore.vec=c('Cytoplasmic Translation machinery' , 
    #              'Actin myosin motors','Early ring transcripts')
    # if (temp.categories.name %in% ignore.vec){ 
    #   next
    # }
    temp.markers.vec=rownames(markers.list[[temp.categories.name]])
    title.str=paste(temp.categories.name,' (n=',length(temp.markers.vec),')',sep='')
    plot.markers.expression.barplot.per.sub.population(rpkm.df = rpkm.df,meta.df = meta.df,
                                                       markers.vec =temp.markers.vec,
                                                       title.str =  title.str)
    }
}
plot.markers.expression.boxplot.per.sub.population.per.category=function(rpkm.df,meta.df,markers.df,
                                                                         title.str='Test'){
  markers.list=split(markers.df,f=markers.df$category)
  categories.names=names(markers.list)
  len.categories.names=length(categories.names)
  for(m in 1:len.categories.names){
    temp.categories.name=categories.names[m]
    # ignore.vec=c('Cytoplasmic Translation machinery' , 
    #              'Actin myosin motors','Early ring transcripts')
    # if (temp.categories.name %in% ignore.vec){ 
    #   next
    # }
    temp.markers.vec=rownames(markers.list[[temp.categories.name]])
    title.str=paste(temp.categories.name,' (n=',length(temp.markers.vec),')',sep='')
    plot.specific.genes.expression.boxplot.per.sub.population(rpkm.df = rpkm.df,meta.df = meta.df,
                                                              gene.vec  = temp.markers.vec,
                                                              title.str = temp.categories.name)
  }
}
plot.specific.genes.expression.boxplot.per.sub.population=
  function(rpkm.df,meta.df,gene.vec,title.str='Test boxplot'){
    all.genes=rownames(rpkm.df)
    markers.vec=intersect(gene.vec,y =all.genes)
    in.rpkm.df=rpkm.df[markers.vec,]
    #in.rpkm.df=filter.none.expressed.genes(in.rpkm.df)
    #in.rpkm.df=filter.none.expressed.samples(input.data  = filter.rpkm.less.than.cutoff(df =  in.rpkm.df,rpkm.cutoff = 1,no.samples = 2))
    meta.df=meta.df[colnames(in.rpkm.df),]
    expressed.markers.vec=rownames(in.rpkm.df)
    len.expressed.markers.vec=length(expressed.markers.vec)
    #meta.list=split(meta.df,f = meta.df$markers.cluster.groups)
    #meta.list=split(meta.df,f = meta.df$ordered.SP) 
    meta.list=split(meta.df,f = meta.df$subpopulation)
    development.stages=names(meta.list)
    len.development.stages=length(development.stages)
    expr.list=list()
    for(m in 1:len.development.stages){
      development.stage=development.stages[m]
      temp.meta.df=meta.list[[development.stage]]
      temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
      #temp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data   =temp.rpkm.df))
      if(is.null(dim(temp.rpkm.df))){
        expr.list[[development.stage]]=0
        next
      }
      expressed.marker=rownames(temp.rpkm.df)
      temp.rpkm.df=temp.rpkm.df[expressed.marker,]
      temp.rpkm.vec=as.numeric(melt(temp.rpkm.df)$value)
      median.exp.vec=as.numeric(apply(temp.rpkm.df,1,median))
      mean.exp.vec=as.numeric(apply(temp.rpkm.df,1,mean))
      #expr.list[[development.stage]]=log2(median.exp.vec+1)
      #expr.list[[development.stage]]=log2(mean.exp.vec+1)
      #expr.list[[development.stage]]=mean.exp.vec
      trim.mean.exp.vec=as.numeric(apply(temp.rpkm.df,2,function(gn.expr){
        temp.expr=as.numeric(gn.expr)
        #temp.expr=temp.expr[temp.expr>0]
        #trimmed.mean=mean(robustHD::winsorize(x = temp.expr))
        trimmed.mean=mean(temp.expr,trim=0.1)
        gn.mean=mean(temp.expr)
        return(gn.mean)
      }))
      print(length(trim.mean.exp.vec))
      temp.gn.xpr.vec=as.numeric(melt(temp.rpkm.df)$value)
      #expr.list[[development.stage]]=temp.gn.xpr.vec
      expr.list[[development.stage]]=trim.mean.exp.vec
    }
    total.counts=max(as.numeric(lapply(expr.list,length)))
    expr.df.list=lapply(expr.list,function(lst){
      len.lst=length(lst)
      temp.lst=rep(NA,total.counts-len.lst)
      if(length(temp.lst)>=1){
        lst=c(lst,temp.lst)
        return(lst)
      }
      else{
        return(lst)
      }
    })
    expr.df=as.data.frame(melt(data.frame(expr.df.list)))
    colnames(expr.df)=c('SP','RPKM')
    limit=sp.outlier.limits[title.str]
    if(!is.na(limit)){
      temp.box.stats.list=boxplot(as.numeric(expr.df$RPKM),plot=F)
      #p <- ggplot(expr.df, aes(SP,RPKM))+geom_jitter(na.rm = T)
      #p <- ggplot(expr.df, aes(SP,RPKM,fill=SP))+geom_bar(na.rm = T,stat = 'identity') +theme_classic()
      #p <- ggplot(expr.df, aes(SP,RPKM,fill=SP))+geom_boxplot(na.rm = T,outlier.size = NA,outlier.shape = NA) +theme_classic()
      p <- ggplot(expr.df, aes(SP,RPKM,fill=SP))+geom_boxplot(na.rm = T,outlier.shape=NA)
      #p <- ggplot(expr.df, aes(SP,RPKM,fill=SP))+geom_boxplot(na.rm = T)
      #p <-p +theme_classic()
      p <-p +theme_classic()+scale_y_continuous(limits=c(0,limit))
      p<-p+geom_jitter(width = 0.1)
      p<-p+scale_fill_manual(breaks = names(ordered.sp.cols),values =as.character(ordered.sp.cols))
      #p<-p+scale_x_continuous(limit=(0,2000))
      p<-p + ggtitle(title.str)
      #p<-p+geom_point(position = position_jitter(width = 0.1),size=1.5)
      p<-p+scale_color_manual(breaks = names(ordered.sp.cols),values =as.character(ordered.sp.cols) )
      print(p)
    }
    else{
      temp.box.stats.list=boxplot(as.numeric(expr.df$RPKM),plot=F)
      #p <- ggplot(expr.df, aes(SP,RPKM))+geom_jitter(na.rm = T)
      #p <- ggplot(expr.df, aes(SP,RPKM,fill=SP))+geom_bar(na.rm = T,stat = 'identity') +theme_classic()
      p <- ggplot(expr.df, aes(SP,RPKM,fill=SP))+geom_boxplot(na.rm = T,outlier.size = NA,outlier.shape = NA) +theme_classic()
      #p <- ggplot(expr.df, aes(SP,RPKM,fill=SP))+geom_boxplot(na.rm = T,outlier.shape=NA) +theme_classic()+scale_y_continuous(limits=c(0,min(temp.box.stats.list$out)))
      #p <- ggplot(expr.df, aes(SP,RPKM,fill=SP))+geom_boxplot(na.rm = T,outlier.shape=NA) +theme_classic()
      #p<-p+geom_jitter(width = 0.1)
      p<-p+scale_fill_manual(breaks = names(ordered.sp.cols),values =as.character(ordered.sp.cols))
      #p<-p+scale_x_continuous(limit=(0,2000))
      p<-p + ggtitle(title.str)
      #p<-p+geom_point(position = position_jitter(width = 0.1),size=1.5)
      p<-p+scale_color_manual(breaks = names(ordered.sp.cols),values =as.character(ordered.sp.cols) )
      #print(p)
    }
  }
plot.markers.expression.violin.per.sub.population.per.category=function(rpkm.df,meta.df,markers.df,
                                                                         title.str='Test'){
  markers.list=split(markers.df,f=markers.df$category)
  categories.names=names(markers.list)
  len.categories.names=length(categories.names)
  for(m in 1:len.categories.names){
    temp.categories.name=categories.names[m]
    # ignore.vec=c('Cytoplasmic Translation machinery' , 
    #              'Actin myosin motors','Early ring transcripts')
    # if (temp.categories.name %in% ignore.vec){ 
    #   next
    # }
    temp.markers.vec=rownames(markers.list[[temp.categories.name]])
    title.str=paste(temp.categories.name,' (n=',length(temp.markers.vec),')',sep='')
    plot.specific.genes.expression.violin.per.sub.population(rpkm.df = rpkm.df,meta.df = meta.df,
                                                              gene.vec  = temp.markers.vec,
                                                              title.str = title.str)
  }
}
plot.specific.genes.expression.violin.per.sub.population=
  function(rpkm.df,meta.df,gene.vec,title.str='Test boxplot'){
    all.genes=rownames(rpkm.df)
    markers.vec=intersect(gene.vec,y =all.genes)
    in.rpkm.df=rpkm.df[markers.vec,]
    #in.rpkm.df=filter.none.expressed.genes(in.rpkm.df)
    #in.rpkm.df=filter.none.expressed.samples(input.data  = filter.rpkm.less.than.cutoff(df =  in.rpkm.df,rpkm.cutoff = 1,no.samples = 2))
    meta.df=meta.df[colnames(in.rpkm.df),]
    expressed.markers.vec=rownames(in.rpkm.df)
    len.expressed.markers.vec=length(expressed.markers.vec)
    #meta.list=split(meta.df,f = meta.df$markers.cluster.groups)
    meta.list=split(meta.df,f = meta.df$ordered.SP) 
    development.stages=names(meta.list)
    len.development.stages=length(development.stages)
    expr.list=list()
    for(m in 1:len.development.stages){
      development.stage=development.stages[m]
      temp.meta.df=meta.list[[development.stage]]
      temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
      #temp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data   =temp.rpkm.df))
      if(is.null(dim(temp.rpkm.df))){
        expr.list[[development.stage]]=0
        next
      }
      expressed.marker=rownames(temp.rpkm.df)
      temp.rpkm.df=temp.rpkm.df[expressed.marker,]
      temp.rpkm.vec=as.numeric(melt(temp.rpkm.df)$value)
      median.exp.vec=as.numeric(apply(temp.rpkm.df,1,median))
      mean.exp.vec=as.numeric(apply(temp.rpkm.df,1,mean))
      #expr.list[[development.stage]]=log2(median.exp.vec+1)
      #expr.list[[development.stage]]=log2(mean.exp.vec+1)
      temp.rpkm.vec[is.na(temp.rpkm.vec)]=0
      #expr.list[[development.stage]]=log2(temp.rpkm.vec+1)
      expr.list[[development.stage]]=log2(temp.rpkm.vec[temp.rpkm.vec>0]+1)
    }
    total.counts=max(as.numeric(lapply(expr.list,length)))
    expr.df.list=lapply(expr.list,function(lst){
      len.lst=length(lst)
      temp.lst=rep(NA,total.counts-len.lst)
      if(length(temp.lst)>=1){
        lst=c(lst,temp.lst)
        return(lst)
      }
      else{
        return(lst)
      }
    })
    expr.df=as.data.frame(melt(data.frame(expr.df.list)))
    colnames(expr.df)=c('SP','RPKM')
    #p <- ggplot(expr.df, aes(SP,RPKM))+geom_jitter(na.rm = T)
    #p <- ggplot(expr.df, aes(SP,RPKM,fill=SP))+geom_bar(na.rm = T,stat = 'identity') +theme_classic()
    #p <- ggplot(expr.df, aes(SP,RPKM,fill=SP))+geom_boxplot(na.rm = T) +theme_classic()
    p <- ggplot(expr.df, aes(SP,RPKM,fill=SP))+geom_violin(na.rm = T,trim=T) +theme_classic()
    #p<-p+geom_jitter(width = 0.1)
    p<-p+scale_fill_manual(breaks = names(ordered.sp.cols),values =as.character(ordered.sp.cols))
    p<-p + ggtitle(title.str)
    #p<-p+geom_point(position = position_jitter(width = 0.1),size=1.5)
    p<-p+scale_color_manual(breaks = names(ordered.sp.cols),values =as.character(ordered.sp.cols) )
    print(p)
  }
plot.markers.expression.boxplot.per.sub.population=function(rpkm.df,meta.df,markers.vec,title.str='Test'){
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = rpkm.df))
  all.genes=rownames(in.rpkm.df)
  markers.vec=intersect(x=markers.vec,y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  pseudo.value=1
  meta.df=meta.df[colnames(in.rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  #meta.df=order.by.marker.sub.population.meta.df(meta.df = meta.df,order.vec = sub.populations.grp.order.names.vec)
  #meta.list=split(meta.df,f = meta.df$markers.cluster.groups)
  meta.list=split(meta.df,f = meta.df$ordered.SP)
  sub.populations.names.vec=names(meta.list)
  len.sub.populations.names.vec=length(sub.populations.names.vec)
  marker.expr.list=list()
  values.vec=c()
  mean.vec=c()
  se.vec=c()
  col.vec=c()
  for(n in 1:len.sub.populations.names.vec){
    temp.sub.populations.grp.order.names.vec=sub.populations.names.vec[n]
    temp.meta.df=meta.list[[temp.sub.populations.grp.order.names.vec]]
    temp.in.rpkm.df=subset(in.rpkm.df,select = rownames(temp.meta.df))
    temp.markers.rpkm.df=temp.in.rpkm.df[markers.vec,]
    expr.vec=as.numeric(melt(temp.markers.rpkm.df)$value)
    log.expr.vec=log2(expr.vec+1)
    marker.expr.list[[temp.sub.populations.grp.order.names.vec]]=expr.vec+1
    col.vec=append(x = col.vec,values = ordered.sp.cols[temp.sub.populations.grp.order.names.vec],after = length(col.vec))
  }
  boxplot2(marker.expr.list,main=title.str,log='y',las=2,pch=19,col=col.vec)
}
plot.specific.gene.expression.boxplot.per.sub.population=function(rpkm.df,meta.df,gene.id,title.str='Test'){
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = rpkm.df))
  all.genes=rownames(in.rpkm.df)
  markers.vec=intersect(x=gene.id,y =all.genes)
  if(length(markers.vec)==0){print('Gene not found')
  }
  else{
    in.rpkm.df=rpkm.df[markers.vec,]
    #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
    expressed.markers.vec=rownames(in.rpkm.df)
    len.expressed.markers.vec=length(expressed.markers.vec)
    pseudo.value=1
    meta.df=meta.df[colnames(in.rpkm.df),]
    meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
    meta.list=split(meta.df,f = meta.df$ordered.SP)
    sub.populations.names.vec=names(meta.list)
    len.sub.populations.names.vec=length(sub.populations.names.vec)
    marker.expr.list=list()
    values.vec=c()
    mean.vec=c()
    se.vec=c()
    col.vec=c()
    for(n in 1:len.sub.populations.names.vec){
      temp.sub.populations.grp.order.names.vec=sub.populations.names.vec[n]
      temp.meta.df=meta.list[[temp.sub.populations.grp.order.names.vec]]
      temp.in.rpkm.df=subset(in.rpkm.df,select = rownames(temp.meta.df))
      temp.markers.rpkm.df=temp.in.rpkm.df[markers.vec,]
      expr.vec=as.numeric(melt(temp.markers.rpkm.df)$value)
      log.expr.vec=log2(expr.vec+1)
      marker.expr.list[[temp.sub.populations.grp.order.names.vec]]=log.expr.vec
      col.vec=append(x = col.vec,values = ordered.sp.cols[temp.sub.populations.grp.order.names.vec],after = length(col.vec))
    }
    boxplot(marker.expr.list,main=gene.id,las=2,pch=19,col=col.vec,frame.plot=F)
  }
}
plot.markers.median.expression.barplot.per.sub.population.per.category=function(rpkm.df,meta.df,markers.df,
                                                                         title.str='Test'){
  markers.list=split(markers.df,f=markers.df$category)
  categories.names=names(markers.list)
  len.categories.names=length(categories.names)
  for(m in 1:len.categories.names){
    temp.categories.name=categories.names[m]
    # ignore.vec=c('Cytoplasmic Translation machinery' , 
    #              'Actin myosin motors','Early ring transcripts')
    # if (temp.categories.name %in% ignore.vec){
    #   next
    # }
    temp.markers.vec=rownames(markers.list[[temp.categories.name]])
    title.str=paste(temp.categories.name,' (n=',length(temp.markers.vec),')',sep='')
    plot.markers.median.expression.barplot.per.sub.population(rpkm.df = rpkm.df,meta.df = meta.df,
                                                       markers.vec =temp.markers.vec,
                                                       title.str =  title.str)
  }
}
plot.markers.median.expression.barplot.per.sub.population=function(rpkm.df,meta.df,markers.vec,title.str='Test'){
  markers.sub.pop.col.list=get.color.list.for.pheatmap(in.vec = as.character(meta.df$markers.cluster.groups))
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = rpkm.df))
  all.genes=rownames(in.rpkm.df)
  markers.vec=intersect(x=markers.vec,y =all.genes)
  in.rpkm.df=rpkm.df[markers.vec,]
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  pseudo.value=1
  meta.df=meta.df[colnames(in.rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  #meta.df=order.by.marker.sub.population.meta.df(meta.df = meta.df,order.vec = sub.populations.grp.order.names.vec)
  meta.list=split(meta.df,f = meta.df$markers.cluster.groups)
  sub.populations.names.vec=names(meta.list)
  len.sub.populations.names.vec=length(sub.populations.names.vec)
  marker.expr.list=list()
  values.vec=c()
  mean.vec=c()
  se.vec=c()
  median.vec=c()
  for(n in 1:len.sub.populations.names.vec){
    temp.sub.populations.grp.order.names.vec=sub.populations.names.vec[n]
    temp.meta.df=meta.list[[temp.sub.populations.grp.order.names.vec]]
    temp.in.rpkm.df=subset(in.rpkm.df,select = rownames(temp.meta.df))
    temp.markers.rpkm.df=temp.in.rpkm.df[markers.vec,]
    expr.vec=as.numeric(melt(temp.markers.rpkm.df)$value)
    log.expr.vec=log2(expr.vec+1)
    #log.expr.vec=log2(expr.vec+pseudo.value)
    expr.stats=describe(expr.vec)
    print(expr.stats)
    values.vec=append(x = values.vec,values =n ,after = length(values.vec))
    mean.vec=append(x = mean.vec,values = expr.stats$mean,after = length(mean.vec))
    se.vec=append(x = se.vec,values =expr.stats$se ,after = length(se.vec))
    median.vec=append(x=median.vec,values =expr.stats$median ,after = length(median.vec))
    marker.expr.list[[temp.sub.populations.grp.order.names.vec]]=temp.markers.rpkm.df
  }
  names(median.vec)=sub.populations.names.vec
  print(sub.populations.names.vec)
  barplot2(height = median.vec,col = markers.sub.pop.col.list[sub.populations.names.vec],
           main=title.str)
}
plot.markers.expression.barplot.per.sub.population=function(rpkm.df,meta.df,markers.vec,title.str='Test'){
  markers.sub.pop.col.list=get.color.list.for.pheatmap(in.vec = as.character(meta.df$markers.cluster.groups))
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = rpkm.df))
  all.genes=rownames(in.rpkm.df)
  markers.vec=intersect(x=markers.vec,y =all.genes)
  in.rpkm.df=rpkm.df[markers.vec,]
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,
  #rpkm.cutoff = 1,no.samples = 1))
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  pseudo.value=1
  meta.df=meta.df[colnames(in.rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  #meta.df=order.by.marker.sub.population.meta.df(meta.df = meta.df,order.vec = sub.populations.grp.order.names.vec)
  meta.list=split(meta.df,f = meta.df$markers.cluster.groups)
  sub.populations.names.vec=names(meta.list)
  len.sub.populations.names.vec=length(sub.populations.names.vec)
  marker.expr.list=list()
  values.vec=c()
  mean.vec=c()
  se.vec=c()
  for(n in 1:len.sub.populations.names.vec){
    temp.sub.populations.grp.order.names.vec=sub.populations.names.vec[n]
    temp.meta.df=meta.list[[temp.sub.populations.grp.order.names.vec]]
    temp.in.rpkm.df=subset(in.rpkm.df,select = rownames(temp.meta.df))
    temp.markers.rpkm.df=temp.in.rpkm.df[markers.vec,]
    expr.vec=as.numeric(melt(temp.markers.rpkm.df)$value)
    log.expr.vec=log2(expr.vec+1)
    #log.expr.vec=log2(expr.vec+pseudo.value)
    expr.stats=describe(expr.vec)
    print(expr.stats)
    #order.num=as.numeric(sub.populations.grp.order.names.vec[temp.sub.populations.grp.order.names.vec])
    values.vec=append(x = values.vec,values =n ,after = length(values.vec))
    mean.vec=append(x = mean.vec,values = expr.stats$mean,after = length(mean.vec))
    se.vec=append(x = se.vec,values =expr.stats$se ,after = length(se.vec))
    marker.expr.list[[temp.sub.populations.grp.order.names.vec]]=temp.markers.rpkm.df
  }
  marker.stats.df=data.frame(values=values.vec,mean=mean.vec,se=se.vec)
  ylim = c(0,max(marker.stats.df$mean+marker.stats.df$se)+100)
  rownames(marker.stats.df)=sub.populations.names.vec
  col.vec=markers.sub.pop.col.list[rownames(marker.stats.df)]
  marker.stats.df$col.str=col.vec
  marker.stats.df$idc.order=sub.populations.grp.order.names.vec[rownames(marker.stats.df)]
  marker.stats.df=marker.stats.df[order(marker.stats.df$idc.order),]
  #ordered.sub.populations.vec=paste('SP.',as.character(marker.stats.df$values),sep = '')
  #rownames(marker.stats.df)=ordered.sub.populations.vec
  #marker.stats.df=marker.stats.df[sort(rownames(marker.stats.df)),]
  ordered.col.vec=markers.sub.pop.col.list[rownames(marker.stats.df)]
  editted.error.bars(stats=marker.stats.df,bars=T,main=title.str,ylab = 'rpkm',
                     xlab='',eyes = F,col.str=ordered.col.vec,width.size.vec = 2)
}
plot.markers.expression.barplot.per.timepoint.per.category=function(rpkm.df,meta.df,markers.df,title.str='Test'){
  markers.list=split(markers.df,f=markers.df$category)
  categories.names=names(markers.list)
  len.categories.names=length(categories.names)
  for(m in 1:len.categories.names){
    temp.categories.name=categories.names[m]
    temp.markers.vec=rownames(markers.list[[temp.categories.name]])
    title.str=paste(temp.categories.name,' (n=',length(temp.markers.vec),')',sep='')
    if(is.null((markers.peaking.stages.list[[temp.categories.name]]))){
      next
    }
    peaking.stage.vec=as.character(unlist(lapply(X = strsplit(x = markers.peaking.stages.list[[temp.categories.name]],split = ','),FUN = discard.lead.and.tail)))
    #peaking.stage.vec=as.character(unlist(lapply(peaking.stage.vec,discard.lead.and.tail)))
    plot.markers.expression.barplot.per.timepoint(rpkm.df = rpkm.df,meta.df = meta.df,markers.vec =temp.markers.vec,title.str =  title.str,peaking.stage=peaking.stage.vec)
  }
}
plot.markers.expression.barplot.per.timepoint=function(rpkm.df,meta.df,markers.vec,title.str='Test',peaking.stage=NULL){
  #markers.sub.pop.col.list=get.color.list.for.pheatmap(in.vec = as.character(meta.df$markers.cluster.groups))
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = rpkm.df))
  all.genes=rownames(in.rpkm.df)
  markers.vec=intersect(x=markers.vec,y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  pseudo.value=min(in.rpkm.df[in.rpkm.df!=0])/2
  meta.df=meta.df[colnames(in.rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  ordered.stages=c('R','LR','ET','T','ES','S')
  meta.df=order.by.target.development.stage.meta.df(meta.df = meta.df,order.vec = names(abbrev.std.stages.col.vec))
  #meta.df=order.by.marker.sub.population.meta.df(meta.df = meta.df,order.vec = sub.populations.grp.order.names.vec)
  meta.list=split(meta.df,f = meta.df$development.stage)
  sub.populations.names.vec=intersect(ordered.stages,names(meta.list))
  len.sub.populations.names.vec=length(sub.populations.names.vec)
  marker.expr.list=list()
  values.vec=c()
  mean.vec=c()
  se.vec=c()
  for(n in 1:len.sub.populations.names.vec){
    temp.sub.populations.grp.order.names.vec=sub.populations.names.vec[n]
    temp.meta.df=meta.list[[temp.sub.populations.grp.order.names.vec]]
    temp.in.rpkm.df=subset(in.rpkm.df,select = rownames(temp.meta.df))
    temp.markers.rpkm.df=temp.in.rpkm.df[markers.vec,]
    expr.vec=as.numeric(melt(temp.markers.rpkm.df)$value)
    log.expr.vec=log2(expr.vec+pseudo.value)
    expr.stats=describe(expr.vec)
    #order.num=as.numeric(sub.populations.grp.order.names.vec[temp.sub.populations.grp.order.names.vec])
    values.vec=append(x = values.vec,values =n ,after = length(values.vec))
    mean.vec=append(x = mean.vec,values = expr.stats$mean,after = length(mean.vec))
    se.vec=append(x = se.vec,values =expr.stats$se ,after = length(se.vec))
    marker.expr.list[[temp.sub.populations.grp.order.names.vec]]=temp.markers.rpkm.df
  }
  marker.stats.df=data.frame(values=values.vec,mean=mean.vec,se=se.vec)
  ylim = c(0,max(marker.stats.df$mean+marker.stats.df$se)+100)
  rownames(marker.stats.df)=sub.populations.names.vec
  col.vec=c()
  if(is.null(peaking.stage)){
    col.vec=abbrev.std.stages.col.vec[rownames(marker.stats.df)]
  }
  else{
    col.vec=ifelse(rownames(marker.stats.df) %in% peaking.stage ,'blue','gray')
  }
  marker.stats.df$col.str=col.vec
  #marker.stats.df$idc.order=sub.populations.grp.order.names.vec[rownames(marker.stats.df)]
  #marker.stats.df=marker.stats.df[order(marker.stats.df$idc.order),]
  #ordered.sub.populations.vec=paste('SP.',as.character(marker.stats.df$values),sep = '')
  #rownames(marker.stats.df)=ordered.sub.populations.vec
  #marker.stats.df=marker.stats.df[sort(rownames(marker.stats.df)),]
  #ordered.col.vec=markers.sub.pop.col.list[rownames(marker.stats.df)]
  editted.error.bars(stats=marker.stats.df,bars=T,main=title.str,ylab = 'rpkm',xlab='',eyes = F,col.str=col.vec,width.size.vec = rep(0.2,times = 6))
}
editted.error.bars=function (x, stats = NULL, ylab = "Dependent Variable", xlab = "Independent Variable", 
          main = NULL, eyes = TRUE, ylim = NULL, xlim = NULL, alpha = 0.05, 
          sd = FALSE, labels = NULL, pos = NULL, arrow.len = 0.05, 
          arrow.col = "black", add = FALSE, bars = FALSE, within = FALSE, 
          col.str = c("blue"),width.size.vec, ...) 
{
  SCALE = 0.5
  if (is.null(stats)) {
    x.stats <- describe(x)
    if (within) {
      x.smc <- smc(x, covar = TRUE)
      x.stats$se <- sqrt((x.stats$sd^2 - x.smc)/x.stats$n)
    }
    if (is.null(dim(x))) {
      z <- 1
    }
    else {
      z <- dim(x)[2]
    }
    names <- colnames(x)
  }
  else {
    x.stats <- stats
    z <- dim(x.stats)[1]
    names <- rownames(stats)
  }
  min.x <- min(x.stats$mean, na.rm = TRUE)
  max.x <- max(x.stats$mean, na.rm = TRUE)
  max.se <- max(x.stats$se, na.rm = TRUE)
{
    if (!sd) {
      if (is.null(stats)) {
        ci <- qt(1 - alpha/2, x.stats$n - 1)
      }
      else {
        ci <- rep(1, z)
      }
    }
    else {
      ci <- sqrt(x.stats$n)
      max.se <- max(ci * x.stats$se, na.rm = TRUE)
    }
  }
if (is.null(main)) {
  if (!sd) {
    main = paste((1 - alpha) * 100, "% confidence limits", 
                 sep = "")
  }
  else {
    main = paste("Means and standard deviations")
  }
}
if (is.null(ylim)) {
  if (is.na(max.x) | is.na(max.se) | is.na(min.x) | is.infinite(max.x) | 
        is.infinite(min.x) | is.infinite(max.se)) {
    ylim = c(0, 1)
  }
  else {
    if (bars) {
      ylim = c(min(0, min.x - 2 * max.se), max.x + 
                 2 * max.se)
    }
    else {
      ylim = c(min.x - 2 * max.se, max.x + 2 * max.se)
    }
  }
}
if (bars) {
  mp = barplot(x.stats$mean, ylim = ylim, xlab = xlab, 
               ylab = ylab, main = main,col=col.str,space=0.3,
               width=width.size.vec, ...)
  axis(1, mp[1:z], names)
  axis(2)
  box()
}
else {
  if (!add) {
    if (missing(xlim)) 
      xlim <- c(0.5, z + 0.5)
    if (is.null(x.stats$values)) {
      plot(x.stats$mean, ylim = ylim, xlab = xlab, 
           ylab = ylab, xlim = xlim, axes = FALSE, main = main, 
           ...)
      axis(1, 1:z, names, ...)
      axis(2)
      box()
    }
    else {
      plot(x.stats$values, x.stats$mean, ylim = ylim, 
           xlab = xlab, ylab = ylab, main = main, ...)
    }
  }
  else {
    points(x.stats$mean, ...)
  }
}
if (!is.null(labels)) {
  lab <- labels
}
else {
  lab <- paste("V", 1:z, sep = "")
}
if (length(pos) == 0) {
  locate <- rep(1, z)
}
else {
  locate <- pos
}
if (length(labels) == 0) 
  lab <- rep("", z)
else lab <- labels
s <- c(1:z)
if (bars) {
  arrows(mp[s], x.stats$mean[s] - ci[s] * x.stats$se[s], 
         mp[s], x.stats$mean[s] + ci[s] * x.stats$se[s], length = arrow.len, 
         angle = 90, code = 3, col = par("fg"), lty = NULL, 
         lwd = par("lwd"), xpd = NULL)
}
else {
  if (is.null(x.stats$values)) {
    arrows(s[s], x.stats$mean[s] - ci[s] * x.stats$se[s], 
           s[s], x.stats$mean[s] + ci[s] * x.stats$se[s], 
           length = arrow.len, angle = 90, code = 3, col = arrow.col)
  }
  else {
    arrows(x.stats$values, x.stats$mean[s] - ci[s] * 
             x.stats$se[s], x.stats$values, x.stats$mean[s] + 
             ci[s] * x.stats$se[s], length = arrow.len, angle = 90, 
           code = 3, col = arrow.col)
  }
  if (eyes) {
    if (length(col) == 1) 
      col <- rep(col, z)
    ln <- seq(-3, 3, 0.1)
    rev <- (length(ln):1)
    for (s in 1:z) {
      if (!is.null(x.stats$n[s])) {
        catseyes(x = s, y = x.stats$mean[s], se = x.stats$se[s], 
                 n = x.stats$n[s], alpha = alpha, density = -10, 
                 col = col[s])
      }
    }
  }
}
}
discard.lead.and.tail=function(in.str){
  out.str=gsub(pattern = '^\\s+',replacement = '',x = gsub(pattern = '\\s+$',replacement = '',in.str))
  return(out.str)
}
plot.markers.expression.heatmap.per.sample.in.bulk.controls=function(rpkm.df,meta.df,title.str='Test',cellwidth = 5,cellheight = 8,fontsize_row = 2,fontsize_col = 2,show_rownames = T,show_colnames = T,markers.df ){
  sub.group.col.list=get.color.rainbow.list.for.pheatmap(in.vec = meta.df$sub.group)
  all.genes=rownames(rpkm.df)
  in.rpkm.df=rpkm.df
  in.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(in.rpkm.df))
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  pseudo.value=min(in.rpkm.df[in.rpkm.df!=0])/2
  markers.vec=rownames(markers.df)
  markers.vec=intersect(x=markers.vec,y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  meta.df=meta.df[colnames(in.rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  #meta.df=order.df.by.specific.col(meta.df = meta.df,order.vec = std.stages.col.abbrev.legend.str)
  sub.populations.grp.names.vec=names(sub.group.col.list)
  sub.populations.grp.order.names.vec=sort(unique(sub.populations.grp.names.vec))
  meta.df=order.df.by.specific.col(meta.df = meta.df,order.vec = sub.populations.grp.order.names.vec)
  meta.df$sub.group=factor(x = meta.df$sub.group ,levels = as.factor(sub.populations.grp.order.names.vec))
  in.rpkm.df=subset(in.rpkm.df,select = rownames(meta.df))
  #development.stages.col=get.col.factor(as.character(meta.df[,'development.stage']))
  #development.stages.col=get.col.factor(as.character(meta.df[,'SP']))
  #development.stages.col=get.col.factor(as.character(meta.df[,'markers.cluster.groups']))
  in.rpkm.mat=as.matrix(log2(in.rpkm.df+pseudo.value))
  col.annot.col.df=subset(meta.df,select = c(development.stage,sub.group))
  col.annot.col.list=list(development.stage=abbrev.std.stages.col.vec,sub.group=sub.group.col.list)
  pheatmap(mat = in.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = 4,cellheight = 4,fontsize_row = 2,fontsize_col = 2,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,annotation_legend = F,show_rownames = T,show_colnames = T)
  #pheatmap(mat = in.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = cellwidth,cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,annotation_legend = F,cluster_rows = F,cluster_cols = F,show_rownames = show_rownames,show_colnames = show_colnames,border_color = NA)
  pheatmap(mat = in.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = cellwidth,cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,legend = T,annotation_legend = T,cluster_rows = F,cluster_cols = F,show_rownames = show_rownames,show_colnames = show_colnames,border_color = NA)
  sub.markers.df=markers.df[rownames(in.rpkm.mat),]
  show(dim(in.rpkm.mat))
  #rownames(in.rpkm.mat)=as.character(sub.markers.df$Product_description)
  rownames(in.rpkm.mat)=as.character(sub.markers.df$gene.name)
  pheatmap(mat = in.rpkm.mat,main =title.str,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = cellwidth,cellheight = cellheight,fontsize_row = fontsize_row,fontsize_col = fontsize_col,legend = T,annotation_legend = F,cluster_rows = F,cluster_cols = F,annotation_col = col.annot.col.df,annotation_colors = col.annot.col.list,show_rownames = show_rownames,show_colnames = show_colnames,border_color = NA)
  return(markers.vec)
}
order.by.target.development.stage.meta.df=function(meta.df,order.vec){
  meta.list=split(meta.df,f = meta.df$development.stage)
  len.order.vec=length(order.vec)
  ordered.cols=c()
  for(n in 1:len.order.vec){
    temp.name=order.vec[n]
    temp.meta.df=meta.list[[temp.name]]
    ordered.cols=append(ordered.cols,values =rownames(temp.meta.df) ,after = length(ordered.cols))
  }
  ordered.meta.df=meta.df[ordered.cols,]
  return(ordered.meta.df)
}
order.by.marker.sub.population.meta.df=function(meta.df,order.vec){
  meta.list=split(meta.df,f = meta.df$markers.cluster.groups)
  order.vec=names(order.vec)
  len.order.vec=length(order.vec)
  ordered.cols=c()
  for(n in 1:len.order.vec){
    temp.name=order.vec[n]
    temp.meta.df=meta.list[[temp.name]]
    ordered.cols=append(ordered.cols,values =rownames(temp.meta.df) ,after = length(ordered.cols))
  }
  ordered.meta.df=meta.df[ordered.cols,]
  return(ordered.meta.df)
}
order.df.by.specific.col=function(meta.df,order.vec){
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.list=split(meta.df,f = meta.df$development.stage)
  order.vec=intersect(order.vec,names(meta.list))
  len.order.vec=length(order.vec)
  ordered.cols=c()
  for(n in 1:len.order.vec){
    temp.name=order.vec[n]
    temp.meta.df=meta.list[[temp.name]]
    ordered.cols=append(ordered.cols,values =rownames(temp.meta.df) ,after = length(ordered.cols))
  }
  ordered.meta.df=meta.df[ordered.cols,]
  return(ordered.meta.df)
}
order.any.df.by.specific.col=function(meta.list,order.vec){
  order.vec=intersect(order.vec,names(meta.list))
  len.order.vec=length(order.vec)
  ordered.list=list()
  row.names.vec=c()
  for(n in 1:len.order.vec){
    temp.name=order.vec[n]
    temp.meta.df=meta.list[[temp.name]]
    ordered.list[[temp.name]]=temp.meta.df
    row.names.vec=append(x = row.names.vec,values =rownames(temp.meta.df) ,after = length(row.names.vec))
  }
  out.df=do.call('rbind',ordered.list )
  rownames(out.df)=row.names.vec
  #show(out.df)
  return(out.df) 
}
plot.specific.genes.expression.boxplot=function(rpkm.df,meta.df,gene.vec,parental=F,col.str='blue'){
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(gene.vec,y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.genes(in.rpkm.df)
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  meta.df=meta.df[colnames(in.rpkm.df),]
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=list()
  intersect.development.stages=c()
  if(parental){
    meta.df=create.parental.stage(meta.df =  meta.df)
    meta.list=split(meta.df,f = meta.df$parent.development.stage)
    development.stages=names(meta.list)
    ordered.development.stages=c("R" , "T" ,  "S")
    intersect.development.stages = intersect(ordered.development.stages,development.stages)
    len.intersect.development.stages=length(intersect.development.stages)
    for(n in 1:len.expressed.markers.vec){
      expressed.marker=expressed.markers.vec[n]
      temp.list=list()
      for(m in 1:len.intersect.development.stages){
        intersect.development.stage=intersect.development.stages[m]
        temp.meta.df=meta.list[[intersect.development.stage]]
        temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
        temp.rpkm=as.numeric(as.character(temp.rpkm.df[expressed.marker,]))
        temp.rpkm=ifelse(temp.rpkm==0,0.001,temp.rpkm)
        temp.list[[intersect.development.stage]]=temp.rpkm
      }
      boxplot(temp.list,las=2,main=expressed.marker,log = 'y',pars = list(boxwex = 0.5, staplewex = 0.3, outwex = 0.3),pch=19,col=col.str)
    }
  }
  else{
    meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
    meta.list=split(meta.df,f = meta.df$development.stage)
    development.stages=names(meta.list)
    ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
    intersect.development.stages = intersect(ordered.development.stages,development.stages)
    len.intersect.development.stages=length(intersect.development.stages)
    for(n in 1:len.expressed.markers.vec){
      expressed.marker=expressed.markers.vec[n]
      temp.list=list()
      for(m in 1:len.intersect.development.stages){
        intersect.development.stage=intersect.development.stages[m]
        temp.meta.df=meta.list[[intersect.development.stage]]
        temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
        temp.rpkm=as.numeric(as.character(temp.rpkm.df[expressed.marker,]))
        temp.rpkm=ifelse(temp.rpkm==0,0.001,temp.rpkm)
        temp.list[[intersect.development.stage]]=temp.rpkm
      }
      boxplot(temp.list,las=1,main=expressed.marker,log = 'y',col = col.str,pch=19)
    }
  }
}
filter.reid.data=function(input.data,lib.size=25000,gene.counts=1000,reads.cut.off=1,no.samples=3){
  filt.input.data=filter.none.expressed.samples(input.data = filter.none.expressed.genes(input.data = input.data))
  filt.samples=colnames(filt.input.data)[colSums(filt.input.data)>=lib.size]
  filt.gns=rownames(filt.input.data)[rowSums(filt.input.data>reads.cut.off)>no.samples]
  out.df=input.data[filt.gns,filt.samples]
  return(out.df[,apply(out.df,2,Matrix::nnzero)>gene.counts])
}
scran.calc <- function(countData, spike=F) {
  if(spike==FALSE) {
    message(paste0("Using computeSumFactors, i.e. deconvolution over all cells!"))
    cnts = countData
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(cnts)))
    if(ncol(countData)<=14) {
      sizes <- c(round(seq(from=2, to=trunc(ncol(countData)/2), by = 1)))
      sf <- scran::computeSumFactors(sce,
                                     sizes=sizes,
                                     positive=FALSE, sf.out=TRUE)
    }
    if(ncol(countData)>14 & ncol(countData)<=50) {
      sizes <- c(round(seq(from=2, to=trunc(ncol(countData)/2), length.out=6)))
      sf <- scran::computeSumFactors(sce,
                                     sizes=sizes,
                                     positive=FALSE, sf.out=TRUE)
    }
    if(ncol(countData)>50 & ncol(countData)<=1000) {
      sizes <- c(round(seq(from=10, to=trunc(ncol(countData)/2), length.out=6)))
      sf <- scran::computeSumFactors(sce,sizes=sizes,
                                     positive=FALSE, sf.out=TRUE)
    }
    if(trunc(ncol(countData))>1000) {
      sizes <- c(round(seq(from=20, to=trunc(ncol(countData)/2), length.out=6)))
      qclust <- scran::quickCluster(sce, min.size = 30)
      sf <- scran::computeSumFactors(sce, clusters=qclust,positive=FALSE, sf.out=TRUE)
    }
  }
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  return(res)
}
plot.reid.data.expression.boxplot.per.subpop=function(f.rpkm.df,f.meta.df,s.rpkm.df,s.meta.df,f.col,s.col,gene.vec){
  f.grps=unique(f.meta.df[f.col])
  s.grps=unique(s.meta.df[s.col])
  for (gn in gene.vec){
    for (f.grp in f.grps){
      temp.f.samples.vec=rownames(f.meta.df)[f.meta.df[f.col]==f.grp]
      temp.f.samples.vec=temp.f.samples.vec[!is.na(temp.f.samples.vec)]
      temp.f.samples.vec=intersect(colnames(f.rpkm.df),temp.f.samples.vec)
      temp.f.expr=subset(f.rpkm.df,select=temp.f.samples.vec)
      temp.f.dim=dim(temp.f.expr)
      for (s.grp in s.grps){
        temp.s.samples.vec=rownames(s.meta.df)[s.meta.df[s.col]==s.grp]
        temp.s.samples.vec=temp.s.samples.vec[!is.na(temp.s.samples.vec)]
        temp.s.samples.vec=intersect(colnames(s.rpkm.df),temp.s.samples.vec)
        temp.s.expr=subset(s.rpkm.df,select=temp.s.samples.vec)
        temp.s.dim=dim(temp.s.expr)
        temp.title=paste(c(as.character(f.grp),'_vs_',as.character(s.grp),gn),sep  ='')
        print(temp.title)
        plot.reid.data.expression.boxplot(first.rpkm.df = temp.f.expr,sec.rpkm.df =temp.s.expr,gene.id =gn,title.str = temp.title)
        #print(c(as.character(f.grp),as.character(s.grp)))
      }
    }
  }
}
plot.reid.data.expression.boxplot.per.gene=function(f.rpkm.df,s.rpkm.df,gene.vec,f.name='SP1',s.name='SP2',title.str=NULL){
  for (gn in gene.vec){
    temp.title=gn
    plot.reid.data.expression.boxplot(first.rpkm.df=f.rpkm.df,sec.rpkm.df=s.rpkm.df,
                                      gene.id=gn,f.name=f.name,s.name=s.name,title.str=temp.title)
  }
}
plot.reid.data.expression.boxplot=function(first.rpkm.df,sec.rpkm.df,gene.id,f.name='SP1',s.name='SP2',title.str=NULL){
  first.rpkm.df=filter.none.expressed.genes(first.rpkm.df)
  sec.rpkm.df=filter.none.expressed.genes(sec.rpkm.df)
  if(!gene.id %in% rownames(first.rpkm.df) && !gene.id %in% rownames(sec.rpkm.df)){
    print('Gene not found!!')
  }
  else{
    temp.f.expr.vec=try(as.numeric(first.rpkm.df[gene.id,]))
    if(class(temp.f.expr.vec) == "try-error"){
      temp.f.expr.vec=rep(0.0,10)
    }
    temp.f.expr.vec[is.na(temp.f.expr.vec)]=0
    temp.s.expr.vec=try(as.numeric(sec.rpkm.df[gene.id,]))
    if(class(temp.s.expr.vec)=="try-error"){
      temp.s.expr.vec=rep(0.0,10)
    }
    temp.s.expr.vec[is.na(temp.s.expr.vec)]=0
    f.box=boxplot(temp.f.expr.vec,plot=F)
    temp.f.expr.vec=temp.f.expr.vec[!temp.f.expr.vec %in% f.box$out]
    s.box=boxplot(temp.s.expr.vec,plot=F)
    temp.s.expr.vec=temp.s.expr.vec[!temp.s.expr.vec %in% s.box$out]
    title.str=ifelse(is.null(title.str),gene.id,title.str)
    boxplot(list(temp.f.expr.vec+1,temp.s.expr.vec+1),names=c(f.name,s.name),log='y',main=title.str,pch=19,frame.plot=F,col=c('blue','red'))
    #boxplot(list(temp.f.expr.vec,temp.s.expr.vec),names=c(f.name,s.name),main=title.str,pch=19,frame.plot=F,col=c('blue','red'))
    #boxplot(list(log2(temp.f.expr.vec+1),log2(temp.s.expr.vec+1)),names=c(f.name,s.name),main=title.str,pch=19,frame.plot=F,col=c('blue','red'))
  }
}
get.diff.expression.reid.data=function(first.counts.df,first.meta.df,sec.counts.df,sec.meta.df,
                                       gene.vec,first.col,sec.col){
  first.grps=unique(first.meta.df[,first.col])
  sec.grps=unique(sec.meta.df[,sec.col])
  expr.list=list()
  col.vec=c()
  for (f in first.grps){
    if(is.na(f)){
      next
    }
    if(f=='NA'){
      next
    }
    for (s in sec.grps){
      if(is.na(s)){
        next
      }
      temp.f.samples.vec=rownames(first.meta.df)[first.meta.df[first.col]==f]
      temp.f.samples.vec=temp.f.samples.vec[!is.na(temp.f.samples.vec)]
      temp.f.samples.vec=intersect(colnames(first.counts.df),temp.f.samples.vec)
      temp.f.cnts.df=subset(first.counts.df,select=temp.f.samples.vec)
      temp.f.cnts.df=filter.none.expressed.genes(temp.f.cnts.df)
      temp.s.samples.vec=rownames(sec.meta.df)[sec.meta.df[sec.col]==s]
      temp.s.samples.vec=temp.s.samples.vec[!is.na(temp.s.samples.vec)]
      temp.s.samples.vec=intersect(colnames(sec.counts.df),temp.s.samples.vec)
      temp.s.cnts.df=subset(sec.counts.df,select=temp.s.samples.vec)
      temp.s.cnts.df=filter.none.expressed.genes(temp.s.cnts.df)
      temp.f.dim=dim(temp.f.cnts.df)
      temp.s.dim=dim(temp.s.cnts.df)
      if(temp.f.dim[2]<50 || temp.s.dim[2]<50){
        next
      }
      intersect.gns=intersect(rownames(temp.f.cnts.df),rownames(temp.s.cnts.df))
      temp.cnts.df=data.frame(temp.f.cnts.df[intersect.gns,],temp.s.cnts.df[intersect.gns,])
      temp.meta.df=data.frame(row.names = colnames(temp.cnts.df),samples=colnames(temp.cnts.df),
                              grps=c(rep('asex',temp.f.dim[2]),rep('sex',temp.s.dim[2])))
      #temp.scde.res=run.scde.between.grps(counts.df = temp.cnts.df,meta.df = temp.meta.df,grp = 'grps')
      temp.grp.name=paste(f,'_vs_',s,'f',sep='')
      #expr.list[[temp.grp.name]]=temp.scde.res
      colnames(temp.cnts.df)=c(rep('asex',temp.f.dim[2]),rep('sex',temp.s.dim[2]))
      #outfname=paste(c('~/Downloads/',temp.grp.name,'.tab'),collapse ='')
      write.table(temp.cnts.df,file = outfname,sep='\t',quote = F)
    }
  }
  return(expr.list)
}
plot.specific.genes.expression.boxplot.btw.sc.and.bulk=function(rpkm.df,meta.df,pop.rpkm.df,pop.meta.df,gene.vec,parental=F){
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(gene.vec,y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.genes(in.rpkm.df)
  pop.rpkm.df=filter.none.expressed.genes(input.data = pop.rpkm.df)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  meta.df=meta.df[colnames(in.rpkm.df),]
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=list()
  intersect.development.stages=c()
  if(parental){
    col.str=c()
    meta.df=create.parental.stage(meta.df =  meta.df)
    pop.meta.df=create.parental.stage(meta.df = pop.meta.df)
    meta.list=split(meta.df,f = meta.df$parent.development.stage)
    pop.meta.list=split(pop.meta.df,f = pop.meta.df$parent.development.stage)
    development.stages=names(meta.list)
    ordered.development.stages=c("R" , "T" ,  "S")
    intersect.development.stages = intersect(ordered.development.stages,development.stages)
    len.intersect.development.stages=length(intersect.development.stages)
    for(n in 1:len.expressed.markers.vec){
      expressed.marker=expressed.markers.vec[n]
      temp.list=list()
      for(m in 1:len.intersect.development.stages){
        intersect.development.stage=intersect.development.stages[m]
        temp.meta.df=meta.list[[intersect.development.stage]]
        temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
        temp.rpkm=as.numeric(as.character(temp.rpkm.df[expressed.marker,]))
        temp.rpkm=ifelse(temp.rpkm==0,0.001,temp.rpkm)
        temp.pop.meta.df=pop.meta.list[[intersect.development.stage]]
        temp.pop.rpkm.df=pop.rpkm.df[,rownames(temp.pop.meta.df)]
        temp.pop.rpkm=as.numeric(as.character(temp.pop.rpkm.df[expressed.marker,]))
        temp.pop.rpkm=ifelse(temp.pop.rpkm==0,0.001,temp.pop.rpkm)
        sc.str=paste('SCs:',intersect.development.stage,collapse='')
        pop.str=paste('Pop:',intersect.development.stage,collapse='')
        temp.list[[sc.str]]=temp.rpkm
        temp.list[[pop.str]]=temp.pop.rpkm
        col.str=append(col.str,values = c('red','blue'),length(col.str))
      }
      boxplot(temp.list,las=2,main=expressed.marker,log = 'y',pars = list(boxwex = 0.5, staplewex = 0.3, outwex = 0.3),pch=19,col=col.str,outline = T)
    }
  }
  else{
    meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
    meta.list=split(meta.df,f = meta.df$development.stage)
    development.stages=names(meta.list)
    ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
    intersect.development.stages = intersect(ordered.development.stages,development.stages)
    len.intersect.development.stages=length(intersect.development.stages)
    for(n in 1:len.expressed.markers.vec){
      expressed.marker=expressed.markers.vec[n]
      temp.list=list()
      for(m in 1:len.intersect.development.stages){
        intersect.development.stage=intersect.development.stages[m]
        temp.meta.df=meta.list[[intersect.development.stage]]
        temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
        temp.rpkm=as.numeric(as.character(temp.rpkm.df[expressed.marker,]))
        temp.rpkm=ifelse(temp.rpkm==0,0.001,temp.rpkm)
        temp.list[[intersect.development.stage]]=temp.rpkm
      }
      boxplot(temp.list,las=2,main=expressed.marker,log = 'y')
    }
  }
}
plot.specific.genes.expression.barplot.btw.sc.and.bulk=function(rpkm.df,meta.df,pop.rpkm.df,pop.meta.df,gene.vec,parental=F){
  all.genes=rownames(rpkm.df)
  #gene.symbols.vec=names(gene.vec)
  markers.vec=intersect(gene.vec,y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  #in.rpkm.df=filter.none.expressed.genes(in.rpkm.df)
  #pop.rpkm.df=filter.none.expressed.genes(input.data = pop.rpkm.df)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  meta.df=meta.df[colnames(in.rpkm.df),]
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=list()
  intersect.development.stages=c()
  if(parental){
    col.str=c()
    meta.df=create.parental.stage(meta.df =  meta.df)
    pop.meta.df=create.parental.stage(meta.df = pop.meta.df)
    meta.list=split(meta.df,f = meta.df$parent.development.stage)
    pop.meta.list=split(pop.meta.df,f = pop.meta.df$parent.development.stage)
    development.stages=names(meta.list)
    ordered.development.stages=c("R" , "T" ,  "S")
    intersect.development.stages = intersect(ordered.development.stages,development.stages)
    len.intersect.development.stages=length(intersect.development.stages)
    for(n in 1:len.expressed.markers.vec){
      expressed.marker=expressed.markers.vec[n]
      temp.list=list()
      values.vec=c()
      mean.vec=c()
      se.vec=c()
      col.vec=c()
      row.names.vec=c()
      for(m in 1:len.intersect.development.stages){
        intersect.development.stage=intersect.development.stages[m]
        temp.meta.df=meta.list[[intersect.development.stage]]
        temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
        temp.rpkm=as.numeric(as.character(temp.rpkm.df[expressed.marker,]))
        temp.pop.meta.df=pop.meta.list[[intersect.development.stage]]
        temp.pop.rpkm.df=pop.rpkm.df[,rownames(temp.pop.meta.df)]
        temp.pop.rpkm=as.numeric(as.character(temp.pop.rpkm.df[expressed.marker,]))
        sc.str=paste('SCs:',intersect.development.stage,collapse='')
        pop.str=paste('Pop:',intersect.development.stage,collapse='')
        row.names.vec=append(x = row.names.vec,values =c(sc.str,pop.str) ,after = length(row.names.vec))
        sc.expr.stats=describe(temp.rpkm)
        pop.expr.stats=describe(temp.pop.rpkm)
        values.vec=append(x = values.vec,values =c(m,m+1) ,after = length(values.vec))
        mean.vec=append(x = mean.vec,values = c(sc.expr.stats$mean,pop.expr.stats$mean),after = length(mean.vec))
        se.vec=append(x = se.vec,values =c(sc.expr.stats$se,pop.expr.stats$se) ,after = length(se.vec))
        col.vec=append(col.vec,values = c('red','blue'),length(col.vec))
      }
      marker.stats.df=data.frame(values=values.vec,mean=mean.vec,se=se.vec)
      rownames(marker.stats.df)=row.names.vec
      width.size.vec=c(rep(0.1,times = dim(marker.stats.df)[1]))
      editted.error.bars(stats=marker.stats.df,bars=T,main=expressed.marker,ylab = 'rpkm',xlab='',eyes = F,col.str =col.vec,width.size.vec=width.size.vec)
    }
  }
  else{
    sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
    sc.meta.list=split(sc.meta.df,f = sc.meta.df$development.stage)
    pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = pop.meta.df)
    pop.meta.list=split(pop.meta.df,f=pop.meta.df$development.stage)
    development.stages=as.character(intersect(names(sc.meta.list),names(pop.meta.list)))
    ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
    intersect.development.stages = intersect(ordered.development.stages,development.stages)
    len.intersect.development.stages=length(intersect.development.stages)
    for(n in 1:len.expressed.markers.vec){
      expressed.marker=expressed.markers.vec[n]
      temp.list=list()
      values.vec=c()
      mean.vec=c()
      se.vec=c()
      col.vec=c()
      row.names.vec=c()
      for(m in 1:len.intersect.development.stages){
        intersect.development.stage=intersect.development.stages[m]
        temp.meta.df=sc.meta.list[[intersect.development.stage]]
        temp.rpkm.df=subset(in.rpkm.df,select = rownames(temp.meta.df))
        temp.rpkm=as.numeric(as.character(temp.rpkm.df[expressed.marker,]))
        temp.rpkm=ifelse(temp.rpkm==0,0.001,temp.rpkm)
        temp.pop.meta.df=pop.meta.list[[intersect.development.stage]]
        temp.pop.rpkm.df=subset(pop.rpkm.df,select=rownames(temp.pop.meta.df))
        #show(dim(temp.pop.rpkm.df))
        pop.rpkm.vec=as.numeric(as.character(temp.pop.rpkm.df[expressed.marker,]))
        sc.str=paste('SCs:',intersect.development.stage,collapse='')
        pop.str=paste('Pop:',intersect.development.stage,collapse='')
        row.names.vec=append(x = row.names.vec,values =c(sc.str,pop.str) ,after = length(row.names.vec))
        sc.expr.stats=describe(temp.rpkm)
        pop.expr.stats=describe(pop.rpkm.vec)
        values.vec=append(x = values.vec,values =c(m,m+1) ,after = length(values.vec))
        mean.vec=append(x = mean.vec,values = c(sc.expr.stats$mean,pop.expr.stats$mean),after = length(mean.vec))
        se.vec=append(x = se.vec,values =c(sc.expr.stats$se,pop.expr.stats$se) ,after = length(se.vec))
        col.vec=append(col.vec,values = c('red','blue'),length(col.vec))
        temp.list[[intersect.development.stage]]=temp.rpkm
      }
      marker.stats.df=data.frame(values=values.vec,mean=mean.vec,se=se.vec)
      rownames(marker.stats.df)=row.names.vec
      width.size.vec=c(rep(5,times = dim(marker.stats.df)[1]))
      editted.error.bars(stats=marker.stats.df,bars=T,
                         main=expressed.marker,ylab = 'rpkm',xlab='',
                         eyes = F,col.str =col.vec,width.size.vec=width.size.vec)
      #error.bars(stats=marker.stats.df,bars=T,
                  #main=expressed.marker,ylab = 'rpkm',xlab='',
                  #eyes = F,col =col.vec)
    }
  }
}
plot.surface.markers.genes.expression.boxplot=function(rpkm.df,meta.df,markers.df){
  rpkm.df=filter.none.expressed.genes( df = rpkm.df)
  all.genes=rownames(rpkm.df)
  markers.vec=intersect(rownames(markers.df),y =all.genes  )
  in.rpkm.df=rpkm.df[markers.vec,]
  in.rpkm.df=filter.none.expressed.genes(in.rpkm.df)
  #in.rpkm.df=filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = in.rpkm.df,rpkm.cutoff = 1,no.samples = 1))
  meta.df=meta.df[colnames(in.rpkm.df),]
  expressed.markers.vec=rownames(in.rpkm.df)
  len.expressed.markers.vec=length(expressed.markers.vec)
  meta.list=split(meta.df,f = meta.df$development.stage)
  development.stages=names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "trophozoite" , "early.schizont"  , "schizont")
  intersect.development.stages = intersect(ordered.development.stages,development.stages)
  len.intersect.development.stages=length(intersect.development.stages)
  for(n in 1:len.expressed.markers.vec){
    expressed.marker=expressed.markers.vec[n]
    category=paste('(',as.character(markers.df[expressed.marker,'category']),')',sep='')
    lab.str=paste(expressed.marker,category,sep=' ')
    temp.list=list()
    temp.median=c()
    temp.mean=c()
    for(m in 1:len.intersect.development.stages){
      intersect.development.stage=intersect.development.stages[m]
      temp.meta.df=meta.list[[intersect.development.stage]]
      temp.rpkm.df=in.rpkm.df[,rownames(temp.meta.df)]
      temp.rpkm=as.numeric(as.character(temp.rpkm.df[expressed.marker,]))
      temp.median=append(temp.median,median(temp.rpkm),length(temp.median))
      temp.mean=append(temp.mean,mean(temp.rpkm),length(temp.mean))
      temp.rpkm=ifelse(temp.rpkm==0,0.001,temp.rpkm)
      temp.list[[intersect.development.stage]]=temp.rpkm
    }
    #if (!all(c(len.timepoint.one.samples>=3,len.timepoint.two.samples>=3))){
      #next
    #}
    if(all(temp.median==0.0)){
      next
    }
    #if(all(temp.mean<=10.0)){
      #next
    #}
    stage.names=names(temp.list)
    stage.names=gsub("ring" ,'R',gsub("schizont" ,'S',gsub("early.schizont",'ES',gsub("trophozoite",'LT',gsub("early.trophozoite" ,'ET', gsub("late.ring" ,'LR',stage.names))))))              
    boxplot(temp.list,las=2,main=lab.str,log = 'y',names = stage.names)
    #temp.ggplot=ggplot(temp.markers.dotplot.df, aes(x = ordered.samples.categories, y = names)) + geom_point(aes(size=log.mean.expression,colour=fraction.detected)) + ggtitle(label = markers.dotplot.list.name)+ylab(label = '')+xlab(label = '') +theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
}
encode.development.stage=function(meta.df){
  timepoints=as.character(meta.df$timepoint)
  out.timepoint=c()
  samples.names=rownames(meta.df)
  development.stage=c()
  for(m in 1:length(timepoints)){
    timepoint=timepoints[m]
    if(timepoint=='none'){
      timepoint=samples.names[m]
      if(grepl('10h',timepoint)){
        development.stage=append(development.stage,values = 'ring',after=length(development.stage))
        out.timepoint=append(out.timepoint,'10h',after = length(out.timepoint))
      }
      else if(grepl('16h',timepoint)){
        development.stage=append(development.stage,values = 'late.ring',after=length(development.stage))
        out.timepoint=append(out.timepoint,'16h',after = length(out.timepoint))
      }
      else if(grepl('15h',timepoint)){
        development.stage=append(development.stage,values = 'late.ring',after=length(development.stage))
        out.timepoint=append(out.timepoint,'15h',after = length(out.timepoint))
      }
      else if(grepl('20h',timepoint)){
        development.stage=append(development.stage,values = 'early.trophozoite',after=length(development.stage))
        out.timepoint=append(out.timepoint,'20h',after = length(out.timepoint))
      }
      else if(grepl('24h',timepoint)){
        development.stage=append(development.stage,values = 'late.trophozoite',after=length(development.stage))
        out.timepoint=append(out.timepoint,'24h',after = length(out.timepoint))
      }
      else if(grepl('tp.24h',timepoint)){
        development.stage=append(development.stage,values = 'late.trophozoite',after=length(development.stage))
        out.timepoint=append(out.timepoint,'24h',after = length(out.timepoint))
      }
      else if(grepl('25h',timepoint)){
        development.stage=append(development.stage,values = 'late.trophozoite',after=length(development.stage))
        out.timepoint=append(out.timepoint,'25h',after = length(out.timepoint))
      }
      else if(grepl('30h',timepoint)){
        development.stage=append(development.stage,values = 'late.trophozoite',after=length(development.stage))
        out.timepoint=append(out.timepoint,'30h',after = length(out.timepoint))
      }
      else if(grepl('36h',timepoint)){
        development.stage=append(development.stage,values = 'early.schizont',after=length(development.stage))
        out.timepoint=append(out.timepoint,'36h',after = length(out.timepoint))
      }
      else if(grepl('tp.36h',timepoint)){
        development.stage=append(development.stage,values = 'early.schizont',after=length(development.stage))
        out.timepoint=append(out.timepoint,'36h',after = length(out.timepoint))
      }
      else if(grepl('40h',timepoint)){
        development.stage=append(development.stage,values = 'schizont',after=length(development.stage))
        out.timepoint=append(out.timepoint,'40h',after = length(out.timepoint))
      }
      else if(grepl('tp1',timepoint)){
        development.stage=append(development.stage,values = 'ring',after=length(development.stage))
        out.timepoint=append(out.timepoint,'10h',after = length(out.timepoint))
      }
      else if(grepl('tp2',timepoint)){
        development.stage=append(development.stage,values = 'early.trophozoite',after=length(development.stage))
        out.timepoint=append(out.timepoint,'20h',after = length(out.timepoint))
      }
      else if(grepl('tp3',timepoint)){
        development.stage=append(development.stage,values = 'late.trophozoite',after=length(development.stage))
        out.timepoint=append(out.timepoint,'24h',after = length(out.timepoint))
      }
      else if(grepl('tp4',timepoint)){
        development.stage=append(development.stage,values = 'early.schizont',after=length(development.stage))
        out.timepoint=append(out.timepoint,'30h',after = length(out.timepoint))
      }
      else if(grepl('tp5',timepoint)){
        development.stage=append(development.stage,values = 'schizont',after=length(development.stage))
        out.timepoint=append(out.timepoint,'36h',after = length(out.timepoint))
      }
      else if(grepl('c_troph',timepoint)){
        development.stage=append(development.stage,values = 'late.trophozoite',after=length(development.stage))
        out.timepoint=append(out.timepoint,'24h',after = length(out.timepoint))
      }
      else{
        development.stage=append(development.stage,values = 'none',after=length(development.stage))
        out.timepoint=append(out.timepoint,'0h',after = length(out.timepoint))
      }
    } 
    else{
      timepoint=timepoints[m] 
    if(timepoint=='tp.16h'){
      development.stage=append(development.stage,values = 'late.ring',after=length(development.stage))
      out.timepoint=append(out.timepoint,'16h',after = length(out.timepoint))
    }
    else if(timepoint=='tp.15h'){
      development.stage=append(development.stage,values = 'late.ring',after=length(development.stage))
      out.timepoint=append(out.timepoint,'15h',after = length(out.timepoint))
    }
    else if(timepoint=='tp.24-30h'){
      development.stage=append(development.stage,values = 'late.trophozoite',after=length(development.stage))
      out.timepoint=append(out.timepoint,'24h',after = length(out.timepoint))
    }
    else if(timepoint=='tp.24h'){
      development.stage=append(development.stage,values = 'late.trophozoite',after=length(development.stage))
      out.timepoint=append(out.timepoint,'24h',after = length(out.timepoint))
    }
    else if(timepoint=='tp.36h'){
      development.stage=append(development.stage,values = 'early.schizont',after=length(development.stage))
      out.timepoint=append(out.timepoint,'36h',after = length(out.timepoint))
    }
    else if(timepoint=='tp.40h'){
      development.stage=append(development.stage,values = 'schizont',after=length(development.stage))
      out.timepoint=append(out.timepoint,'40h',after = length(out.timepoint))
    }
    else if(timepoint=='tp1'){
      development.stage=append(development.stage,values = 'ring',after=length(development.stage))
      out.timepoint=append(out.timepoint,'10h',after = length(out.timepoint))
    }
    else if(timepoint=='tp2'){
      development.stage=append(development.stage,values = 'early.trophozoite',after=length(development.stage))
      out.timepoint=append(out.timepoint,'20h',after = length(out.timepoint))
    }
    else if(timepoint=='tp3'){
      development.stage=append(development.stage,values = 'late.trophozoite',after=length(development.stage))
      out.timepoint=append(out.timepoint,'24h',after = length(out.timepoint))
    }
    else if(timepoint=='tp4'){
      development.stage=append(development.stage,values = 'early.schizont',after=length(development.stage))
      out.timepoint=append(out.timepoint,'36h',after = length(out.timepoint))
    }
    else if(timepoint=='tp5'){
      development.stage=append(development.stage,values = 'schizont',after=length(development.stage))
      out.timepoint=append(out.timepoint,'40h',after = length(out.timepoint))
    }
    else{
      development.stage=append(development.stage,values = 'none',after=length(development.stage))
      out.timepoint=append(out.timepoint,'0h',after = length(out.timepoint))
    }
    }
  }
  meta.df$development.stage=development.stage
  meta.df$timepoint=out.timepoint
  return(meta.df)
}
plot.pca.at.different.var.cut.off = function(in.rpkm.df,in.meta.df){
  in.meta.df=in.meta.df[colnames(in.rpkm.df),]
  genes.var=as.numeric(apply(in.rpkm.df,1,sd))
  #var.quantiles=as.numeric(quantile(genes.var))[1:3]
  cut.offs=seq(0,to = 1000,by = 100)
  #cut.offs=var.quantiles
  for(m in 1:length(cut.offs)){
    temp.cut.off=cut.offs[m]
    temp.rpkm.df=filter.non.variable.rows(df = in.rpkm.df,cut.off =temp.cut.off)
    temp.rpkm.df=filter.none.expressed.samples(df = temp.rpkm.df)
    temp.meta.df=in.meta.df[colnames(temp.rpkm.df),]
    plot.new()
    label.str=paste('Gene variance cut off:',temp.cut.off,sep=' ')
    text(x = .5,y = .5,labels = label.str)
    temp.plot.pca(in.rpkm.df=temp.rpkm.df,trans=F,meta.data=temp.meta.df,title.str='SCs pca',log.pca.scores=F,plot.pairs = F)
  }
}
plot.sc.sub.pop.vs.pop.per.development.stage.cor.boxplots=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df,down.sampling.replicates=10){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df)
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$markers.cluster.groups)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = pop.meta.df)
  pop.meta.df$development.stage=factor(pop.meta.df$development.stage,factor(std.stages.col.abbrev.legend.str))
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$development.stage)
  pop.timepoints=names(pop.meta.list)
  len.pop.timepoints=length(pop.timepoints)
  intersect.timepoints=sub.populations.grp.names
  stages.col.list=get.color.list.for.pheatmap(in.vec  =intersect.timepoints )
  len.stages.col.list=length(stages.col.list)
  col.list=list()
  for(n in 1:length(intersect.timepoints)){
    temp.cor.list=list()
    out.list=list()
    p.values=c()
    p.values=append( p.values,'NA',length( p.values))
    intersect.timepoint=intersect.timepoints[n]
    temp.sc.meta.df=sc.meta.list[[intersect.timepoint]]
    temp.sc.rpkm.df=filter.none.expressed.genes(input.data = sc.rpkm.df[,rownames(temp.sc.meta.df)])
    temp.col.vec=c()
    temp.all.cor.vec=c()
    cor.fact=c()
    for(m in 1:len.pop.timepoints){
      temp.pop.timepoints=pop.timepoints[m]
      temp.pop.meta.df=pop.meta.list[[temp.pop.timepoints]]
      temp.pop.rpkm.df=filter.none.expressed.genes(input.data = pop.rpkm.df[,rownames(temp.pop.meta.df)])
      temp.cor.score.vec=get.cor.between.two.df.cor(first.df = temp.sc.rpkm.df,second.df = temp.pop.rpkm.df)
      temp.cor.list[[temp.pop.timepoints]]=temp.cor.score.vec
      temp.all.cor.vec=append(temp.all.cor.vec,values = temp.cor.score.vec,after = length(temp.all.cor.vec))
      cor.fact=append(cor.fact,values =rep(temp.pop.timepoints,times = length(temp.cor.score.vec)) ,after = length(cor.fact))
      temp.col.vec=append(temp.col.vec,abbrev.std.stages.col.list[[temp.pop.timepoints]],length(temp.col.vec))
    }
    temp.all.downsampled.cor.vec=c()
    temp.all.downsampled.cor.fact.vec=c()
    temp.cor.df=data.frame(cor.coeff=temp.all.cor.vec)
    median.cor.mat=matrix(ncol = length(std.stages.col.abbrev.legend.str),nrow = down.sampling.replicates)
    down.sampled.list.names=intersect(std.stages.col.abbrev.legend.str,unique(cor.fact))
    len.down.sampled.list.names=length(down.sampled.list.names)
    for(d in 1:down.sampling.replicates){
      down.sampled.list=downSample(x = temp.cor.df,y =factor(cor.fact),list = T )
      down.sampled.df=down.sampled.list$x
      down.sampled.df$class=as.character(down.sampled.list$y)
      temp.all.downsampled.cor.vec=append(x = temp.all.downsampled.cor.vec,as.numeric(down.sampled.df$cor.coeff),after = length(temp.all.downsampled.cor.vec))
      down.sampled.list=split(down.sampled.df,f=down.sampled.df$class)
      temp.all.downsampled.cor.fact.vec=append(x = temp.all.downsampled.cor.fact.vec,values = as.character(down.sampled.df$class),after = length(temp.all.downsampled.cor.fact.vec))
    }
    temp.down.sampled.df=data.frame(cor.coeff=temp.all.downsampled.cor.vec,category=temp.all.downsampled.cor.fact.vec)
    temp.down.sampled.list=split(temp.down.sampled.df,f = temp.down.sampled.df$category)
    down.sampled.cor.list=lapply(temp.down.sampled.list,function(list.item){
      temp.df=list.item
      cor.coeff=as.numeric(temp.df$cor.coeff)
      return(cor.coeff)
    })
    down.sampled.ordered.cor.list=list()
    for(b in 1:len.down.sampled.list.names){
      temp.stage=down.sampled.list.names[b]
      down.sampled.ordered.cor.list[[temp.stage]]=down.sampled.cor.list[[temp.stage]]
    }
    title.str=gsub(pattern = '\\.',replacement = ': ', x = intersect.timepoint)
    boxplot(down.sampled.ordered.cor.list,main=title.str,cex=.8,col = temp.col.vec,pch=19)
  }
}
plot.sc.sub.pop.average.vs.pop.per.development.stage.cor.boxplots=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df,down.sampling.replicates=10){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df)
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$markers.cluster.groups)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = pop.meta.df)
  pop.meta.df$development.stage=factor(pop.meta.df$development.stage,factor(std.stages.col.abbrev.legend.str))
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$development.stage)
  pop.timepoints=names(pop.meta.list)
  len.pop.timepoints=length(pop.timepoints)
  intersect.timepoints=sub.populations.grp.names
  stages.col.list=get.color.list.for.pheatmap(in.vec  =intersect.timepoints )
  len.stages.col.list=length(stages.col.list)
  col.list=list()
  for(n in 1:length(intersect.timepoints)){
    temp.cor.list=list()
    out.list=list()
    p.values=c()
    p.values=append( p.values,'NA',length( p.values))
    intersect.timepoint=intersect.timepoints[n]
    temp.sc.meta.df=sc.meta.list[[intersect.timepoint]]
    temp.sc.rpkm.df=filter.none.expressed.genes(input.data = sc.rpkm.df[,rownames(temp.sc.meta.df)])
    temp.average.sc.rpkm.df=data.frame('average.sc'=as.numeric(apply(temp.sc.rpkm.df,1,mean)))
    rownames(temp.average.sc.rpkm.df) = rownames(temp.sc.rpkm.df)
    temp.col.vec=c()
    temp.all.cor.vec=c()
    cor.fact=c()
    for(m in 1:len.pop.timepoints){
      temp.pop.timepoints=pop.timepoints[m]
      temp.pop.meta.df=pop.meta.list[[temp.pop.timepoints]]
      temp.pop.rpkm.df=filter.none.expressed.genes(input.data = pop.rpkm.df[,rownames(temp.pop.meta.df)])
      intersect.genes=intersect(rownames(temp.average.sc.rpkm.df),rownames(temp.pop.rpkm.df))
      intersect.temp.sc.average.rpkm.df=data.frame(t(subset(t(temp.average.sc.rpkm.df),select = intersect.genes)))
      sc.expr.vec=as.numeric(as.character(intersect.temp.sc.average.rpkm.df[,'average.sc']))
      temp.pop.rpkm.df=temp.pop.rpkm.df[intersect.genes,]
      temp.cor.score.vec=as.numeric(unlist(apply(temp.pop.rpkm.df,2,function(pop.col.vec){
        pop.col.vec=as.numeric(pop.col.vec)
        temp.cor.coeff=cor(x=pop.col.vec,y = sc.expr.vec,method = 'spearman')
        return(temp.cor.coeff)
      })))
      temp.cor.list[[temp.pop.timepoints]]=temp.cor.score.vec
      temp.col.vec=append(temp.col.vec,abbrev.std.stages.col.list[[temp.pop.timepoints]],length(temp.col.vec))
    }
    stages.ordered.cor.list=list()
    for(b in 1:len.pop.timepoints){
      temp.stage=pop.timepoints[b]
      stages.ordered.cor.list[[temp.stage]]=temp.cor.list[[temp.stage]]
    }
    title.str=gsub(pattern = '\\.',replacement = ': ', x = intersect.timepoint)
    boxplot(stages.ordered.cor.list,main=title.str,cex=.8,col = temp.col.vec,pch=19)
  }
}
plot.pop.and.sc.cor.boxplots=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df)
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$development.stage)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = pop.meta.df)
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$development.stage)
  pop.timepoints=names(pop.meta.list)
  intersect.timepoints=intersect(sc.timepoints,pop.timepoints)
  input.timepoints=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  intersect.timepoints=intersect(input.timepoints,intersect.timepoints)
  stages.col.list=get.col.factor(col.factor =intersect.timepoints )
  #stages.col.str=stages.col.list$legend.col
  stages.col.str=names(abbrev.std.stages.col.list)
  #stages.name.str=stages.col.list$legend.str
  stages.name.str=as.character(unlist(abbrev.std.stages.col.list))
  stages.name.str=get.abbrev.std.stage.cols(in.stages.vec =intersect.timepoints )
  col.list=list()
  for(d in 1:length(stages.col.str)){
    col.list[[stages.name.str[d]]]=stages.col.str[d]
  }
  for(n in 1:length(intersect.timepoints)){
    temp.cor.list=list()
    out.list=list()
    temp.col.vec=c()
    p.values=c()
    p.values=append( p.values,'NA',length( p.values))
    intersect.timepoint=intersect.timepoints[n]
    #     if(intersect.timepoint=="ring" |intersect.timepoint== "late.ring"){
    #       
    #       next
    #     }
    #temp.col.vec=append(temp.col.vec,col.list[[intersect.timepoint]],length(temp.col.vec))
    temp.col.vec=append(temp.col.vec,abbrev.std.stages.col.list[[intersect.timepoint]],length(temp.col.vec))
    other.timepoints=setdiff(intersect.timepoints,intersect.timepoint)
    temp.sc.meta.df=sc.meta.list[[intersect.timepoint]]
    temp.sc.samples=rownames(temp.sc.meta.df)
    temp.pop.meta.df=pop.meta.list[[intersect.timepoint]]
    temp.pop.samples=rownames(temp.pop.meta.df)
    temp.sc.rpkm.df=sc.rpkm.df[,temp.sc.samples]
    temp.sc.rpkm.df=filter.none.expressed.genes(input.data = temp.sc.rpkm.df)
    temp.pop.rpkm.df=pop.rpkm.df[,temp.pop.samples]
    temp.pop.rpkm.df=filter.none.expressed.genes(temp.pop.rpkm.df)
    temp.same.stage.cor=get.cor.between.two.df.cor(first.df = temp.sc.rpkm.df,second.df = temp.pop.rpkm.df)
    temp.cor.list[[intersect.timepoint]]=temp.same.stage.cor
    len.other.timepoints=length(other.timepoints)
    all.timepoints.cor.vec=temp.same.stage.cor
    all.timepoints.cor.fact=rep(intersect.timepoint,times = length(temp.same.stage.cor))
    for(m in 1:len.other.timepoints){
      temp.other.timepoints=other.timepoints[m]
      temp.other.sc.meta.df=sc.meta.list[[temp.other.timepoints]]
      temp.sc.samples=rownames(temp.other.sc.meta.df)
      other.sc.rpkm.df=sc.rpkm.df[,temp.sc.samples]
      other.sc.rpkm.df=filter.none.expressed.genes(other.sc.rpkm.df)
      #get.cor.between.two.df.cor(first.df = other.sc.rpkm.df,second.df = temp.pop.rpkm.df)
      other.timepoints.cor=get.cor.between.two.df.cor(first.df = other.sc.rpkm.df,second.df = temp.pop.rpkm.df)
      all.timepoints.cor.vec=append(all.timepoints.cor.vec,other.timepoints.cor,length(all.timepoints.cor.vec))
      all.timepoints.cor.fact=append(all.timepoints.cor.fact,rep(temp.other.timepoints,times=length(other.timepoints.cor)),length(all.timepoints.cor.fact))
      temp.all.cor.vec=c(temp.same.stage.cor,other.timepoints.cor)
      temp.all.cor.df=data.frame(cor.coeff=temp.all.cor.vec)
      cor.fact=factor(x = c(rep('s',times = length(temp.same.stage.cor)),rep('d',times = length(other.timepoints.cor))))
      down.sampled.list=downSample(x = temp.all.cor.df,y =cor.fact,list = T )
      down.sampled.df=down.sampled.list$x
      down.sampled.df$class=as.character(down.sampled.list$y)
      down.sampled.list=split(down.sampled.df,f=down.sampled.df$class)
      down.sampled.same.vec=as.numeric(as.character(down.sampled.list$s$cor.coeff))
      down.sampled.different.vec=as.numeric(as.character(down.sampled.list$d$cor.coeff))
      temp.cor.list[[temp.other.timepoints]]=down.sampled.different.vec
      #temp.cor.list[[temp.other.timepoints]]=other.timepoints.cor
      #temp.col.vec=append(temp.col.vec,col.list[[temp.other.timepoints]],length(temp.col.vec))
      temp.col.vec=append(temp.col.vec,abbrev.std.stages.col.list[[temp.other.timepoints]],length(temp.col.vec))
      #wilcox.test.p.value=format(wilcox.test(x=down.sampled.same.vec,y=down.sampled.different.vec,paired=F,alternative = 'greater')$p.value,digits=3)
      wilcox.test.p.value=format(wilcox.test(x=temp.same.stage.cor,y=other.timepoints.cor,paired=F,alternative = 'greater')$p.value,digits=3)
      p.values=append( p.values,wilcox.test.p.value,length( p.values))
    }
    all.timepoints.df=data.frame(cor.coeff=all.timepoints.cor.vec)
    all.timepoints.cor.fact=factor(all.timepoints.cor.fact)
    all.timepoints.downsampled.list=downSample(all.timepoints.df,y = all.timepoints.cor.fact)
    all.timepoints.downsampled.df=all.timepoints.downsampled.list$x
    all.timepoints.downsampled.df$class=all.timepoints.downsampled.list$y
    #intersect.timepoint=convert.to.title.case(paste(sub(pattern ='\\.',replacement =' ',x =   intersect.timepoint),'bulk-RNAseq',sep=' '))
    intersect.timepoint=paste(sub(pattern ='\\.',replacement =' ',x =   intersect.timepoint),'bulk-RNAseq',sep=' ')
    re.ordered.cor.list=list()
    #input.timepoints=c("R" ,"LR", "ET", "T" , "ES"  , "S")
    len.input.timepoints=length(input.timepoints)
    temp.re.ordered.col.vec=c()
    for(g in 1:len.input.timepoints){
      temp.tm.point=input.timepoints[g]
      re.ordered.cor.list[[temp.tm.point]]=temp.cor.list[[temp.tm.point]]
      temp.re.ordered.col.vec=append(temp.re.ordered.col.vec,abbrev.std.stages.col.list[[temp.tm.point]],length(temp.re.ordered.col.vec))
    }
    #boxplot2(temp.cor.list,main=intersect.timepoint,col=temp.col.vec,las=2,names = p.values,cex=.4)
    #boxplot2(re.ordered.cor.list,main=intersect.timepoint,col=temp.re.ordered.col.vec,las=2,names = p.values,cex=.2)
    #boxplot(re.ordered.cor.list,main=intersect.timepoint,col=temp.re.ordered.col.vec,las=2,names = p.values,cex=.2)
    boxplot(re.ordered.cor.list,main='',col=temp.re.ordered.col.vec,las=2,names = p.values,cex=.2)
    #boxplot2(temp.cor.list,main=intersect.timepoint,col=temp.col.vec,las=2,cex=.4)
    #legend('topright',legend=stages.name.str,fill=stages.col.str,box.lty = 0,cex = 1.5)
    #plot.new()
    #legend('center',legend=stages.name.str,fill=stages.col.str,box.lty = 0,cex = 1.5)
    #boxplot(temp.cor.list,main='',col=temp.col.vec,las=2,names = rep('',time=length(p.values)))
  }
  #plot.new()
  #legend('center',legend=rep('',length(std.stages.col.abbrev.legend.str)),fill=std.stages.col.legend.fill,box.lty = 0,cex = 1.5)
}
plot.pop.and.sc.cor.boxplots=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df)
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$development.stage)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = pop.meta.df)
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$development.stage)
  pop.timepoints=names(pop.meta.list)
  intersect.timepoints=intersect(sc.timepoints,pop.timepoints)
  input.timepoints=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  intersect.timepoints=intersect(input.timepoints,intersect.timepoints)
  stages.col.list=get.col.factor(col.factor =intersect.timepoints )
  #stages.col.str=stages.col.list$legend.col
  stages.col.str=names(abbrev.std.stages.col.list)
  #stages.name.str=stages.col.list$legend.str
  stages.name.str=as.character(unlist(abbrev.std.stages.col.list))
  stages.name.str=get.abbrev.std.stage.cols(in.stages.vec =intersect.timepoints )
  col.list=list()
  for(d in 1:length(stages.col.str)){
    col.list[[stages.name.str[d]]]=stages.col.str[d]
  }
  for(n in 1:length(intersect.timepoints)){
    temp.cor.list=list()
    out.list=list()
    temp.col.vec=c()
    p.values=c()
    p.values=append( p.values,'NA',length( p.values))
    intersect.timepoint=intersect.timepoints[n]
    #     if(intersect.timepoint=="ring" |intersect.timepoint== "late.ring"){
    #       
    #       next
    #     }
    #temp.col.vec=append(temp.col.vec,col.list[[intersect.timepoint]],length(temp.col.vec))
    temp.col.vec=append(temp.col.vec,abbrev.std.stages.col.list[[intersect.timepoint]],length(temp.col.vec))
    other.timepoints=setdiff(intersect.timepoints,intersect.timepoint)
    temp.sc.meta.df=sc.meta.list[[intersect.timepoint]]
    temp.sc.samples=rownames(temp.sc.meta.df)
    temp.pop.meta.df=pop.meta.list[[intersect.timepoint]]
    temp.pop.samples=rownames(temp.pop.meta.df)
    temp.sc.rpkm.df=sc.rpkm.df[,temp.sc.samples]
    temp.sc.rpkm.df=filter.none.expressed.genes(input.data = temp.sc.rpkm.df)
    temp.pop.rpkm.df=pop.rpkm.df[,temp.pop.samples]
    temp.pop.rpkm.df=filter.none.expressed.genes(temp.pop.rpkm.df)
    temp.same.stage.cor=get.cor.between.two.df.cor(first.df = temp.sc.rpkm.df,second.df = temp.pop.rpkm.df)
    temp.cor.list[[intersect.timepoint]]=temp.same.stage.cor
    len.other.timepoints=length(other.timepoints)
    all.timepoints.cor.vec=temp.same.stage.cor
    all.timepoints.cor.fact=rep(intersect.timepoint,times = length(temp.same.stage.cor))
    for(m in 1:len.other.timepoints){
      temp.other.timepoints=other.timepoints[m]
      temp.other.sc.meta.df=sc.meta.list[[temp.other.timepoints]]
      temp.sc.samples=rownames(temp.other.sc.meta.df)
      other.sc.rpkm.df=sc.rpkm.df[,temp.sc.samples]
      other.sc.rpkm.df=filter.none.expressed.genes(other.sc.rpkm.df)
      #get.cor.between.two.df.cor(first.df = other.sc.rpkm.df,second.df = temp.pop.rpkm.df)
      other.timepoints.cor=get.cor.between.two.df.cor(first.df = other.sc.rpkm.df,second.df = temp.pop.rpkm.df)
      all.timepoints.cor.vec=append(all.timepoints.cor.vec,other.timepoints.cor,length(all.timepoints.cor.vec))
      all.timepoints.cor.fact=append(all.timepoints.cor.fact,rep(temp.other.timepoints,times=length(other.timepoints.cor)),length(all.timepoints.cor.fact))
      temp.all.cor.vec=c(temp.same.stage.cor,other.timepoints.cor)
      temp.all.cor.df=data.frame(cor.coeff=temp.all.cor.vec)
      cor.fact=factor(x = c(rep('s',times = length(temp.same.stage.cor)),rep('d',times = length(other.timepoints.cor))))
      down.sampled.list=downSample(x = temp.all.cor.df,y =cor.fact,list = T )
      down.sampled.df=down.sampled.list$x
      down.sampled.df$class=as.character(down.sampled.list$y)
      down.sampled.list=split(down.sampled.df,f=down.sampled.df$class)
      down.sampled.same.vec=as.numeric(as.character(down.sampled.list$s$cor.coeff))
      down.sampled.different.vec=as.numeric(as.character(down.sampled.list$d$cor.coeff))
      temp.cor.list[[temp.other.timepoints]]=down.sampled.different.vec
      #temp.cor.list[[temp.other.timepoints]]=other.timepoints.cor
      #temp.col.vec=append(temp.col.vec,col.list[[temp.other.timepoints]],length(temp.col.vec))
      temp.col.vec=append(temp.col.vec,abbrev.std.stages.col.list[[temp.other.timepoints]],length(temp.col.vec))
      #wilcox.test.p.value=format(wilcox.test(x=down.sampled.same.vec,y=down.sampled.different.vec,paired=F,alternative = 'greater')$p.value,digits=3)
      wilcox.test.p.value=format(wilcox.test(x=temp.same.stage.cor,y=other.timepoints.cor,paired=F,alternative = 'greater')$p.value,digits=3)
      p.values=append( p.values,wilcox.test.p.value,length( p.values))
    }
    all.timepoints.df=data.frame(cor.coeff=all.timepoints.cor.vec)
    all.timepoints.cor.fact=factor(all.timepoints.cor.fact)
    all.timepoints.downsampled.list=downSample(all.timepoints.df,y = all.timepoints.cor.fact)
    all.timepoints.downsampled.df=all.timepoints.downsampled.list$x
    all.timepoints.downsampled.df$class=all.timepoints.downsampled.list$y
    #intersect.timepoint=convert.to.title.case(paste(sub(pattern ='\\.',replacement =' ',x =   intersect.timepoint),'bulk-RNAseq',sep=' '))
    intersect.timepoint=paste(sub(pattern ='\\.',replacement =' ',x =   intersect.timepoint),'bulk-RNAseq',sep=' ')
    re.ordered.cor.list=list()
    #input.timepoints=c("R" ,"LR", "ET", "T" , "ES"  , "S")
    len.input.timepoints=length(input.timepoints)
    temp.re.ordered.col.vec=c()
    for(g in 1:len.input.timepoints){
      temp.tm.point=input.timepoints[g]
      re.ordered.cor.list[[temp.tm.point]]=temp.cor.list[[temp.tm.point]]
      temp.re.ordered.col.vec=append(temp.re.ordered.col.vec,abbrev.std.stages.col.list[[temp.tm.point]],length(temp.re.ordered.col.vec))
    }
    #boxplot2(temp.cor.list,main=intersect.timepoint,col=temp.col.vec,las=2,names = p.values,cex=.4)
    #boxplot2(re.ordered.cor.list,main=intersect.timepoint,col=temp.re.ordered.col.vec,las=2,names = p.values,cex=.2)
    #boxplot(re.ordered.cor.list,main=intersect.timepoint,col=temp.re.ordered.col.vec,las=2,names = p.values,cex=.2)
    boxplot(re.ordered.cor.list,main='',col=temp.re.ordered.col.vec,las=2,names = p.values,cex=.2)
    #boxplot2(temp.cor.list,main=intersect.timepoint,col=temp.col.vec,las=2,cex=.4)
    #legend('topright',legend=stages.name.str,fill=stages.col.str,box.lty = 0,cex = 1.5)
    #plot.new()
    #legend('center',legend=stages.name.str,fill=stages.col.str,box.lty = 0,cex = 1.5)
    #boxplot(temp.cor.list,main='',col=temp.col.vec,las=2,names = rep('',time=length(p.values)))
  }
  #plot.new()
  #legend('center',legend=rep('',length(std.stages.col.abbrev.legend.str)),fill=std.stages.col.legend.fill,box.lty = 0,cex = 1.5)
}
plot.sc.sub.states.and.pop.cor.boxplots=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$development.stage)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$development.stage)
  pop.timepoints=names(pop.meta.list)
  intersect.timepoints=pop.timepoints
  input.timepoints=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  intersect.timepoints=intersect(input.timepoints,intersect.timepoints)
  stages.col.list=get.col.factor(col.factor =intersect.timepoints )
  stages.col.str=stages.col.list$legend.col
  stages.name.str=stages.col.list$legend.str
  col.list=list()
  for(d in 1:length(stages.col.str)){
    col.list[[stages.name.str[d]]]=stages.col.str[d]
  }
  len.sc.timepoints=length(sc.timepoints)
  len.pop.timepoints=length(pop.timepoints)
  for(n in 1:len.sc.timepoints){
    sc.timepoint=sc.timepoints[n]
    temp.sc.meta.df=sc.meta.list[[sc.timepoint]]
    temp.sc.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(sc.rpkm.df[,rownames(temp.sc.meta.df)]))
    temp.sc.meta.df= temp.sc.meta.df[colnames(temp.sc.rpkm.df),]
    temp.sc.meta.list=split(temp.sc.meta.df,f=temp.sc.meta.df$sub.group)
    temp.names=names(temp.sc.meta.list)
    len.temp.names=length(temp.names)
    for(m in 1:len.temp.names){
      temp.name=temp.names[m]
      temp.sc.sub.grp.meta.df=temp.sc.meta.list[[temp.name]]
      temp.sc.sub.grp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(temp.sc.rpkm.df[,rownames(temp.sc.sub.grp.meta.df)]))
      temp.sc.sub.grp.meta.df= temp.sc.sub.grp.meta.df[colnames(temp.sc.sub.grp.rpkm.df),]
      temp.dist.list=list()
      for(s in 1:len.pop.timepoints){
        pop.timepoint=pop.timepoints[s]
        temp.pop.meta.df=pop.meta.list[[pop.timepoint]]
        #temp.pop.rpkm.df=pop.rpkm.df[,rownames(temp.pop.meta.df)]
        temp.pop.meta.list=split(temp.pop.meta.df,f=temp.pop.meta.df$sub.group)
        pop.temp.names=names(temp.pop.meta.list)
        len.pop.temp.names=length(pop.temp.names)
        for(f in 1:len.pop.temp.names){
          temp.pop.sub.grp.name=pop.temp.names[f]
          temp.pop.sub.grp.meta.df=temp.pop.meta.list[[temp.pop.sub.grp.name]]
          temp.pop.sub.grp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(pop.rpkm.df[,rownames(temp.pop.sub.grp.meta.df)]))
          temp.dist.list[[paste(pop.timepoint,temp.pop.sub.grp.name,collapse = ':')]]=get.cor.between.two.df.cor(first.df = temp.sc.sub.grp.rpkm.df,second.df =temp.pop.sub.grp.rpkm.df,in.method = 'spearman' )
          #temp.dist.list[[paste(pop.timepoint,temp.pop.sub.grp.name,collapse = ':')]]=get.euclidian.dist.between.two.df.cor(first.df  = temp.sc.sub.grp.rpkm.df,second.df = temp.pop.sub.grp.rpkm.df )
        }
      }
      title.str=gsub('\\.',replacement = ' ',x = convert.to.title.case(temp.name))
      stages.col.list=get.col.factor(col.factor =names(temp.dist.list) )
      stages.col.str=stages.col.list$legend.col
      stages.name.str=stages.col.list$legend.str
      #boxplot2(temp.dist.list,las=2,cex=.4,main=title.str,col=stages.col.str)
      row.labels=names(temp.dist.list)
      boxplot2(temp.dist.list,names=rep('',length(row.labels)),las=2,cex=.4,main='',col=stages.col.str)
      abline(h = as.numeric(quantile(as.numeric(unlist(temp.dist.list)),.75)),b = 1,col='lightgray')
      #legend('topright',legend=stages.name.str,fill=stages.col.str)
    }
  }
}
plot.sc.sub.states.average.and.pop.cor.boxplots=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$development.stage)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$development.stage)
  pop.timepoints=names(pop.meta.list)
  intersect.timepoints=pop.timepoints
  input.timepoints=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  intersect.timepoints=intersect(input.timepoints,intersect.timepoints)
  stages.col.list=get.col.factor(col.factor =intersect.timepoints )
  stages.col.str=stages.col.list$legend.col
  stages.name.str=stages.col.list$legend.str
  col.list=list()
  for(d in 1:length(stages.col.str)){
    col.list[[stages.name.str[d]]]=stages.col.str[d]
  }
  len.sc.timepoints=length(sc.timepoints)
  len.pop.timepoints=length(pop.timepoints)
  for(n in 1:len.sc.timepoints){
    sc.timepoint=sc.timepoints[n]
    temp.sc.meta.df=sc.meta.list[[sc.timepoint]]
    temp.sc.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(sc.rpkm.df[,rownames(temp.sc.meta.df)]))
    temp.sc.meta.df= temp.sc.meta.df[colnames(temp.sc.rpkm.df),]
    temp.sc.meta.list=split(temp.sc.meta.df,f=temp.sc.meta.df$sub.group)
    temp.names=names(temp.sc.meta.list)
    len.temp.names=length(temp.names)
    for(m in 1:len.temp.names){
      temp.name=temp.names[m]
      temp.sc.sub.grp.meta.df=temp.sc.meta.list[[temp.name]]
      temp.sc.sub.grp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(temp.sc.rpkm.df[,rownames(temp.sc.sub.grp.meta.df)]))
      temp.sc.average.rpkm.df=data.frame(average.sc=apply(temp.sc.sub.grp.rpkm.df,1,mean),row.names =rownames(temp.sc.sub.grp.rpkm.df ))
      temp.sc.sub.grp.meta.df= temp.sc.sub.grp.meta.df[colnames(temp.sc.sub.grp.rpkm.df),]
      #temp.sc.sub.grp.rpkm.df=temp.sc.average.rpkm.df
      temp.dist.list=list()
      for(s in 1:len.pop.timepoints){
        pop.timepoint=pop.timepoints[s]
        temp.pop.meta.df=pop.meta.list[[pop.timepoint]]
        #temp.pop.rpkm.df=pop.rpkm.df[,rownames(temp.pop.meta.df)]
        temp.pop.meta.list=split(temp.pop.meta.df,f=temp.pop.meta.df$sub.group)
        pop.temp.names=names(temp.pop.meta.list)
        len.pop.temp.names=length(pop.temp.names)
        for(f in 1:len.pop.temp.names){
          temp.pop.sub.grp.name=pop.temp.names[f]
          temp.pop.sub.grp.meta.df=temp.pop.meta.list[[temp.pop.sub.grp.name]]
          temp.pop.sub.grp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(pop.rpkm.df[,rownames(temp.pop.sub.grp.meta.df)]))
          intersect.genes=intersect(rownames(temp.sc.average.rpkm.df),rownames(temp.pop.sub.grp.rpkm.df))
          intersect.temp.sc.average.rpkm.df=temp.sc.average.rpkm.df[intersect.genes,]
          temp.pop.sub.grp.rpkm.df=temp.pop.sub.grp.rpkm.df[intersect.genes,]
          combn.df=data.frame('average.sc'=intersect.temp.sc.average.rpkm.df,temp.pop.sub.grp.rpkm.df)
          pseudo.value=min(as.numeric(combn.df[combn.df!=0]))/2
          cor.dist.mat=cor(log10(combn.df+pseudo.value),method = 'spearman')
          combn.dist.df=data.frame(as.matrix(dist(t(combn.df),method = 'euclidean')))
          #temp.dist.list[[paste(pop.timepoint,temp.pop.sub.grp.name,collapse = ':')]]=log10(as.numeric(combn.dist.df['average.sc',colnames(temp.pop.sub.grp.rpkm.df)]))
          temp.dist.list[[paste(pop.timepoint,temp.pop.sub.grp.name,collapse = ':')]]=cor.dist.mat['average.sc',colnames(temp.pop.sub.grp.rpkm.df)]
        }
      }
      title.str=gsub('\\.',replacement = ' ',x = convert.to.title.case(temp.name))
      stages.col.list=get.col.factor(col.factor =names(temp.dist.list) )
      stages.col.str=stages.col.list$legend.col
      stages.name.str=stages.col.list$legend.str
      boxplot2(temp.dist.list,las=2,cex=.4,main=title.str,col=stages.col.str)
      abline(h = as.numeric(quantile(as.numeric(unlist(temp.dist.list)),.75)),b = 1,col='lightgray')
      legend('topright',legend=stages.name.str,fill=stages.col.str)
    }
  }
}
plot.sc.and.pop.sub.groups.cor.boxplots=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$development.stage)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$development.stage)
  pop.timepoints=names(pop.meta.list)
  intersect.timepoints=pop.timepoints
  input.timepoints=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  intersect.timepoints=intersect(input.timepoints,intersect.timepoints)
  stages.col.list=get.col.factor(col.factor =intersect.timepoints )
  stages.col.str=stages.col.list$legend.col
  stages.name.str=stages.col.list$legend.str
  col.list=list()
  for(d in 1:length(stages.col.str)){
    col.list[[stages.name.str[d]]]=stages.col.str[d]
  }
  len.sc.timepoints=length(sc.timepoints)
  len.pop.timepoints=length(pop.timepoints)
  for(n in 1:len.sc.timepoints){
    sc.timepoint=sc.timepoints[n]
    temp.sc.meta.df=sc.meta.list[[sc.timepoint]]
    temp.sc.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(sc.rpkm.df[,rownames(temp.sc.meta.df)]))
    temp.sc.meta.df= temp.sc.meta.df[colnames(temp.sc.rpkm.df),]
    #temp.sc.meta.list=split(temp.sc.meta.df,f=temp.sc.meta.df$sub.group)
    temp.names=colnames(temp.sc.rpkm.df)
    len.temp.names=length(temp.names)
    for(m in 1:len.temp.names){
      temp.dist.list=list()
      temp.name=temp.names[m]
      #temp.sc.sub.grp.meta.df=temp.sc.meta.list[[temp.name]]
      temp.rpkm.vec=as.numeric(as.character(temp.sc.rpkm.df[,temp.name]))
      sc.sample.rpkm.df=data.frame(sc.rpkm=temp.rpkm.vec,row.names =rownames(temp.sc.rpkm.df))
      for(s in 1:len.pop.timepoints){
        pop.timepoint=pop.timepoints[s]
        temp.pop.meta.df=pop.meta.list[[pop.timepoint]]
        #temp.pop.rpkm.df=pop.rpkm.df[,rownames(temp.pop.meta.df)]
        temp.pop.meta.list=split(temp.pop.meta.df,f=temp.pop.meta.df$sub.group)
        pop.temp.names=names(temp.pop.meta.list)
        len.pop.temp.names=length(pop.temp.names)
        for(f in 1:len.pop.temp.names){
          temp.pop.sub.grp.name=pop.temp.names[f]
          temp.pop.sub.grp.meta.df=temp.pop.meta.list[[temp.pop.sub.grp.name]]
          temp.pop.sub.grp.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(pop.rpkm.df[,rownames(temp.pop.sub.grp.meta.df)]))
          intersect.genes=intersect(rownames(sc.sample.rpkm.df),rownames(temp.pop.sub.grp.rpkm.df))
          intersect.temp.sc.average.rpkm.df=sc.sample.rpkm.df[intersect.genes,]
          temp.pop.sub.grp.rpkm.df=temp.pop.sub.grp.rpkm.df[intersect.genes,]
          combn.df=data.frame('sample.sc'=intersect.temp.sc.average.rpkm.df,temp.pop.sub.grp.rpkm.df)
          pseudo.value=min(as.numeric(combn.df[combn.df!=0]))/2
          cor.dist.mat=cor(log10(combn.df+pseudo.value),method = 'spearman')
          #combn.dist.df=data.frame(as.matrix(dist(t(combn.df),method = 'euclidean')))
          #temp.dist.list[[paste(pop.timepoint,temp.pop.sub.grp.name,collapse = ':')]]=log10(as.numeric(combn.dist.df['average.sc',colnames(temp.pop.sub.grp.rpkm.df)]))
          temp.dist.list[[paste(pop.timepoint,temp.pop.sub.grp.name,collapse = ':')]]=cor.dist.mat['sample.sc',colnames(temp.pop.sub.grp.rpkm.df)]
        }
      }
      title.str=gsub('\\.',replacement = ' ',x = convert.to.title.case(temp.name))
      stages.col.list=get.col.factor(col.factor =names(temp.dist.list) )
      stages.col.str=stages.col.list$legend.col
      stages.name.str=stages.col.list$legend.str
      boxplot2(temp.dist.list,las=2,cex=.4,main=title.str,col=stages.col.str)
      abline(h = as.numeric(quantile(as.numeric(unlist(temp.dist.list)),.75)),b = 1,col='lightgray')
      legend('topright',legend=stages.name.str,fill=stages.col.str)
    }
  }
}
plot.sc.sub.pop.and.pop.sub.groups.cor.boxplots=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$markers.cluster.groups)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$sub.group)
  pop.timepoints=names(pop.meta.list)
  intersect.timepoints=pop.timepoints
  #input.timepoints=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  #intersect.timepoints=intersect(input.timepoints,intersect.timepoints)
  stages.col.list=get.col.factor(col.factor =intersect.timepoints )
  stages.col.str=stages.col.list$legend.col
  stages.name.str=stages.col.list$legend.str
  col.list=list()
  for(d in 1:length(stages.col.str)){
    col.list[[stages.name.str[d]]]=stages.col.str[d]
  }
  len.sc.timepoints=length(sc.timepoints)
  len.pop.timepoints=length(pop.timepoints)
  for(n in 1:len.sc.timepoints){
    sc.timepoint=sc.timepoints[n]
    temp.sc.meta.df=sc.meta.list[[sc.timepoint]]
    temp.sc.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(sc.rpkm.df[,rownames(temp.sc.meta.df)]))
    temp.sc.meta.df= temp.sc.meta.df[colnames(temp.sc.rpkm.df),]
    #temp.sc.meta.list=split(temp.sc.meta.df,f=temp.sc.meta.df$sub.group)
    temp.names=colnames(temp.sc.rpkm.df)
    len.temp.names=length(temp.names)
    temp.dist.list=list()
    for(s in 1:len.pop.timepoints){
      pop.timepoint=pop.timepoints[s]
      temp.pop.meta.df=pop.meta.list[[pop.timepoint]]
      temp.pop.rpkm.df=filter.none.expressed.genes(input.data = pop.rpkm.df[,rownames(temp.pop.meta.df)])
      cor.score.vec=get.cor.between.two.df.cor(first.df = temp.sc.rpkm.df,second.df =temp.pop.rpkm.df,in.method = 'spearman' )
      temp.dist.list[[pop.timepoint]]=cor.score.vec
      #temp.pop.meta.list=split(temp.pop.meta.df,f=temp.pop.meta.df$sub.group)
    }
    title.str=gsub('\\.',replacement = ' ',x = convert.to.title.case(sc.timepoint))
    stages.col.list=get.col.factor(col.factor =names(temp.dist.list) )
    stages.col.str=stages.col.list$legend.col
    stages.name.str=stages.col.list$legend.str
    unlabel.str=rep('',times = length(temp.dist.list))
    boxplot(temp.dist.list,las=2,cex=.8,main=title.str,col=stages.col.str,pch=19,names = unlabel.str)
    #boxplot(temp.dist.list,las=2,cex=.8,main=title.str,col=stages.col.str,pch=19)
    #abline(h = as.numeric(quantile(as.numeric(unlist(temp.dist.list)),.75)),b = 1,col='lightgray')
    #plot.new()
    #legend('topright',legend=stages.name.str,fill=stages.col.str,box.lty = 0,cex = 1.5,border = F)
  }
}
plot.sc.sub.pop.and.pop.sub.groups.cor.boxplots.based.on.idc.markers=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df,markers.vec){
  sc.rpkm.df=sc.rpkm.df[intersect(rownames(sc.rpkm.df),markers.vec),]
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$markers.cluster.groups)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$sub.group)
  pop.timepoints=names(pop.meta.list)
  intersect.timepoints=pop.timepoints
  #input.timepoints=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  #intersect.timepoints=intersect(input.timepoints,intersect.timepoints)
  stages.col.list=get.col.factor(col.factor =intersect.timepoints )
  stages.col.str=stages.col.list$legend.col
  stages.name.str=stages.col.list$legend.str
  col.list=list()
  for(d in 1:length(stages.col.str)){
    col.list[[stages.name.str[d]]]=stages.col.str[d]
  }
  len.sc.timepoints=length(sc.timepoints)
  len.pop.timepoints=length(pop.timepoints)
  for(n in 1:len.sc.timepoints){
    sc.timepoint=sc.timepoints[n]
    temp.sc.meta.df=sc.meta.list[[sc.timepoint]]
    temp.sc.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(sc.rpkm.df[,rownames(temp.sc.meta.df)]))
    temp.sc.meta.df= temp.sc.meta.df[colnames(temp.sc.rpkm.df),]
    #temp.sc.meta.list=split(temp.sc.meta.df,f=temp.sc.meta.df$sub.group)
    temp.names=colnames(temp.sc.rpkm.df)
    len.temp.names=length(temp.names)
    temp.dist.list=list()
    for(s in 1:len.pop.timepoints){
      pop.timepoint=pop.timepoints[s]
      temp.pop.meta.df=pop.meta.list[[pop.timepoint]]
      temp.pop.rpkm.df=filter.none.expressed.genes(input.data = pop.rpkm.df[,rownames(temp.pop.meta.df)])
      cor.score.vec=get.cor.between.two.df.cor(first.df = temp.sc.rpkm.df,second.df =temp.pop.rpkm.df,in.method = 'spearman' )
      temp.dist.list[[pop.timepoint]]=cor.score.vec
      #temp.pop.meta.list=split(temp.pop.meta.df,f=temp.pop.meta.df$sub.group)
    }
    title.str=gsub('\\.',replacement = ': ',x = toupper(sc.timepoint))
    stages.col.list=get.col.factor(col.factor =names(temp.dist.list) )
    stages.col.str=stages.col.list$legend.col
    stages.name.str=stages.col.list$legend.str
    unlabel.str=rep('',times = length(temp.dist.list))
    #boxplot(temp.dist.list,las=2,cex=.8,main=title.str,col=stages.col.str,pch=19,names = unlabel.str)
    boxplot(temp.dist.list,las=2,cex=.8,main=title.str,col=stages.col.str,pch=19)
    #abline(h = as.numeric(quantile(as.numeric(unlist(temp.dist.list)),.75)),b = 1,col='lightgray')
    #plot.new()
    #legend('topright',legend=stages.name.str,fill=stages.col.str,box.lty = 0,cex = 1.5,border = F)
  }
}
plot.sc.sub.pop.average.and.pop.sub.groups.cor.boxplots=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$markers.cluster.groups)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$sub.group)
  pop.timepoints=names(pop.meta.list)
  intersect.timepoints=pop.timepoints
  #input.timepoints=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  #intersect.timepoints=intersect(input.timepoints,intersect.timepoints)
  stages.col.list=get.col.factor(col.factor =intersect.timepoints )
  stages.col.str=stages.col.list$legend.col
  stages.name.str=stages.col.list$legend.str
  col.list=list()
  for(d in 1:length(stages.col.str)){
    col.list[[stages.name.str[d]]]=stages.col.str[d]
  }
  len.sc.timepoints=length(sc.timepoints)
  len.pop.timepoints=length(pop.timepoints)
  average.sc.pop.cor.mat=matrix(nrow =len.sc.timepoints,ncol =  len.pop.timepoints)
  for(n in 1:len.sc.timepoints){
    sc.timepoint=sc.timepoints[n]
    temp.sc.meta.df=sc.meta.list[[sc.timepoint]]
    temp.sc.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(sc.rpkm.df[,rownames(temp.sc.meta.df)]))
    temp.sc.meta.df= temp.sc.meta.df[colnames(temp.sc.rpkm.df),]
    #temp.sc.meta.list=split(temp.sc.meta.df,f=temp.sc.meta.df$sub.group)
    temp.names=colnames(temp.sc.rpkm.df)
    len.temp.names=length(temp.names)
    temp.dist.list=list()
    temp.cor.mean=c()
    for(s in 1:len.pop.timepoints){
      pop.timepoint=pop.timepoints[s]
      temp.pop.meta.df=pop.meta.list[[pop.timepoint]]
      temp.pop.rpkm.df=filter.none.expressed.genes(input.data = pop.rpkm.df[,rownames(temp.pop.meta.df)])
      intersect.genes=intersect(rownames(temp.sc.rpkm.df),rownames(temp.pop.rpkm.df))
      temp.pop.rpkm.df=temp.pop.rpkm.df[intersect.genes,]
      temp.sc.rpkm.df=temp.sc.rpkm.df[intersect.genes,]
      temp.average.rpkm.vec=as.numeric(apply(temp.sc.rpkm.df,1,mean))
      #cor.score.vec=get.cor.between.two.df.cor(first.df = temp.sc.rpkm.df,second.df =temp.pop.rpkm.df,in.method = 'spearman' )
      cor.score.vec=as.numeric(unlist(apply(temp.pop.rpkm.df,2,function(temp.rpkm){
        temp.rpkm.vec=as.numeric(temp.rpkm)
        cor.coeff=cor(x =temp.average.rpkm.vec,y = temp.rpkm.vec,method = 'spearman' )
        return(as.numeric(cor.coeff))
      })))
      temp.dist.list[[pop.timepoint]]=cor.score.vec
      temp.cor.mean=append(x = temp.cor.mean,values = mean(cor.score.vec),after = length(temp.cor.mean))
      #temp.pop.meta.list=split(temp.pop.meta.df,f=temp.pop.meta.df$sub.group)
    }
    average.sc.pop.cor.mat[n,]=as.numeric(temp.cor.mean)
    title.str=gsub('\\.',replacement = ' ',x = convert.to.title.case(sc.timepoint))
    stages.col.list=get.col.factor(col.factor =names(temp.dist.list) )
    stages.col.str=stages.col.list$legend.col
    stages.name.str=stages.col.list$legend.str
    unlabel.str=rep('',times = length(temp.dist.list))
    #boxplot(temp.dist.list,las=2,cex=.8,main=title.str,col=stages.col.str,pch=19,names = unlabel.str)
    boxplot(temp.dist.list,las=2,cex=.8,main=title.str,col=stages.col.str,pch=19,border = NA)
    #abline(h = as.numeric(quantile(as.numeric(unlist(temp.dist.list)),.75)),b = 1,col='lightgray')
    #plot.new()
    #legend('topright',legend=stages.name.str,fill=stages.col.str,box.lty = 0,cex = 1.5,border = F)
  }
}
plot.sc.sub.pop.average.and.pop.sub.groups.cor.heatmap=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$markers.cluster.groups)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$sub.group)
  pop.timepoints=names(pop.meta.list)
  intersect.timepoints=pop.timepoints
  #input.timepoints=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  #intersect.timepoints=intersect(input.timepoints,intersect.timepoints)
  stages.col.list=get.col.factor(col.factor =intersect.timepoints )
  stages.col.str=stages.col.list$legend.col
  stages.name.str=stages.col.list$legend.str
  col.list=list()
  for(d in 1:length(stages.col.str)){
    col.list[[stages.name.str[d]]]=stages.col.str[d]
  }
  len.sc.timepoints=length(sc.timepoints)
  len.pop.timepoints=length(pop.timepoints)
  average.sc.pop.cor.mat=matrix(nrow =len.sc.timepoints,ncol =  len.pop.timepoints)
  for(n in 1:len.sc.timepoints){
    sc.timepoint=sc.timepoints[n]
    temp.sc.meta.df=sc.meta.list[[sc.timepoint]]
    temp.sc.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(sc.rpkm.df[,rownames(temp.sc.meta.df)]))
    temp.sc.meta.df= temp.sc.meta.df[colnames(temp.sc.rpkm.df),]
    #temp.sc.meta.list=split(temp.sc.meta.df,f=temp.sc.meta.df$sub.group)
    temp.names=colnames(temp.sc.rpkm.df)
    len.temp.names=length(temp.names)
    temp.dist.list=list()
    temp.cor.mean=c()
    for(s in 1:len.pop.timepoints){
      pop.timepoint=pop.timepoints[s]
      temp.pop.meta.df=pop.meta.list[[pop.timepoint]]
      temp.pop.rpkm.df=filter.none.expressed.genes(input.data = pop.rpkm.df[,rownames(temp.pop.meta.df)])
      intersect.genes=intersect(rownames(temp.sc.rpkm.df),rownames(temp.pop.rpkm.df))
      temp.pop.rpkm.df=temp.pop.rpkm.df[intersect.genes,]
      temp.sc.rpkm.df=temp.sc.rpkm.df[intersect.genes,]
      temp.average.rpkm.vec=as.numeric(apply(temp.sc.rpkm.df,1,mean))
      #cor.score.vec=get.cor.between.two.df.cor(first.df = temp.sc.rpkm.df,second.df =temp.pop.rpkm.df,in.method = 'spearman' )
      cor.score.vec=as.numeric(unlist(apply(temp.pop.rpkm.df,2,function(temp.rpkm){
        temp.rpkm.vec=as.numeric(temp.rpkm)
        cor.coeff=cor(x =temp.average.rpkm.vec,y = temp.rpkm.vec,method = 'spearman' )
        return(as.numeric(cor.coeff))
      })))
      temp.dist.list[[pop.timepoint]]=cor.score.vec
      temp.cor.mean=append(x = temp.cor.mean,values = mean(cor.score.vec),after = length(temp.cor.mean))
      #temp.pop.meta.list=split(temp.pop.meta.df,f=temp.pop.meta.df$sub.group)
    }
    average.sc.pop.cor.mat[n,]=as.numeric(temp.cor.mean)
    title.str=gsub('\\.',replacement = ' ',x = convert.to.title.case(sc.timepoint))
    stages.col.list=get.col.factor(col.factor =names(temp.dist.list) )
    stages.col.str=stages.col.list$legend.col
    stages.name.str=stages.col.list$legend.str
    unlabel.str=rep('',times = length(temp.dist.list))
  }
  rownames(average.sc.pop.cor.mat)=sc.timepoints
  colnames(average.sc.pop.cor.mat)= pop.timepoints
  sorted.average.sc.pop.cor.mat=average.sc.pop.cor.mat[names(sort(sub.populations.grp.order.vec,decreasing = T)),names(sort(bulk.sub.grps.order.vec,decreasing = T))]
  #colnames(average.sc.pop.cor.mat)=gsub(pattern = 'SP',replacement = 'SB',x = pop.timepoints)
  pheatmap(mat = average.sc.pop.cor.mat,main = '',color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),cellwidth = 30,cellheight = 30,fontsize_row = 10,fontsize_col = 10,border_color = NA,cluster_rows = F,cluster_cols = F)
  pheatmap(mat = sorted.average.sc.pop.cor.mat,main = '',color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),fontsize_row = 20,fontsize_col = 20,cellwidth = 40,cellheight = 40,border_color = NA,cluster_rows = F,cluster_cols = F)
}
plot.sc.sub.pop.average.and.pop.sub.groups.cor.using.idc.markers.boxplots=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df,markers.vec){
  sc.rpkm.df=sc.rpkm.df[intersect(rownames(sc.rpkm.df),markers.vec),]
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$markers.cluster.groups)
  sc.timepoints=names(sc.meta.list)
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$sub.group)
  pop.timepoints=names(pop.meta.list)
  intersect.timepoints=pop.timepoints
  #input.timepoints=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  #intersect.timepoints=intersect(input.timepoints,intersect.timepoints)
  stages.col.list=get.col.factor(col.factor =intersect.timepoints )
  stages.col.str=stages.col.list$legend.col
  stages.name.str=stages.col.list$legend.str
  col.list=list()
  for(d in 1:length(stages.col.str)){
    col.list[[stages.name.str[d]]]=stages.col.str[d]
  }
  len.sc.timepoints=length(sc.timepoints)
  len.pop.timepoints=length(pop.timepoints)
  for(n in 1:len.sc.timepoints){
    sc.timepoint=sc.timepoints[n]
    temp.sc.meta.df=sc.meta.list[[sc.timepoint]]
    temp.sc.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(sc.rpkm.df[,rownames(temp.sc.meta.df)]))
    show(dim(temp.sc.rpkm.df))
    temp.sc.meta.df= temp.sc.meta.df[colnames(temp.sc.rpkm.df),]
    #temp.sc.meta.list=split(temp.sc.meta.df,f=temp.sc.meta.df$sub.group)
    temp.names=colnames(temp.sc.rpkm.df)
    len.temp.names=length(temp.names)
    temp.dist.list=list()
    for(s in 1:len.pop.timepoints){
      pop.timepoint=pop.timepoints[s]
      temp.pop.meta.df=pop.meta.list[[pop.timepoint]]
      temp.pop.rpkm.df=filter.none.expressed.genes(input.data = pop.rpkm.df[,rownames(temp.pop.meta.df)])
      intersect.genes=intersect(rownames(temp.sc.rpkm.df),rownames(temp.pop.rpkm.df))
      temp.pop.rpkm.df=temp.pop.rpkm.df[intersect.genes,]
      temp.sc.rpkm.df=temp.sc.rpkm.df[intersect.genes,]
      temp.average.rpkm.vec=as.numeric(apply(temp.sc.rpkm.df,1,mean))
      #cor.score.vec=get.cor.between.two.df.cor(first.df = temp.sc.rpkm.df,second.df =temp.pop.rpkm.df,in.method = 'spearman' )
      cor.score.vec=as.numeric(unlist(apply(temp.pop.rpkm.df,2,function(temp.rpkm){
        temp.rpkm.vec=as.numeric(temp.rpkm)
        cor.coeff=cor(x =temp.average.rpkm.vec,y = temp.rpkm.vec,method = 'spearman' )
        return(as.numeric(cor.coeff))
      })))
      temp.dist.list[[pop.timepoint]]=cor.score.vec
      #temp.pop.meta.list=split(temp.pop.meta.df,f=temp.pop.meta.df$sub.group)
    }
    title.str=gsub('\\.',replacement = ': ',x = toupper(sc.timepoint))
    stages.col.list=get.col.factor(col.factor =names(temp.dist.list) )
    stages.col.str=stages.col.list$legend.col
    stages.name.str=stages.col.list$legend.str
    #unlabel.str=rep('',times = length(temp.dist.list))
    #boxplot(temp.dist.list,las=2,cex=.8,main=title.str,col=stages.col.str,pch=19,names = unlabel.str)
    boxplot(temp.dist.list,las=2,cex=.8,main=title.str,col=stages.col.str,pch=19)
    #abline(h = as.numeric(quantile(as.numeric(unlist(temp.dist.list)),.75)),b = 1,col='lightgray')
    #plot.new()
    #legend('topright',legend=stages.name.str,fill=stages.col.str,box.lty = 0,cex = 1.5,border = F)
  }
}
temp.plot.pca=function(in.rpkm.df,trans=F,meta.data,title.str='Test',first.pc='PC1',second.pc='PC2',
                       log.pca.scores=F,plot.pairs=T,three.dim.scatter=F,log.rpkm=T){
  rpkm.df=in.rpkm.df
  meta.data=meta.data[as.character(colnames(rpkm.df)),]
  meta.data=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.data)
  rpkm.df=if(trans){t(rpkm.df)} else {rpkm.df}
  pseudo.value=min(rpkm.df[rpkm.df>0])/2
  pseudo.value=1.0
  if(log.rpkm){
    rpkm.df=log2(rpkm.df+pseudo.value)
  }
  no.samples=dim(rpkm.df)[2]
  no.genes=dim(rpkm.df)[1]
  out.results.list=list()
  samples.pca=prcomp(t(rpkm.df),retx=T,center=T,scale.=T)
  samples.pca.summary=summary(samples.pca)$importance
  out.results.list[['pca.summary']]=samples.pca.summary
  #barplot(100*samples.pca.summary[2,],main='Components explained variance',ylab='%')
  #abline(h=100/length(colnames(samples.pca.summary)),b=1)
  first.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,first.pc]))*100,2)
  second.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,second.pc]))*100,2)
  out.results.list[['pca.obj']]= samples.pca
  samples.pca.scores=samples.pca$x
  pc.scores.meta.df=data.frame(samples.pca.scores,meta.data)
  x.points=as.numeric(as.character(pc.scores.meta.df[,first.pc]))
  y.points=as.numeric(as.character(pc.scores.meta.df[,second.pc]))
  if(log.pca.scores){
    x.points=log2(x.points+1)
    y.points=log2(y.points+1)
  }
  pca.col.factor=as.character(pc.scores.meta.df$development.stage)
  #pca.col.factor=as.character(pc.scores.meta.df$markers.cluster.groups)
  #pca.col.factor=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  #pca.col.factor=as.character(pc.scores.meta.df$group)
  #pca.col.factor=as.character(pc.scores.meta.df$batch)
  #pca.col.factor=as.character(pc.scores.meta.df[,'blood.group'])
  pca.col.list=get.col.factor(col.factor=pca.col.factor)
  pca.col.str =abbrev.std.stages.col.vec[pca.col.factor]
  #pca.col.str=pca.col.list[['col.str']]
  first.pc.explained.var=paste('(Explained var: ',paste(first.pc.explained.var,'%)',sep=''),sep='')
  second.pc.explained.var=paste('(Explained var: ',paste(second.pc.explained.var,'%)',sep=''),sep='')
  first.pc.lab=paste(first.pc,first.pc.explained.var,'')
  second.pc.lab=paste(second.pc,second.pc.explained.var,'')
  x.lim=c(min(x.points,y.points),max(x.points,y.points))
  y.lim=c(min(x.points,y.points),max(x.points,y.points))
  temp.pc.scores.meta.df=pc.scores.meta.df
  points.size=as.numeric(lapply(strsplit(x=as.character(temp.pc.scores.meta.df[,'no.of.detected.genes']),split=':'),function(counts.category){
    counts.category=as.character(counts.category)
    return(counts.category[1])
  }))
  point.size.fact=points.size/100
  point.size.fact=point.size.fact*.2
  plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,ylab=second.pc.lab,main=title.str,
       pch=19,cex=point.size.fact)
  #abline(h=0,b=1,col = "lightgray")
  #abline(v=0,b=1,col = "lightgray")
  #plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,ylab=second.pc.lab,main=title.str,pch=19,cex=point.size.fact)
  #abline(h=0,b=1,col = "lightgray")
  #abline(v=0,b=1,col = "lightgray")
  #legend('bottomright',legend=pca.col.list[['legend.str']],fill=pca.col.list[['legend.cols']],cex=1.5)
  plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,
       ylab=second.pc.lab,main=title.str,pch=19,cex=2.0)
  #legend(x = 3,y = -5,legend=pca.col.list[['legend.str']],fill=pca.col.list[['legend.cols']],cex=1)
  #abline(h=0,b=1,col = "lightgray")
  #abline(v=0,b=1,col = "lightgray")
  #plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,ylab=second.pc.lab,main=title.str,pch=19,cex=1.5)
  #abline(h=0,b=1,col = "lightgray")
  #abline(v=0,b=1,col = "lightgray")
  #legend('topright',legend=pca.col.list[['legend.str']],fill=pca.col.list[['legend.cols']])
  #smoothScatter(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,ylab=second.pc.lab,main=title.str,pch=19,cex=1.5)
  #abline(h=0,b=1,col = "lightgray")
  #abline(v=0,b=1,col = "lightgray")
  #plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,ylab=second.pc.lab,main=title.str,pch=19,cex=1.5)
  #abline(h=0,b=1,col = "lightgray")
  #abline(v=0,b=1,col = "lightgray")
  total.no.pc=dim(samples.pca.scores)[2]
  pca.pair.mat=samples.pca.scores
  if(total.no.pc>=6){
    pca.pair.mat=samples.pca.scores[,1:6]
  }
  if(plot.pairs){
    pairs(pca.pair.mat,cex=1.5,pch=19,col=pca.col.str,main=title.str)
    #legend('topright',legend=pca.col.list[['legend.str']],fill=pca.col.list[['legend.cols']])
  }
}
plot.gene.clusters=function(temp.rpkm.df,temp.meta.df){
  timepoints=as.character(temp.meta.df$development.stage)
  timepoints.col.list=get.col.factor(timepoints)
  timepoints.col=timepoints.col.list$col.str
  spearman.cor=cor(temp.rpkm.df,method = 'spearman')
  spearman.cor[is.na(spearman.cor)]=0
  spearman.cor.dist=as.dist(1-spearman.cor)
  dist.mat=dist(x=log10(t(temp.rpkm.df)+1))
  dist.hclust=hclust(spearman.cor.dist, method = "complete")
  dist.dend=as.dendrogram(dist.hclust)
  plot(dist.hclust,cex=.3)
  test.dist.heat.mat=as.matrix(dist.mat)
  col.factor.list=timepoints.col.list
  col.str=col.factor.list[['col.str']]
  legend.str=col.factor.list[['legend.str']]
  legend.cols=col.factor.list[['legend.cols']]
  heatmap.2(spearman.cor, trace="none",margin=c(6, 6),srtCol=45,srtRow=45,cexRow=0.5,cexCol=0.8,dendrogram='col',RowSideColors=col.str,ColSideColors=col.str,Colv = dist.dend)
  legend('topright',legend=legend.str,fill=legend.cols)
  temp.rpkm.mat=as.matrix(temp.rpkm.df)
  samples.spearman.cor=cor(temp.rpkm.mat,method = 'spearman')
  samples.spearman.cor[is.na(samples.spearman.cor)]=0
  samples.spearman.cor=as.dist(1-samples.spearman.cor)
  samples.spearman.cor.hclust=hclust(samples.spearman.cor,method = 'average')
  samples.dendrogram=as.dendrogram(samples.spearman.cor.hclust)
  heatmap.2(log10(temp.rpkm.mat+1),ColSideColors =col.str,trace='none' ,col = greenred(5),Colv =samples.dendrogram)
  rect.hclust(tree = dist.hclust,k = 4)
  heatmap.2(temp.rpkm.mat,ColSideColors =col.str,trace='none' ,col = greenred(5),Colv =samples.dendrogram)
}
select.the.right.genes.for.pca=function(rpkm.df,in.trans=F,in.meta.data,in.title.str,in.var.cut.off,in.first.pc='PC1',in.second.pc='PC2',log.pca.scores=F,plot.pairs=T,three.dim.scatter=F){
  rpkm.df=filter.non.variable.rows(rpkm.df,cut.off=1)
  temp.pca=prcomp(x=t(rpkm.df),retx=T,center=T,scale.=T)
  important.components=select.gene.with.sig.component.weights(pca.obj=temp.pca)
  genes.list=important.components[['pc.loadings.list']]
  important.pc.cut.off=important.components[['no.pc.cut.off']]
  genes=unique(unlist(important.components))
  out.list=list()
  for(n in 1:5){
    loop.no=paste('loop_number:',n,sep='')
    show(length(genes))
    temp.rpkm.df=rpkm.df[genes,]
    out.list[[loop.no]]=genes
    temp.title.str=paste(in.title.str,loop.no,sep=': ')
    temp.out=plot.pca(in.rpkm.df=temp.rpkm.df,trans=in.trans,meta.data=in.meta.data,title.str=temp.title.str,plot.pairs=T,var.cut.off=in.var.cut.off)
    temp.pca=temp.out[['pca.obj']]
    temp.component.genes=temp.out[['component.genes']]
    temp.genes=unique(c(unlist(temp.component.genes$PC1)),unlist(temp.component.genes$PC2))
    genes=temp.genes
  }
  return(out.list)
}
select.gene.with.sig.component.weights=function(pca.obj){
  out.list=list()
  samples.pca=pca.obj
  pca.importance.df=data.frame(t(summary(samples.pca)$importance))
  important.pcs.cut.off =1/dim(pca.importance.df)[1]
  important.pcs.df=pca.importance.df[which(as.character(as.numeric(pca.importance.df$'Proportion.of.Variance'))>important.pcs.cut.off),]
  no.of.filt.pcs=dim(important.pcs.df)[1]
  pca.rotation=samples.pca$rotation
  pca.rotation.cut.off=1/nrow(pca.rotation)
  pca.rotation.squared=pca.rotation^2
  #pca.rotation.squared=abs(pca.rotation)
  genes.with.high.weight.test=pca.rotation.squared>pca.rotation.cut.off
  pc.names=colnames(pca.rotation.squared)
  #pc.high.weight.genes=list(names<-pc.names)
  pc.sig.gene.loadings.list=list()
  for(m in 1:no.of.filt.pcs){
    pc=pc.names[m]
    temp.df=data.frame(genes=rownames(pca.rotation.squared),row.names=rownames(pca.rotation.squared),weight.test=genes.with.high.weight.test[,pc])
    temp.df=temp.df[which(temp.df[,'weight.test']==T),]
    pc.sig.gene.loadings.list[[pc]]=rownames(temp.df)
  }
  out.list[['pc.loadings.list']]=pc.sig.gene.loadings.list
  out.list[['no.pc.cut.off']]=no.of.filt.pcs 
  return(out.list)
}
filter.genes.not.expressed.in.all.samples=function(df){
  df=data.frame(df)
  show(dim(df))
  no.samples=dim(df)[2]
  df$gene.counts=as.numeric(apply(df,1,nnzero))
  out.df=subset(df,gene.counts==no.samples)
  out.df=subset(out.df,select=-gene.counts)
  return(out.df)
}
filter.genes.not.expressed.in.a.number.of.samples=function(df,sample.count=2){
  df=data.frame(df)
  df$no.samples=as.numeric(as.character(apply(df,1,nnzero)))
  out.df=subset(df,no.samples>=sample.count)
  out.df=subset(out.df,select = -no.samples)
  return(out.df)
}
encode.stage.to.markers=function(plos.one.plasmodb.df,plos.df){
  genes=rownames(plos.one.plasmodb.df)
  category=c()
  for (m in 1: length(genes)){
    gene =genes[m]
    old.id=as.character(plos.one.plasmodb.df[gene,'X.Input.ID.'])
    #show(old.id)
    stage=unique(as.character(plos.df[which(plos.df$ID==old.id),'stage']))
    if(length(stage)!=1){
      category=append(category,'none',length(category))
    }
    else{
      category=append(category,stage,length(category))
    }
  }
  plos.one.plasmodb.df$category=category
  return(plos.one.plasmodb.df)
}
run.scde.between.devt.stages=function(counts.df,meta.df,correct.batch.effect=F,parental.development=F){
  all.samples=colnames(counts.df)
  counts.df=counts.df[,all.samples]
  meta.df=meta.df[all.samples,]
  meta.list=split(meta.df,f = meta.df$development.stage)
  if(parental.development){
    meta.list=split(meta.df,f = meta.df$parent.development.stage)
  }
  stages.names=names(meta.list)
  show(table(stages.names))
  #stages.names=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  #stages.names=c("ring","trophozoite" ,"schizont" )
  stages.combn=combn(stages.names,m = 2)
  no.combn=dim(stages.combn)
  out.list=list()
  #for(n in 1:2){
  for(n in 1:no.combn[2]){
    temp.out.list=list()
    temp.combn=stages.combn[,n]
    first.stage.name=temp.combn[1]
    sec.stage.name=temp.combn[2]
    temp.name=paste(first.stage.name,sec.stage.name,sep = '_vs_')
    first.meta.df=meta.list[[first.stage.name]]
    sec.meta.df=meta.list[[sec.stage.name]]
    groups=factor(c(rep(x =first.stage.name ,times = dim(first.meta.df)[1]),rep(x =sec.stage.name ,times = dim(sec.meta.df)[1])))
    batch=factor(c(first.meta.df[,'batch'],sec.meta.df[,'batch']))
    first.samples=rownames(first.meta.df)
    sec.samples=rownames(sec.meta.df)
    temp.counts.df=counts.df[,c(first.samples,sec.samples)]
    temp.counts.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = temp.counts.df))
    show(dim(temp.counts.df))
    if(is.null(dim(temp.counts.df))){
      next
    }
    temp.meta.df=meta.df[colnames(temp.counts.df),]
    groups=factor(as.character(temp.meta.df$development.stage))
    if(parental.development){
      groups=factor(as.character(temp.meta.df$parent.development.stage))
    }
    names(groups)=as.character(colnames(temp.counts.df))
    no.cores=2
    o.ifm <- try(scde.error.models(counts =temp.counts.df,groups=groups,n.cores=no.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1));
    if(class(o.ifm)=='try-error'){
      next
    }
    valid.cells <- o.ifm$corr.a >0;
    o.ifm <- o.ifm[valid.cells,];
    temp.counts.df=counts.df[,rownames(o.ifm)]
    temp.meta.df=meta.df[rownames(o.ifm),]
    temp.sample.development.stage=as.character(temp.meta.df$development.stage)
    if(parental.development){
      temp.sample.development.stage=as.character(temp.meta.df$parent.development.stage)
    }
    if (length(unique(temp.sample.development.stage))<2){
      next
    }
    no.samples.per.group=as.numeric(as.character(table(temp.sample.development.stage)))
    groups=factor(as.character(temp.meta.df$development.stage))
    if(parental.development){
      groups=factor(as.character(temp.meta.df$parent.development.stage))
    }
    names(groups)=as.character(rownames(temp.meta.df))
    batch=factor(as.character(temp.meta.df$batch))
    temp.out.list[['groups']]=groups
    temp.out.list[['batch']]=batch
    o.prior <- try(scde.expression.prior(models=o.ifm,counts=temp.counts.df,length.out=400,show.plot=F))
    if(class(o.prior)=='try-error'){
      next
    }
    no.batches=length(unique(batch))
    ediff=data.frame()
    if(correct.batch.effect){
      ediff <- scde.expression.difference(models = o.ifm,counts = temp.counts.df,prior =o.prior,groups=groups,n.randomizations=100,n.cores=no.cores,verbose=1,batch = batch,return.posteriors = T)
    }
    else{
      ediff <- scde.expression.difference(models = o.ifm,counts = temp.counts.df,prior =o.prior,groups=groups,n.randomizations=100,n.cores=no.cores,verbose=1,return.posteriors = T)
    }
    ediff.res=ediff$results
    ediff.res$p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$Z))))
    ediff.res$adj.p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$cZ))))
    temp.out.list[['o.ifm']]=o.ifm
    temp.out.list[['o.prior']]=o.prior
    temp.out.list[['ediff']]=ediff
    temp.out.list[['results']]=ediff.res
    temp.out.list[['count']]=temp.counts.df
    temp.out.list[['meta']]=temp.meta.df
    out.list[[temp.name]]=temp.out.list
  }
  return(out.list)
}
run.scde.between.devt.stages.in.order=function(counts.df,meta.df,correct.batch.effect=F,parental.development=F){
  all.samples=colnames(counts.df)
  meta.df=meta.df[all.samples,]
  meta.list=split(meta.df,f = meta.df$development.stage)
  if(parental.development){
    meta.list=split(meta.df,f = meta.df$parent.development.stage)
  }
  stages.names=names(meta.list)
  ordered.stages.names=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  if(parental.development){
    ordered.stages.names=c("R" , "T", "S")
  }
  stages.names=intersect(ordered.stages.names,stages.names)
  len.stages=length(stages.names)
  out.list=list()
  for(n in 1:len.stages){
    if(n==len.stages){
      next
    }
    first.stage.name=stages.names[n]
    sec.stage.name=stages.names[n+1]
    temp.name=paste(first.stage.name,sec.stage.name,sep = '_vs_')
    first.meta.df=meta.list[[first.stage.name]]
    sec.meta.df=meta.list[[sec.stage.name]]
    groups=factor(c(rep(x =first.stage.name ,times = dim(first.meta.df)[1]),rep(x =sec.stage.name ,times = dim(sec.meta.df)[1])))
    batch=factor(c(first.meta.df[,'batch'],sec.meta.df[,'batch']))
    first.samples=rownames(first.meta.df)
    sec.samples=rownames(sec.meta.df)
    temp.counts.df=counts.df[,c(first.samples,sec.samples)]
    temp.counts.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = temp.counts.df))
    if(is.null(dim(temp.counts.df))){
      next
    }
    temp.meta.df=meta.df[colnames(temp.counts.df),]
    groups=factor(as.character(temp.meta.df$development.stage))
    if(parental.development){
      groups=factor(as.character(temp.meta.df$parent.development.stage))
    }
    temp.out.list=list()
    names(groups)=as.character(colnames(temp.counts.df))
    no.cores=2
    o.ifm <- try(scde.error.models(counts =temp.counts.df,groups=groups,n.cores=no.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1));
    if(class(o.ifm)=='try-error'){
      next
    }
    valid.cells <- o.ifm$corr.a >0;
    o.ifm <- o.ifm[valid.cells,];
    temp.counts.df=counts.df[,rownames(o.ifm)]
    temp.meta.df=meta.df[rownames(o.ifm),]
    temp.sample.development.stage=as.character(temp.meta.df$development.stage)
    if(parental.development){
      temp.sample.development.stage=as.character(temp.meta.df$parent.development.stage)
    }
    if (length(unique(temp.sample.development.stage))<2){
      next
    }
    no.samples.per.group=as.numeric(as.character(table(temp.sample.development.stage)))
    show(no.samples.per.group)
    groups=factor(as.character(temp.meta.df$development.stage))
    if(parental.development){
      groups=factor(as.character(temp.meta.df$parent.development.stage))
    }
    names(groups)=as.character(rownames(temp.meta.df))
    batch=factor(as.character(temp.meta.df$batch))
    temp.out.list[['groups']]=groups
    temp.out.list[['batch']]=batch
    o.prior <- try(scde.expression.prior(models=o.ifm,counts=temp.counts.df,length.out=400,show.plot=F))
    if(class(o.prior)=='try-error'){
      next
    }
    no.batches=length(unique(batch))
    ediff=data.frame()
    if(correct.batch.effect){
      ediff <- scde.expression.difference(models = o.ifm,counts = temp.counts.df,prior =o.prior,groups=groups,n.randomizations=100,n.cores=no.cores,verbose=1,batch = batch,return.posteriors = T)
    }
    else{
      ediff <- scde.expression.difference(models = o.ifm,counts = temp.counts.df,prior =o.prior,groups=groups,n.randomizations=100,n.cores=no.cores,verbose=1,return.posteriors = T)
    }
    ediff.res=ediff$results
    ediff.res$p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$Z))))
    ediff.res$adj.p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$cZ))))
    temp.out.list[['o.ifm']]=o.ifm
    temp.out.list[['o.prior']]=o.prior
    temp.out.list[['ediff']]=ediff
    temp.out.list[['results']]=ediff.res
    temp.out.list[['count']]=temp.counts.df
    temp.out.list[['meta']]=temp.meta.df
    out.list[[temp.name]]=temp.out.list
  }
  return(out.list)
}
run.scde.between.sub.types=function(counts.df,meta.df){
  meta.df=meta.df[colnames(counts.df),]
  all.samples=rownames(meta.df)
  meta.list=split(meta.df,f = meta.df$development.stage)
  stages.names=names(meta.list)
  #stages.names=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  len.stages=length(stages.names)
  out.list=list()
  for(n in 1:len.stages){
    development.stage=stages.names[n]
    temp.meta.df=meta.list[[development.stage]]
    temp.counts.df=filter.none.expressed.genes(counts.df[,rownames(temp.meta.df)])
    groups=as.factor(as.character(temp.meta.df$sub.group))
    names(groups)=as.character(colnames(temp.counts.df))
    no.cores=2
    groups.categories=unique(groups)
    len.groups.categories=length(groups.categories)
    if(len.groups.categories==2){
      temp.out.list=list()
      o.ifm <- try(scde.error.models(counts =temp.counts.df,groups=groups,n.cores=no.cores,threshold.segmentation=T,save.crossfit.plots=T,save.model.plots=T,verbose=1));
      if(class(o.ifm)=='try-error'){
        next
      }
      valid.cells <- o.ifm$corr.a >0;
      o.ifm <- o.ifm[valid.cells,];
      temp.counts.df=temp.counts.df[,rownames(o.ifm)]
      temp.meta.df=meta.df[rownames(o.ifm),]
      names(groups)=as.character(colnames(temp.counts.df))
      o.prior <- try(scde.expression.prior(models=o.ifm,counts=temp.counts.df,length.out=400,show.plot=T))
      if(class(o.prior)=='try-error'){
        next
      }
      ediff <- scde.expression.difference(models = o.ifm,counts = temp.counts.df,prior =o.prior,groups=groups,n.randomizations=100,n.cores=no.cores,verbose=1,return.posteriors = T)
      ediff.res=ediff$results
      ediff.res$p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$Z))))
      ediff.res$adj.p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$cZ))))
      temp.out.list[['o.ifm']]=o.ifm
      temp.out.list[['o.prior']]=o.prior
      temp.out.list[['ediff']]=ediff
      temp.out.list[['results']]=ediff.res
      temp.out.list[['count']]=temp.counts.df
      temp.out.list[['meta']]=temp.meta.df
      out.list[[development.stage]]=temp.out.list
    }
    else{
      temp.meta.list=split(temp.meta.df,f=temp.meta.df$sub.group)
      groups.combn=combn(names(temp.meta.list),2)
      len.groups.combn=dim(groups.combn)[2]
      for(m in 1:len.groups.combn){
        temp.out.list=list()
        sub.groups.pair=as.character(groups.combn[,m])
        first.temp.meta.df=temp.meta.list[[sub.groups.pair[1]]]
        sec.temp.meta.df=temp.meta.list[[sub.groups.pair[2]]]
        sub.temp.meta.df=temp.meta.df[as.character(c(rownames(first.temp.meta.df),rownames(sec.temp.meta.df))),]
        sub.temp.meta.df$sub.group=as.factor(as.character(sub.temp.meta.df$sub.group))
        sub.groups=as.factor(as.character(sub.temp.meta.df$sub.group))
        sub.temp.counts.df=filter.none.expressed.genes(input.data = temp.counts.df[,rownames(sub.temp.meta.df)])
        names(sub.groups)=as.character(colnames(sub.temp.counts.df))
        o.ifm <- try(scde.error.models(counts =sub.temp.counts.df,groups=sub.groups,n.cores=no.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1));
        if(class(o.ifm)=='try-error'){
          next
        }
        valid.cells <- o.ifm$corr.a >0;
        o.ifm <- o.ifm[valid.cells,];
        sub.temp.counts.df=sub.temp.counts.df[,rownames(o.ifm)]
        sub.temp.meta.df=sub.temp.meta.df[rownames(o.ifm),]
        sub.groups=as.factor(as.character(sub.temp.meta.df$sub.group))
        names(sub.groups)=as.character(colnames(sub.temp.counts.df))
        o.prior <- try(scde.expression.prior(models=o.ifm,counts=sub.temp.counts.df,length.out=400,show.plot=T))
        if(class(o.prior)=='try-error'){
          next
        }
        ediff <- scde.expression.difference(models = o.ifm,counts = sub.temp.counts.df,prior =o.prior,groups=sub.groups,n.randomizations=100,n.cores=no.cores,verbose=1,return.posteriors = T)
        ediff.res=ediff$results
        ediff.res$p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$Z))))
        ediff.res$adj.p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$cZ))))
        temp.out.list[['o.ifm']]=o.ifm
        temp.out.list[['o.prior']]=o.prior
        temp.out.list[['ediff']]=ediff
        temp.out.list[['results']]=ediff.res
        temp.out.list[['count']]=sub.temp.counts.df
        temp.out.list[['meta']]=sub.temp.meta.df
        out.list[[paste(development.stage,m,sep='_')]]=temp.out.list
      }
    }
  }
  return(out.list)
}
run.scde.between.committed.sc=function(counts.df,meta.df,correct.batch.effect=F){
  all.samples=colnames(counts.df)
  counts.df=counts.df[,all.samples]
  meta.df=meta.df[all.samples,]
  meta.list=split(meta.df,f = meta.df$sex.comit)
  stages.names=names(meta.list)
  stages.combn=combn(stages.names,m = 2)
  no.combn=dim(stages.combn)
  out.list=list()
  #for(n in 1:2){
  for(n in 1:no.combn[2]){
    temp.out.list=list()
    temp.combn=stages.combn[,n]
    first.stage.name=temp.combn[1]
    sec.stage.name=temp.combn[2]
    temp.name=paste(first.stage.name,sec.stage.name,sep = '_vs_')
    first.meta.df=meta.list[[first.stage.name]]
    sec.meta.df=meta.list[[sec.stage.name]]
    groups=factor(c(rep(x =first.stage.name ,times = dim(first.meta.df)[1]),rep(x =sec.stage.name ,times = dim(sec.meta.df)[1])))
    batch=factor(c(first.meta.df[,'batch'],sec.meta.df[,'batch']))
    first.samples=rownames(first.meta.df)
    sec.samples=rownames(sec.meta.df)
    temp.counts.df=counts.df[,c(first.samples,sec.samples)]
    temp.counts.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = temp.counts.df))
    if(is.null(dim(temp.counts.df))){
      next
    }
    temp.meta.df=meta.df[colnames(temp.counts.df),]
    groups=factor(as.character(temp.meta.df$sex.comit))
    names(groups)=as.character(colnames(temp.counts.df))
    no.cores=2
    o.ifm <- try(scde.error.models(counts =temp.counts.df,groups=groups,n.cores=no.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1));
    if(class(o.ifm)=='try-error'){
      next
    }
    valid.cells <- o.ifm$corr.a >0;
    o.ifm <- o.ifm[valid.cells,];
    temp.counts.df=counts.df[,rownames(o.ifm)]
    temp.meta.df=meta.df[rownames(o.ifm),]
    temp.sample.development.stage=as.character(temp.meta.df$sex.comit)
    if (length(unique(temp.sample.development.stage))<2){
      next
    }
    no.samples.per.group=as.numeric(as.character(table(temp.sample.development.stage)))
    groups=factor(as.character(temp.meta.df$sex.comit))
    names(groups)=as.character(rownames(temp.meta.df))
    batch=factor(as.character(temp.meta.df$batch))
    temp.out.list[['groups']]=groups
    temp.out.list[['batch']]=batch
    o.prior <- try(scde.expression.prior(models=o.ifm,counts=temp.counts.df,length.out=400,show.plot=F))
    if(class(o.prior)=='try-error'){
      next
    }
    no.batches=length(unique(batch))
    ediff=data.frame()
    if(correct.batch.effect){
      ediff <- scde.expression.difference(models = o.ifm,counts = temp.counts.df,prior =o.prior,groups=groups,n.randomizations=100,n.cores=no.cores,verbose=1,batch = batch,return.posteriors = T)
    }
    else{
      ediff <- try(scde.expression.difference(models = o.ifm,counts = temp.counts.df,prior =o.prior,groups=groups,n.randomizations=100,n.cores=no.cores,verbose=1,return.posteriors = T))
      if(class(ediff)=='try-error'){
        next
      }
    }
    ediff.res=ediff$results
    ediff.res$p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$Z))))
    ediff.res$adj.p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$cZ))))
    temp.out.list[['o.ifm']]=o.ifm
    temp.out.list[['o.prior']]=o.prior
    temp.out.list[['ediff']]=ediff
    temp.out.list[['results']]=ediff.res
    temp.out.list[['count']]=temp.counts.df
    temp.out.list[['meta']]=temp.meta.df
    out.list[[temp.name]]=temp.out.list
  }
  return(out.list)
}
run.scde.between.grps=function(counts.df,meta.df,grp,correct.batch.effect=F){
  counts.mat=as.matrix(counts.df)
  counts.mat=scde::clean.counts(counts = counts.mat,min.lib.size = 2000,min.reads = 2,
                                min.detected = 3)
  all.samples=colnames(counts.mat)
  meta.df=meta.df[all.samples,]
  groups=meta.df$grp 
  names(groups)=colnames(counts.mat)
  out.list=list()
  no.cores=2
  o.ifm <- scde::scde.error.models(counts =counts.mat,n.cores=2,min.nonfailed = 5,
                                   threshold.segmentation=T,save.crossfit.plots=F,
                                   linear.fit = T,save.model.plots=F,verbose=1,
                                   min.count.threshold = 2);
  valid.cells <- o.ifm$corr.a >0;
  o.ifm <- o.ifm[valid.cells,];
  counts.mat=counts.mat[,rownames(o.ifm)]
  meta.df=meta.df[rownames(o.ifm),]
  groups=factor(as.character(meta.df[,grp])) 
  names(groups)=as.character(colnames(counts.mat))
  out.list[['groups']]=groups
  o.prior <- scde::scde.expression.prior(models=o.ifm,counts=counts.mat,length.out=400,show.plot=F)
  ediff <- scde::scde.expression.difference(models = o.ifm,counts = counts.mat,prior =o.prior,
                                            groups=groups,n.randomizations=100,n.cores=2,
                                            verbose=1,return.posteriors = T)
  ediff.res=ediff$results
  ediff.res$p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$Z))))
  ediff.res$adj.p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$cZ))))
  out.list[['o.ifm']]=o.ifm
  out.list[['o.prior']]=o.prior
  out.list[['ediff']]=ediff
  out.list[['results']]=ediff.res
  out.list[['count']]=counts.df
  out.list[['meta']]=meta.df
  return(out.list)
}
run.scde.between.committed.sc.test=function(counts.df,meta.df,correct.batch.effect=F){
  counts.df=filter.none.expressed.samples(input.data = filter.none.expressed.genes(input.data = counts.df))
  all.samples=colnames(counts.df)
  counts.df=counts.df[,all.samples]
  meta.df=meta.df[all.samples,]
  batch=factor(meta.df[,'batch'])
  groups=factor(as.character(meta.df$sex.comit)) 
  names(groups)=as.character(colnames(counts.df))
  out.list=list()
  no.cores=2
  o.ifm <- try(scde::scde.error.models(counts =counts.df,groups=groups,n.cores=no.cores,threshold.segmentation=T,
                                       save.crossfit.plots=F,save.model.plots=F,verbose=1));
  valid.cells <- o.ifm$corr.a >0;
  o.ifm <- o.ifm[valid.cells,];
  counts.df=counts.df[,rownames(o.ifm)]
  meta.df=meta.df[rownames(o.ifm),]
  temp.sample.development.stage=as.character(meta.df$sex.comit)
  no.samples.per.group=as.numeric(as.character(table(temp.sample.development.stage)))
  groups=factor(as.character(meta.df$sex.comit))
  names(groups)=as.character(rownames(meta.df))
  batch=factor(as.character(meta.df$batch))
  out.list[['groups']]=groups
  out.list[['batch']]=batch
  o.prior <- try(scde.expression.prior(models=o.ifm,counts=counts.df,length.out=400,show.plot=T))
  no.batches=length(unique(batch))
  ediff=data.frame()
  if(correct.batch.effect){
    ediff <- scde.expression.difference(models = o.ifm,counts = counts.df,prior =o.prior,groups=groups,
                                        n.randomizations=100,n.cores=no.cores,verbose=1,batch = batch,
                                        return.posteriors = T)
  }
  else{
    ediff <- try(scde.expression.difference(models = o.ifm,counts = counts.df,prior =o.prior,groups=groups,
                                            n.randomizations=100,n.cores=no.cores,verbose=1,
                                            return.posteriors = T))
  }
  ediff.res=ediff$results
  ediff.res$p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$Z))))
  ediff.res$adj.p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$cZ))))
  out.list[['o.ifm']]=o.ifm
  out.list[['o.prior']]=o.prior
  out.list[['ediff']]=ediff
  out.list[['results']]=ediff.res
  out.list[['count']]=counts.df
  out.list[['meta']]=meta.df
  return(out.list)
}
run.scde.between.sex.committed.with.sample.balancing.test=function(counts.df,meta.df,correct.batch.effect=F){
  counts.df=filter.none.expressed.samples(input.data  = counts.df)
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  lr.meta.df=subset(meta.df,development.stage=='LR')
  counts.df=subset(counts.df,select = rownames(lr.meta.df))
  all.samples=colnames(counts.df)
  meta.df=meta.df[all.samples,]
  meta.list=split(meta.df,f = meta.df$mRFP1.expr)
  stages.names=names(meta.list)
  stages.combn=combn(stages.names,m = 2)
  no.combn=dim(stages.combn)
  out.list=run.scde.between.committed.sc(counts.df=counts.df,meta.df=meta.df,correct.batch.effect=F)
  return(out.list)
}
run.scde.between.schizont.sub.grps=function(counts.df,meta.df,correct.batch.effect=F){
  all.samples=colnames(counts.df)
  counts.df=counts.df[,all.samples]
  meta.df=meta.df[all.samples,]
  meta.list=split(meta.df,f = meta.df$in.stage.sub.grp)
  stages.names=names(meta.list)
  stages.combn=combn(stages.names,m = 2)
  no.combn=dim(stages.combn)
  out.list=list()
  #for(n in 1:2){
  for(n in 1:no.combn[2]){
    temp.out.list=list()
    temp.combn=stages.combn[,n]
    first.stage.name=temp.combn[1]
    sec.stage.name=temp.combn[2]
    temp.name=paste(first.stage.name,sec.stage.name,sep = '_vs_')
    first.meta.df=meta.list[[first.stage.name]]
    sec.meta.df=meta.list[[sec.stage.name]]
    groups=factor(c(rep(x =first.stage.name ,times = dim(first.meta.df)[1]),rep(x =sec.stage.name ,times = dim(sec.meta.df)[1])))
    batch=factor(c(first.meta.df[,'batch'],sec.meta.df[,'batch']))
    first.samples=rownames(first.meta.df)
    sec.samples=rownames(sec.meta.df)
    temp.counts.df=counts.df[,c(first.samples,sec.samples)]
    temp.counts.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = temp.counts.df))
    if(is.null(dim(temp.counts.df))){
      next
    }
    temp.meta.df=meta.df[colnames(temp.counts.df),]
    groups=factor(as.character(temp.meta.df$in.stage.sub.grp))
    names(groups)=as.character(colnames(temp.counts.df))
    no.cores=2
    o.ifm <- try(scde.error.models(counts =temp.counts.df,groups=groups,n.cores=no.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1));
    if(class(o.ifm)=='try-error'){
      next
    }
    valid.cells <- o.ifm$corr.a >0;
    o.ifm <- o.ifm[valid.cells,];
    temp.counts.df=counts.df[,rownames(o.ifm)]
    temp.meta.df=meta.df[rownames(o.ifm),]
    temp.sample.development.stage=as.character(temp.meta.df$in.stage.sub.grp)
    if (length(unique(temp.sample.development.stage))<2){
      next
    }
    no.samples.per.group=as.numeric(as.character(table(temp.sample.development.stage)))
    groups=factor(as.character(temp.meta.df$in.stage.sub.grp))
    names(groups)=as.character(rownames(temp.meta.df))
    batch=factor(as.character(temp.meta.df$batch))
    temp.out.list[['groups']]=groups
    temp.out.list[['batch']]=batch
    o.prior <- try(scde.expression.prior(models=o.ifm,counts=temp.counts.df,length.out=400,show.plot=F))
    if(class(o.prior)=='try-error'){
      next
    }
    no.batches=length(unique(batch))
    ediff=data.frame()
    if(correct.batch.effect){
      ediff <- scde.expression.difference(models = o.ifm,counts = temp.counts.df,prior =o.prior,groups=groups,n.randomizations=100,n.cores=no.cores,verbose=1,batch = batch,return.posteriors = T)
    }
    else{
      ediff <- try(scde.expression.difference(models = o.ifm,counts = temp.counts.df,prior =o.prior,groups=groups,n.randomizations=100,n.cores=no.cores,verbose=1,return.posteriors = T))
      if(class(ediff)=='try-error'){
        next
      }
    }
    ediff.res=ediff$results
    ediff.res$p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$Z))))
    ediff.res$adj.p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$cZ))))
    temp.out.list[['o.ifm']]=o.ifm
    temp.out.list[['o.prior']]=o.prior
    temp.out.list[['ediff']]=ediff
    temp.out.list[['results']]=ediff.res
    temp.out.list[['count']]=temp.counts.df
    temp.out.list[['meta']]=temp.meta.df
    out.list[[temp.name]]=temp.out.list
  }
  return(out.list)
}
run.scde.between.marker.based.clusters=function(counts.df,meta.df,correct.batch.effect=F){
  all.samples=colnames(counts.df)
  counts.df=counts.df[,all.samples]
  meta.df=meta.df[all.samples,]
  meta.list=split(meta.df,f = meta.df$markers.cluster.groups)
  stages.names=names(meta.list)
  stages.combn=combn(stages.names,m = 2)
  no.combn=dim(stages.combn)
  out.list=list()
  #for(n in 1:2){
  for(n in 1:no.combn[2]){
    temp.out.list=list()
    temp.combn=stages.combn[,n]
    first.stage.name=temp.combn[1]
    sec.stage.name=temp.combn[2]
    temp.name=paste(first.stage.name,sec.stage.name,sep = '_vs_')
    first.meta.df=meta.list[[first.stage.name]]
    sec.meta.df=meta.list[[sec.stage.name]]
    groups=factor(c(rep(x =first.stage.name ,times = dim(first.meta.df)[1]),rep(x =sec.stage.name ,times = dim(sec.meta.df)[1])))
    batch=factor(c(first.meta.df[,'batch'],sec.meta.df[,'batch']))
    first.samples=rownames(first.meta.df)
    sec.samples=rownames(sec.meta.df)
    temp.counts.df=counts.df[,c(first.samples,sec.samples)]
    temp.counts.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = temp.counts.df))
    if(is.null(dim(temp.counts.df))){
      next
    }
    temp.meta.df=meta.df[colnames(temp.counts.df),]
    groups=factor(as.character(temp.meta.df$markers.cluster.groups))
    names(groups)=as.character(colnames(temp.counts.df))
    no.cores=2
    o.ifm <- try(scde.error.models(counts =temp.counts.df,groups=groups,n.cores=no.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1));
    if(class(o.ifm)=='try-error'){
      next
    }
    valid.cells <- o.ifm$corr.a >0;
    o.ifm <- o.ifm[valid.cells,];
    temp.counts.df=counts.df[,rownames(o.ifm)]
    temp.meta.df=meta.df[rownames(o.ifm),]
    temp.sample.development.stage=as.character(temp.meta.df$markers.cluster.groups)
    if (length(unique(temp.sample.development.stage))<2){
      next
    }
    no.samples.per.group=as.numeric(as.character(table(temp.sample.development.stage)))
    groups=factor(as.character(temp.meta.df$markers.cluster.groups))
    names(groups)=as.character(rownames(temp.meta.df))
    batch=factor(as.character(temp.meta.df$batch))
    temp.out.list[['groups']]=groups
    temp.out.list[['batch']]=batch
    o.prior <- try(scde.expression.prior(models=o.ifm,counts=temp.counts.df,length.out=400,show.plot=F))
    if(class(o.prior)=='try-error'){
      next
    }
    no.batches=length(unique(batch))
    ediff=data.frame()
    if(correct.batch.effect){
      ediff <- scde.expression.difference(models = o.ifm,counts = temp.counts.df,prior =o.prior,groups=groups,n.randomizations=100,n.cores=no.cores,verbose=1,batch = batch,return.posteriors = T)
    }
    else{
      ediff <- try(scde.expression.difference(models = o.ifm,counts = temp.counts.df,prior =o.prior,groups=groups,n.randomizations=100,n.cores=no.cores,verbose=1,return.posteriors = T))
      if(class(ediff)=='try-error'){
        next
      }
    }
    ediff.res=ediff$results
    ediff.res$p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$Z))))
    ediff.res$adj.p.value=2*pnorm(-abs(as.numeric(as.character(ediff.res$cZ))))
    temp.out.list[['o.ifm']]=o.ifm
    temp.out.list[['o.prior']]=o.prior
    temp.out.list[['ediff']]=ediff
    temp.out.list[['results']]=ediff.res
    temp.out.list[['count']]=temp.counts.df
    temp.out.list[['meta']]=temp.meta.df
    out.list[[temp.name]]=temp.out.list
  }
  return(out.list)
}
run.scde.between.cluster.at.diff.intervals=function(counts.df,marker.intervals.res.list){
  intervals.names=names(marker.intervals.res.list)
  len.intervals.names=length(intervals.names)
  out.scde.res.list=list()
  for(m in 1:len.intervals.names){
    intervals.name=intervals.names[m]
    temp.interval.clusters.list=marker.intervals.res.list[[intervals.name]]
    temp.meta.df=temp.interval.clusters.list$meta.df
    temp.counts.df=filter.none.expressed.genes(input.data = counts.df[,rownames(temp.meta.df)])
    temp.scde.res.list=run.scde.between.marker.based.clusters(counts.df = temp.counts.df,meta.df =temp.meta.df,correct.batch.effect = T )
    out.scde.res.list[[intervals.name]]=temp.scde.res.list
  }
  return(out.scde.res.list)
}
split.marker.based.sub.grps.to.balance=function(meta.df){
  meta.list=split(meta.df,f=meta.df$markers.cluster.groups)
  new.sub.groups=c()
  sub.grp.names=names(meta.list)
  len.sub.grp.names=length(sub.grp.names)
  for(m in 1:len.sub.grp.names){
    sub.grp.name=sub.grp.names[m]
    temp.meta.df=meta.list[[sub.grp.name]]
    if(sub.grp.name=='SC.grp.1'){
      sample.names=rownames(temp.meta.df)
      first.set=sample(x = sample.names,size = 20)
      sec.set=setdiff(sample.names,first.set)
    }
    else{
    }
  }
}
run.scde.between.marker.based.clusters.with.sample.balancing=function(counts.df,meta.df,correct.batch.effect=F){
  all.samples=colnames(counts.df)
  meta.df=meta.df[all.samples,]
  meta.list=split(meta.df,f = meta.df$markers.cluster.groups)
  stages.names=names(meta.list)
  stages.combn=combn(stages.names,m = 2)
  no.combn=dim(stages.combn)
  out.list=list()
  #for(n in 1:2){
  for(n in 1:no.combn[2]){
    temp.out.list=list()
    temp.combn=stages.combn[,n]
    first.stage.name=temp.combn[1]
    sec.stage.name=temp.combn[2]
    temp.name=paste(first.stage.name,sec.stage.name,sep = '_vs_')
    first.meta.df=meta.list[[first.stage.name]]
    sec.meta.df=meta.list[[sec.stage.name]]
    first.samples.counts=dim(first.meta.df)[1]
    sec.samples.counts=dim(sec.meta.df)[1]
    temp.name=paste(first.stage.name,sec.stage.name,sep = '_vs_')
    leading.sample.count=max(first.samples.counts,sec.samples.counts)
    least.sample.count=min(first.samples.counts,sec.samples.counts)
    fold.diff=trunc(leading.sample.count / least.sample.count)
    remainder.sample.count=leading.sample.count %% least.sample.count
    if(fold.diff>=2){
      leading.meta.df=data.frame()
      least.meta.df=data.frame()
      if(dim(first.meta.df)[1]==leading.sample.count){
        leading.meta.df=first.meta.df
        least.meta.df=sec.meta.df
      }
      else{
        leading.meta.df=sec.meta.df
        least.meta.df=first.meta.df
      }
      leading.meta.samples=rownames(leading.meta.df)
      least.meta.samples=rownames(least.meta.df)
      subset.samples.list=list()
      sample.tracter=c()
      for( m in 1:fold.diff){
        temp.samples.vec=c()
        if(m!=fold.diff){
          temp.leading.samples.updated.vec=setdiff(leading.meta.samples,sample.tracter)
          temp.leading.samples.vec=sample(temp.leading.samples.updated.vec,size = least.sample.count)
          sample.tracter=append(sample.tracter,temp.leading.samples.vec,length(sample.tracter))
          #show(length(temp.leading.samples.vec))
          subset.samples.list[[m]]=temp.leading.samples.vec
          temp.samples.vec=c(temp.leading.samples.vec,least.meta.samples)
        }
        else{
          if(remainder.sample.count<3){
            temp.leading.samples.updated.vec=setdiff(leading.meta.samples,sample.tracter)
            temp.leading.samples.vec=temp.leading.samples.updated.vec
            sample.tracter=append(sample.tracter,temp.leading.samples.vec,length(sample.tracter))
            subset.samples.list[[m]]=temp.leading.samples.vec
            temp.samples.vec=c(temp.leading.samples.vec,least.meta.samples)
          }
          else{
            for(s in 1:2){
              if(s==1){
                temp.leading.samples.updated.vec=setdiff(leading.meta.samples,sample.tracter)
                temp.leading.samples.vec=sample(temp.leading.samples.updated.vec,size = least.sample.count)
                sample.tracter=append(sample.tracter,temp.leading.samples.vec,length(sample.tracter))
                subset.samples.list[[m]]=temp.leading.samples.vec
                temp.samples.vec=c(temp.leading.samples.vec,least.meta.samples)
              }
              else{
                temp.leading.samples.updated.vec=setdiff(leading.meta.samples,sample.tracter)
                temp.leading.samples.vec=sample(temp.leading.samples.updated.vec,size = remainder.sample.count)
                sample.tracter=append(sample.tracter,temp.leading.samples.vec,length(sample.tracter))
                subset.samples.list[[m]]=temp.leading.samples.vec
                last.samples.vec=sample(sample.tracter,size = least.sample.count-remainder.sample.count)
                temp.leading.samples.vec=append(temp.leading.samples.vec,last.samples.vec,length(temp.leading.samples.vec))
                temp.samples.vec=c(temp.leading.samples.vec,least.meta.samples)
              }
            }
          }
        }
        temp.counts.df=filter.none.expressed.genes( df = counts.df[,temp.samples.vec])
        temp.meta.df=meta.df[colnames(temp.counts.df),]
        temp.scde.res=run.scde.between.marker.based.clusters(counts.df = temp.counts.df,meta.df =temp.meta.df,correct.batch.effect = F )
        temp.out.list[[paste('scde.res.',m,sep = '')]]=temp.scde.res
      }
    }
    else{
        if(remainder.sample.count<3){
          temp.sample.vec=c(rownames( first.meta.df),rownames( sec.meta.df))
          temp.counts.df=filter.none.expressed.genes( df = counts.df[,temp.samples.vec])
          temp.meta.df=meta.df[colnames(temp.counts.df),]
          temp.scde.res=run.scde.between.marker.based.clusters(counts.df = temp.counts.df,meta.df =temp.meta.df,correct.batch.effect = F )
          temp.out.list[['scde.res.']]=temp.scde.res
        }
      else{
        sample.tracter=c()
        for(s in 1:2){
          if(s==1){
            temp.leading.samples.updated.vec=setdiff(leading.meta.samples,sample.tracter)
            temp.leading.samples.vec=sample(temp.leading.samples.updated.vec,size = least.sample.count)
            sample.tracter=append(sample.tracter,temp.leading.samples.vec,length(sample.tracter))
            subset.samples.list[[m]]=temp.leading.samples.vec
            temp.samples.vec=c(temp.leading.samples.vec,least.meta.samples)
            temp.counts.df=filter.none.expressed.genes( df = counts.df[,temp.samples.vec])
            temp.meta.df=meta.df[colnames(temp.counts.df),]
            temp.scde.res=run.scde.between.marker.based.clusters(counts.df = temp.counts.df,meta.df =temp.meta.df,correct.batch.effect = F )
            temp.out.list[[paste('scde.res.',s,sep='')]]=temp.scde.res
          }
          else{
            temp.leading.samples.updated.vec=setdiff(leading.meta.samples,sample.tracter)
            temp.leading.samples.vec=sample(temp.leading.samples.updated.vec,size = remainder.sample.count)
            sample.tracter=append(sample.tracter,temp.leading.samples.vec,length(sample.tracter))
            subset.samples.list[[m]]=temp.leading.samples.vec
            last.samples.vec=sample(sample.tracter,size = least.sample.count-remainder.sample.count)
            temp.leading.samples.vec=append(temp.leading.samples.vec,last.samples.vec,length(temp.leading.samples.vec))
            temp.samples.vec=c(temp.leading.samples.vec,least.meta.samples)
            temp.counts.df=filter.none.expressed.genes( df = counts.df[,temp.samples.vec])
            temp.meta.df=meta.df[colnames(temp.counts.df),]
            temp.scde.res=run.scde.between.marker.based.clusters(counts.df = temp.counts.df,meta.df =temp.meta.df,correct.batch.effect = F )
            temp.out.list[[paste('scde.res.',s,sep='')]]=temp.scde.res
          }
        }
      }
  }
  out.list[[temp.name]]=temp.out.list
  }
  return(out.list)
}
run.scde.between.sex.committed.with.sample.balancing=function(counts.df,meta.df,correct.batch.effect=F){
  counts.df=filter.none.expressed.samples(df = counts.df)
  all.samples=colnames(counts.df)
  meta.df=meta.df[all.samples,]
  meta.list=split(meta.df,f = meta.df$mRFP1.expr)
  stages.names=names(meta.list)
  stages.combn=combn(stages.names,m = 2)
  no.combn=dim(stages.combn)
  out.list=list()
  #for(n in 1:2){
  for(n in 1:no.combn[2]){
    temp.out.list=list()
    temp.combn=stages.combn[,n]
    first.stage.name=temp.combn[1]
    sec.stage.name=temp.combn[2]
    first.meta.df=meta.list[[first.stage.name]]
    sec.meta.df=meta.list[[sec.stage.name]]
    first.samples.counts=dim(first.meta.df)[1]
    sec.samples.counts=dim(sec.meta.df)[1]
    temp.name=paste(first.stage.name,sec.stage.name,sep = '_vs_')
    leading.sample.count=max(first.samples.counts,sec.samples.counts)
    least.sample.count=min(first.samples.counts,sec.samples.counts)
    fold.diff=trunc(leading.sample.count / least.sample.count)
    remainder.sample.count=leading.sample.count %% least.sample.count
    if(fold.diff>=2){
      leading.meta.df=data.frame()
      least.meta.df=data.frame()
      if(dim(first.meta.df)[1]==leading.sample.count){
        leading.meta.df=first.meta.df
        least.meta.df=sec.meta.df
      }
      else{
        leading.meta.df=sec.meta.df
        least.meta.df=first.meta.df
      }
      leading.meta.samples=rownames(leading.meta.df)
      least.meta.samples=rownames(least.meta.df)
      subset.samples.list=list()
      sample.tracter=c()
      for( m in 1:fold.diff){
        temp.samples.vec=c()
        if(m!=fold.diff){
          temp.leading.samples.updated.vec=setdiff(leading.meta.samples,sample.tracter)
          temp.leading.samples.vec=sample(temp.leading.samples.updated.vec,size = least.sample.count)
          sample.tracter=append(sample.tracter,temp.leading.samples.vec,length(sample.tracter))
          subset.samples.list[[m]]=temp.leading.samples.vec
          temp.samples.vec=c(temp.leading.samples.vec,least.meta.samples)
        }
        else{
          if(remainder.sample.count<3){
            temp.leading.samples.updated.vec=setdiff(leading.meta.samples,sample.tracter)
            temp.leading.samples.vec=temp.leading.samples.updated.vec
            sample.tracter=append(sample.tracter,temp.leading.samples.vec,length(sample.tracter))
            subset.samples.list[[m]]=temp.leading.samples.vec
            temp.samples.vec=c(temp.leading.samples.vec,least.meta.samples)
          }
          else{
            for(s in 1:2){
              if(s==1){
                temp.leading.samples.updated.vec=setdiff(leading.meta.samples,sample.tracter)
                temp.leading.samples.vec=sample(temp.leading.samples.updated.vec,size = least.sample.count)
                sample.tracter=append(sample.tracter,temp.leading.samples.vec,length(sample.tracter))
                subset.samples.list[[m]]=temp.leading.samples.vec
                temp.samples.vec=c(temp.leading.samples.vec,least.meta.samples)
              }
              else{
                temp.leading.samples.updated.vec=setdiff(leading.meta.samples,sample.tracter)
                temp.leading.samples.vec=sample(temp.leading.samples.updated.vec,size = remainder.sample.count)
                sample.tracter=append(sample.tracter,temp.leading.samples.vec,length(sample.tracter))
                subset.samples.list[[m]]=temp.leading.samples.vec
                last.samples.vec=sample(sample.tracter,size = least.sample.count-remainder.sample.count)
                temp.leading.samples.vec=append(temp.leading.samples.vec,last.samples.vec,length(temp.leading.samples.vec))
                temp.samples.vec=c(temp.leading.samples.vec,least.meta.samples)
              }
            }
          }
        }
        temp.counts.df=filter.none.expressed.genes( df = counts.df[,temp.samples.vec])
        temp.meta.df=meta.df[colnames(temp.counts.df),]
        temp.scde.res=run.scde.between.committed.sc(counts.df = temp.counts.df,meta.df =temp.meta.df,correct.batch.effect = F )
        temp.out.list[[paste('scde.res.',m,sep = '')]]=temp.scde.res
      }
    }
    out.list[[temp.name]]=temp.out.list
  }
  return(out.list)
}
cluster.by.markers=function(rpkm.df,meta.df,markers.df,title.str){
  genes=intersect(rownames(rpkm.df),rownames(markers.df))
  markers.rpkm.df=rpkm.df[genes,]
  markers.rpkm.df=reset.rpkm.less.than.cut.off(df =markers.rpkm.df,rpkm.cutoff = 1 )
  markers.rpkm.df= filter.none.expressed.genes(input.data = markers.rpkm.df)
  markers.rpkm.df=filter.non.variable.rows(df = markers.rpkm.df,cut.off = 1)
  genes=rownames(markers.rpkm.df)
  #log.markers.rpkm.df=log10(markers.rpkm.df+1)
  #log.markers.rpkm.df=filter.non.variable.rows(df = log.markers.rpkm.df,cut.off = 1)
  #markers.rpkm.df=log.markers.rpkm.df
  samples.col.fact.list=get.col.factor(col.factor = as.character(meta.df$development.stage))
  genes.col.factor.str=as.character(markers.df[genes,'category'])
  genes.col.factor.list= get.col.factor(col.factor = genes.col.factor.str)
  genes.col.str=genes.col.factor.list$col.str
  gene.cor.dist=cor(x = t(markers.rpkm.df),method = 'spearman')
  gene.cor.dist=as.dist(1-gene.cor.dist)
  gene.cor.hclust=hclust(d = gene.cor.dist,method = 'average')
  gene.dist.dend=as.dendrogram(gene.cor.hclust)
  plot(gene.cor.hclust)
  meta.df=meta.df[colnames(markers.rpkm.df),]
  meta.df=meta.df[order(meta.df$development.stage),]
  heatmap.2(x = as.matrix(log10(markers.rpkm.df+1)),ColSideColors = samples.col.fact.list$col.str,RowSideColors = genes.col.str,dendrogram = 'col',trace='none',margins = c(10,10),cexRow = .5,labCol = '',col=greenred(3))
  heatmap.2(x = as.matrix(log10(markers.rpkm.df+1)),ColSideColors = samples.col.fact.list$col.str,RowSideColors = genes.col.str,dendrogram = 'col',trace='none',margins = c(10,10),cexRow = .5,labCol = '',col=greenred(3))
  legend('topright',legend=samples.col.fact.list$legend.str,fill=samples.col.fact.list$legend.col)
  legend('bottomleft',legend=genes.col.factor.list$legend.str,fill=genes.col.factor.list$legend.col)
  heatmap.2(x = as.matrix(gene.cor.dist),RowSideColors = genes.col.str,dendrogram = 'none',trace='none',margins = c(10,10),cexRow = .5,labCol = '')
  legend('left',legend=genes.col.factor.list$legend.str,fill=genes.col.factor.list$legend.col)
  plot.samples.correlations.heatmap(in.rpkm.df =markers.rpkm.df,title.str = title.str,in.meta.df =meta.df,filter.non.var.rows =1 )
}
create.df.from.file.list=function(files.list){
  no.files=length(files.list)
  samples.names=c()
  genes=c()
  out.list=list()
  for (m in 1:no.files){
    file.name=files.list[m]
    temp.df=read.table(file.name,header = T,sep = '\t',row.names=1,as.is = T)
    out.list=append(out.list,temp.df,after=length(out.list))
    if(m ==1){
      genes=append(genes,rownames(temp.df),after = length(genes))
    }
  }
  out.df=data.frame(out.list)
  rownames(out.df)=genes
  out.df=edit.samples.df.names(df =out.df )
  #out.df=edit.gene.ids(rpkm.df = out.df)
  return(out.df)
}
create.gene.body.cov.bias.df=function(files.list){
  no.files=length(files.list)
  samples.names=c()
  genes=c()
  out.list=list()
  for (m in 1:no.files){
    file.name=files.list[m]
    temp.df=read.table(file.name,header = T,sep = '\t',row.names=1,as.is = T)
    out.list=append(out.list,temp.df,after=length(out.list))
    if(m ==1){
      genes=append(genes,rownames(temp.df),after = length(genes))
    }
  }
  out.df=data.frame(out.list)
  rownames(out.df)=genes
  out.df=edit.samples.df.names(df =out.df )
  out.df=edit.gene.ids(rpkm.df = out.df)
  return(out.df)
}
encode.gene.names.to.refseq.id=function(ref.refseq.gene.symbol,query.refseq.id.vec){
  len.query.refseq.id.vec=length(query.refseq.id.vec)
  for(n in 1:len.query.refseq.id.vec){
    query.gene=query.refseq.id.vec[n]
    query.parts=strsplit(x = query.gene,split = '\\.')
    len.query.parts=length(query.parts)
    if(len.query.parts==1){
    }
    else{
    }
  }
}
filter.samples.based.on.gene.body.cov=function(gene.body.cov.df,rpkm.df,meta.df,cov=.3){
  genes=intersect(rownames(gene.body.cov.df),rownames(rpkm.df))
  gene.body.cov.df=gene.body.cov.df[genes,]
  rpkm.df=rpkm.df[genes,]
  rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  gene.body.cov.df=gene.body.cov.df[rownames(rpkm.df),]
  samples.rpkm=colnames(rpkm.df)
  rpkm.df=rpkm.df[,samples.rpkm]
  meta.df=meta.df[samples.rpkm,]
  temp.names=rownames(meta.df)
  intersect.samples=intersect(colnames(gene.body.cov.df),temp.names)
  gene.body.cov.df=gene.body.cov.df[,intersect.samples]
  meta.df=meta.df[intersect.samples,]
  rpkm.df=rpkm.df[,intersect.samples]
  rpkm.df[gene.body.cov.df<as.numeric(cov)]=0
  rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  rpkm.df=filter.none.expressed.samples(df = rpkm.df)
  meta.df=meta.df[colnames(rpkm.df),]
  gene.body.cov.df=gene.body.cov.df[rownames(rpkm.df),colnames(rpkm.df)]
  out.list=list(coverage=gene.body.cov.df,rpkm=rpkm.df,meta=meta.df)
  return(out.list)
}
plot.average.cell.to.pop.cor.per.timepoint=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df,no.scs =2){
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df[colnames(sc.rpkm.df),])
  pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = pop.meta.df[colnames(pop.rpkm.df),])
  sc.meta.df$development.stage=factor(sc.meta.df$development.stage,levels = c("R" ,"LR", "ET", "T" , "ES"  , "S"))
  pop.meta.df$development.stage=factor(pop.meta.df$development.stage,levels = c("R" ,"LR", "ET", "T" , "ES"  , "S"))
  sc.timepoints=c(sc.meta.df$development.stage)
  pop.timepoints=c(pop.meta.df$development.stage)
  sc.meta.list=split(x=sc.meta.df,f=sc.meta.df$development.stage)
  pop.meta.list=split(x=pop.meta.df,f=pop.meta.df$development.stage)
  #development.stages=intersect(names(sc.meta.list),names(pop.meta.list))
  #development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  stages.col.list=get.abbrev.std.stage.cols(in.stages.vec = development.stages )
  #stages.col.str=stages.col.list$legend.col
  #stages.col.str=stages.col.list$legend.col
  #stages.name.str=stages.col.list$legend.str
  stages.name.str=development.stages
  col.list=abbrev.std.stages.col.list
  for(m in 1:length(development.stages)){
    development.stage=development.stages[m]
    temp.col.vec=c()
    all.timepoints.cor.vec=c()
    all.timepoints.cor.fact.vec=c()
    temp.col.vec=append(temp.col.vec,col.list[[development.stage]],length(temp.col.vec))
    temp.sc.meta.df=sc.meta.list[[development.stage]]
    temp.sc.rpkm.df=sc.rpkm.df[,rownames(temp.sc.meta.df)]
    temp.pop.meta.df=pop.meta.list[[development.stage]]
    temp.pop.rpkm.df=pop.rpkm.df[,rownames(temp.pop.meta.df)]
    other.pop.samples=setdiff(colnames(pop.rpkm.df),colnames(temp.pop.rpkm.df))
    other.pop.rpkms.df=pop.rpkm.df[,other.pop.samples]
    other.pop.rpkms.df=filter.none.expressed.genes(input.data=other.pop.rpkms.df)
    no.sc.samples=dim(temp.sc.rpkm.df)[2]
    sc.pop.cor.list=list()
    temp.aggregate.sc.rpkm.df=create.aggregate.cells(sc.df =temp.sc.rpkm.df,no.scs =no.scs  )
    if(dim(temp.aggregate.sc.rpkm.df)[2]<2){
      next
    }
    cor.same.pop.vec=get.cor.between.two.df.cor(first.df =temp.aggregate.sc.rpkm.df,second.df = temp.pop.rpkm.df )
    sc.pop.cor.list[[development.stage]]=cor.same.pop.vec
    other.stages=setdiff(development.stages,development.stage)
    len.other.stages=length(other.stages)
    p.values=c('NA')
    for(n in 1:len.other.stages){
      temp.other.stages=other.stages[n]
      temp.col.vec=append(temp.col.vec,col.list[[temp.other.stages]],length(temp.col.vec))
      temp.other.stages.meta.df=pop.meta.list[[temp.other.stages]]
      temp.other.stages.samples=rownames(temp.other.stages.meta.df)
      temp.other.stages.rpkm.df=pop.rpkm.df[,temp.other.stages.samples]
      temp.other.stages.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(temp.other.stages.rpkm.df))
      temp.pop.other.stages.cor=get.cor.between.two.df.cor(first.df =temp.aggregate.sc.rpkm.df ,second.df = temp.other.stages.rpkm.df)
      both.timepoints.cor.vec=as.numeric(c(cor.same.pop.vec,temp.pop.other.stages.cor))
      cor.fact.vec=c(rep(development.stage,length(cor.same.pop.vec)),rep(temp.other.stages,length(temp.pop.other.stages.cor)))
      temp.cor.down.df=data.frame(cor.coeff=both.timepoints.cor.vec)
      cor.fact=factor(x = cor.fact.vec)
      down.sampled.list=downSample(x = temp.cor.down.df,y =cor.fact,list = T )
      down.sampled.df=down.sampled.list$x
      down.sampled.df$class=as.character(down.sampled.list$y)
      down.sampled.list=split(down.sampled.df,f=down.sampled.df$class)
      down.sampled.same.df=subset(down.sampled.list[[development.stage]],select = cor.coeff)
      down.sampled.same.vec=as.numeric(down.sampled.same.df[,1])
      down.sampled.diff.df=subset(down.sampled.list[[temp.other.stages]],select = cor.coeff)
      down.sampled.diff.vec=as.numeric(down.sampled.diff.df[,1])
      #sc.pop.cor.list[[temp.other.stages]]=temp.pop.other.stages.cor
      sc.pop.cor.list[[temp.other.stages]]=down.sampled.diff.vec
      #wilcox.test.p.value=format(wilcox.test(x=cor.same.pop.vec,y=temp.pop.other.stages.cor,paired=F,alternative = 'greater')$p.value,digits=3)
      wilcox.test.p.value=format(wilcox.test(x=down.sampled.same.vec,y=down.sampled.diff.vec,paired=F,alternative = 'greater')$p.value,digits=3)
      p.values=append(p.values,wilcox.test.p.value,length( p.values))
    }
    names(p.values)=c(development.stage,other.stages)
    #cor.diff.pop.vec=get.cor.between.two.df.cor(first.df =temp.aggregate.sc.rpkm.df,second.df = other.pop.rpkms.df)
    #t.test.p.value=format(t.test(cor.same.pop.vec,cor.diff.pop.vec, alternative = 'greater')$p.value,digits=3)
    #wilcox.test.p.value=format(wilcox.test(x=cor.same.pop.vec,y=cor.diff.pop.vec,paired=F,alternative = 'greater')$p.value,digits=3)
    #title.str=paste('P-value: ',wilcox.test.p.value,sep='')
    #temp.boxplot=boxplot(x=list(cor.same.pop.vec,cor.diff.pop.vec),col=c('blue','red'),main=development.stage,names=c('',''),ylab='Spearman cor.')
    #text(x=1.5,y=max(temp.boxplot$stats),labels=c(title.str),cex=.7)
    re.ordered.cor.list=list()
    input.timepoints=c("R" ,"LR", "ET", "T" , "ES"  , "S")
    len.input.timepoints=length(input.timepoints)
    temp.re.ordered.col.vec=c()
    re.ordered.p.values=c()
    for(g in 1:len.input.timepoints){
      temp.tm.point=input.timepoints[g]
      re.ordered.cor.list[[temp.tm.point]]=sc.pop.cor.list[[temp.tm.point]]
      temp.re.ordered.col.vec=append(temp.re.ordered.col.vec,abbrev.std.stages.col.list[[temp.tm.point]],length(temp.re.ordered.col.vec))
      re.ordered.p.values=append(re.ordered.p.values,p.values[temp.tm.point],length(re.ordered.p.values))
    }
    #temp.boxplot=boxplot(x=sc.pop.cor.list,las=2,main=development.stage,col=temp.col.vec,cex=.3,names=p.values)
    #temp.boxplot=boxplot2(x=re.ordered.cor.list,las=2,main=development.stage,col=temp.re.ordered.col.vec,names = re.ordered.p.values,cex=.3)
    temp.boxplot=boxplot(x=re.ordered.cor.list,las=2,main='',col=temp.re.ordered.col.vec,names = re.ordered.p.values,cex=.8)
    #plot.new()
    #legend('center',legend=rep('',length(stages.name.str)),fill=stages.col.str,cex=1.5,box.lty=0)
  }
}
plot.all.cells.average.cell.to.pop.cor.per.timepoint=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,
                                                              pop.meta.df,
                                                              in.col.ramp=colorRampPalette(c("gray",'red'))(1000)){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df)
  pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = pop.meta.df)
  sc.timepoints=c(sc.meta.df$development.stage)
  pop.timepoints=c(pop.meta.df$development.stage)
  sc.meta.list=split(x=sc.meta.df,f=sc.meta.df$development.stage)
  pop.meta.list=split(x=pop.meta.df,f=pop.meta.df$development.stage)
  development.stages=intersect(names(sc.meta.list),names(pop.meta.list))
  #ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  development.stages=intersect(ordered.development.stages,development.stages)
  stages.col.list=get.col.factor(col.factor =development.stages )
  stages.col.str=stages.col.list$legend.col
  stages.name.str=stages.col.list$legend.str
  len.development.stages=length(development.stages)
  out.mat=matrix(nrow = len.development.stages,ncol = len.development.stages)
  for(m in 1:len.development.stages){
    development.stage=development.stages[m]
    temp.sc.meta.df=sc.meta.list[[development.stage]]
    temp.sc.rpkm.df=filter.none.expressed.genes(sc.rpkm.df[,rownames(temp.sc.meta.df)])
    no.sc.samples=dim(temp.sc.rpkm.df)[2]
    if(dim(temp.sc.rpkm.df)[2]<2){
      out.mat[m,]=rep(0,times = len.development.stages)
      next
    }
    temp.aggregate.sc.rpkm.df=data.frame(average.cell=as.numeric(apply(temp.sc.rpkm.df,1,mean)),
                                         row.names = rownames(temp.sc.rpkm.df))
    for(n in 1:len.development.stages){
      temp.pop.stage=development.stages[n]
      temp.pop.meta.df=pop.meta.list[[temp.pop.stage]]
      temp.pop.rpkm.df=subset(pop.rpkm.df,select=rownames(temp.pop.meta.df))
      intersect.genes.vec=intersect(rownames(temp.pop.rpkm.df),rownames(temp.aggregate.sc.rpkm.df))
      temp.aggregate.sc.rpkm.df=temp.aggregate.sc.rpkm.df[intersect.genes.vec,]
      temp.pop.rpkm.df=temp.pop.rpkm.df[intersect.genes.vec,]
      temp.aggregate.sc.rpkm.df=data.frame(temp.aggregate.sc.rpkm.df,row.names =intersect.genes.vec )
      aggregate.rpkm.vec=as.numeric(temp.aggregate.sc.rpkm.df[,1])
      cor.vec.list=as.numeric(apply(temp.pop.rpkm.df,2,function(pop.rpkm.vec){
        cor.p.values=c()
        pop.rpkm.vec=as.numeric(pop.rpkm.vec)
        cor.coeff=as.numeric(cor(aggregate.rpkm.vec,pop.rpkm.vec,method = 'spearman'))
        cor.coeff.list=cor.test(aggregate.rpkm.vec,pop.rpkm.vec,method = 'spearman')
        temp.cor.p.value=cor.coeff.list$p.value
        return(temp.cor.p.value)
        #return(cor.coeff.list$p.value)
      }))
      title.str=paste(development.stage,'aggregate',sep=' ')
      xlab.str=paste(temp.pop.stage,'bulk',sep=' ')
      plot(cor.vec.list,pch=19,cex.points = .5,ylab = 'p.values',main=title.str,xlab=xlab.str)
      cor.vec=as.numeric(apply(temp.pop.rpkm.df,2,function(pop.rpkm.vec){
        pop.rpkm.vec=as.numeric(pop.rpkm.vec)
        cor.coeff=as.numeric(cor(aggregate.rpkm.vec,pop.rpkm.vec,method = 'spearman'))
        return(cor.coeff)
      }))
      out.mat[m,n]=mean(cor.vec)
    }
  }
  colnames(out.mat)=development.stages
  rownames(out.mat)=development.stages
  reverse.ordered.development.stages=rev(c("S" ,"ES", "T", "ET" , "LR"  , "R"))
  out.mat=out.mat[reverse.ordered.development.stages,reverse.ordered.development.stages]
  #heatmap.2(x = out.mat,trace='none',main='Average cell vs. bulk',col=bluered(10000),margins=c(13,13),dendrogram = 'none',Rowv = development.stages,Colv = development.stages,key.title = 'Spearman cor.',key.xlab = '',key.ylab = '')
  #heatmap.2(x = out.mat,trace='none',main='',col=bluered(10000),margins=c(13,13),dendrogram = 'none',Rowv = development.stages,Colv = development.stages,key.title = '',key.xlab = '',key.ylab = '',labRow =rep('',times=length(development.stages)),labCol  =rep('',times=length(development.stages)))
  col.pelette=colorRampPalette(c("darkgray","gray","white",
                                 "yellow","orange",'red'))(5)
  col.pelette=in.col.ramp
  pheatmap(mat =out.mat,color = col.pelette, 
           cluster_rows = F,cluster_cols = F,
           main='Average cell vs. bulk',cellwidth = 50,
           cellheight = 50,fontsize_row = 10,fontsize_col = 10,
           border_color = 'lightgray')
  pheatmap(mat =out.mat,color = col.pelette, 
           cluster_rows = F,cluster_cols = F,
           main='',cellwidth = 50, cellheight = 50,legend=T,show_colnames = F,
           fontsize_row = 50,fontsize_col = 10,border_color = 'lightgray')
}
pretty.heatmap=function (mat, color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                       name = "RdYlBu")))(100), kmeans_k = NA, breaks = NA, border_color = "grey60", 
          cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE, 
          cluster_cols = TRUE, clustering_distance_rows = "euclidean", 
          clustering_distance_cols = "euclidean", clustering_method = "complete", 
          clustering_callback = identity2, cutree_rows = NA, cutree_cols = NA, 
          treeheight_row = ifelse(cluster_rows, 50, 0), treeheight_col = ifelse(cluster_cols, 
                                                                                50, 0), legend = TRUE, legend_breaks = NA, legend_labels = NA, 
          annotation_row = NA, annotation_col = NA, annotation = NA, 
          annotation_colors = NA, annotation_legend = TRUE, drop_levels = TRUE, 
          show_rownames = T, show_colnames = T, main = NA, fontsize = 10, 
          fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = F, 
          number_format = "%.2f", number_color = "grey30", fontsize_number = 0.8 * 
            fontsize, gaps_row = NULL, gaps_col = NULL, labels_row = NULL, 
          labels_col = NULL, filename = NA, width = NA, height = NA, 
          silent = FALSE, ...) 
{
  if (is.null(labels_row)) {
    labels_row = rownames(mat)
  }
  if (is.null(labels_col)) {
    labels_col = colnames(mat)
  }
  mat = as.matrix(mat)
  if (scale != "none") {
    mat = scale_mat(mat, scale)
    if (is.na2(breaks)) {
      breaks = generate_breaks(mat, length(color), center = T)
    }
  }
  if (!is.na(kmeans_k)) {
    km = kmeans(mat, kmeans_k, iter.max = 100)
    mat = km$centers
    t = table(km$cluster)
    labels_row = sprintf("Cluster: %s Size: %d", names(t), 
                         t)
  }
  else {
    km = NA
  }
  if (is.matrix(display_numbers) | is.data.frame(display_numbers)) {
    if (nrow(display_numbers) != nrow(mat) | ncol(display_numbers) != 
          ncol(mat)) {
      stop("If display_numbers provided as matrix, its dimensions have to match with mat")
    }
    display_numbers = as.matrix(display_numbers)
    fmat = matrix(as.character(display_numbers), nrow = nrow(display_numbers), 
                  ncol = ncol(display_numbers))
    fmat_draw = TRUE
  }
  else {
    if (display_numbers) {
      fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), 
                    ncol = ncol(mat))
      fmat_draw = TRUE
    }
    else {
      fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
      fmat_draw = FALSE
    }
  }
  if (cluster_rows) {
    tree_row = cluster_mat(mat, distance = clustering_distance_rows, 
                           method = clustering_method)
    tree_row = clustering_callback(tree_row, mat)
    mat = mat[tree_row$order, , drop = FALSE]
    fmat = fmat[tree_row$order, , drop = FALSE]
    labels_row = labels_row[tree_row$order]
    if (!is.na(cutree_rows)) {
      gaps_row = find_gaps(tree_row, cutree_rows)
    }
    else {
      gaps_row = NULL
    }
  }
  else {
    tree_row = NA
    treeheight_row = 0
  }
  if (cluster_cols) {
    tree_col = cluster_mat(t(mat), distance = clustering_distance_cols, 
                           method = clustering_method)
    tree_col = clustering_callback(tree_col, t(mat))
    mat = mat[, tree_col$order, drop = FALSE]
    fmat = fmat[, tree_col$order, drop = FALSE]
    labels_col = labels_col[tree_col$order]
    if (!is.na(cutree_cols)) {
      gaps_col = find_gaps(tree_col, cutree_cols)
    }
    else {
      gaps_col = NULL
    }
  }
  else {
    tree_col = NA
    treeheight_col = 0
  }
  attr(fmat, "draw") = fmat_draw
  if (!is.na2(legend_breaks) & !is.na2(legend_labels)) {
    if (length(legend_breaks) != length(legend_labels)) {
      stop("Lengths of legend_breaks and legend_labels must be the same")
    }
  }
  if (is.na2(breaks)) {
    breaks = generate_breaks(as.vector(mat), length(color))
  }
  if (legend & is.na2(legend_breaks)) {
    legend = grid.pretty(range(as.vector(breaks)))
    names(legend) = legend
  }
  else if (legend & !is.na2(legend_breaks)) {
    legend = legend_breaks[legend_breaks >= min(breaks) & 
                             legend_breaks <= max(breaks)]
    if (!is.na2(legend_labels)) {
      legend_labels = legend_labels[legend_breaks >= min(breaks) & 
                                      legend_breaks <= max(breaks)]
      names(legend) = legend_labels
    }
    else {
      names(legend) = legend
    }
  }
  else {
    legend = NA
  }
  mat = scale_colours(mat, col = color, breaks = breaks)
  if (is.na2(annotation_col) & !is.na2(annotation)) {
    annotation_col = annotation
  }
  if (!is.na2(annotation_col)) {
    annotation_col = annotation_col[colnames(mat), , drop = F]
  }
  if (!is.na2(annotation_row)) {
    annotation_row = annotation_row[rownames(mat), , drop = F]
  }
  annotation = c(annotation_row, annotation_col)
  annotation = annotation[unlist(lapply(annotation, function(x) !is.na2(x)))]
  if (length(annotation) != 0) {
    annotation_colors = generate_annotation_colours(annotation, 
                                                    annotation_colors, drop = drop_levels)
  }
  else {
    annotation_colors = NA
  }
  if (!show_rownames) {
    labels_row = NULL
  }
  if (!show_colnames) {
    labels_col = NULL
  }
  gt = heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, 
                     cellheight = cellheight, treeheight_col = treeheight_col, 
                     treeheight_row = treeheight_row, tree_col = tree_col, 
                     tree_row = tree_row, filename = filename, width = width, 
                     height = height, breaks = breaks, color = color, legend = legend, 
                     annotation_row = annotation_row, annotation_col = annotation_col, 
                     annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
                     main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                     fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, 
                     number_color = number_color, gaps_row = gaps_row, gaps_col = gaps_col, 
                     labels_row = labels_row, labels_col = labels_col, ...)
  if (is.na(filename) & !silent) {
    grid.newpage()
    grid.draw(gt)
  }
  invisible(list(tree_row = tree_row, tree_col = tree_col, 
                 kmeans = km, gtable = gt))
}
plot.all.cells.average.cell.to.pop.euclidean.per.timepoint=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df)
  pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = pop.meta.df)
  sc.timepoints=c(sc.meta.df$development.stage)
  pop.timepoints=c(pop.meta.df$development.stage)
  sc.meta.list=split(x=sc.meta.df,f=sc.meta.df$development.stage)
  pop.meta.list=split(x=pop.meta.df,f=pop.meta.df$development.stage)
  development.stages=intersect(names(sc.meta.list),names(pop.meta.list))
  #ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  development.stages=intersect(ordered.development.stages,development.stages)
  stages.col.list=get.col.factor(col.factor =development.stages )
  stages.col.str=stages.col.list$legend.col
  stages.name.str=stages.col.list$legend.str
  len.development.stages=length(development.stages)
  out.mat=matrix(nrow = len.development.stages,ncol = len.development.stages)
  for(m in 1:len.development.stages){
    development.stage=development.stages[m]
    temp.sc.meta.df=sc.meta.list[[development.stage]]
    temp.sc.rpkm.df=filter.none.expressed.genes(sc.rpkm.df[,rownames(temp.sc.meta.df)])
    no.sc.samples=dim(temp.sc.rpkm.df)[2]
    if(dim(temp.sc.rpkm.df)[2]<2){
      out.mat[m,]=rep(0,times = len.development.stages)
      next
    }
    temp.aggregate.sc.rpkm.df=data.frame(average.cell=as.numeric(apply(temp.sc.rpkm.df,1,mean)),row.names = rownames(temp.sc.rpkm.df))
    for(n in 1:len.development.stages){
      temp.pop.stage=development.stages[n]
      temp.pop.meta.df=pop.meta.list[[temp.pop.stage]]
      temp.pop.rpkm.df=subset(pop.rpkm.df,select=rownames(temp.pop.meta.df))
      intersect.genes.vec=intersect(rownames(temp.pop.rpkm.df),rownames(temp.aggregate.sc.rpkm.df))
      temp.aggregate.sc.rpkm.df=temp.aggregate.sc.rpkm.df[intersect.genes.vec,]
      temp.pop.rpkm.df=temp.pop.rpkm.df[intersect.genes.vec,]
      temp.aggregate.sc.rpkm.df=data.frame(temp.aggregate.sc.rpkm.df,row.names =intersect.genes.vec )
      aggregate.rpkm.vec=as.numeric(temp.aggregate.sc.rpkm.df[,1])
      euclidean.vec=as.numeric(apply(temp.pop.rpkm.df,2,function(pop.rpkm.vec){
        pop.rpkm.vec=as.numeric(pop.rpkm.vec)
        euclead.dist=log10(sqrt(sum((aggregate.rpkm.vec-pop.rpkm.vec)^2)))
        return(euclead.dist)
      }))
      out.mat[m,n]=mean(euclidean.vec)
    }
  }
  colnames(out.mat)=development.stages
  rownames(out.mat)=development.stages
  heatmap.2(x = out.mat,trace='none',main='Average cell vs. bulk',col=bluered(10000),margins=c(13,13),dendrogram = 'none',Rowv = development.stages,Colv = development.stages,key.title = 'Spearman cor.',key.xlab = '',key.ylab = '')
  heatmap.2(x = out.mat,trace='none',main='',col=bluered(10000),margins=c(13,13),dendrogram = 'none',Rowv = development.stages,Colv = development.stages,key.title = '',key.xlab = '',key.ylab = '',labRow =rep('',times=length(development.stages)),labCol  =rep('',times=length(development.stages)))
}
get.cor.between.two.df.cor=function(first.df,second.df,in.method='spearman'){
  intersect.genes=intersect(rownames(first.df),rownames(second.df))
  first.df=data.frame(t(subset(t(first.df),select=intersect.genes)))
  first.samples=colnames(first.df)
  second.df=second.df[intersect.genes,]
  second.samples=colnames(second.df)
  combn.df=data.frame(first.df,second.df)
  colnames(combn.df)=c(first.samples,second.samples)
  cor.mat=data.frame(cor(combn.df,method = in.method))
  if(in.method=='pearson'){
    pseudo.rpkm=min(combn.df[combn.df!=0])/2
    cor.mat=data.frame(cor(log10(combn.df+pseudo.rpkm),method = in.method))
  }
  cor.out=as.numeric(as.character(melt(cor.mat[first.samples,second.samples])$value))
  return(cor.out)
}
get.euclidian.dist.between.two.df.cor=function(first.df,second.df){
  intersect.genes=intersect(rownames(first.df),rownames(second.df))
  first.df=data.frame(t(subset(t(first.df),select=intersect.genes)))
  first.samples=colnames(first.df)
  second.df=data.frame(t(subset(t(second.df),select=intersect.genes)))
  second.samples=colnames(second.df)
  combn.df=data.frame(first.df,second.df)
  colnames(combn.df)=c(first.samples,second.samples)
  pseudo.value=min(as.numeric(combn.df[combn.df!=0]))/2
  euclidian.dist.mat =as.matrix(dist(x = log10(t(combn.df)+pseudo.value),method = 'euclidean'))
  dist.out=as.numeric(as.character(melt(euclidian.dist.mat[first.samples,second.samples])$value))
  return(dist.out)
}
compare.average.cell.to.pop.cor.per.timepoint.at.diff.scs.cut.off=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  for(m in 2:5){
    par(mfrow=c(1,1))
    plot.new()
    label.str=paste('Number of SCs averaged: ',m,sep='')
    text(x = .5,y=.5,labels = label.str,cex = 2.0)
    plot.new()
    legend('center',legend=c('',''),fill = c('blue','red'),cex = 1.5,box.lty = 0)
    compare.average.cell.to.pop.cor.per.timepoint(sc.rpkm.df = sc.rpkm.df,sc.meta.df =sc.meta.df,pop.rpkm.df =pop.rpkm.df,pop.meta.df = pop.meta.df,no.scs = m  )
  }
}
run.deseq=function(counts.mat,conditions.names,pvalue.cut.off,in.deseq2=T){
  counts.mat=as.matrix(counts.mat)
  col.df=data.frame(condition=factor(conditions.names))
  rownames(col.df)=colnames(counts.mat)
  deseq.results.sorted=data.frame()
  if(in.deseq2){
    timepoints.factor=sort(unique(conditions.names))
    col.df$condition =factor(col.df$condition,levels=timepoints.factor)
    counts.dds=DESeqDataSetFromMatrix(counts.mat,colData=col.df,design=~condition)
    counts.dds.out=DESeq(counts.dds)
    #plotDispEsts(counts.dds.out)
    #counts.dds.out.rld <- rlogTransformation(counts.dds.out,)
    #   rld=rlogTransformation(object=counts.dds,blind=F)
    #counts.dds.out.vsd <- varianceStabilizingTransformation(counts.dds.out)
    #   vsd <-varianceStabilizingTransformation(object=counts.dds,blind=F) 
    #counts.dds.out.rlogMat <- assay(counts.dds.out.rld)
    #counts.dds.out.vstMat <- assay(counts.dds.out.vsd)
    #select <- order(rowMeans(counts(counts.dds.out,normalized=TRUE)),decreasing=TRUE)[1:30]
    #hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    #heatmap.2(counts(counts.dds.out,normalized=TRUE)[select,],main='Raw counts(top 30 genes)', col = hmcol,
              #Rowv = FALSE, Colv = FALSE, scale="none",
              #dendrogram="none", trace="none", margin=c(10,6))
    #heatmap.2(assay(counts.dds.out.rld)[select,],main='Log-transformed counts', col = hmcol,
              #Rowv = FALSE, Colv = FALSE, scale="none",
              #dendrogram="none", trace="none", margin=c(10, 6))
    #heatmap.2(assay(counts.dds.out.vsd)[select,],main='Variance-transformed counts', col = hmcol,
    #Rowv = FALSE, Colv = FALSE, scale="none",
    #dendrogram="none", trace="none", margin=c(10, 6))
    #distsRL <- dist(t(counts.dds.out.rlogMat))
    #mat <- as.matrix(distsRL)
    #rownames(mat) <- colnames(mat) <- with(colData(counts.dds.out),paste(condition, sep=" : "))
    #heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13),main='Euclidian dist(rlog)')
    #print(plotPCA(counts.dds.out.rld, intgroup=c("condition")))
    #print(plotPCA(counts.dds.out.vsd, intgroup=c("condition")))
    deseq.results=results(counts.dds.out)
    #deseq.results=deseq.results[which(deseq.results$padj<=pvalue.cut.off),]
    deseq.results.sorted=deseq.results[order(deseq.results[,'padj']),]
    #deseq.results.sorted=deseq.results.sorted
    deseq.results.sorted=data.frame(deseq.results.sorted)
  }
  else{
    deseq.timepoints.conditions=conditions.names
    deseq.conditions.names.sorted=sort(unique(deseq.timepoints.conditions))
    counts.dds=newCountDataSet(countData=counts.mat,conditions=deseq.timepoints.conditions)
    counts.dds <- estimateSizeFactors(counts.dds)
    counts.dds <- estimateDispersions(counts.dds,method='blind',sharingMode='fit-only',fitType='local')
    #counts.dds <- estimateDispersions(counts.dds,fitType='local')
    counts.dds.out=nbinomTest(counts.dds,deseq.conditions.names.sorted[1],deseq.conditions.names.sorted[2])
    deseq.results=counts.dds.out
    rownames(deseq.results)=as.character(counts.dds.out[,'id'])
    deseq.results=deseq.results[deseq.results$padj<=pvalue.cut.off,]
    deseq.results.sorted=deseq.results[order(deseq.results[,'padj']),]
    deseq.results.sorted=deseq.results.sorted
    deseq.results.sorted=data.frame(deseq.results.sorted)
  }
  return(deseq.results.sorted)
}
run.pop.deseq=function(pop.counts.df,meta.df,others=F,deseq2=T,return.sig.genes.only=F,parental.development=F,log2FoldChange=2){
  pop.counts.df=filter.none.expressed.genes(pop.counts.df)
  samples=as.character(colnames(pop.counts.df))
  meta.df=meta.df[samples,]
  meta.df.list =split(meta.df,f=meta.df$development.stage)
  if(parental.development){
    meta.df.list =split(meta.df,f=meta.df$parent.development.stage)
  }
  timepoints=as.character(names(meta.df.list))
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  #ordered.development.stages=c("late.ring",  "late.trophozoite" ,  "schizont")
  timepoints=intersect(timepoints,ordered.development.stages)
  out.list=list()
  no.of.timepoints=length(timepoints)
  out.df.list=list()
  if(others){
    for(m in 1:no.of.timepoints){
      timepoint=timepoints[m]
      timepoint.samples=rownames(meta.df.list[[timepoint]])
      other.samples=setdiff(samples,timepoint.samples)
      timepoint.samples.df=filter.none.expressed.samples(filter.none.expressed.genes(subset(pop.counts.df,select=timepoint.samples)))
      other.samples.df=filter.none.expressed.samples(filter.none.expressed.genes(subset(pop.counts.df,select=other.samples)))
      len.timepoint.one.samples=dim(timepoint.samples.df)[2]
      len.timepoint.two.samples=dim(other.samples.df)[2]
      if (!all(c(len.timepoint.one.samples>=3,len.timepoint.two.samples>=3))){
        next
      }
      intersect.genes=intersect(rownames(timepoint.samples.df),rownames(other.samples.df))
      timepoint.samples.df=timepoint.samples.df[intersect.genes,]
      other.samples.df=other.samples.df[intersect.genes,]
      diff.df=data.frame(timepoint.samples.df,other.samples.df)
      timepoint.samples=colnames(timepoint.samples.df)
      other.samples=colnames(other.samples.df)
      #diff.df=filter.none.expressed.genes(diff.df)
      deseq.factor=c(rep(timepoint,times=length(timepoint.samples)),rep('control',times=length(other.samples)))
      timepoints.comp=paste(c(timepoint,'others'),collapse='_vs_')
      if(dim(diff.df)[1]>10 & dim(diff.df)[2]>2){
        deseq.res=try(run.deseq(counts.mat=diff.df,conditions.names=deseq.factor,pvalue.cut.off=0.05,in.deseq2 = T))
        if(class(deseq.res)=="try-error") {
          #out.list[[timepoints.comp]]=='Error'
          next; 
        }
        if(return.sig.genes.only){
          deseq.res=data.frame(deseq.res[which(deseq.res$padj<=.05 & abs(deseq.res$log2FoldChange)>=log2FoldChange),])
          regulation=as.numeric(as.character(deseq.res$log2FoldChange))
          regulation=ifelse(regulation>=0,'up','down')
          deseq.res$regulation=regulation
        }
        out.list[[timepoints.comp]]=deseq.res
      }
    }
  }
  else{
    timepoints.combn=combn(x=timepoints,m=2)
    no.timepoints.combn=dim(timepoints.combn)[2]
    for (n in 1:no.timepoints.combn){
      timepoint.combn=timepoints.combn[,n]
      temp.timepoint.one.meta.df=meta.df.list[[timepoint.combn[1]]]
      temp.timepoint.two.meta.df=meta.df.list[[timepoint.combn[2]]]
      timepoint.one.samples=rownames(temp.timepoint.one.meta.df)
      timepoint.two.samples=rownames(temp.timepoint.two.meta.df)
      len.timepoint.one.samples=length(timepoint.one.samples)
      len.timepoint.two.samples=length(timepoint.two.samples)
      #if (!all(c(len.timepoint.one.samples>=2,len.timepoint.two.samples>=2))){
        #next
      #}
      if (!all(c(len.timepoint.one.samples>=3,len.timepoint.two.samples>=3))){
        next
      }
      timepoint.one.counts.df=subset(pop.counts.df,select=timepoint.one.samples)
      timepoint.two.counts.df=subset(pop.counts.df,select=timepoint.two.samples)
      temp.diff.df=data.frame(timepoint.one.counts.df,timepoint.two.counts.df)
      temp.diff.df=filter.none.expressed.genes(temp.diff.df)
      deseq.factor=c(rep(timepoint.combn[1],times=length(timepoint.one.samples)),rep(timepoint.combn[2],times=length(timepoint.two.samples)))
      if(dim(temp.diff.df)[1]>10 & dim(temp.diff.df)[2]>=2){
        timepoints.comp=paste(c(timepoint.combn[1],timepoint.combn[2]),collapse='_vs_')
        deseq.res=try(run.deseq(counts.mat=temp.diff.df,conditions.names=deseq.factor,pvalue.cut.off=0.05))
        if(class(deseq.res)=="try-error") {
          #out.list[[timepoints.comp]]=='Error'
          next; 
        }
        if(return.sig.genes.only){
          deseq.res=data.frame(deseq.res[which(deseq.res$padj<=.05 & abs(deseq.res$log2FoldChange)>=log2FoldChange),])
          regulation=as.numeric(as.character(deseq.res$log2FoldChange))
          regulation=ifelse(regulation>=0,'up','down')
          deseq.res$regulation=regulation
        }
        out.list[[timepoints.comp]]=deseq.res
        #temp.list=list(counts=temp.diff.df,conditions.names=deseq.factor)
        #out.df.list[[timepoints.comp]]=temp.list
      }
    }
  }
  return(out.list)
}
run.pop.deseq.ordered=function(pop.counts.df,meta.df,others=F,deseq2=T,return.sig.genes.only=F,parental.development=F,log2FoldChange=2){
  pop.counts.df=filter.none.expressed.genes(pop.counts.df)
  samples=as.character(colnames(pop.counts.df))
  meta.df=meta.df[samples,]
  meta.df.list =split(meta.df,f=meta.df$development.stage)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  ordered.development.stages=intersect(ordered.development.stages,names(meta.df.list))
  if(parental.development){
    meta.df.list =split(meta.df,f=meta.df$parent.development.stage)
  }
  timepoints=as.character(names(meta.df.list))
  out.list=list()
  no.of.timepoints=length(timepoints)
  out.df.list=list()
  if(others){
    for(m in 1:no.of.timepoints){
      timepoint=timepoints[m]
      timepoint.samples=rownames(meta.df.list[[timepoint]])
      other.samples=setdiff(samples,timepoint.samples)
      timepoint.samples.df=filter.none.expressed.samples(filter.none.expressed.genes(subset(pop.counts.df,select=timepoint.samples)))
      other.samples.df=filter.none.expressed.samples(filter.none.expressed.genes(subset(pop.counts.df,select=other.samples)))
      len.timepoint.one.samples=dim(timepoint.samples.df)[2]
      len.timepoint.two.samples=dim(other.samples.df)[2]
      if (!all(c(len.timepoint.one.samples>=3,len.timepoint.two.samples>=3))){
        next
      }
      intersect.genes=intersect(rownames(timepoint.samples.df),rownames(other.samples.df))
      timepoint.samples.df=timepoint.samples.df[intersect.genes,]
      other.samples.df=other.samples.df[intersect.genes,]
      diff.df=data.frame(timepoint.samples.df,other.samples.df)
      timepoint.samples=colnames(timepoint.samples.df)
      other.samples=colnames(other.samples.df)
      #diff.df=filter.none.expressed.genes(diff.df)
      deseq.factor=c(rep(timepoint,times=length(timepoint.samples)),rep('control',times=length(other.samples)))
      timepoints.comp=paste(c(timepoint,'others'),collapse='_vs_')
      if(dim(diff.df)[1]>10 & dim(diff.df)[2]>2){
        deseq.res=try(run.deseq(counts.mat=diff.df,conditions.names=deseq.factor,pvalue.cut.off=0.05,in.deseq2 = T))
        if(class(deseq.res)=="try-error") {
          #out.list[[timepoints.comp]]=='Error'
          next; 
        }
        if(return.sig.genes.only){
          deseq.res=data.frame(deseq.res[which(deseq.res$padj<=.05 & abs(deseq.res$log2FoldChange)>=log2FoldChange),])
          regulation=as.numeric(as.character(deseq.res$log2FoldChange))
          regulation=ifelse(regulation>=0,'up','down')
          deseq.res$regulation=regulation
        }
        out.list[[timepoints.comp]]=deseq.res
      }
    }
  }
  else{
    timepoints.combn=combn(x=timepoints,m=2)
    no.timepoints.combn=dim(timepoints.combn)[2]
    for (n in 1:no.timepoints.combn){
      timepoint.combn=timepoints.combn[,n]
      temp.timepoint.one.meta.df=meta.df.list[[timepoint.combn[1]]]
      temp.timepoint.two.meta.df=meta.df.list[[timepoint.combn[2]]]
      timepoint.one.samples=rownames(temp.timepoint.one.meta.df)
      timepoint.two.samples=rownames(temp.timepoint.two.meta.df)
      len.timepoint.one.samples=length(timepoint.one.samples)
      len.timepoint.two.samples=length(timepoint.two.samples)
      #if (!all(c(len.timepoint.one.samples>=2,len.timepoint.two.samples>=2))){
      #next
      #}
      if (!all(c(len.timepoint.one.samples>=3,len.timepoint.two.samples>=3))){
        next
      }
      timepoint.one.counts.df=subset(pop.counts.df,select=timepoint.one.samples)
      timepoint.two.counts.df=subset(pop.counts.df,select=timepoint.two.samples)
      temp.diff.df=data.frame(timepoint.one.counts.df,timepoint.two.counts.df)
      temp.diff.df=filter.none.expressed.genes(temp.diff.df)
      deseq.factor=c(rep(timepoint.combn[1],times=length(timepoint.one.samples)),rep(timepoint.combn[2],times=length(timepoint.two.samples)))
      if(dim(temp.diff.df)[1]>10 & dim(temp.diff.df)[2]>=2){
        timepoints.comp=paste(c(timepoint.combn[1],timepoint.combn[2]),collapse='_vs_')
        deseq.res=try(run.deseq(counts.mat=temp.diff.df,conditions.names=deseq.factor,pvalue.cut.off=0.05))
        if(class(deseq.res)=="try-error") {
          #out.list[[timepoints.comp]]=='Error'
          next; 
        }
        if(return.sig.genes.only){
          deseq.res=data.frame(deseq.res[which(deseq.res$padj<=.05 & abs(deseq.res$log2FoldChange)>=log2FoldChange),])
          regulation=as.numeric(as.character(deseq.res$log2FoldChange))
          regulation=ifelse(regulation>=0,'up','down')
          deseq.res$regulation=regulation
        }
        out.list[[timepoints.comp]]=deseq.res
        #temp.list=list(counts=temp.diff.df,conditions.names=deseq.factor)
        #out.df.list[[timepoints.comp]]=temp.list
      }
    }
  }
  return(out.list)
}
get.sig.genes.from.deseq.res=function(deseq.res.list,padj.cut.off=0.05,log.fold.cut.off=2,return.genes.only=F){
  comparison.names=names(deseq.res.list)
  len.comparison.names=length(comparison.names)
  out.list=list()
  for(n in 1:len.comparison.names){
    comparison.name=comparison.names[n]
    temp.deseq.df=deseq.res.list[[comparison.name]]
    #temp.sig.genes.df=subset(temp.deseq.df, padj<=padj.cut.off)
    temp.sig.genes.df=subset(temp.deseq.df, padj<=padj.cut.off & abs(log2FoldChange) >=log.fold.cut.off)
    out.list[[comparison.name]]=temp.sig.genes.df
    if(return.genes.only){
      out.list[[comparison.name]]=rownames(temp.sig.genes.df)
    }
  }
  return(out.list)
}
get.stage.specific.genes.from.deseq.sig.genes.res=function(deseq.res.list){
  comparison.names=names(deseq.res.list)
  len.comparison.names=length(comparison.names)
  out.list=list()
  development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  out.list=list()
  len.development.stages=length(development.stages)
  for(m in 1:len.development.stages){
    temp.stage=development.stages[m]
    other.stages=setdiff(development.stages,temp.stage)
    len.other.stages=length(other.stages)
    temp.list=list()
    for(n in 1:len.other.stages){
      other.stage=other.stages[n]
      pair.name=paste(sort(c(temp.stage,other.stage)),collapse = '_vs_')
      temp.list[[n]]=deseq.res.list[[pair.name]]
    }
    out.list[[temp.stage]]=Reduce(intersect,temp.list)
  }
  return(out.list)
}
cluster.sc.based.on.pop.diff.expressed.genes =function(sc.rpkm.df,sc.meta.df,pop.count.df,pop.meta.df,foldch=0,pvalue=.05,parental.development=F){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  pop.meta.df=pop.meta.df[colnames(pop.count.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$development.stage)
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$development.stage)
  intersect.stages=intersect(c(sc.meta.df$development.stage),c(pop.meta.df$development.stage))
  if(parental.development){
    sc.meta.list=split(sc.meta.df,f=sc.meta.df$parental.development.stage)
    pop.meta.list=split(pop.meta.df,f=pop.meta.df$parental.development.stage)
    intersect.stages=intersect(c(sc.meta.df$parental.development.stage),c(pop.meta.df$parental.development.stage))
  }
  intersect.stages.combn=combn(intersect.stages,2)
  no.pairs=dim(intersect.stages.combn)[2]
  for (m in 1:no.pairs){
    temp.pair=intersect.stages.combn[,m]
    first.stage=temp.pair[1]
    sec.stage=temp.pair[2]
    first.pop.count.df=filter.none.expressed.samples(filter.none.expressed.genes(pop.count.df[,rownames(pop.meta.list[[first.stage]])]))
    sec.pop.count.df=filter.none.expressed.samples(filter.none.expressed.genes(pop.count.df[,rownames(pop.meta.list[[sec.stage]])]))
    category.names=c(rep(first.stage,times = dim(first.pop.count.df)[2]),rep(sec.stage,times = dim(sec.pop.count.df)[2]))
    first.pop.count.df=filter.none.expressed.samples(filter.none.expressed.genes(first.pop.count.df))
    sec.pop.count.df=filter.none.expressed.samples(filter.none.expressed.genes(sec.pop.count.df))
    intersect.genes=intersect(rownames(first.pop.count.df),rownames(sec.pop.count.df))
    first.pop.count.df=first.pop.count.df[intersect.genes,]
    sec.pop.count.df=sec.pop.count.df[intersect.genes,]
    diff.df=data.frame(first.pop.count.df, sec.pop.count.df)
    deseq.out.df=try(run.deseq(counts.mat =diff.df,conditions.names = category.names,in.deseq2 = T  ))
    if(class(deseq.out.df)=='try-error'){
      next
    }
    sig.genes.df=deseq.out.df[which(deseq.out.df$pvalue<=pvalue & abs(deseq.out.df$log2FoldChange)>=foldch),]
    #sig.genes.df=deseq.out.df[which(deseq.out.df$pvalue<=pvalue),]
    #show(sort(as.numeric(sig.genes.df$log2FoldChange)))
    sc.temp.meta.df=rbind(sc.meta.list[[first.stage]],sc.meta.list[[sec.stage]])
    sc.rpkm.temp.df=sc.rpkm.df[,rownames(sc.temp.meta.df)]
    intersect.genes=intersect(rownames(sig.genes.df),rownames(sc.rpkm.temp.df))
    sc.rpkm.temp.df=sc.rpkm.temp.df[intersect.genes,]
    sc.rpkm.temp.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = sc.rpkm.temp.df))
    #sc.rpkm.temp.df=filter.genes.not.expressed.in.all.samples(df = sc.rpkm.temp.df)
    sc.temp.meta.df=sc.temp.meta.df[colnames(sc.rpkm.temp.df),]
    title.str=paste(c(first.stage,sec.stage),collapse = ' vs ')
    title.str=gsub(pattern = '\\.',replacement = ' ',x = title.str)
    #Get different clusters
    #Spearman
    samples.spearman.cor=cor(sc.rpkm.temp.df,method = 'spearman')
    samples.spearman.cor.dist=as.dist(1-samples.spearman.cor)
    samples.spearman.cor.dist.hclust=hclust(d = samples.spearman.cor.dist,method = 'complete')
    samples.spearman.cor.dist.hclust.dendo=as.dendrogram(object = samples.spearman.cor.dist.hclust)
    #Pearson
    samples.pearson.cor=cor(sc.rpkm.temp.df,method = 'pearson')
    samples.pearson.cor.dist=as.dist(1-samples.pearson.cor)
    samples.pearson.cor.dist.hclust=hclust(d = samples.pearson.cor.dist,method = 'complete')
    samples.pearson.cor.dist.hclust.dendo=as.dendrogram(object = samples.pearson.cor.dist.hclust)
    #Genes cluster
    #Spearman
    genes.spearman.cor=cor(t(sc.rpkm.temp.df),method = 'spearman')
    genes.spearman.cor.dist=as.dist(1-genes.spearman.cor)
    genes.spearman.cor.dist.hclust=hclust(d = genes.spearman.cor.dist,method = 'average')
    genes.spearman.cor.dist.hclust.dendro=as.dendrogram(object =genes.spearman.cor.dist.hclust )
    #pearson
    genes.pearson.cor=cor(t(sc.rpkm.temp.df),method = 'pearson')
    genes.pearson.cor.dist=as.dist(1-genes.pearson.cor)
    genes.pearson.cor.dist.hclust=hclust(d = genes.pearson.cor.dist,method = 'average')
    genes.pearson.cor.dist.hclust.dendro=as.dendrogram(object =genes.pearson.cor.dist.hclust )
    #temp.heatmap=try(plot.samples.correlations.heatmap(in.rpkm.df =  sc.rpkm.temp.df,title.str =title.str,in.meta.df = sc.temp.meta.df,filter.non.var.rows = T ))
    #if(class(temp.heatmap)=='try-error'){
      #next
    #}
    samples.col.factor=as.character(sc.temp.meta.df$development.stage)
    if(parental.development){
      samples.col.factor=as.character(sc.temp.meta.df$parental.development.stage)
    }
    samples.col.factor.list=get.col.factor(col.factor = samples.col.factor)
    samples.col.factor.str=samples.col.factor.list$col.str
    genes.col.factor=as.numeric(as.character(sig.genes.df[rownames(sc.rpkm.temp.df),]$log2FoldChange))
    genes.col.factor=ifelse(genes.col.factor>0,'darkred','darkgreen')
    rpkm.mat=as.matrix(sc.rpkm.temp.df)
    log.rpkm.mat=log10(rpkm.mat+1)
    temp.heatmap=heatmap.2(x =log.rpkm.mat,RowSideColors = genes.col.factor,ColSideColors = samples.col.factor.str,main = title.str,trace='none',Colv = samples.spearman.cor.dist.hclust.dendo,Rowv=genes.spearman.cor.dist.hclust.dendro,labCol = '',labRow = '')
    legend('topright',legend=samples.col.factor.list$legend.str,fill=samples.col.factor.list$legend.col)
  }
}
#Encodes the parental development stage
encode.parental.development.stage=function(meta.df){
  development.stages=as.character(meta.df$development.stage)
  len.development.stages=length(development.stages)
  parental.stages=c()
  for(m in 1:len.development.stages){
    development.stage=development.stages[m]
    if(grepl('ring',development.stage)){
      parental.stages=append(parental.stages,'ring',length(parental.stages))
    }
    if(grepl('trophozoite',development.stage)){
      parental.stages=append(parental.stages,'trophozoite',length(parental.stages))
    }
    if(grepl('schizont',development.stage)){
      parental.stages=append(parental.stages,'schizont',length(parental.stages))
    }
  }
  meta.df$parental.development.stage=parental.stages
  return(meta.df)
}
plot.sc.diff.expressed.genes.heatmap=function(sc.rpkm.df,sc.meta.df,deseq.res.list){
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$development.stage)
  sc.categories=names(sc.meta.list)
  deseq.res.list.names=names(deseq.res.list)
  #deseq.res.list.names=deseq.res.list.names[1:2]
  len.deseq.res.list.names=length(deseq.res.list.names)
  show(deseq.res.list.names)
  for(m in 1:len.deseq.res.list.names){
    deseq.res.list.name=deseq.res.list.names[m]
    temp.deseq.res.df=deseq.res.list[[deseq.res.list.name]]
    categories.parts=strsplit(deseq.res.list.name,split = '_vs_')
    first.grp=categories.parts[[1]][1]
    sec.grp=categories.parts[[1]][2]
    if(!all(first.grp %in% sc.categories & sec.grp %in% sc.categories)){
      next
    }
    first.sc.meta.df=sc.meta.list[[first.grp]]
    sec.sc.meta.df=sc.meta.list[[sec.grp]]
    first.sc.rpkm.df=sc.rpkm.df[,rownames(first.sc.meta.df)]
    sec.sc.rpkm.df=sc.rpkm.df[,rownames(sec.sc.meta.df)]
    combined.sc.rpkm.df=data.frame(first.sc.rpkm.df,sec.sc.rpkm.df)
    combined.sc.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(combined.sc.rpkm.df))
    combined.sc.meta.df=sc.meta.df[colnames(combined.sc.rpkm.df),]
    samples.col.factor=as.character(combined.sc.meta.df$development.stage)
    samples.col.factor.list=get.col.factor(col.factor = samples.col.factor)
    intersect.genes=intersect(rownames(temp.deseq.res.df),rownames(combined.sc.rpkm.df))
    temp.deseq.res.df=temp.deseq.res.df[intersect.genes,]
    temp.deseq.res.df=temp.deseq.res.df[order(temp.deseq.res.df$regulation),]
    combined.sc.rpkm.df=combined.sc.rpkm.df[intersect.genes,]
    if(dim(combined.sc.rpkm.df)[1]<5){
      next
    }
    #show(head(temp.deseq.res.df,100))
    temp.deseq.res.df=temp.deseq.res.df[intersect.genes,]
    temp.deseq.res.df=temp.deseq.res.df[order(temp.deseq.res.df$regulation),]
    genes.regulation=as.character(temp.deseq.res.df[,'regulation'])
    #genes.regulation=as.character(temp.deseq.res.df$regulation)
    gene.cols=ifelse(genes.regulation=='up','darkred','darkgreen')
    combined.sc.rpkm.mat=as.matrix(log2(combined.sc.rpkm.df+.00000001))
    #combined.sc.rpkm.mat=as.matrix(combined.sc.rpkm.df)
    #combined.sc.rpkm.df=filter.non.variable.rows(df = combined.sc.rpkm.df,cut.off = 1)
    temp.no.samples=dim(combined.sc.rpkm.df)[2]
    samples.cor.mat=cor(combined.sc.rpkm.df,method = 'spearman')
    title.str=gsub(pattern = '\\.',' ', x = deseq.res.list.name)
    title.str=gsub(pattern = '_vs_',' vs ', x = title.str)
    #heatmap.2(x = combined.sc.rpkm.mat,trace = 'none',col=greenred(100),margins = c(10,10),ColSideColors = samples.col.factor.list$col.str,RowSideColors = gene.cols,labCol = '',labRow = '',cexRow = .3,cexCol = .3,cex=.2,dendrogram = 'col',main=title.str)
    samples.col.factor.legend.str=samples.col.factor.list$legend.str
    samples.col.factor.legend.str=gsub(pattern = '\\.',replacement = ' ',samples.col.factor.legend.str)
    samples.col.factor.legend.col=samples.col.factor.list$legend.col
    samples.cor.dendro=as.dendrogram(hclust(as.dist(samples.cor.mat)))
    #heatmap.2(x = samples.cor.mat,trace = 'none',margins = c(10,10),ColSideColors = samples.col.factor.list$col.str,labCol = '',labRow = '',cexRow = .3,cexCol = .3,cex=.2,dendrogram = 'col',main=title.str,key.xlab = '',key.ylab = '',key.title = 'Correlation key',col=bluered(1000),Colv=samples.cor.dendro)
    heatmap.2(x = combined.sc.rpkm.mat,trace = 'none',col=bluered(1000),margins = c(13,13),ColSideColors = samples.col.factor.list$col.str,RowSideColors = gene.cols,cexRow = .3,cexCol = .3, main=deseq.res.list.name,cex=.2,dendrogram = 'col' ,Rowv = F)
    legend('topright',legend=samples.col.factor.legend.str,fill=samples.col.factor.legend.col,box.lty = 0,cex=1.5)
    heatmap.2(x = combined.sc.rpkm.mat,trace = 'none',col=bluered(1000),margins = c(13,13),ColSideColors = samples.col.factor.list$col.str,RowSideColors = gene.cols,cexRow = .3,cexCol = .3, main='',cex=.2,key.title = '',key.xlab = '',key.ylab = '',labRow = 'none')
    #plot.new()
    #plot.new()
    #legend('center',legend=samples.col.factor.legend.str,fill=samples.col.factor.legend.col,box.lty = 0,cex=1.5)
    #plot.new()
    #legend('center',legend=rep('',length(samples.col.factor.legend.str)),fill=samples.col.factor.legend.col,box.lty = 0,cex=1.5)
    #plot.new()
    #legend('top',legend=c('up-regulated','down-regulated'),fill=c('darkred','darkgreen'),box.lty = 0,cex=1.5)
    #plot.new()
    #legend('bottom',legend=rep('',2),fill=c('darkred','darkgreen'),box.lty = 0,cex=1.5)
    #samples.cor.dist=as.dist(1- samples.cor.mat)
    #samples.cor.hclust=hclust(samples.cor.dist,method = 'average')
    #samples.cor.dendrogram=as.dendrogram(object = samples.cor.hclust)
    #genes.cor.mat=cor(t(combined.sc.rpkm.df),method='spearman')
    #genes.cor.mat=cor(t(log10(combined.sc.rpkm.df+1)))
    #genes.cor.dist=as.dist(1-genes.cor.mat)
    #genes.cor.dist=as.dist(genes.cor.mat)
    #genes.cor.hclust=hclust(genes.cor.dist,method = 'complete')
    #genes.cor.dendrogram =as.dendrogram(genes.cor.hclust)
  }
}
temp.plot.aggregated.samples.pca=function(aggregates.list,first.pc='PC1',second.pc='PC2',title.str='Test',log.samples.pca.scores=F){
  average.categories=names(aggregates.list)
  len.average.categories=length(average.categories)
  for(m in 1:len.average.categories){
    average.category=average.categories[m]
    if(average.category=='average.scs.2'){
      first.pc='PC2'
      second.pc='PC3'
    }
    if(average.category=='average.scs.3'){
      first.pc='PC2'
      second.pc='PC3'
    }
    if(average.category=='average.scs.4'){
      first.pc='PC2'
      second.pc='PC3'
    }
    if(average.category=='average.scs.5'){
      first.pc='PC1'
      second.pc='PC2'
    }
    aggregate.rpkm.df=aggregates.list[[average.category]]
    aggregate.rpkm.df=filter.non.variable.rows(df =aggregate.rpkm.df,cut.off = 1 )
    stages.label=colnames(aggregate.rpkm.df)
    stages.label=gsub(pattern = 'sample_.*.',replacement = '',stages.label)
    stages.label=gsub(pattern = '\\.$',replacement = '',stages.label)
    col.factor=stages.label
    col.factor.list=get.col.factor(col.factor = col.factor)
    samples.pca=prcomp(t(aggregate.rpkm.df),retx=T,center=T,scale.=T)
    samples.pca.summary=summary(samples.pca)$importance
    first.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,first.pc]))*100,2)
    second.pc.explained.var=round(as.numeric(as.character(samples.pca.summary[2,second.pc]))*100,2)
    samples.pca.scores=samples.pca$x
    x.points=as.numeric(as.character(samples.pca.scores[,first.pc]))
    y.points=as.numeric(as.character(samples.pca.scores[,second.pc]))
    if(log.samples.pca.scores){
      x.points=log2(x.points+1)
      y.points=log2(y.points+1)
    }
    pca.col.list=get.col.factor(col.factor = col.factor)
    pca.col.str=pca.col.list[['col.str']]
    first.pc.explained.var=paste('(Explained var: ',paste(first.pc.explained.var,'%)',sep=''),sep='')
    second.pc.explained.var=paste('(Explained var: ',paste(second.pc.explained.var,'%)',sep=''),sep='')
    first.pc.lab=paste(first.pc,first.pc.explained.var,'')
    second.pc.lab=paste(second.pc,second.pc.explained.var,'')
    x.lim=c(min(x.points,y.points),max(x.points,y.points))
    y.lim=c(min(x.points,y.points),max(x.points,y.points))
    plot(x=x.points,y=y.points,col=pca.col.str,xlab=first.pc.lab,ylab=second.pc.lab,main='',pch=19,cex=2)
    #legend('topright',legend = pca.col.list$legend.str,fill=pca.col.list$legend.cols,box.lty = 0,cex=1.5)
    #pairs(samples.pca.scores[,1:5],main=average.category,col=pca.col.str,pch=19,cex=.8)
  }
}
add.ercc.row.to.df=function(genes.count.df,ercc.counts.df){
  samples=colnames(genes.count.df)
  len.samples=length(samples)
  ercc.samples=colnames(ercc.counts.df)
  ercc.ids=rownames(ercc.counts.df)
  len.ercc.ids=length(ercc.ids)
  out.genes.count.df=genes.count.df
  for(m in 1:len.ercc.ids){
    ercc.id=ercc.ids[m]
    ercc.counts.vec=c()
    for(n in 1:len.samples){
      sample=samples[n]
      temp.ercc.count=as.numeric(ercc.counts.df[ercc.id,sample])
      ercc.count=ifelse(length(temp.ercc.count)>0,temp.ercc.count,0)
      ercc.counts.vec=append(ercc.counts.vec,ercc.count,length(ercc.counts.vec))
    }
    out.genes.count.df[ercc.id,]=ercc.counts.vec
  }
  return(out.genes.count.df)
}
get.stage.specific.deseq.normalized.counts=function(counts.df,meta.df,method=median){
  counts.df=filter.none.expressed.genes(input.data = counts.df)
  all.genes=rownames(counts.df)
  meta.df=meta.df[colnames(counts.df),]
  meta.list=split(meta.df,f =meta.df$development.stage)
  development.stages = names(meta.list)
  ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "late.trophozoite" , "early.schizont"  , "schizont")
  development.stages = intersect(ordered.development.stages,development.stages)
  no.development.stage=length(development.stages)
  out.list=list()
  colnames.vec=c()
  for (i in 1:no.development.stage){
    development.stage=development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    if(dim(temp.meta.df)[1]<2){
      next
    }
    temp.counts.df=counts.df[,rownames(temp.meta.df)]
    temp.counts.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = temp.counts.df))
    show(dim(temp.counts.df))
    temp.meta.df=temp.meta.df[colnames(temp.counts.df),]
    size.factors=estimateSizeFactorsForMatrix(counts =as.matrix(temp.counts.df),locfunc = method)
    temp.rpkm.df=data.frame(t(t(temp.counts.df)/size.factors))
    none.expressed.genes=setdiff(all.genes,rownames(temp.rpkm.df))
    temp.df=data.frame(matrix(0,nrow =length(none.expressed.genes),ncol = dim(temp.rpkm.df)[2] ),row.names = none.expressed.genes)
    colnames(temp.df)=colnames(temp.rpkm.df)
    temp.rpkm.df=rbind(temp.rpkm.df,temp.df)
    out.list[[development.stage]]=temp.rpkm.df
    colnames.vec=append(colnames.vec,colnames(temp.rpkm.df),after = length(colnames.vec))
  }
  out.df=as.data.frame(out.list)
  colnames(out.df)=colnames.vec
  return(out.df)
}
get.deseq.normalized.counts=function(counts.df,method=median){
  counts.df=filter.none.expressed.genes(input.data = counts.df)
  size.factors=estimateSizeFactorsForMatrix(counts =as.matrix(counts.df),locfunc = method)
  temp.rpkm.df=data.frame(t(t(counts.df)/size.factors))
  return(temp.rpkm.df)
}
annotate.monocle.genes =function(rpkm.df){
  bio.type.char=as.character(unlist(apply(rpkm.df,1,function(rpkm.row){
    rpkm.row='protein_coding'
  })))
  out.df=data.frame(gene.id=rownames(rpkm.df),bio.type=bio.type.char,gene_short_name=rownames(rpkm.df))
  rownames(out.df)=rownames(rpkm.df)
  return(out.df)
}
run.monocle=function(rpkm.df,meta.df){
  rpkm.mat=as.matrix(rpkm.df)
  sample.sheet=meta.df[colnames(rpkm.df),]
  gene.ann=annotate.monocle.genes(rpkm.df=rpkm.df)
  pd <- new("AnnotatedDataFrame", data = sample.sheet)
  fd <- new("AnnotatedDataFrame", data = gene.ann)
  monocle.dataset <- new("CellDataSet", exprs = rpkm.mat, phenoData = pd, featureData = fd)
  diff.test.res <- differentialGeneTest(monocle.dataset,fullModelFormulaStr="expression~development.stage")
  ordering.genes <- row.names(subset(diff.test.res, qval < 0.01))
  diff.test.res <- setOrderingFilter(diff.test.res, ordering.genes)
  diff.test.res <- reduceDimension(diff.test.res, use_irlba = F)
  diff.test.res <- orderCells(diff.test.res, num_paths = 1, reverse = F)
  return(diff.test.res)
}
run.monocle.at.different.param=function(rpkm.df,meta.df){
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df$development.stage=factor(meta.df$development.stage,levels = factor(c('R','LR','ET','LT','ES','S')))
  rpkm.cut.off=c(1,3,5,10)
  rpkm.cut.off=c(1)
  #no.samples.cut.off=c(2,3,5,10)
  no.samples.cut.off=c(2)
  len.rpkm.cut.off=length(rpkm.cut.off)
  len.no.samples.cut.off=length(no.samples.cut.off)
  for(m in 1:len.rpkm.cut.off){
    for(n in 1:len.no.samples.cut.off){
      samples.count=as.numeric(no.samples.cut.off[n])
      temp.rpkm.cut.off=as.numeric(rpkm.cut.off[m])
      pd=new('AnnotatedDataFrame',meta.df)
      fd=new('AnnotatedDataFrame',annotate.monocle.genes(rpkm.df = rpkm.df))
      HSMM <- newCellDataSet(as.matrix(rpkm.df), phenoData = pd, featureData = fd)
      #show(dim(exprs(HSMM)))
      HSMM <- detectGenes(HSMM, min_expr = temp.rpkm.cut.off)
      expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= samples.count))
      HSMM <- HSMM[expressed_genes,]
      ordering.genes=expressed_genes
      HSMM <- setOrderingFilter(HSMM, ordering.genes)
      HSMM <- reduceDimension(HSMM, use_irlba = F)
      HSMM <- orderCells(HSMM, num_paths = 1, reverse = T)
      plot.spanning.tree(HSMM,color_by = 'development.stage',path.size =1.5,point.size = 6,show_cell_names = T,cell_name_size = .3 )
    }
  }
}
create.monocle.CellDataSet=function(rpkm.df,meta.df){
  rpkm.mat=as.matrix(rpkm.df)
  sample.sheet=meta.df[colnames(rpkm.df),]
  gene.ann=annotate.monocle.genes(rpkm.df=rpkm.df)
  pd <- new("AnnotatedDataFrame", data = sample.sheet)
  fd <- new("AnnotatedDataFrame", data = gene.ann)
  monocle.dataset <- new("CellDataSet", exprs = rpkm.mat, phenoData = pd, featureData = fd)
  return(monocle.dataset)
}
run.msp.at.different.params=function(HSMM,pdf.file){
  rpkm.intervals=c(0.1,1,5,10)
  rpkm.intervals=c(10)
  no.of.cells=c(2,4,6,8,10)
  no.of.cells=c(40)
  len.rpkm.intervals=length(rpkm.intervals)
  len.no.of.cells=length(no.of.cells)
  pdf(pdf.file)
  for(m in 1:len.rpkm.intervals){
    for(n in 1:len.no.of.cells){
      #plot.new()
      lab.str=paste(c('rpkm cut.off: ',rpkm.intervals[m],'no.cells.cut.off: ',no.of.cells[n]),collapse = ' ')
      show(lab.str)
      #text(x = .5,y = .5,labels =lab.str )
      HSMM <- detectGenes(HSMM, min_expr =rpkm.intervals[m])
      expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= no.of.cells[n]))
      ordering.genes=expressed_genes
      show(length(ordering.genes))
      HSMM <- setOrderingFilter(HSMM, ordering.genes)
      HSMM <- reduceDimension(HSMM, use_irlba = F)
      HSMM <- orderCells(HSMM, num_paths = 2, reverse = T)
      plot.spanning.tree(HSMM,color_by = 'development.stage')
    }
  }
  dev.off()
}
plot.spanning.tree=function (cds, x = 1, y = 2, color_by = "State", show_tree = TRUE, 
                             show_backbone = TRUE, backbone_color = "black", markers = NULL, 
                             show_cell_names = FALSE, cell_name_size = 1,path.size=1.0,point.size=2.5) 
{
  lib_info_with_pseudo <- pData(cds)
  S_matrix <- reducedDimS(cds)
  if (is.null(S_matrix)) {
    stop("You must first call reduceDimension() before using this function")
  }
  ica_space_df <- data.frame(t(S_matrix[c(x, y), ]))
  colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
  ica_space_df$sample_name <- row.names(ica_space_df)
  ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, 
                                   by.x = "sample_name", by.y = "row.names")
  dp_mst <- minSpanningTree(cds)
  if (is.null(dp_mst)) {
    stop("You must first call orderCells() before using this function")
  }
  edge_list <- as.data.frame(get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")
  edge_df <- merge(ica_space_with_state_df, edge_list, by.x = "sample_name", 
                   by.y = "source", all = TRUE)
  edge_df <- rename(edge_df, c(ICA_dim_1 = "source_ICA_dim_1", 
                               ICA_dim_2 = "source_ICA_dim_2"))
  edge_df <- merge(edge_df, ica_space_with_state_df[, c("sample_name", 
                                                        "ICA_dim_1", "ICA_dim_2")], by.x = "target", by.y = "sample_name", 
                   all = TRUE)
  edge_df <- rename(edge_df, c(ICA_dim_1 = "target_ICA_dim_1", 
                               ICA_dim_2 = "target_ICA_dim_2"))
  diam <- as.data.frame(as.vector(V(dp_mst)[get.diameter(dp_mst, 
                                                         weights = NA)]$name))
  colnames(diam) <- c("sample_name")
  diam <- arrange(merge(ica_space_with_state_df, diam, by.x = "sample_name", 
                        by.y = "sample_name"), Pseudotime)
  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% 
                              markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- melt(exprs(cds[row.names(markers_fData), 
                                      ]))
      markers_exprs <- merge(markers_exprs, markers_fData, 
                             by.x = "Var1", by.y = "row.names")
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
        0) {
    edge_df <- merge(edge_df, markers_exprs, by.x = "sample_name", 
                     by.y = "Var2")
    g <- ggplot(data = edge_df, aes(x = source_ICA_dim_1, 
                                    y = source_ICA_dim_2, size = log10(value + 0.1))) + 
      facet_wrap(~feature_label)
  }
  else {
    g <- ggplot(data = edge_df, aes(x = source_ICA_dim_1, 
                                    y = source_ICA_dim_2))
  }
  if (show_tree) {
    g <- g + geom_segment(aes_string(xend = "target_ICA_dim_1", 
                                     yend = "target_ICA_dim_2", color = color_by), size = 0.3, 
                          linetype = "solid", na.rm = TRUE)
  }
  #abbrev.std.stages.col=c('R'='green','LR'='green4','ET'='blue','LT'='blue4','ES'='deeppink3','S'='red')
  abbrev.std.stages.col=c('R'='cyan','LR'='green4','ET'='blue4','T'='orange','ES'='purple','S'='firebrick')
  g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE,size=point.size)
  if (show_backbone) {
    g <- g + geom_path(aes(x = ICA_dim_1, y = ICA_dim_2), 
                       color = I(backbone_color), size = path.size, data = diam, 
                       na.rm = TRUE) + geom_point(aes_string(x = "ICA_dim_1", 
                                                             y = "ICA_dim_2", color = color_by), size = I(1.5), 
                                                  data = diam, na.rm = TRUE)
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  g <- g + theme(panel.border = element_blank(), axis.line = element_line()) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    ylab("Component 1") + xlab("Component 2")  + theme(panel.background = element_rect(fill = "white")) +
    theme(legend.position = "top", legend.key.height = unit(0.35, "in"))+
    theme(legend.key = element_blank()) +
    scale_colour_manual(name='Development stage',values = abbrev.std.stages.col)
  g
}
plot.spanning.tree.without.std.cols=function (cds, x = 1, y = 2, color_by = "State", show_tree = TRUE, 
                             show_backbone = TRUE, backbone_color = "black", markers = NULL, 
                             show_cell_names = FALSE, cell_name_size = 1,path.size=1.0,point.size=2.5) 
{
  lib_info_with_pseudo <- pData(cds)
  S_matrix <- reducedDimS(cds)
  if (is.null(S_matrix)) {
    stop("You must first call reduceDimension() before using this function")
  }
  ica_space_df <- data.frame(t(S_matrix[c(x, y), ]))
  colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
  ica_space_df$sample_name <- row.names(ica_space_df)
  ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, 
                                   by.x = "sample_name", by.y = "row.names")
  dp_mst <- minSpanningTree(cds)
  if (is.null(dp_mst)) {
    stop("You must first call orderCells() before using this function")
  }
  edge_list <- as.data.frame(get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")
  edge_df <- merge(ica_space_with_state_df, edge_list, by.x = "sample_name", 
                   by.y = "source", all = TRUE)
  edge_df <- rename(edge_df, c(ICA_dim_1 = "source_ICA_dim_1", 
                               ICA_dim_2 = "source_ICA_dim_2"))
  edge_df <- merge(edge_df, ica_space_with_state_df[, c("sample_name", 
                                                        "ICA_dim_1", "ICA_dim_2")], by.x = "target", by.y = "sample_name", 
                   all = TRUE)
  edge_df <- rename(edge_df, c(ICA_dim_1 = "target_ICA_dim_1", 
                               ICA_dim_2 = "target_ICA_dim_2"))
  diam <- as.data.frame(as.vector(V(dp_mst)[get.diameter(dp_mst, 
                                                         weights = NA)]$name))
  colnames(diam) <- c("sample_name")
  diam <- arrange(merge(ica_space_with_state_df, diam, by.x = "sample_name", 
                        by.y = "sample_name"), Pseudotime)
  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% 
                              markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- melt(exprs(cds[row.names(markers_fData), 
                                      ]))
      markers_exprs <- merge(markers_exprs, markers_fData, 
                             by.x = "Var1", by.y = "row.names")
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
        0) {
    edge_df <- merge(edge_df, markers_exprs, by.x = "sample_name", 
                     by.y = "Var2")
    g <- ggplot(data = edge_df, aes(x = source_ICA_dim_1, 
                                    y = source_ICA_dim_2, size = log10(value + 0.1))) + 
      facet_wrap(~feature_label)
  }
  else {
    g <- ggplot(data = edge_df, aes(x = source_ICA_dim_1, 
                                    y = source_ICA_dim_2))
  }
  if (show_tree) {
    g <- g + geom_segment(aes_string(xend = "target_ICA_dim_1", 
                                     yend = "target_ICA_dim_2", color = color_by), size = 0.3, 
                          linetype = "solid", na.rm = TRUE)
  }
  g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE,size=point.size)
  if (show_backbone) {
    g <- g + geom_path(aes(x = ICA_dim_1, y = ICA_dim_2), 
                       color = I(backbone_color), size = path.size, data = diam, 
                       na.rm = TRUE) + geom_point(aes_string(x = "ICA_dim_1", 
                                                             y = "ICA_dim_2", color = color_by), size = I(1.5), 
                                                  data = diam, na.rm = TRUE)
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  g <- g + theme(panel.border = element_blank(), axis.line = element_line()) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    ylab("Component 1") + xlab("Component 2")  + theme(panel.background = element_rect(fill = "white")) +
    theme(legend.position = "top", legend.key.height = unit(0.35, "in"))+
    theme(legend.key = element_blank()) 
  g
}
geneScatterplot <- function( x, y, xlab, ylab, col ,main.title='Test title') {
  plot( NULL, xlim=c( -.1, 6.2 ), ylim=c( -1, 6.2 ),
        xaxt="n", yaxt="n", xaxs="i", yaxs="i", asp=1,
        xlab=xlab, ylab=ylab,main= main.title)
  abline( a=-1, b=1, col = "lightgray", lwd=2 )
  abline( a=0, b=1, col = "lightgray", lwd=2 )
  abline( a=1, b=1, col = "lightgray", lwd=2 )
  abline( h=c(0,2,4,6), v=c(0,2,4,6), col = "lightgray", lwd=2 )
  points(
    ifelse( x > 0, log10(x), -.7 ),
    ifelse( y > 0, log10(y), -.7 ),
    pch=19, cex=.2, col = col )
  axis( 1, c( -.7, 0:6 ), c( "0", "1", "10", "100", expression(10^3), expression(10^4),
                             expression(10^5), expression(10^6) ) )
  axis( 2, c( -.7, 0:6 ),
        c( "0", "1", "10", "100", expression(10^3), expression(10^4),
           expression(10^5), expression(10^6) ), las=2 )
  axis( 1, -.35, "//", tick=FALSE, line=-.7 )
  axis( 2, -.35, "\\\\", tick=FALSE, line=-.7 )
}
filter.genes.with.cov.bias=function(rpkm.df,cov.bias.test.df){
  samples=colnames(rpkm.df)
  cov.bias.test.df=cov.bias.test.df[,samples]
  intersect.genes=intersect(rownames(rpkm.df),rownames(cov.bias.test.df))
  rpkm.df=rpkm.df[intersect.genes,]
  cov.bias.test.df=cov.bias.test.df[intersect.genes,]
  len.intersect.genes=length(intersect.genes)
  len.samples=length(samples)
  out.mat=matrix(nrow = len.intersect.genes,ncol =len.samples )
  bias.genes.list=c()
  bias.samples.list=c()
  bias.counts=c()
  for(m in 1:len.intersect.genes){
    intersect.gene=intersect.genes[m]
    out.rpkm=ifelse(cov.bias.test.df[intersect.gene,],as.numeric(as.character(rpkm.df[intersect.gene,])),0.0)
    out.mat[m,]=out.rpkm
    test.outcome=cov.bias.test.df[intersect.gene,]
    out.test=ifelse(test.outcome,NA,as.character(colnames(rpkm.df[intersect.gene,])))
    out.test=as.character(out.test[!is.na(out.test)])
    len.out.test=length(out.test)
    if(len.out.test>0){
      bias.genes.list=append(bias.genes.list,values =intersect.gene ,length(bias.genes.list))
      bias.samples.list=append(bias.samples.list,paste(out.test,collapse = ':'),length(bias.samples.list))
      bias.counts=append(bias.counts,len.out.test,length(bias.counts))
    }
  }
  out.df=data.frame(out.mat)
  colnames(out.df)=samples
  rownames(out.df)=intersect.genes
  out.bias.genes.df=data.frame(genes=bias.genes.list,samples.counts=bias.counts,samples=bias.samples.list)
  rownames(out.bias.genes.df)=bias.genes.list
  out.list=list(rpkm=out.df,biased.genes=out.bias.genes.df)
  return(out.list)
}
plot.gene.jitters=function (cds_subset, grouping = "State", min_expr = 0.1, cell_size = 0.75, 
          nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL, 
          plot_trend = F, labelled=T) 
{
  cds_exprs <- melt(exprs(cds_subset))
  colnames(cds_exprs) <- c("gene_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "gene_id", 
                     by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  ordered.development.stages=intersect(c("R" ,"LR", "ET", "T" , "ES"  , "S"),unique(as.character(cds_exprs$development.stage)))
  cds_exprs$development.stage=factor(x =cds_exprs$development.stage,levels = ordered.development.stages)
  col.factor.list=get.col.factor( col.factor = as.character(cds_exprs$development.stage))
  col.str=get.abbrev.std.stage.cols(in.stages.vec =as.character(cds_exprs$development.stage))
  #col.str=col.factor.list$col.str
  col.legend=col.factor.list$legend.cols
  col.legend=std.stages.col.legend.fill
  str.legend=col.factor.list$legend.str
  cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  if (is.null(cds_exprs$gene_short_name) == FALSE) {
    cds_exprs$gene_label <- cds_exprs$gene_short_name
    cds_exprs$gene_label[is.na(cds_exprs$gene_label)] <- cds_exprs$gene_id
  }
  else {
    cds_exprs$gene_label <- cds_exprs$gene_id
  }
  if (is.null(panel_order) == FALSE) {
    cds_subset$gene_label <- factor(cds_subset$gene_label, 
                                    levels = panel_order)
  }
  if(labelled){
    #q <- ggplot(aes_string(x = grouping, y = "expression"), data = cds_exprs)
    q <- ggplot(aes_string(x = grouping, y = "expression",colour='development.stage'), data = cds_exprs)+scale_colour_manual(names=,values = col.legend)
  }
  else{
    #show(head(cds_exprs))
    #q <- ggplot(aes_string(x = grouping, y = "expression",colour='development.stage'), data = cds_exprs)+theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.text =element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",panel.border=element_blank(),strip.text.x = element_blank(),panel.margin = unit(c(1,2), "lines"))+scale_colour_manual(values = col.legend)
    q <- ggplot(aes_string(x = grouping, y = "expression",colour='development.stage'), data = cds_exprs)+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",panel.border=element_blank(),strip.text.x = element_blank(),panel.margin = unit(c(1,2), "lines"))+scale_colour_manual(values = col.legend)
  }
  if (is.null(color_by) == FALSE) {
    q <- q + geom_jitter(size = I(cell_size)) +scale_colour_manual(values = col.legend)
  }
  else {
    q <- q + geom_jitter(size = I(cell_size))
  }
  if (plot_trend == TRUE) {
    q <- q + stat_summary(aes_string(color = color_by), fun.data = "mean_cl_boot", 
                          size = .8)
    #q <- q + stat_summary(fun.data = "mean_cl_boot", size = .8)+scale_colour_manual(values = col.str)
    q <- q + stat_summary(aes_string(x = grouping, y = "expression", 
                                     color = color_by, group = color_by), fun.data = "mean_cl_boot", 
                          size = .5, geom = "line")
    q <- q + scale_y_log10() + facet_wrap(~gene_label, nrow = nrow, 
                                          ncol = ncol, scales = "free_y")
    q <- q + ylab("") + xlab("")
  }
  print(q)
  #plot.new()
  #legend('center',legend=str.legend,fill=col.legend)
  #plot.new()
  #legend('center',legend=rep('',times = length(str.legend)),fill=col.legend,box.lty=0,cex=1.5)
}
plot.gene.jitters.without.std.cols=function (cds_subset, grouping = "State", min_expr = 0.1, cell_size = 0.75, 
                            nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL, 
                            plot_trend = F, labelled=T) 
{
  cds_exprs <- melt(exprs(cds_subset))
  colnames(cds_exprs) <- c("gene_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "gene_id", 
                     by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  if (is.null(cds_exprs$gene_short_name) == FALSE) {
    cds_exprs$gene_label <- cds_exprs$gene_short_name
    cds_exprs$gene_label[is.na(cds_exprs$gene_label)] <- cds_exprs$gene_id
  }
  else {
    cds_exprs$gene_label <- cds_exprs$gene_id
  }
  if (is.null(panel_order) == FALSE) {
    cds_subset$gene_label <- factor(cds_subset$gene_label, 
                                    levels = panel_order)
  }
  if(labelled){
    #q <- ggplot(aes_string(x = grouping, y = "expression"), data = cds_exprs)
    q <- ggplot(aes_string(x = grouping, y = "expression",colour=grouping), data = cds_exprs)
  }
  else{
    q <- ggplot(aes_string(x = grouping, y = "expression",colour=grouping), data = cds_exprs)+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",panel.border=element_blank(),strip.text.x = element_blank(),panel.margin = unit(c(1,2), "lines"))
  }
  if (is.null(color_by) == FALSE) {
    q <- q + geom_jitter(size = I(cell_size))
  }
  else {
    q <- q + geom_jitter(size = I(cell_size))
  }
  if (plot_trend == TRUE) {
    q <- q + stat_summary(aes_string(color = color_by), fun.data = "mean_cl_boot", 
                          size = .8)
    #q <- q + stat_summary(fun.data = "mean_cl_boot", size = .8)+scale_colour_manual(values = col.str)
    q <- q + stat_summary(aes_string(x = grouping, y = "expression", 
                                     color = color_by, group = color_by), fun.data = "mean_cl_boot", 
                          size = .5, geom = "line")
    q <- q + scale_y_log10() + facet_wrap(~gene_label, nrow = nrow, 
                                          ncol = ncol, scales = "free_y")
    q <- q + ylab("") + xlab("")
  }
  print(q)
  #plot.new()
  #legend('center',legend=str.legend,fill=col.legend)
  #plot.new()
  #legend('center',legend=rep('',times = length(str.legend)),fill=col.legend,box.lty=0,cex=1.5)
}
plot.jitters.at.gene.intervals=function(HSMM, grouping = "State", min_expr = 0.1, cell_size = 0.75, 
                                        nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL, 
                                        plot_trend = F, labelled=T,intervals=4,genes.vec){
  genes.subsets=list()
  in.vec=genes.vec
  len.in.vec=length(in.vec)
  chunk.sizes=intervals
  chunk.index=seq(1,len.in.vec,by =chunk.sizes)
  len.chunk.index=length(chunk.index)
  out.list=list()
  for(m in 1:len.chunk.index){
    d=chunk.index[m]
    start.id=d
    end.in=d+chunk.sizes-1
    if(end.in>len.in.vec){
      end.in=len.in.vec
    }
    sub.vec=in.vec[start.id:end.in]
    plot.gene.jitters(HSMM[sub.vec,], grouping=grouping, ncol=ncol, color_by=color_by,plot_trend =plot_trend)
  }
}
plot.jitters.at.gene.intervals.without.std.col=function(HSMM, grouping = "State", min_expr = 0.1, cell_size = 0.75, 
                                        nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL, 
                                        plot_trend = F, labelled=T,intervals=4,genes.vec){
  genes.subsets=list()
  in.vec=genes.vec
  len.in.vec=length(in.vec)
  chunk.sizes=intervals
  chunk.index=seq(1,len.in.vec,by =chunk.sizes)
  len.chunk.index=length(chunk.index)
  out.list=list()
  for(m in 1:len.chunk.index){
    d=chunk.index[m]
    start.id=d
    end.in=d+chunk.sizes-1
    if(end.in>len.in.vec){
      end.in=len.in.vec
    }
    sub.vec=in.vec[start.id:end.in]
    plot.gene.jitters.without.std.cols(HSMM[sub.vec,], grouping=grouping, ncol=ncol, color_by=color_by,plot_trend =plot_trend)
  }
}
#Filters artefacts
filter.artefacts.from.df=function(artefact.df,counts.df){
  artefact.genes=rownames(artefact.df)
  count.genes=rownames(counts.df)
  intersect.genes=intersect(artefact.genes,count.genes)
  samples=colnames(counts.df)
  len.genes=length(count.genes)
  len.samples=length(samples)
  out.mat=matrix(nrow = len.genes,ncol = len.samples)
  for (n in 1:len.genes){
    gene=count.genes[n]
    if(gene %in% intersect.genes ){
      for (m in 1:len.samples){
        sample=samples[m]
        gene.count=as.numeric(as.character(counts.df[gene,sample]))
        if(gene.count==0){
          out.mat[n,m]=0
          next
        }
        gene.artefact.test=as.character(artefact.df[gene,sample])
        if(gene.artefact.test=='fail'){
          out.mat[n,m]=0
        }
        else{
          out.mat[n,m]=gene.count
        }
      }
    }
    else{
      gene.parts=strsplit(x = gene,split = '\\+')[[1]]
      #gene.parts.artefacts.test=c()
      for (m in 1:len.samples){
        sample=samples[m]
        gene.count=as.numeric(as.character(counts.df[gene,sample]))
        gene.parts.artefacts.test=as.character(artefact.df[gene.parts,sample])
        if(all(gene.parts.artefacts.test!='fail')){
          out.mat[n,m]=as.numeric(as.character(counts.df[gene,sample]))
        }
        else{
          out.mat[n,m]=0
        }
      }
    }
  }
  out.df=data.frame(out.mat)
  rownames(out.df)=count.genes
  colnames(out.df)=samples
  return(out.df)
}
get.artefact.genes=function(artefact.df){
  samples=colnames(artefact.df)
  len.samples=length(samples)
  genes=rownames(artefact.df)
  len.genes=length(genes)
  samples.names=c()
  gene.names=c()
  for(m in 1:len.genes){
    gene=genes[m]
    for(n in 1:len.samples){
      sample=samples[n]
      artefact.test=as.character(artefact.df[gene,sample])
      if(artefact.test=='fail'){
        samples.names=append(samples.names,sample,length(samples.names))
        gene.names=append(gene.names,gene,length(gene.names))
      }
      else{
        next
      }
    }
  }
  out.df=data.frame(sample=samples.names,gene=gene.names)
  return(out.df)
}
classify.sc.using.bulk.subgroups=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$sub.group)
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sub.grps.names=names(pop.meta.list)
  len.pop.meta.names=length(sub.grps.names)
  sc.samples=colnames(sc.rpkm.df)
  len.sc.samples=length(sc.samples)
  sc.to.keep=c()
  sub.group.vec=c()
  for(m in 1:len.sc.samples){
    sc.sample=sc.samples[m]
    temp.sc.rpkm.df=subset(sc.rpkm.df,select=sc.sample)
    cor.list=list()
    for(n in 1:len.pop.meta.names){
      sub.grps.name=sub.grps.names[n]
      temp.pop.meta.df=pop.meta.list[[sub.grps.name]]
      temp.pop.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(subset(pop.rpkm.df,select = rownames(temp.pop.meta.df))))
      intersect.genes.vec=intersect(rownames(temp.pop.rpkm.df),rownames(temp.sc.rpkm.df))
      temp.cor.vec=get.cor.between.two.df.cor(first.df =temp.sc.rpkm.df,second.df =  temp.pop.rpkm.df)
      cor.list[[sub.grps.name]]=temp.cor.vec
    }
    compared.stages=names(cor.list)
    len.compared.stages=length(compared.stages)
    for(l in 1:len.compared.stages){
      compared.stage=compared.stages[l]
      compared.stage.cor.vec=cor.list[[compared.stage]]
      other.stages=setdiff(compared.stages,compared.stage)
      len.other.stages=length(other.stages)
      test.vec=c()
      for(d in 1:len.other.stages){
        other.stage=other.stages[d]
        other.stage.cor.vec=cor.list[[other.stage]]
        p.value.score=as.numeric(wilcox.test(x =compared.stage.cor.vec,y = other.stage.cor.vec ,alternative = 'greater')$p.value)
        test.vec=append(test.vec,ifelse(p.value.score<=0.05,T,F),after = length(test.vec))
      }
      median.cor.score=median(compared.stage.cor.vec)
      if(all(test.vec)){
        if(median.cor.score >=.3){
          sc.to.keep=append(sc.to.keep,sc.sample,after = length(sc.to.keep))
          sub.group.vec=append( sub.group.vec,compared.stage,after = length( sub.group.vec))
          boxplot2(cor.list,las=2, main =sc.sample)
        }
      }
    }
  }
  sc.out.rpkm.df=subset(sc.rpkm.df,select=sc.to.keep)
  sc.out.meta.df=data.frame(t(subset(t(sc.meta.df),select = colnames(sc.out.rpkm.df))))
  sc.out.meta.df$top.sub.grp=sub.group.vec
  sc.out.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(sc.out.rpkm.df))
  sc.out.meta.df=data.frame(t(subset(t(sc.out.meta.df),select=colnames(sc.out.rpkm.df))))
  out.list=list(rpkm=sc.out.rpkm.df,meta=sc.out.meta.df)
  return(out.list)
}
classify.sc.using.bulk.samples=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df,cor.cut.off=.3){
  sc.rpkm.df=filter.none.expressed.samples(df =sc.rpkm.df )
  pop.rpkm.df=filter.none.expressed.samples(df =pop.rpkm.df )
  pop.samples=colnames(pop.rpkm.df)
  pop.meta.df=pop.meta.df[pop.samples,]
  pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df =pop.meta.df )
  sc.samples=colnames(sc.rpkm.df)
  sc.meta.df=sc.meta.df[sc.samples,]
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df =sc.meta.df )
  len.sc.samples=length(sc.samples)
  len.pop.samples=length(pop.samples)
  sc.to.keep=c()
  sub.group.vec=c()
  intersect.genes.vec=intersect(rownames(sc.rpkm.df),rownames(pop.rpkm.df))
  combn.df=data.frame(pop.rpkm.df[intersect.genes.vec,],sc.rpkm.df[intersect.genes.vec,],row.names = intersect.genes.vec)
  samples.cor.score.mat=data.frame(cor(combn.df,method = 'spearman'))
  samples.cor.score.mat=samples.cor.score.mat[sc.samples,pop.samples]
  matching.stages.vec=c()
  sc.samples.df.list=list()
  for(m in 1:len.sc.samples){
    sc.sample=sc.samples[m]
    stage.cor.vec=c()
    cor.stages.names.vec=c()
    bulk.samples.vec=c()
    for(n in 1: len.pop.samples){
      pop.sample=pop.samples[n]
      temp.cor=as.numeric(samples.cor.score.mat[sc.sample,pop.sample])
      temp.cor=ifelse(is.na(temp.cor),0,temp.cor)
      if(temp.cor>=cor.cut.off){
        cor.detected.stage=as.character(pop.meta.df[pop.sample,'development.stage'])
        stage.cor.vec=append(stage.cor.vec,temp.cor,after = length(stage.cor.vec))
        cor.stages.names.vec=append(cor.stages.names.vec,cor.detected.stage,after = length(cor.stages.names.vec))
        bulk.samples.vec=append(bulk.samples.vec,pop.sample,after = length(bulk.samples.vec))
      }
      else{
        next
      }
    }
    if(length(stage.cor.vec)>=1){
      sc.to.keep=append(sc.to.keep,sc.sample,after = length(sc.to.keep))
      matching.stages=paste(as.character(names(table(cor.stages.names.vec))),collapse = '.')
      matching.stages.vec=append(matching.stages.vec,matching.stages,after = length(matching.stages.vec))
      temp.df=data.frame(bulk.sample=bulk.samples.vec,bulk.stage=cor.stages.names.vec,cor.score=stage.cor.vec)
      sc.samples.df.list[[sc.sample]]=temp.df
    }
  }
  sc.out.meta.df=sc.meta.df[sc.to.keep,]
  sc.out.meta.df$matching.bulks=sort(matching.stages.vec)
  sc.out.rpkm.df=sc.rpkm.df[,sc.to.keep]
  sc.out.rpkm.df=filter.none.expressed.genes(sc.out.rpkm.df)
  sc.out.meta.df=sc.out.meta.df[colnames(sc.out.rpkm.df),]
  out.list=list(sc.meta=sc.out.meta.df,sc.rpkm=sc.out.rpkm.df,sc.bulk.list=sc.samples.df.list)
  return(out.list)
}
get.sc.variable.genes.brenneck=function(rpkm.df,
                                        counts.df,
                                        mean.fit.genes.cut.off=.3,
                                        mean.exp.quantile.cut.off=.5,
                                        title.str='Test'){
  counts.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = counts.df))
  samples=colnames(counts.df)
  #rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(rpkm.df[,samples]))
  deseq.sf=estimateSizeFactorsForMatrix(counts = counts.df)
  temp.counts.mat=as.matrix(counts.df)
  temp.rpkm.df=t(t(counts.df)/deseq.sf)
  #temp.rpkm.df=rpkm.df
  col.points <- "#00207040"
  temp.rpkm.mat=as.matrix(temp.rpkm.df)
  means.rpkm <- rowMeans(temp.rpkm.mat )
  vars.rpkm <- rowVars( temp.rpkm.mat )
  cv2.rpkm <- vars.rpkm / means.rpkm^2
  minMeanForFit <- unname(quantile(means.rpkm[ which( cv2.rpkm > mean.fit.genes.cut.off ) ],  mean.exp.quantile.cut.off) )
  useForFit <- means.rpkm >= minMeanForFit
  fit <- glmgam.fit(cbind( a0 = 1, a1tilde = 1/means.rpkm[useForFit] ),cv2.rpkm[useForFit] )
  xi <- mean(1 / deseq.sf )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"] - xi )
  col.points <- "#00207040"
  # Prepare the plot (scales, grid, labels, etc.)
  psia1theta <- mean( 1 / deseq.sf ) + a1 * mean( deseq.sf )
  #NB : The original approach for psia1theta
  #psia1theta <- mean( 1 / sfAt ) + a1 * mean( sfHeLa / sfAt )
  minBiolDisp <- .5^2
  m <- ncol(temp.counts.mat)
  cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
  testDenom <- ( means.rpkm * psia1theta + means.rpkm^2 * cv2th ) / ( 1 + cv2th/m )
  p <- 1 - pchisq( vars.rpkm * (m-1) / testDenom, m-1 )
  padj <- p.adjust( p, "BH" )
  sig <- padj < .1
  sig[is.na(sig)] <- FALSE
  genes.test.res=data.frame(genes=names(sig),pval=as.character(sig))
  rownames(genes.test.res)=names(sig)
  variable.genes.vec=rownames(subset(genes.test.res,pval==T))
  # Add the data points
  plot(means.rpkm, cv2.rpkm, pch=20, cex=1.0,log="xy", 
       col=ifelse(padj < .1, "#C0007090",col.points),
       xlab = "average normalized read count", 
       ylab = "squared coefficient of variation (CV^2)",main=title.str )
  # Plot the fitted curve
  xg <- 10^seq( -2, 6, length.out=1000 )
  lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
  # Plot quantile lines around the fit
  df <- ncol(temp.counts.mat) - 1
  lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .975, df ) / df,
         col="#FF000080", lwd=2, lty="dashed" )
  lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .025, df ) / df,
         col="#FF000080", lwd=2, lty="dashed" )
  plot(means.rpkm, cv2.rpkm, pch=20, cex=1.0,log="xy", col=ifelse(padj < .1, "#C0007090",col.points),xlab = "", ylab = "",main='' )
  # Plot the fitted curve
  xg <- 10^seq( -2, 6, length.out=1000 )
  lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
  # Plot quantile lines around the fit
  df <- ncol(temp.counts.mat) - 1
  lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .975, df ) / df,
         col="#FF000080", lwd=2, lty="dashed" )
  lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .025, df ) / df,
         col="#FF000080", lwd=2, lty="dashed" )
  return(variable.genes.vec)
}
get.rpkm=function(counts.df,meta.df,gene.lengths.df){
  meta.df=meta.df[colnames(counts.df),]
  no.reads=as.numeric(meta.df$Uniquely_mapped_reads_number_)
  normalization.fact=1000000000/no.reads
  out.df=apply(counts.df,1,function(count.row){
    norm.expr=as.numeric(count.row)*normalization.fact
    return(norm.expr)
  })
  out.df=data.frame(t(out.df))
  colnames(out.df)=colnames(counts.df)
  gene.len.vec=get.gene.lengths(genes.vec = rownames(out.df),gene.lengths.df = gene.lengths.df)
  temp.out.df=apply(out.df,2,function(col.vec){
    col.vec=as.numeric(col.vec)/as.numeric(gene.len.vec)
    return(col.vec)
  })
  out.df=data.frame(temp.out.df)
  rownames(out.df)=names(gene.len.vec)
  return(out.df)
}
get.gene.lengths=function(genes.vec,gene.lengths.df){
  gene.len.vec=as.numeric(gene.lengths.df$Size)
  names(gene.len.vec)=rownames(gene.lengths.df)
  no.genes=length(genes.vec)
  out.gene.len.vec=c()
  for(m in 1:no.genes){
    temp.gene=genes.vec[m]
    if(grepl(pattern = ',',x = temp.gene)){
      temp.genes.vec=as.character(unlist(strsplit(temp.gene,split = ',')[[1]]))
      gene.len=as.numeric(max(gene.len.vec[temp.genes.vec]))
      out.gene.len.vec=append(x = out.gene.len.vec,values = gene.len,after = length(out.gene.len.vec))
    }
    else{
      gene.len=as.numeric(gene.len.vec[temp.gene])
      out.gene.len.vec=append(x = out.gene.len.vec,values = gene.len,after = length(out.gene.len.vec))
    }
  }
  names(out.gene.len.vec)=genes.vec
  return(out.gene.len.vec)
}
get.mean.normalize.exprs=function(rpkm.df){
  col.means.vec=as.numeric(apply(rpkm.df,2,mean))
  out.df=data.frame(t(t(rpkm.df)-col.means.vec))
  return(out.df)
}
get.distance.btw.average.reclassified.sc.and.pop=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$markers.cluster.groups)
  pop.meta.list=split(pop.meta.df,f=pop.meta.df$markers.cluster.groups)
  sc.grps=names(sc.meta.list)
  pop.grps=names(pop.meta.list)
  len.pop.grps=length(pop.grps)
  len.sc.grps=length(sc.grps)
  sc.average.mat=matrix(nrow = dim(sc.rpkm.df)[1],ncol =len.sc.grps )
  sc.colnames.vec=c()
  for(n in 1:len.sc.grps){
    sc.grp=sc.grps[n]
    temp.sc.meta.df=sc.meta.list[[sc.grp]]
    temp.sc.rpkm.df=subset(sc.rpkm.df,select=rownames(temp.sc.meta.df))
    average.rpkm=as.numeric(apply(temp.sc.rpkm.df,1,mean))
    sc.average.mat[,n]=average.rpkm
    sc.grp.lab=gsub(pattern = ':',replacement = '',x = gsub(pattern = ' ',replacement = '',x = sc.grp))
    sc.colnames.vec=append(sc.colnames.vec,sc.grp.lab,after = length(sc.colnames.vec))
  }
  sc.average.df=data.frame(sc.average.mat)
  colnames(sc.average.df)=sc.colnames.vec
  rownames(sc.average.df)=rownames(sc.rpkm.df)
  intersect.genes=intersect(rownames(sc.average.df),rownames(pop.rpkm.df))
  combn.df=data.frame(sc.average.df[intersect.genes,],pop.rpkm.df[intersect.genes,])
  cor.dist.mat=data.frame(cor(combn.df,method = 'spearman'))
  euclidian.dist.mat=data.frame(as.matrix(dist(t(combn.df))))
  plot.euclidian.dist.mat=matrix(nrow = dim(sc.average.df)[2],ncol =dim(pop.rpkm.df)[2] )
  plot.cor.dist.mat=matrix(nrow = dim(sc.average.df)[2],ncol =dim(pop.rpkm.df)[2] )
  len.sc.colnames.vec=length(sc.colnames.vec)
  pop.samples.vec=colnames(pop.rpkm.df)
  len.pop.samples.vec=length(pop.samples.vec)
  for(m in 1:len.sc.colnames.vec){
    for(l in 1:len.pop.samples.vec){
      sc.grp=sc.colnames.vec[m]
      pop.sample=pop.samples.vec[l]
      cor.dist=cor.dist.mat[sc.grp,pop.sample]
      euclidian.dist=euclidian.dist.mat[sc.grp,pop.sample]
      plot.cor.dist.mat[m,l]=cor.dist
      plot.euclidian.dist.mat[m,l]=euclidian.dist
    }
  }
  plot.cor.dist.df=data.frame(plot.cor.dist.mat)
  plot.euclidian.dist.df=data.frame(plot.euclidian.dist.mat)
  rownames(plot.cor.dist.df)=sc.colnames.vec
  rownames(plot.euclidian.dist.df)=sc.colnames.vec
  colnames(plot.cor.dist.df)=pop.samples.vec
  colnames(plot.euclidian.dist.df)=pop.samples.vec
  pop.grp.vec=pop.meta.df[colnames(plot.cor.dist.df),]
  pop.grp.vec=as.character(pop.grp.vec$markers.cluster.groups)
  pop.grp.col.list=get.col.factor(col.factor = pop.grp.vec)
  heatmap.2(as.matrix(plot.cor.dist.df),trace='none',margins = c(10,10),ColSideColors = pop.grp.col.list$col.str,cexRow = .3,cexCol = .3,main = 'Correlation distance',col=bluered(10000))
  legend('topright',legend=pop.grp.col.list$legend.str,fill=pop.grp.col.list$legend.col,cex=0.5,box.lty = 0)
  pop.grp.vec=pop.meta.df[colnames(plot.euclidian.dist.df),]
  pop.grp.vec=as.character(pop.grp.vec$markers.cluster.groups)
  pop.grp.col.list=get.col.factor(col.factor = pop.grp.vec)
  heatmap.2(as.matrix(log10(plot.euclidian.dist.df)),trace='none',margins = c(10,10),ColSideColors = pop.grp.col.list$col.str,cexRow = .3,cexCol = .3,main = 'Euclidean distance',col=bluered(10000))
  legend('topright',legend=pop.grp.col.list$legend.str,fill=pop.grp.col.list$legend.col,cex=0.5,box.lty = 0)
}
get.distance.btw.reclassified.sc.and.timepoints.groups=function(sc.rpkm.df,sc.meta.df, title.str='Test'){
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.list=split(sc.meta.df,f=sc.meta.df$markers.cluster.groups)
  sc.grps=names(sc.meta.list)
  len.sc.grps=length(sc.grps)
  reclassified.cor.dist=c()
  reclassified.euclidean.dist=c()
  dist.list=list()
  eucliean.dist.list=list()
  for(n in 1:len.sc.grps){
    sc.grp=sc.grps[n]
    temp.sc.meta.df=sc.meta.list[[sc.grp]]
    temp.sc.rpkm.df=subset(sc.rpkm.df,select=rownames(temp.sc.meta.df))
    temp.cor.dist=data.frame(cor(temp.sc.rpkm.df,method = 'spearman'))
    temp.samples=colnames(temp.sc.rpkm.df)
    len.temp.samples=length(temp.samples)
    reclassified.cor.dist=append(x = reclassified.cor.dist,values = get.cor.list(df = temp.sc.rpkm.df),after = length(reclassified.cor.dist))
    euclidean.dist=get.euclean.dist.list( df = temp.sc.rpkm.df,log.rpkm = T)
    reclassified.euclidean.dist=append(reclassified.euclidean.dist,euclidean.dist,length(reclassified.euclidean.dist))
  }
  dist.list[['regrouped']]=reclassified.cor.dist
  eucliean.dist.list[['regrouped']]=reclassified.euclidean.dist
  sc.meta.timepoint.list=split(sc.meta.df,f=sc.meta.df$development.stage)
  sc.time.point.grps=names(sc.meta.timepoint.list)
  len.sc.time.point.grps=length(sc.time.point.grps)
  timepoint.based.cor.dist=c()
  timepoint.euclidean.dist=c()
  for(m in 1:len.sc.time.point.grps){
    sc.grp=sc.time.point.grps[m]
    temp.sc.meta.df=sc.meta.timepoint.list[[sc.grp]]
    temp.sc.rpkm.df=subset(sc.rpkm.df,select=rownames(temp.sc.meta.df))
    if(dim(temp.sc.rpkm.df)[2]<2){
      next
    }
    timepoint.based.cor.dist=append(x = timepoint.based.cor.dist,values = get.cor.list(df = temp.sc.rpkm.df),after = length(timepoint.based.cor.dist))
    euclidean.dist=get.euclean.dist.list( df = temp.sc.rpkm.df,log.rpkm = T)
    timepoint.euclidean.dist=append(timepoint.euclidean.dist,euclidean.dist,length(timepoint.euclidean.dist))
  }
  dist.list[['timepoint']]=timepoint.based.cor.dist
  eucliean.dist.list[['timepoint']]=timepoint.euclidean.dist
  cor.pvalue=as.numeric(wilcox.test(x = reclassified.cor.dist, y = timepoint.based.cor.dist,alternative = 'greater')$p.value)
  cor.pvalue=format(cor.pvalue,digits = 3)
  cor.title.str=paste('P.value = ',cor.pvalue,sep='')
  temp.boxplot=boxplot2(dist.list,las=2,ylab='Correlation coefficient')
  text(x=1.5,y=max(temp.boxplot$stats),labels=cor.title.str,cex=.7)
  euclidean.pvalue=as.numeric(wilcox.test(x = reclassified.euclidean.dist, y = timepoint.euclidean.dist,alternative = 'less')$p.value)
  euclidean.pvalue=format(euclidean.pvalue,digits = 3)
  euclidean.pvalue.lab=paste('P.value = ',euclidean.pvalue,sep='')
  temp.boxplot.2=boxplot2(eucliean.dist.list,las=2,ylab='Euclidean distance')
  text(x=1.5,y=max(temp.boxplot.2$stats),labels=euclidean.pvalue.lab,cex=.7)
}
get.euclean.dist.list=function(df,log.rpkm=F){
  dist.mat=matrix(nrow = dim(df)[2],ncol =dim(df)[2] )
  if(log.rpkm){
    pseudo.value=min(df[df>0])/2
    dist.mat=as.matrix(dist(log10(t(df+pseudo.value))))
  }
  else{
    dist.mat=as.matrix(dist(t(df)))
  }
  samples=colnames(df)
  len.samples=length(samples)
  dist.vec=c()
  samples.combn=combn(samples,2)
  no.pairs=dim(samples.combn)[2]
  if(len.samples==2){
    dist.vec=append(dist.vec,dist.mat[1,2],length(dist.vec))
  }
  else if(len.samples>2){
    for(n in 1:no.pairs){
      temp.pair=samples.combn[,n]
      first.sample=temp.pair[1]
      sec.sample=temp.pair[2]
      temp.dist=as.numeric(dist.mat[first.sample,sec.sample])
      dist.vec=append(dist.vec,temp.dist,length(dist.vec))
    }
  }
  else{
    dist.vec=append(dist.vec,NA,length(dist.vec))
  }
  return(dist.vec)
}
bootstrap.clusters.from.different.intervals=function(cluster.res.list,in.nboot=10){
  cluster.res.list.names=names(cluster.res.list)
  len.cluster.res.list.names=length(cluster.res.list.names)
  out.list=list()
  for(n in 1:len.cluster.res.list.names){
  #for(n in 1:1){
    cluster.res.list.name=cluster.res.list.names[n]
    rec.cut.off=as.numeric(gsub(pattern = 'interval.', replacement = '',x =  cluster.res.list.name))
    temp.interval.res.list=cluster.res.list[[cluster.res.list.name]]
    temp.markers.score.df=temp.interval.res.list[['markers.score']]
    temp.pvclust=pvclust(data = temp.markers.score.df,method.hclust = 'complete',
                         nboot = in.nboot,store = T,method.dist = 'euclidean')
    out.list[[cluster.res.list.name]]=temp.pvclust
    temp.order.df=data.frame(order=as.numeric(temp.pvclust$hclust$order),label=as.character(temp.pvclust$hclust$labels))
    rownames(temp.order.df)=as.character(temp.pvclust$hclust$labels)
    plot(temp.pvclust,cex=.2,cex.pv=.3)
    #text(temp.pvclust, col=c(2,3,8), print.num=TRUE, float=0.01, cex=.2, font=NULL)
    temp.markers.meta.df=temp.interval.res.list[['meta.df']]
    #temp.pvclust$hclust$labels=as.character(temp.markers.meta.df[as.character(temp.pvclust$hclust$labels),'markers.cluster.groups'])
    temp.order.df=temp.order.df[order(temp.order.df$order),]
    hclust.obj=temp.pvclust$hclust
    plot(hclust.obj,cex=.3,hang=-1)
    hcd = as.dendrogram(hclust.obj)
    # cut dendrogram in 4 clusters
    sp.col.vec=get.color.list.for.pheatmap(as.character(temp.markers.meta.df[as.character(temp.pvclust$hclust$labels),'markers.cluster.groups']))
    #labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951")
    labelColors =as.character(get.color.list.for.pheatmap(as.character(temp.markers.meta.df[as.character(temp.pvclust$hclust$labels),'markers.cluster.groups'])))
    clusMember = cutree(hclust.obj, 8)
    #clusMember = cutree(hclust.obj, 4)
    # function to get color labels
    colLab <- function(n) {
      if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
        #attr(n, "nodePar") <- c(a$nodePar, lab.cex = .2)
      }
      n
    }
    resize.lab <- function(n) {
      if (is.leaf(n)) {
        a <- attributes(n)
        #labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.cex = .2)
      }
      n
    }
    # using dendrapply
    clusDendro = dendrapply(hcd, colLab)
    #clusDendro = dendrapply(hcd,resize.lab)
    # make plot
    plot(clusDendro,cex=.5, 
         main = "Sub-populations")
    #pvrect(temp.pvclust,border = rec.cut.off,alpha =.50,pv = 'au',max.only = T)
    #seplot(object =temp.pvclust,type = 'au',identify = F,main = cluster.res.list.name,pch=19 )
    #final.out.mat=
    #pheatmap(mat=final.out.mat,cellwidth = 2,cellheight = 2,main='',fontsize_row  = 1,fontsize_col = 1,treeheight_row = 1,cutree_cols = in.sc.cut.off,clustering_distance_cols = 'euclidean',clustering_distance_rows ='euclidean',annotation_col = final.out.annotation.df,annotation_colors = final.out.annotation.col.list,annotation_row = row.annotation.df,annotation_legend = F,color = bluered(10000),show_rownames = F,show_colnames = F)
  }
  return(out.list)
}
pvclust.plot.editted=function (object, type = c("au", "bp"), identify = FALSE, main = NULL, 
          xlab = NULL, ylab = NULL, ...) 
{
  if (!is.na(pm <- pmatch(type[1], c("au", "bp")))) {
    wh <- c("au", "bp")[pm]
    if (is.null(main)) 
      main <- "p-value vs standard error plot"
    if (is.null(xlab)) 
      xlab <- c("AU p-value", "BP value")[pm]
    if (is.null(ylab)) 
      ylab <- "Standard Error"
    plot(object$edges[, wh], object$edges[, paste("se", wh, 
                                                  sep = ".")], main = main, xlab = xlab, ylab = ylab, 
         ...)
    if (identify) 
      identify(x = object$edges[, wh], y = object$edges[, 
                                                        paste("se", wh, sep = ".")], labels = row.names(object$edges))
  }
  else stop("'type' should be \"au\" or \"bp\".")
}
plot.sig.idc.genes.heatmap.btw.sub.populations=function(rpkm.df,meta.df,scde.sub.pop.scde.res.list,
                                                        idc.markers.df,title.str='Test',
                                                        show_colnames    = F,annotation_legend = F,
                                                        in.col.ramp=colorRampPalette(c("gray",'red'))(1000)){
  #genes.col.vec=get.color.rainbow.list.for.pheatmap(in.vec =  as.character(idc.markers.df$category))
  genes.col.vec=gene.cat.col.names
  col.pelette=in.col.ramp
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  sub.pop.col.vec=get.color.list.for.pheatmap(in.vec =as.character(meta.df$markers.cluster.groups))
  idc.markers.vec=rownames(idc.markers.df)
  scde.sub.pop.scde.res.list.names=names(scde.sub.pop.scde.res.list)
  len.scde.sub.pop.scde.res.list.names=length(scde.sub.pop.scde.res.list.names)
  for(n in 1:len.scde.sub.pop.scde.res.list.names){
    scde.sub.pop.scde.res.list.name=scde.sub.pop.scde.res.list.names[n]
    scde.sub.pop.scde.res=scde.sub.pop.scde.res.list[[scde.sub.pop.scde.res.list.name]]
    scde.sub.pop.scde.meta.df=scde.sub.pop.scde.res$meta
    scde.sub.pop.scde.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = scde.sub.pop.scde.meta.df)
    scde.sub.pop.scde.meta.df=encode.timepoint(meta.df = scde.sub.pop.scde.meta.df)
    scde.sub.pop.scde.results.df=scde.sub.pop.scde.res$results
    sig.scde.genes.sub.pop.scde.results.df=subset(scde.sub.pop.scde.results.df,adj.p.value<=.05 & abs(cZ) >=2)
    temp.samples.vec=rownames(scde.sub.pop.scde.meta.df)
    temp.rpkm.df=filter.none.expressed.genes(subset(rpkm.df,select = temp.samples.vec))
    sig.diff.idc.genes=intersect(idc.markers.vec,
                                 unique(rownames(sig.scde.genes.sub.pop.scde.results.df)))
    if(length(sig.diff.idc.genes)<2){
      next
    }
    idc.genes.rpkm.df=temp.rpkm.df[sig.diff.idc.genes,]
    idc.genes.rpkm.mat=as.matrix(idc.genes.rpkm.df)
    pseudo.val=min(idc.genes.rpkm.mat[idc.genes.rpkm.mat!=0])/2
    log.idc.genes.rpkm.mat=log2(idc.genes.rpkm.mat+pseudo.val)
    row.annotation.df=idc.markers.df[sig.diff.idc.genes,]
    row.annotation.df=subset(x = row.annotation.df,select = category)
    #col.annotation.df=subset(scde.sub.pop.scde.meta.df,select=c(markers.cluster.groups,development.stage,
    #                                  Rosetting,timepoint))
    col.annotation.df=subset(scde.sub.pop.scde.meta.df,
                             select=c(markers.cluster.groups,development.stage))
    col.annotation.df$markers.cluster.groups=gsub(pattern = 'SC.grp',replacement ='SP' ,
                                                  x = as.character(col.annotation.df$markers.cluster.groups))
    temp.sub.pop.col.vec=sub.pop.col.vec[unique(as.character(col.annotation.df$markers.cluster.groups))]
    annotation.col.list=list(category=genes.col.vec[unique(row.annotation.df$category)],
                             markers.cluster.groups=temp.sub.pop.col.vec,
                             development.stage=abbrev.std.stages.col.vec)
    title.str=scde.sub.pop.scde.res.list.name
    gene.symbols.vec=as.character(idc.markers.df[rownames(log.idc.genes.rpkm.mat),
                                                 'X.Gene.Name.or.Symbol.'])
    gene.symbols.vec=ifelse(gene.symbols.vec=='null',rownames(log.idc.genes.rpkm.mat),
                            gene.symbols.vec)
    mod.title.str=unlist(lapply(strsplit(scde.sub.pop.scde.res.list.name,'_vs_'),function(sp){
      sp.one=paste('SP.',sub.populations.grp.order.names.vec[sp][1],sep = '')
      sp.two=paste('SP.',sub.populations.grp.order.names.vec[sp][2],sep = '')
      return(paste(sp.one, '_vs_',sp.two,sep=''))
    }))
    #title.str=scde.sub.pop.scde.res.list.name
    title.str=mod.title.str
    pheatmap(mat = log.idc.genes.rpkm.mat,treeheight_row = 20,fontsize_row = 6,fontsize_col = .5,
             cellheight = 4,cellwidth = 4,main = title.str,annotation_row = row.annotation.df,
             annotation_col = col.annotation.df,annotation_colors = annotation.col.list,
             show_colnames=show_colnames,annotation_legend = annotation_legend,cutree_cols = 2,
             border_color = NA,labels_row = gene.symbols.vec,color=col.pelette)
  }
}
write.sig.genes.btw.sub.populations=function(scde.sub.pop.scde.res.list,idc.markers.df,output_dir){
  dir.create(output_dir,recursive = T,mode = '0777')
  idc.markers.vec=rownames(idc.markers.df)
  scde.sub.pop.scde.res.list.names=names(scde.sub.pop.scde.res.list)
  len.scde.sub.pop.scde.res.list.names=length(scde.sub.pop.scde.res.list.names)
  for(n in 1:len.scde.sub.pop.scde.res.list.names){
    scde.sub.pop.scde.res.list.name=scde.sub.pop.scde.res.list.names[n]
    scde.sub.pop.scde.res=scde.sub.pop.scde.res.list[[scde.sub.pop.scde.res.list.name]]
    scde.sub.pop.scde.results.df=scde.sub.pop.scde.res$results
    sig.scde.genes.vec=rownames(subset(scde.sub.pop.scde.results.df,adj.p.value<=.05 & abs(cZ) >=2))
    idc.genes.vec=c()
    idc.grp.vec=c()
    sig.genes.vec=c()
    for (gn in rownames(scde.sub.pop.scde.results.df)){
      idc.gn=ifelse(gn %in% idc.markers.vec,'True','False')
      idc.genes.vec=append(x = idc.genes.vec,values = idc.gn,after = length(idc.genes.vec))
      idc.grp=ifelse(gn %in% idc.markers.vec,idc.markers.df[gn,'category'],'NA')
      idc.grp.vec=append(x = idc.grp.vec,values = idc.grp,after = length(idc.grp.vec))
      sig.gn=ifelse(gn %in% sig.scde.genes.vec,'True','False')
      sig.genes.vec=append(sig.genes.vec,values =sig.gn ,length(sig.genes.vec))
      }
    scde.sub.pop.scde.results.df[,'idc.gene']=idc.genes.vec
    scde.sub.pop.scde.results.df[,'idc.category']=idc.grp.vec
    scde.sub.pop.scde.results.df[,'sig.gene']=sig.genes.vec
    mod.title.str=unlist(lapply(strsplit(scde.sub.pop.scde.res.list.name,'_vs_'),function(sp){
      sp.one=paste('SP.',sub.populations.grp.order.names.vec[sp][1],sep = '')
      sp.two=paste('SP.',sub.populations.grp.order.names.vec[sp][2],sep = '')
      return(paste(sp.one, '_vs_',sp.two,'.tab',sep=''))
    }))
    output_fname=paste(c(output_dir, mod.title.str),collapse = '/')
    write.table(scde.sub.pop.scde.results.df[sig.scde.genes.vec,],file =output_fname,
                sep='\t',col.names = NA,quote = F)
  }
}
plot.sig.idc.genes.per.interval.heatmap.btw.sub.populations=function(rpkm.df,meta.df,scde.sub.pop.interval.scde.res.list,idc.markers.df,title.str='Test'){
  scde.sub.pop.interval.scde.res.list.names=names(scde.sub.pop.interval.scde.res.list)
  len.scde.sub.pop.interval.scde.res.list.names=length(scde.sub.pop.interval.scde.res.list.names)
  for(n in 1:len.scde.sub.pop.interval.scde.res.list.names){
    scde.sub.pop.interval.scde.res.list.name=scde.sub.pop.interval.scde.res.list.names[n]
    scde.sub.pop.scde.res.list=scde.sub.pop.interval.scde.res.list[[scde.sub.pop.interval.scde.res.list.name]]
    plot.sig.idc.genes.heatmap.btw.sub.populations(rpkm.df =rpkm.df,meta.df = meta.df,scde.sub.pop.scde.res.list =scde.sub.pop.scde.res.list,idc.markers.df = idc.markers.df,title.str =  scde.sub.pop.interval.scde.res.list.name )
  }
}
plot.api.tf.genes.heatmap.btw.sub.populations=function(rpkm.df,meta.df,tf.df,title.str='Test',show_colnames    = F,annotation_legend = F){
  sub.pop.col.vec=get.color.list.for.pheatmap(in.vec =as.character(meta.df$markers.cluster.groups))
  tf.markers.vec=rownames(tf.df)
  #rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df = rpkm.df,sample.count = 2)
  meta.df=meta.df[colnames(rpkm.df),]
  #ordered.grps=c(sort(unique(meta.df$markers.cluster.groups)))
  ordered.grps=c("SC.grp.1" ,"SC.grp.2", "SC.grp.3" ,"SC.grp.4" ,"SC.grp.5", "SC.grp.6" ,"SC.grp.7" ,"SC.grp.8")
  meta.df$markers.cluster.groups=factor(meta.df$markers.cluster.groups,levels =factor(ordered.grps))
  meta.df=meta.df[with(meta.df, order(markers.cluster.groups)), ]
  expr.tf.markers.vec=intersect(tf.markers.vec,rownames(rpkm.df))
  if(length(expr.tf.markers.vec)>=2){
    tf.rpkm.df=rpkm.df[expr.tf.markers.vec,]
    tf.rpkm.mat=as.matrix(tf.rpkm.df)
    pseudo.val=min(tf.rpkm.mat[tf.rpkm.mat!=0])/2
    log.tf.rpkm.mat=log2(tf.rpkm.mat+pseudo.val)
    log.tf.rpkm.mat=log.tf.rpkm.mat[,rownames(meta.df)]
    col.annotation.df=subset(meta.df,select=markers.cluster.groups)
    temp.sub.pop.col.vec=sub.pop.col.vec[unique(as.character(col.annotation.df$markers.cluster.groups))]
    annotation.col.list=list(markers.cluster.groups=temp.sub.pop.col.vec)
    pheatmap(mat = log.tf.rpkm.mat,fontsize_row = .5,fontsize_col = .5,cellheight = 7,cellwidth = 3,main = title.str,annotation_colors = annotation.col.list,annotation_col = col.annotation.df,show_colnames=show_colnames,annotation_legend = F,cluster_cols = F,cluster_rows = F)
    pheatmap(mat = log.tf.rpkm.mat,fontsize_row = .5,fontsize_col = .5,cellheight = 7,cellwidth = 3,main = title.str,annotation_colors = annotation.col.list,annotation_col = col.annotation.df,show_colnames=show_colnames,annotation_legend = annotation_legend,cluster_cols = F,cluster_rows = F)
  }
}
plot.variant.surface.marker.genes.heatmap.btw.development.stages=function(rpkm.df,meta.df,markers.df,title.str='Test',show_rownames=F,show_colnames    = F,annotation_legend = F){
  #sub.pop.col.vec=get.color.list.for.pheatmap(in.vec =as.character(meta.df$markers.cluster.groups))
  #sub.pop.col.vec=get.color.rainbow.list.for.pheatmap(in.vec =as.character(meta.df$markers.cluster.groups))
  markers.col.vec=get.color.rainbow.list.for.pheatmap( in.vec =markers.df$category )
  markers.vec=rownames(markers.df)
  rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  #rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df = rpkm.df,sample.count = 2)
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=order.df.by.specific.col(meta.df = meta.df,order.vec =std.stages.ordered.str )
  rpkm.df=subset(rpkm.df,select = rownames(meta.df))
  expr.markers.vec=intersect(markers.vec,rownames(rpkm.df))
  if(length(expr.markers.vec)>=2){
    vsa.rpkm.df=rpkm.df[expr.markers.vec,]
    vsa.rpkm.mat=as.matrix(vsa.rpkm.df)
    pseudo.val=min(vsa.rpkm.mat[vsa.rpkm.mat!=0])/2
    log.vsa.rpkm.mat=log2(vsa.rpkm.mat+pseudo.val)
    log.vsa.rpkm.mat=log.vsa.rpkm.mat[,rownames(meta.df)]
    col.annotation.df=subset(meta.df,select=development.stage)
    #temp.sub.pop.col.vec=sub.pop.col.vec[unique(as.character(col.annotation.df$markers.cluster.groups))]
    row.annotation.df=markers.df[expr.markers.vec,]
    show(table(row.annotation.df$category))
    expr.categories.vec=unique(as.character(markers.df[expr.markers.vec,'category']))
    expr.categories.vec=markers.col.vec[expr.categories.vec]
    row.annotation.df=row.annotation.df[rownames(log.vsa.rpkm.mat),]
    row.annotation.df=subset(markers.df,select=category)
    temp.df=data.frame(markers.df[rownames(log.vsa.rpkm.mat),])
    annotation.col.list=list(category=expr.categories.vec,development.stage=abbrev.std.stages.col.vec)
    pheatmap(mat = log.vsa.rpkm.mat,fontsize_row = 3,fontsize_col = 1,cellheight = 1.5,cellwidth = 1,main = title.str,annotation_colors = annotation.col.list,annotation_row = row.annotation.df,annotation_col = col.annotation.df,show_rownames = show_rownames,show_colnames=show_colnames,annotation_legend = T,cluster_cols = F,cluster_rows = F,border_color = NA)
    pheatmap(mat = log.vsa.rpkm.mat,fontsize_row = 3,fontsize_col = .5,cellheight = 1.5,cellwidth = 1,main = title.str,annotation_colors = annotation.col.list,annotation_row = row.annotation.df,annotation_col = col.annotation.df,show_colnames=show_colnames,annotation_legend = annotation_legend,cluster_cols = F,cluster_rows = F,border_color = NA)
    barplot.col.vec=markers.col.vec[as.character(markers.df[rownames(vsa.rpkm.df),'category'])]
    #barplot(as.numeric(apply(vsa.rpkm.df,1,nnzero))/dim(vsa.rpkm.df)[2],names.arg = rownames(vsa.rpkm.df),las=2,cex.names = .2,col=barplot.col.vec,border = NA)
    #abline(h = 0.02,b = 1)
    top.two.vsa.list=apply(log.vsa.rpkm.mat,2,function(log.vsa.expr.vec){
      top.two.vec=sort(log.vsa.expr.vec,decreasing = T)[1:2]
      return(top.two.vec)
    })
    temp.vsa.df=data.frame(top.two.vsa.list)
    show(dim(temp.vsa.df))
    subset(temp.vsa.df)
    top.vsa.mat=as.matrix(data.frame(top.two.vsa.list))
    starts.vec=seq(1,304,76)
    for(n in starts.vec){
      col.indices=seq(n,n+75)
      barplot(height = top.vsa.mat[,col.indices],col=c('red','blue'))
    }
    temp.df=data.frame(freq=as.numeric(apply(vsa.rpkm.df,1,nnzero))/dim(vsa.rpkm.df)[2],gene=rownames(vsa.rpkm.df),row.names = rownames(vsa.rpkm.df))
    temp.df=subset(temp.df,freq>=0.02)
    temp.barplot.col.vec=markers.col.vec[as.character(markers.df[rownames(temp.df),'category'])]
    barplot(as.numeric(temp.df$freq),names.arg = rownames(temp.df),las=2,cex.names = .2,col=temp.barplot.col.vec,border = NA)
    abline(h = 0.02,b = 1)
  }
  return(expr.markers.vec)
}
plot.variant.surface.marker.genes.heatmap.btw.single.and.bulk=function(rpkm.df,
                                                                       meta.df,
                                                                       markers.df,
                                                                       pop.rpkm.df,
                                                                       pop.meta.df,
                                                                       title.str='Test',
                                                                       show_rownames=F,
                                                                       show_colnames    = F,
                                                                       annotation_legend = F){
  #sub.pop.col.vec=get.color.list.for.pheatmap(in.vec =as.character(meta.df$markers.cluster.groups))
  #sub.pop.col.vec=get.color.rainbow.list.for.pheatmap(in.vec =as.character(meta.df$markers.cluster.groups))
  markers.col.vec=get.color.rainbow.list.for.pheatmap( in.vec =
                                                         markers.df$category )
  markers.vec=rownames(markers.df)
  #rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  sc.rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df = rpkm.df,
                                                               sample.count = 1)
  sc.meta.df=meta.df[colnames(rpkm.df),]
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df)
  #sc.meta.df=order.df.by.specific.col(meta.df = sc.meta.df,
                                      #order.vec =std.stages.ordered.str )
  #sc.meta.df$sp.order=sub.populations.grp.order.names.vec[as.character((sc.meta.df$markers.cluster.groups))]
  #sc.meta.df=sc.meta.df[order(sc.meta.df$sp.order),]
  sc.rpkm.df=subset(sc.rpkm.df,select = rownames(sc.meta.df))
  sc.meta.df=subset(sc.meta.df,development.stage %in% c('R','LR','ET'))
  sc.rpkm.df=subset(sc.rpkm.df,select = rownames(sc.meta.df))
  #pop.5.rpkm.df=filter.rpkm.less.than.cutoff(df = pop.rpkm.df,rpkm.cutoff = 5,no.samples = 2)
  rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df = 
                                                              pop.rpkm.df,sample.count = 2)
  meta.df=pop.meta.df[colnames(rpkm.df),]
  meta.df=subset(meta.df,
                 development.stage %in% c('R','LR','ET'))
  rpkm.df=subset(rpkm.df,select = rownames(meta.df))
  #meta.df=order.df.by.specific.col(meta.df = meta.df,order.vec =std.stages.ordered.str )
  rpkm.df=subset(rpkm.df,select = rownames(meta.df))
  expr.markers.vec=intersect(markers.vec,rownames(sc.rpkm.df))
  if(length(expr.markers.vec)>=2){
    var.and.rifins.vec=intersect(expr.markers.vec,
                                 rownames(subset(markers.df,category %in% c('var'))))
    expr.markers.vec=var.and.rifins.vec
    sc.vsa.rpkm.df=sc.rpkm.df[expr.markers.vec,]
    #sc.vsa.rpkm.df[sc.vsa.rpkm.df<1]=0.0
    sc.vsa.rpkm.df=filter.none.expressed.samples(df = 
                                                   filter.none.expressed.genes(input.data = sc.vsa.rpkm.df))
    sc.vsa.rpkm.mat=as.matrix(sc.vsa.rpkm.df)
    #sc.pseudo=min(sc.vsa.rpkm.mat[sc.vsa.rpkm.mat>0])/2
    sc.pseudo=1
    log.sc.vsa.rpkm.mat=log2(sc.vsa.rpkm.mat+1)
    vsa.rpkm.df=rpkm.df[rownames(log.sc.vsa.rpkm.mat),]
    vsa.rpkm.df[vsa.rpkm.df<5]=0.0
    vsa.rpkm.df[is.na(vsa.rpkm.df)]=0.0
    vsa.rpkm.df=filter.none.expressed.genes(input.data = vsa.rpkm.df)
    vsa.rpkm.mat=as.matrix(vsa.rpkm.df)
    #pseudo.val=min(vsa.rpkm.mat[vsa.rpkm.mat!=0])/2
    pseudo.val=1
    log.vsa.rpkm.mat=log2(vsa.rpkm.mat+1)
    log.vsa.rpkm.mat=log.vsa.rpkm.mat[,rownames(meta.df)]
    col.annotation.df=subset(meta.df,select=development.stage)
    sc.col.annotation.df=subset(sc.meta.df,select=
                                  c(development.stage))
    #temp.sub.pop.col.vec=sub.pop.col.vec[unique(as.character(col.annotation.df$markers.cluster.groups))]
    row.annotation.df=markers.df[expr.markers.vec,]
    expr.categories.vec=unique(as.character(markers.df[expr.markers.vec,'category']))
    expr.categories.vec=markers.col.vec[expr.categories.vec]
    row.annotation.df=row.annotation.df[rownames(log.vsa.rpkm.mat),]
    row.annotation.df=subset(markers.df,select=category)
    pop.row.border.gaps.vec=as.character(markers.df[rownames(log.vsa.rpkm.mat),
                                                    'category'])
    pop.row.border.gaps.indices=get.gap.indices(categories.vec =
                                                  pop.row.border.gaps.vec )
    sc.row.border.gaps.vec=as.character(markers.df[rownames(log.sc.vsa.rpkm.mat),
                                                   'category'])
    sc.row.border.gaps.indices=get.gap.indices(categories.vec =
                                                 sc.row.border.gaps.vec )
    temp.df=data.frame(markers.df[rownames(log.vsa.rpkm.mat),])
    #markers.cluster.groups.col=
      #get.color.list.for.pheatmap(in.vec = as.character(sc.meta.df$markers.cluster.groups))
    annotation.col.list=list(category=expr.categories.vec,
                             development.stage=abbrev.std.stages.col.vec)
    sc.col.color.indices=get.gap.indices(categories.vec = 
                                        as.character(sc.meta.df[
                                          colnames(log.sc.vsa.rpkm.mat),
                                          'development.stage']))
    #color.ramp = colorRampPalette(c("darkgray","gray","orange",'red'))(4)
    color.ramp=colorRampPalette(c("darkgray","gray","yellow",'orange','red'))(5)
    #log.sc.vsa.rpkm.mat=as.matrix(filter.none.expressed.samples(df = filter.rpkm.less.than.cutoff(df = data.frame(log.sc.vsa.rpkm.mat),rpkm.cutoff = 9,no.samples = 1)))
    print(dim(log.vsa.rpkm.mat))
    pheatmap(mat = log.vsa.rpkm.mat,fontsize_row = 3,fontsize_col = 1,cellheight = 2,
             cellwidth = 5,main = title.str,color = color.ramp,
             annotation_colors = annotation.col.list,annotation_row = row.annotation.df,
             annotation_col = col.annotation.df,show_rownames = show_rownames,
             show_colnames=show_colnames,gaps_row = pop.row.border.gaps.indices,
             annotation_legend = T,cluster_cols = F,cluster_rows = F,
             border_color = NA)
    pheatmap(mat = log.vsa.rpkm.mat,fontsize_row = 2,fontsize_col = .5,cellheight = 2,
             color=color.ramp,cellwidth = 5,main = title.str,annotation_colors = annotation.col.list,
             annotation_row = row.annotation.df,annotation_col = col.annotation.df,
             show_colnames=show_colnames,annotation_legend = annotation_legend,
             cluster_cols = F,cluster_rows = F,border_color = NA)
    print(dim(log.sc.vsa.rpkm.mat))
    pheatmap(mat = log.sc.vsa.rpkm.mat,fontsize_row = 3,fontsize_col = 1,cellheight = 2,
             color = color.ramp,cellwidth = 3,main = title.str,
             annotation_colors = annotation.col.list,
             annotation_row = row.annotation.df,annotation_col = sc.col.annotation.df,
             show_rownames = show_rownames,show_colnames=show_colnames,
             annotation_legend = T,cluster_cols = F,
             cluster_rows = F,border_color = NA)
    pheatmap(mat = log.sc.vsa.rpkm.mat,fontsize_row = 2,fontsize_col = .5,cellheight = 4,
             color = color.ramp,cellwidth = 6,main = title.str,
             annotation_colors = annotation.col.list,
             annotation_row = row.annotation.df,annotation_col = sc.col.annotation.df,
             show_colnames=show_colnames,annotation_legend = annotation_legend,
             cluster_cols = F,cluster_rows = F,border_color = NA)
    #sc.vsa.rpkm.df=sc.rpkm.df[expr.markers.vec,]
    # sc.vsa.rpkm.df=sc.rpkm.df[rownames(log.sc.vsa.rpkm.mat),]
    # 
    # sc.vsa.rpkm.mat=as.matrix(sc.vsa.rpkm.df)
    # 
    # sc.pseudo.val=min(sc.vsa.rpkm.mat[sc.vsa.rpkm.mat!=0])/2
    # 
    # sc.log.vsa.rpkm.mat=log2(sc.vsa.rpkm.mat+sc.pseudo.val)
    # 
    # sc.log.vsa.rpkm.mat=sc.log.vsa.rpkm.mat[,rownames(sc.meta.df)]
    sc.barplot.col.vec=markers.col.vec[as.character(markers.df[rownames(log.sc.vsa.rpkm.mat),'category'])]
    barplot(as.numeric(apply(sc.vsa.rpkm.df,1,nnzero))/dim(sc.vsa.rpkm.df)[2],
            names.arg = rownames(sc.vsa.rpkm.df),las=2,cex.names = .2,
            col=sc.barplot.col.vec,border = NA)
    abline(h = 0.02,b = 1)
    temp.df=data.frame(freq=as.numeric(apply(vsa.rpkm.df,1,nnzero))/dim(vsa.rpkm.df)[2],
                       gene=rownames(vsa.rpkm.df),row.names = rownames(vsa.rpkm.df))
    temp.df=subset(temp.df,freq>=0.02)
    #temp.barplot.col.vec=markers.col.vec[as.character(markers.df[rownames(temp.df),'category'])]
    #barplot(as.numeric(temp.df$freq),names.arg = rownames(temp.df),las=2,cex.names = .2,col=temp.barplot.col.vec,border = NA)
    #abline(h = 0.02,b = 1)
  }
  return(log.sc.vsa.rpkm.mat)
}
plot.variant.surface.marker.genes.heatmap.in.single=function(sc.rpkm.df,sc.meta.df,markers.df,
                                                             title.str='Test',show_rownames=F,
                                                             show_colnames    = F,cellheight = 5,
                                                             cellwidth = 5,annotation_legend = F,
                                                             in.col.ramp=colorRampPalette(c("darkgray","gray","yellow",'orange','red'))(5),
                                                             log.expr.cut.off=0){
  markers.col.vec=get.color.rainbow.list.for.pheatmap( in.vec =markers.df$category )
  sc.rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df =  sc.rpkm.df,sample.count = 2)
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df)
  markers.vec=rownames(markers.df)
  sc.rpkm.df=subset(sc.rpkm.df,select = rownames(sc.meta.df))
  sc.meta.df=subset(sc.meta.df,development.stage %in% c('R','LR','ET'))
  sc.rpkm.df=subset(sc.rpkm.df,select = rownames(sc.meta.df))
  expr.markers.vec=intersect(markers.vec,rownames(sc.rpkm.df))
  if(length(expr.markers.vec)>=2){
    #var.and.rifins.vec=intersect(expr.markers.vec,rownames(subset(markers.df,category %in% c('var','rifin'))))
    var.and.rifins.vec=intersect(expr.markers.vec,rownames(subset(markers.df,category %in% c('var'))))
    expr.markers.vec=var.and.rifins.vec
    sc.vsa.rpkm.df=sc.rpkm.df[expr.markers.vec,]
    sc.vsa.rpkm.df=filter.none.expressed.samples(filter.none.expressed.genes(input.data = sc.vsa.rpkm.df))
    sc.vsa.rpkm.mat=as.matrix(sc.vsa.rpkm.df)
    sc.pseudo=1
    #sc.pseudo=min(sc.vsa.rpkm.mat[sc.vsa.rpkm.mat>0])/2
    log.sc.vsa.rpkm.mat=log2(sc.vsa.rpkm.mat+sc.pseudo)
    log.sc.vsa.rpkm.mat=as.matrix(filter.none.expressed.samples(filter.rpkm.less.than.cutoff(
      df = data.frame(log.sc.vsa.rpkm.mat),rpkm.cutoff =log.expr.cut.off ,no.samples = 1)))
    log.sc.vsa.rpkm.mat=log.sc.vsa.rpkm.mat[,colnames(log.sc.vsa.rpkm.mat)[apply(log.sc.vsa.rpkm.mat>=log.expr.cut.off,2,sum)>=1]]
    sc.col.annotation.df=subset(sc.meta.df,select=c(development.stage))
    row.annotation.df=markers.df[rownames(log.sc.vsa.rpkm.mat),]
    row.annotation.df=subset(row.annotation.df,select=category)
    expr.categories.vec=unique(as.character(markers.df[expr.markers.vec,'category']))
    expr.categories.vec=markers.col.vec[expr.categories.vec]
    annotation.col.list=list(category=expr.categories.vec,development.stage=abbrev.std.stages.col.vec)
    sc.col.color.indices=get.gap.indices(categories.vec = 
                                           as.character(sc.meta.df[
                                             colnames(log.sc.vsa.rpkm.mat),
                                             'development.stage']))
    sc.row.col.indices=get.gap.indices(categories.vec = as.character(row.annotation.df$category))
    color.ramp=in.col.ramp
    pheatmap(mat = log.sc.vsa.rpkm.mat,fontsize_row = 3,fontsize_col = 3,cellheight = cellheight,
             color = color.ramp,cellwidth = cellwidth,main = title.str,
             annotation_colors = annotation.col.list,annotation_row = row.annotation.df,
             annotation_col = sc.col.annotation.df,gaps_row = sc.row.col.indices,
             show_colnames=show_colnames,annotation_legend = annotation_legend,
             cluster_cols = F,cluster_rows = F,border_color = NA)
    return(log.sc.vsa.rpkm.mat)
  }
}
get.gap.indices=function(categories.vec){
  row.border.gaps=categories.vec
  len.row.border.gaps=length(row.border.gaps)
  row.border.gaps.indices=c()
  for(m in 1:len.row.border.gaps){
    if(m!=len.row.border.gaps & row.border.gaps[m]!=row.border.gaps[m+1]){
      row.border.gaps.indices=append(row.border.gaps.indices,values = m,
                                     after = length(row.border.gaps.indices))
    }
  }
  return(row.border.gaps.indices)
}
plot.variant.surface.marker.genes.heatmap.btw.sub.populations=function(rpkm.df,meta.df,markers.df,title.str='Test',show_rownames=F,show_colnames    = F,annotation_legend = F){
  sub.pop.col.vec=get.color.list.for.pheatmap(in.vec =as.character(meta.df$markers.cluster.groups))
  #sub.pop.col.vec=get.color.rainbow.list.for.pheatmap(in.vec =as.character(meta.df$markers.cluster.groups))
  markers.col.vec=get.color.rainbow.list.for.pheatmap( in.vec =markers.df$category )
  markers.vec=rownames(markers.df)
  #rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df = rpkm.df,sample.count = 2)
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  #meta.df=order.by.marker.sub.population.meta.df(meta.df = meta.df,order.vec = sub.populations.grp.order.names.vec)
  meta.df=order.df.by.specific.col(meta.df = meta.df,order.vec =std.stages.ordered.str )
  rpkm.df=subset(rpkm.df,select = rownames(meta.df))
  expr.markers.vec=intersect(markers.vec,rownames(rpkm.df))
  if(length(expr.markers.vec)>=2){
    vsa.rpkm.df=rpkm.df[expr.markers.vec,]
    vsa.rpkm.mat=as.matrix(vsa.rpkm.df)
    pseudo.val=min(vsa.rpkm.mat[vsa.rpkm.mat!=0])/2
    log.vsa.rpkm.mat=log2(vsa.rpkm.mat+pseudo.val)
    log.vsa.rpkm.mat=log2(vsa.rpkm.mat+1)
    log.vsa.rpkm.mat=log.vsa.rpkm.mat[,rownames(meta.df)]
    col.annotation.df=subset(meta.df,select=c(markers.cluster.groups,development.stage))
    temp.sub.pop.col.vec=sub.pop.col.vec[unique(as.character(col.annotation.df$markers.cluster.groups))]
    row.annotation.df=markers.df[expr.markers.vec,]
    expr.categories.vec=unique(as.character(markers.df[expr.markers.vec,'category']))
    expr.categories.vec=markers.col.vec[expr.categories.vec]
    row.annotation.df=row.annotation.df[rownames(log.vsa.rpkm.mat),]
    row.annotation.df=subset(markers.df,select=category)
    temp.df=data.frame(markers.df[rownames(log.vsa.rpkm.mat),])
    annotation.col.list=list(markers.cluster.groups=temp.sub.pop.col.vec,category=expr.categories.vec,development.stage=abbrev.std.stages.col.vec)
    pheatmap(mat = log.vsa.rpkm.mat,fontsize_row = 4,fontsize_col = 1,cellheight = 4.5,cellwidth = 4,main = title.str,annotation_colors = annotation.col.list,annotation_row = row.annotation.df,annotation_col = col.annotation.df,show_rownames = show_rownames,show_colnames=show_colnames,annotation_legend = F,cluster_cols = F,cluster_rows = F,border_color = NA)
    pheatmap(mat = log.vsa.rpkm.mat,fontsize_row = .5,fontsize_col = .5,cellheight = 3.5,cellwidth = 3.5,main = title.str,annotation_colors = annotation.col.list,annotation_row = row.annotation.df,annotation_col = col.annotation.df,show_colnames=show_colnames,annotation_legend = annotation_legend,cluster_cols = F,cluster_rows = F,border_color = NA)
  }
  return(expr.markers.vec)
}
plot.variant.surface.marker.genes.heatmap.btw.sub.groups.in.bulk.controls=function(rpkm.df,meta.df,markers.df,title.str='Test',show_rownames=F,show_colnames    = F,annotation_legend = F){
  sub.pop.col.vec=get.color.rainbow.list.for.pheatmap(in.vec =as.character(meta.df$sub.group))
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  #development.stage.col.vec=get.std.stage.cols(in.vec =as.character(meta.df$development.stage))
  #sub.pop.col.vec=get.color.rainbow.list.for.pheatmap(in.vec =as.character(meta.df$markers.cluster.groups))
  markers.col.vec=get.color.rainbow.list.for.pheatmap( in.vec =markers.df$category )
  markers.vec=rownames(markers.df)
  #rpkm.df=filter.none.expressed.genes(input.data = rpkm.df)
  rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df = rpkm.df,sample.count = 2)
  meta.df=meta.df[colnames(rpkm.df),]
  #meta.df=order.by.marker.sub.population.meta.df(meta.df = meta.df,order.vec = sub.populations.grp.order.names.vec)
  rpkm.df=subset(rpkm.df,select = rownames(meta.df))
  expr.markers.vec=intersect(markers.vec,rownames(rpkm.df))
  if(length(expr.markers.vec)>=2){
    vsa.rpkm.df=rpkm.df[expr.markers.vec,]
    vsa.rpkm.mat=as.matrix(vsa.rpkm.df)
    pseudo.val=min(vsa.rpkm.mat[vsa.rpkm.mat!=0])/2
    log.vsa.rpkm.mat=log2(vsa.rpkm.mat+pseudo.val)
    log.vsa.rpkm.mat=log.vsa.rpkm.mat[,rownames(meta.df)]
    sub.grp.order.vec=sort(unique(as.character(meta.df$sub.group)))
    #meta.df=order.df.by.specific.col(meta.df =meta.df,order.vec = sub.grp.order.vec )
    meta.df=order.df.by.specific.col(meta.df =meta.df,order.vec = std.stages.ordered.str )
    log.vsa.rpkm.mat=log.vsa.rpkm.mat[,rownames(meta.df)]
    col.annotation.df=subset(meta.df,select=c(sub.group,development.stage))
    temp.sub.pop.col.vec=sub.pop.col.vec[unique(as.character(col.annotation.df$sub.group))]
    row.annotation.df=markers.df[expr.markers.vec,]
    expr.categories.vec=unique(as.character(markers.df[expr.markers.vec,'category']))
    expr.categories.vec=markers.col.vec[expr.categories.vec]
    row.annotation.df=row.annotation.df[rownames(log.vsa.rpkm.mat),]
    row.annotation.df=subset(markers.df,select=category)
    temp.df=data.frame(markers.df[rownames(log.vsa.rpkm.mat),])
    annotation.col.list=list(sub.group=temp.sub.pop.col.vec,category=expr.categories.vec,development.stage=abbrev.std.stages.col.vec)
    pheatmap(mat = log.vsa.rpkm.mat,fontsize_row = 4,fontsize_col = 4,cellheight = 3,cellwidth = 5,main = title.str,annotation_colors = annotation.col.list,annotation_row = row.annotation.df,annotation_col = col.annotation.df,show_rownames = show_rownames,show_colnames=show_colnames,annotation_legend = F,cluster_cols = F,cluster_rows = F,border_color = NA)
    pheatmap(mat = log.vsa.rpkm.mat,fontsize_row = .5,fontsize_col = 4,cellheight = 3,cellwidth = 5,main = title.str,annotation_colors = annotation.col.list,annotation_row = row.annotation.df,annotation_col = col.annotation.df,show_colnames=show_colnames,annotation_legend = annotation_legend,cluster_cols = F,cluster_rows = F,border_color = NA)
  }
}
plot.human.vs.pf.detected.genes=function(human.counts.df,pf.counts.df,title.str='Detected genes',log.counts=F){
  human.counts.df=filter.none.expressed.samples(human.counts.df)
  hg.samples.counts.vec=apply(human.counts.df,2,nnzero)
  pf.samples.counts.vec=apply(pf.counts.df,2,nnzero)
  intersect.samples=intersect(names(hg.samples.counts.vec),names(pf.samples.counts.vec))
  len.samples=length(intersect.samples)
  counts.mat=matrix(nrow = len.samples,ncol = 2)
  for(n in 1:len.samples){
    sample=intersect.samples[n]
    counts.mat[n,]=c(hg.samples.counts.vec[sample],pf.samples.counts.vec[sample])
  }
  rownames(counts.mat)=intersect.samples
  colnames(counts.mat)=c('hg','pfa')
  cor.score=round(cor(counts.mat,method = 'spearman'),digits = 2)
  cor.score=unique(melt(cor.score[cor.score!=1]))
  show(cor.score)
  if(log.counts){
    plot(x=log10(as.numeric(counts.mat[,1])),y=log10(as.numeric(counts.mat[,2])),xlab='Host gene counts (log-transformed)',ylab='Parasite gene counts (log-transformed)',main=title.str,pch=19)
  }
  else{
    cor.str=paste('(Cor. coeff. = ',cor.score,')',collapse = '')
    temp.title.str=paste(title.str,cor.str,sep=' ')
    plot(x=as.numeric(counts.mat[,1]),y=as.numeric(counts.mat[,2]),xlab='Host gene counts',ylab='Parasite gene counts',main=temp.title.str,pch=19)
  }
}
plot.human.vs.pf.detected.genes.per.development.stage=function(human.counts.df,pf.counts.df,meta.df){
  pf.counts.df=filter.none.expressed.genes(filter.none.expressed.samples(df = pf.counts.df))
  meta.df=meta.df[colnames(pf.counts.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.list=split(meta.df,f =meta.df$development.stage)
  development.stages = names(meta.list)
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  development.stages = intersect(ordered.development.stages,development.stages)
  no.development.stage=length(development.stages)
  for (i in 1:no.development.stage){
    development.stage= development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.counts.df=pf.counts.df[,rownames(temp.meta.df)]
    temp.title.str=gsub(pattern = '\\.',replacement = ' ',x = convert.to.title.case(in.str = development.stage))
    plot.human.vs.pf.detected.genes(human.counts.df =human.counts.df,pf.counts.df =temp.counts.df ,title.str = temp.title.str )
  }
}
plot.human.vs.pf.read.counts.per.development.stage=function(human.counts.df,pf.counts.df,meta.df,log.counts=T){
  #pf.counts.df=filter.none.expressed.genes(filter.none.expressed.samples(df = pf.counts.df))
  meta.df=meta.df[colnames(pf.counts.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.list=split(meta.df,f =meta.df$development.stage)
  development.stages = names(meta.list)
  ordered.development.stages=c("R" ,"LR", "ET", "T" , "ES"  , "S")
  development.stages = intersect(ordered.development.stages,development.stages)
  no.development.stage=length(development.stages)
  out.list=list()
  for (i in 1:no.development.stage){
    development.stage= development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.pfa.counts.df=pf.counts.df[,intersect(colnames(pf.counts.df),rownames(temp.meta.df))]
    temp.human.counts.df=human.counts.df[,intersect(colnames(human.counts.df),rownames(temp.meta.df))]
    all.samples.vec=unique(c(colnames(temp.pfa.counts.df),colnames(temp.human.counts.df)))
    pfa.read.counts.vec=apply(temp.pfa.counts.df,2,sum)
    hg.read.counts.vec=apply(temp.human.counts.df,2,sum)
    all.samples.vec.len=length(all.samples.vec)
    out.mat=matrix(nrow = all.samples.vec.len,ncol = 2)
    for(m in 1:all.samples.vec.len){
      sample=all.samples.vec[m]
      pfa.count=ifelse(is.na(pfa.read.counts.vec[sample]),0,as.numeric(pfa.read.counts.vec[sample]))
      human.count=ifelse(is.na(hg.read.counts.vec[sample]),0,as.numeric(hg.read.counts.vec[sample]))
      out.mat[m,]=c(pfa.count,human.count)
    }
    rownames(out.mat)=all.samples.vec
    colnames(out.mat)=c('Pfa','hg')
    out.list[[development.stage]]=out.mat
    out.mat=out.mat+1
    temp.title.str= paste(abbrev.std.stages.timepoints.vec[development.stage],'hr',sep='')
    if(log.counts){
      barplot(height = t(out.mat),beside = T,main =temp.title.str,log = 'y',col=c('red','green'),border = NA,las=2,cex.names = .1 )
    }
    else{
      barplot(height = out.mat,beside = T,main =temp.title.str,col=c('red','green'),border = NA )
    }
  }
  return(out.list)
}
plot.human.vs.pf.read.counts=function(human.counts.df,pf.counts.df,title.str='Reads',log.counts=F){
  human.counts.df=filter.none.expressed.samples(human.counts.df)
  hg.samples.counts.vec=apply(human.counts.df,2,sum)
  pf.samples.counts.vec=apply(pf.counts.df,2,sum)
  intersect.samples=unique(c(colnames(human.counts.df),colnames(pf.counts.df)))
  len.samples=length(intersect.samples)
  counts.mat=matrix(nrow = len.samples,ncol = 2)
  for(n in 1:len.samples){
    sample=intersect.samples[n]
    hg.count=ifelse(is.na(hg.samples.counts.vec[sample]),hg.samples.counts.vec[sample],0)
    pfa.count=ifelse(is.na(pf.samples.counts.vec[sample]),pf.samples.counts.vec[sample],0)
    counts.mat[n,]=c(hg.count, pfa.count)
  }
  rownames(counts.mat)=intersect.samples
  colnames(counts.mat)=c('hg','pfa')
  cor.score=round(cor(counts.mat,method = 'spearman'),digits = 2)
  cor.score=unique(melt(cor.score[cor.score!=1]))
  show(cor.score)
  if(log.counts){
    plot(x=log10(as.numeric(counts.mat[,1])),y=log10(as.numeric(counts.mat[,2])),xlab='Host gene counts (log-transformed)',ylab='Parasite gene counts (log-transformed)',main=title.str,pch=19)
  }
  else{
    cor.str=paste('(Cor. coeff. = ',cor.score,')',collapse = '')
    temp.title.str=paste(title.str,cor.str,sep=' ')
    plot(x=as.numeric(counts.mat[,1]),y=as.numeric(counts.mat[,2]),xlab='Host gene counts',ylab='Parasite gene counts',main=temp.title.str,pch=19)
  }
}
plot.human.vs.pf.detected.genes.per.sub.population=function(human.counts.df,pf.counts.df,meta.df){
  pf.counts.df= pf.counts.df[,rownames(meta.df)]
  meta.list=split(meta.df,f =meta.df$markers.cluster.groups)
  development.stages = names(meta.list)
  #ordered.development.stages=c("ring" ,"late.ring", "early.trophozoite", "trophozoite" , "early.schizont"  , "schizont")
  development.stages = intersect(names(sub.populations.grp.order.vec),development.stages)
  no.development.stage=length(development.stages)
  for (i in 1:no.development.stage){
    development.stage= development.stages[i]
    temp.meta.df=meta.list[[development.stage]]
    temp.counts.df=pf.counts.df[,rownames(temp.meta.df)]
    temp.title.str=gsub(pattern = '\\.',replacement = ' ',x = convert.to.title.case(in.str = development.stage))
    plot.human.vs.pf.detected.genes(human.counts.df =human.counts.df,pf.counts.df =temp.counts.df ,title.str = temp.title.str )
  }
}
get.highly.expressed.genes=function(rpkm.df,log.rpkm=F,top.cut.off=100){
  mean.expr=as.numeric(apply(rpkm.df,1,median))
  no.samples=dim(rpkm.df)[2]
  fraction.detected.in=as.numeric(apply(rpkm.df,1,nnzero))/no.samples
  if(log.rpkm){
    pseudo.val=min(rpkm.df[rpkm.df>0])/2
    rpkm.df=log2(rpkm.df+pseudo.val)
    mean.expr=as.numeric(apply(rpkm.df,1,mean))
  }
  detection.score=mean.expr*fraction.detected.in
  out.df=data.frame(mean.expr=mean.expr ,fraction.detected.in=fraction.detected.in,detection.score=detection.score)
  rownames(out.df)=rownames(rpkm.df)
  #detection.score.cut.off=as.numeric(quantile(x = detection.score,probs = quantile.cut.off))
  out.df=out.df[with(out.df,order(detection.score,decreasing = T)),]
  out.df=head(out.df,top.cut.off)
  return(out.df)
}
#
plot.heam.vs.parasite.read.counts=function(heam.counts.df,pfa.counts.df,meta.df,per.stage=T,title.str='Test plot'){
  samples=intersect(intersect(colnames(heam.counts.df),colnames(pfa.counts.df)),rownames(meta.df))
  heam.counts.df=subset(heam.counts.df,select = samples)
  pfa.counts.df=subset(pfa.counts.df,select=samples)
  meta.df=meta.df[samples,]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  meta.df$development.stage=factor(meta.df$development.stage,level=c('R','LR','ET','T','ES','S'))
  stages.order=c('R'=1,'LR'=2,'ET'=3,'T'=4,'ES'=5,'S'=6)
  meta.df$development.stage.order=stages.order[as.character(meta.df$development.stage)]
  #meta.df$development.stage=factor(meta.df$development.stage,levels = factor(c('R','LR','ET','T','ES','S')))
  meta.df=meta.df[with(meta.df,order(development.stage.order)),]
  meta.list=split(meta.df,f=meta.df$development.stage)
  development.stages.vec=c('R','LR','ET','T','ES','S')
  development.stages.vec.len=length(development.stages.vec)
  if(per.stage){
    for(n in 1:development.stages.vec.len){
      temp.stage=development.stages.vec[n]
      temp.meta.df=meta.list[[temp.stage]]
      temp.human.counts.df=subset(heam.counts.df,select=rownames(temp.meta.df))
      temp.pfa.counts.df=subset(pfa.counts.df,select=rownames(temp.meta.df))
      #col.str=abbrev.std.stages.col.vec[as.character(meta.df$development.stage)]
      heam.vec=rownames(temp.human.counts.df)
      heam.vec.len=length(heam.vec)
      expr.vec.list=list()
      for(m in 1:heam.vec.len){
        heam=heam.vec[m]
        temp.title.str=paste(abbrev.std.stages.timepoints.vec[temp.stage],'hr',' (',heam,')',sep = '')
        x.points=as.numeric(temp.human.counts.df[heam,])
        y.points=as.numeric(apply(temp.pfa.counts.df,2,sum))
        temp.df=data.frame(heam=x.points,pfa=y.points)
        heam.expr.df=subset(temp.df,heam>0)
        heam.not.expr.df=subset(temp.df,heam==0)
        heam.expr.vec=as.numeric(heam.expr.df$pfa)
        heam.expr.vec=ifelse(length(heam.expr.vec)==0,1,heam.expr.vec)
        heam.not.expr.vec=as.numeric(heam.not.expr.df$pfa)
        heam.not.expr.vec=ifelse(length(heam.not.expr.vec)==0,1,heam.not.expr.vec)
        #plot(y.points~jitter(x.points,2),pch=19,xlab='Heamoglobin [log2(counts+1)]',ylab='Parasite [log2(counts+1)]',main=temp.title.str,cex=1)
        plot(x=x.points,y=y.points,pch=19,xlab='Heamoglobin [log2(counts+1)]',ylab='Parasite [log2(counts+1)]',main=temp.title.str,cex=1)
        #plot(x=x.points,y=y.points,pch=19,xlab='Heamoglobin [log2(counts+1)]',ylab='Parasite [log2(counts+1)]',main=temp.title.str,cex=1)
      }
      #boxplot(x = expr.vec.list,log = 'y')
    }
  }
  else{
      temp.human.counts.df=subset(heam.counts.df,select=rownames(meta.df))
      temp.pfa.counts.df=subset(pfa.counts.df,select=rownames(meta.df))
      col.str=abbrev.std.stages.col.vec[as.character(meta.df$development.stage)]
      heam.vec=rownames(temp.human.counts.df)
      heam.vec.len=length(heam.vec)
      for(m in 1:heam.vec.len){
        heam=heam.vec[m]
        temp.title.str=paste(title.str,' (',heam,')',sep = '')
        x.points=log2(as.numeric(temp.human.counts.df[heam,])+1)
        y.points=log2(as.numeric(apply(temp.pfa.counts.df,2,sum))+1)
        #plot(x=x.points,y=y.points,pch=19,xlab='Heamoglobin [log2(counts+1)]',ylab='Parasite [log2(counts+1)]',main=temp.title.str,cex=1)
        plot(x=x.points,y=y.points,col=col.str,pch=19,xlab='Heamoglobin [log2(counts+1)]',ylab='Parasite [log2(counts+1)]',main=temp.title.str,cex=1)
        legend('bottomright',legend = abbrev.std.stages.timepoints.vec[names(abbrev.std.stages.col.vec)],fill = abbrev.std.stages.col.vec,box.lty = 0,cex=1.0,border = NA)
      }
    }
}
plot.total.heam.vs.parasite.read.counts=function(heam.counts.df,pfa.counts.df,meta.df){
  samples=intersect(intersect(colnames(heam.counts.df),colnames(pfa.counts.df)),rownames(meta.df))
  heam.counts.df=subset(heam.counts.df,select = samples)
  pfa.counts.df=subset(pfa.counts.df,select=samples)
  meta.df=meta.df[samples,]
  show(dim(meta.df))
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  stages.order=c('R'=1,'LR'=2,'ET'=3,'T'=4,'ES'=5,'S'=6)
  meta.df$development.stage.order=stages.order[as.character(meta.df$development.stage)]
  meta.df$development.stage=factor(meta.df$development.stage,levels = factor(c('R','LR','ET','T','ES','S')))
  meta.df=meta.df[with(meta.df,order(development.stage.order)),]
  show(as.character(meta.df$development.stage))
  col.str=abbrev.std.stages.col.vec[as.character(meta.df$development.stage)]
  heam.vec=rownames(heam.counts.df)
  heam.vec.len=length(heam.vec)
  x.points=log2(as.numeric(apply(heam.counts.df,2,sum))+1)
  y.points=log2(as.numeric(apply(pfa.counts.df,2,sum))+1)
  plot(x=x.points,y=y.points,col=col.str,pch=19,xlab='Heamoglobin [log2(counts+1)]',ylab='Parasite [log2(counts+1)]',main='All heamoglobin variants vs. Parasite')
}
get.optimal.clusters.from.expr.df=function(rpkm.df,title.str='Optimal clusters',log.rpkm=T){
  psuedo.value=min(rpkm.df[rpkm.df>0]/2)
  #rpkm.df=filter.none.expressed.genes(filter.none.expressed.samples(df = filter.genes.with.zero.variance(df = rpkm.df)))
  rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df = rpkm.df,sample.count = 2)
  var.cut.off=as.numeric(quantile(as.numeric(apply(rpkm.df,1,sd)),probs = .3))
  rpkm.df=filter.none.expressed.genes(filter.none.expressed.samples(df = filter.non.variable.rows(df = rpkm.df,cut.off = var.cut.off)))
  #rpkm.df=filter.none.expressed.genes(filter.none.expressed.samples(df = filter.genes.with.zero.variance(df = rpkm.df)))
  if(log.rpkm){
    rpkm.df=log2(rpkm.df+psuedo.value)
  }
  pca.obj=prcomp(x = t(rpkm.df),retx = T,scale. = T,center = T)
  mydata=pca.obj$x[,1:3]
  mydata=t(rpkm.df)
  wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  for (i in 2:20) wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
  plot(1:20, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares",main=title.str,pch=19) 
}
plot.kmeans.clusters=function(rpkm.df,no.clusters=2,no.components=2,log.rpkm=T,first.pc='PC1',second.pc='PC2'){
  psuedo.value=min(rpkm.df[rpkm.df>0]/2)
  rpkm.df=filter.none.expressed.genes(filter.none.expressed.samples(df = filter.genes.with.zero.variance(df = rpkm.df)))
  if(log.rpkm){
    rpkm.df=log2(rpkm.df+psuedo.value)
  }
  pca.obj=prcomp(x = t(rpkm.df),retx = T,scale. = T,center = T)
  mydata=pca.obj$x[,1:no.components]
  mydata.cl <- kmeans(mydata, no.clusters)
  plot(mydata[,c(first.pc,second.pc)], col = mydata.cl$cluster,pch=19,cex=1.5)
  #points(mydata.cl$centers, col = 'orange', pch = 19, cex = 2)
  pairs(mydata,col = mydata.cl$cluster,pch=19)
  return(mydata.cl)
}
filter.genes.with.zero.variance=function(df){
  df=data.frame(df)
  df[,'sd']=as.numeric(as.character(apply(df,1,sd)))
  out.df=df[which(df[,'sd']>0),]
  out.df=data.frame(subset(out.df,select=-sd))
  return(out.df)
}
temp_function=function(rpkm.df,count.df,meta.df,gene.symbols.list){
  gene.names.vec=intersect(rownames(rpkm.df),rownames(count.df))
  rpkm.df=rpkm.df[gene.names.vec,]
  count.df=count.df[gene.names.vec,]
  samples.names=intersect(colnames(rpkm.df),colnames(count.df))
  count.df=subset(count.df,select = samples.names)
  rpkm.df=subset(rpkm.df,select = samples.names)
  x.vec=log2(as.numeric(apply(count.df,1,mean))+1)
  y.vec=log2(as.numeric(apply(rpkm.df,1,mean))+1)
  detection.rate.vec=(as.numeric(apply(rpkm.df,1,nnzero)))/dim(rpkm.df)[2]
  out.df=data.frame(mean.count=x.vec,mean.rpkm=y.vec,detection=detection.rate.vec)
  rownames(out.df)=rownames(rpkm.df)
  gene.symbols=as.character(unlist(lapply(rownames(out.df),function(gene.id){
    gene.id=as.character(gene.id)
    temp.str=strsplit(gene.id,split = '\\.')
    temp.list=lapply(temp.str,function(temp.symbol){
      temp.symbol=as.character(temp.symbol)
      return(gene.symbols.list[temp.symbol])
      })
    return(paste(unique(as.character(unlist(temp.list))),collapse = '+'))
  })))
  out.df$gene.name=gene.symbols
  plot(x=x.vec,y=y.vec,pch=19,cex=.3,xlab='counts',ylab='rpkm',main='Human transcriptome')
  text(x=x.vec,y=y.vec,labels =gene.names.vec,cex = .3 )
  return(out.df)
}
plot.highly.expressed.human.genes=function(rpkm.df,meta.df,gene.symbols.list,sample.cut.off=10,title.str='Test plot'){
  rpkm.df=filter.genes.not.expressed.in.a.number.of.samples(df = rpkm.df,sample.count = sample.cut.off)
  samples=intersect(rownames(meta.df),colnames(rpkm.df))
  rpkm.df=subset(rpkm.df,select=samples)
  meta.df=data.frame(t(subset(t(meta.df),select=samples)))
  samples.names=colnames(rpkm.df)
  gene.names.vec=rownames(rpkm.df)
  y.vec=log2(as.numeric(apply(rpkm.df,1,mean))+1)
  detection.rate.vec=(as.numeric(apply(rpkm.df,1,nnzero)))/dim(rpkm.df)[2]
  out.df=data.frame(mean.rpkm=y.vec,detection=detection.rate.vec)
  rownames(out.df)=rownames(rpkm.df)
  gene.symbols=as.character(unlist(lapply(gene.names.vec,function(gene.id){
    gene.id=as.character(gene.id)
    temp.str=strsplit(gene.id,split = '\\.')
    temp.list=lapply(temp.str,function(temp.symbol){
      temp.symbol=as.character(temp.symbol)
      return(gene.symbols.list[temp.symbol])
    })
    return(paste(unique(as.character(unlist(temp.list))),collapse = '+'))
  })))
  out.df$gene.name=gene.symbols
  temp.gene.list=gene.symbols
  names(temp.gene.list)=rownames(rpkm.df)
  out.df=out.df[order(out.df$mean.rpkm,decreasing = T),]
  top.out.df=head(out.df,n=10)
  top.genes.vec=rownames(top.out.df)
  top.rpkm.df=rpkm.df[top.genes.vec,]
  top.genes.symbols.vec=as.character(temp.gene.list[top.genes.vec])
  meta.df=meta.df[colnames(rpkm.df),]
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  temp.rpkm.mat=as.matrix(rpkm.df)
  pseudo.value=min(as.numeric(temp.rpkm.mat[temp.rpkm.mat!=0]))/2
  temp.rpkm.mat=log2(temp.rpkm.mat+pseudo.value)
  rownames(top.rpkm.df)=top.genes.symbols.vec
  values.vec=c()
  mean.vec=c()
  se.vec=c()
  col.vec=c()
  len.top.genes.symbols.vec=length(top.genes.symbols.vec)
  for (m in 1:len.top.genes.symbols.vec){
    temp.gene=top.genes.symbols.vec[m]
    expr.vec=as.numeric(top.rpkm.df[temp.gene,])
    expr.stats=describe(expr.vec)
    values.vec=append(x = values.vec,values =m ,after = length(values.vec))
    mean.vec=append(x = mean.vec,values = expr.stats$mean,after = length(mean.vec))
    se.vec=append(x = se.vec,values =expr.stats$se ,after = length(se.vec))
  }
  marker.stats.df=data.frame(values=values.vec,mean=mean.vec,se=se.vec)
  rownames(marker.stats.df)=top.genes.symbols.vec
  error.bars(stats = marker.stats.df,bars=T,labels = top.genes.symbols.vec,main = 'Hg19 genes expression',pos = 'above', ylab = 'RPKM',xlab='Genes')
  plot.specific.genes.expression.boxplot(rpkm.df = top.rpkm.df,meta.df = meta.df,gene.vec = top.genes.symbols.vec)
  return(out.df)
}
get.highest.correlated.samples.to.bulk=
  function(seq.bulk.rpkm.df,seq.bulk.meta.df,plos.one.exp.df,plos.one.meta.df){
    seq.bulk.rpkm.df=filter.none.expressed.genes(input.data=seq.bulk.rpkm.df)
    seq.bulk.rpkm.df=log2(seq.bulk.rpkm.df+min(seq.bulk.rpkm.df[seq.bulk.rpkm.df>0])/2)
    intersect.genes=intersect(rownames(seq.bulk.rpkm.df),rownames(plos.one.exp.df))
    seq.bulk.rpkm.df=seq.bulk.rpkm.df[intersect.genes,]
    plos.one.exp.df=plos.one.exp.df[intersect.genes,]
    combined.df=data.frame(seq.bulk.rpkm.df,plos.one.exp.df)
    seq.bulk.samples=colnames(seq.bulk.rpkm.df)
    plos.one.samples=colnames(plos.one.exp.df)
    corr.res=c()
    for (m in seq.bulk.samples){
      temp.cor.list=list()
      temp.samples=plos.one.samples
      for (n in plos.one.samples){
        temp.df=filter.none.expressed.genes(input.data = subset(combined.df,select=c(m,n)))
        temp.cor.list=append(temp.cor.list,values = 
                               round(cor(x=temp.df[,m],y=temp.df[,n],method = 'spearman'),digits = 3),
                             after = length(temp.cor.list))
      }
      names(temp.samples)=temp.cor.list
      top.cor=sort(names(temp.samples),decreasing = T)[1:5]
      top.cor.cts=table(plos.one.meta.df[as.character(temp.samples[top.cor]),'grp'])
      top.cor.cts=top.cor.cts[top.cor.cts>0]
      s.names=paste(names(top.cor.cts),collapse=',')
      s.counts=paste(as.character(top.cor.cts),collapse=',')
      s.corr=paste(top.cor,collapse=',')
      s.name.cts=paste(s.names,s.counts,s.corr,sep = ':')
      corr.res=append(corr.res,values = s.name.cts,after = length(corr.res))
    }
    seq.bulk.meta.df=seq.bulk.meta.df[seq.bulk.samples,]
    seq.bulk.meta.df$plos.one.class=corr.res
    return(seq.bulk.meta.df)
}
encode.timepoint=function(meta.df){
  meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = meta.df)
  temp.vec=as.character(meta.df$development.stage)
  temp.vec=gsub(pattern = 'LR',replacement = '12-16hr',x = temp.vec)
  temp.vec=gsub(pattern = 'R',replacement = '8-12hr',x = temp.vec)
  temp.vec=gsub(pattern = 'ET',replacement = '20-24hr',x = temp.vec)
  temp.vec=gsub(pattern = 'T',replacement = '30-34hr',x = temp.vec)
  temp.vec=gsub(pattern = 'ES',replacement = '36-40hr',x = temp.vec)
  temp.vec=gsub(pattern = 'S',replacement = '44-48hr',x = temp.vec)
  meta.df$timepoint=temp.vec
  return(meta.df)
}
calculate.reads.mapped.rRNA=function(expression.df,rRNA.genes.vec){
  rRNA.genes.vec=intersect(rownames(expression.df),rRNA.genes.vec)
  other.genes.vec=setdiff(rownames(expression.df),rRNA.genes.vec)
  mapping.df=data.frame(apply(expression.df,2,function(expr.vec){
    ribosomal.genes=sum(expr.vec[rRNA.genes.vec])
    other.genes=sum(expr.vec[other.genes.vec])
    total=ribosomal.genes+other.genes
    out.vec=c(ribosomal.genes/total,other.genes/total)
    return(out.vec)
  }),row.names = c('ribosomal','others'))
  mapping.mat=as.matrix(mapping.df*100)
  barplot2(mapping.mat,col=c('red','blue'),las=2,border = NA,
           main = 'Ribosomal mapping',cex.names = .3)
  plot.new()
  legend('topright',legend=c('Ribosomal','Others'),cex = 1.5,
         fill=c('blue','red'),box.lty = 0)
}
plot.qpcr.expression=function(qpcr.df){
  qpcr.expt.list=split(qpcr.df,f = qpcr.df$Experiment)
  expt.list=names(qpcr.expt.list)
  expt.list=c("First" ,"Second" ,"third", "fourth")
  ref.gene="Seryl tRNA Synt"
  target.genes.vec=setdiff(qpcr.df$Gene,ref.gene)
  out.list=list()
  for (m in 1:length(expt.list)){
    temp.qpcr.expt.df=qpcr.expt.list[[expt.list[m]]]
    temp.qpcr.sample.list=split(temp.qpcr.expt.df,temp.qpcr.expt.df$Sample)
    temp.qpcr.sample.names=names(temp.qpcr.sample.list)
    for (n in 1:length(temp.qpcr.sample.names)){
      temp.qpcr.expt.grp.df=temp.qpcr.sample.list[[temp.qpcr.sample.names[n]]]
      temp.ref.gene.df=subset(temp.qpcr.expt.grp.df,Gene ==ref.gene)
      ref.mean.cq=mean(temp.ref.gene.df$Cq)
      for (g in 1:length(target.genes.vec)){
        temp.target.gene=target.genes.vec[g]
        temp.target.gene.df=subset(temp.qpcr.expt.grp.df,Gene ==temp.target.gene)
        if(isEmpty(temp.target.gene.df)){
          next
        }
        temp.target.gene.df$Cq.change=temp.target.gene.df$Cq-ref.mean.cq
        temp.replicate.mean.delta.cq=mean(temp.target.gene.df$Cq.change)
        temp.replicate.sd.delta.cq=sd(temp.target.gene.df$Cq.change)
        temp.replicate.delta.Cq=2^-temp.replicate.mean.delta.cq
        list.name=paste(c(temp.target.gene,expt.list[m],temp.qpcr.sample.names[n]),collapse = ',')
        out.list[[list.name]]=c(temp.target.gene,expt.list[m],temp.qpcr.sample.names[n],
                                temp.replicate.mean.delta.cq,temp.replicate.sd.delta.cq,
                                temp.replicate.delta.Cq)
      }
    }
  }
  cq.res.df=data.frame(t(data.frame(out.list)))
  colnames(cq.res.df)=c('Gene','Expriment','Sample.group','mean.Cq','SD.Cq',
                        'delta.Cq')
  print(head(cq.res.df))
  cq.res.list=split(cq.res.df,cq.res.df$Expriment)
  for (m in 1:length(expt.list)){
    temp.cq.res.df=cq.res.list[[expt.list[m]]]
    temp.cq.res.list=split(temp.cq.res.df,temp.cq.res.df$Sample.group)
    temp.cq.res.list.names=names(temp.cq.res.list)
    temp.grps.combn.mat=combn(temp.cq.res.list.names,m = 2)
    temp.grps.combn.col.ctns=dim(temp.grps.combn.mat)[2]
    for(d in 1:temp.grps.combn.col.ctns){
      temp.grps.category=temp.grps.combn.mat[,d]
      temp.grps.one=temp.grps.category[1]
      temp.grps.two=temp.grps.category[2]
      temp.grps.one.df=temp.cq.res.list[[temp.grps.one]]
      temp.grps.two.df=temp.cq.res.list[[temp.grps.two]]
      if (isEmpty(temp.grps.one.df)|isEmpty(temp.grps.two.df)){
        next
      }
      temp.fold.changes.list=c()
      temp.gene.names=c()
      for (r in 1:length(target.genes.vec)){
        temp.target.gene=target.genes.vec[r]
        temp.grps.one.target.gene.df=subset(temp.grps.one.df,Gene ==temp.target.gene)
        temp.grps.two.target.gene.df=subset(temp.grps.two.df,Gene ==temp.target.gene)
        if(isEmpty(temp.grps.one.target.gene.df) | isEmpty(temp.grps.two.target.gene.df)){
          next
        }
        temp.gene.names=append(temp.gene.names,values = temp.target.gene,after=length(temp.gene.names))
        fold.change=log2(mean(as.numeric(as.character(temp.grps.two.target.gene.df[,6])))/mean(as.numeric(as.character(temp.grps.one.target.gene.df[,6]))))
        temp.fold.changes.list=append(temp.fold.changes.list,values = fold.change,after = length(temp.fold.changes.list))
      }
      title.str=paste(c(temp.grps.two,'vs',temp.grps.one),sep=' ')
      barplot2(temp.fold.changes.list,main =title.str,ylab='log10(fold change)',
               border = NA ,xlab =expt.list[m] ,names.arg = temp.gene.names,las=2,
               cex.names = .5,col = 'lightgreen')
    }
    }
}
#' Convert counts to transcripts per million (TPM).
#' 
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#' 
#'    Lior Pachter. Models for transcript quantification from RNA-Seq.
#'    arXiv:1104.3889v2 
#'    
#'    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
#'    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
#'    doi:10.1007/s12064-012-0162-3
#'    
#' @param counts A numeric matrix of raw feature counts i.e.
#'  fragments assigned to each gene.
#' @param featureLength A numeric vector with feature lengths.
#' @param meanFragmentLength A numeric vector with mean fragment lengths.
#' @return tpm A numeric matrix normalized by library size and feature length.
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}
edit.gene.names=function(expr.df,plasmodb.gene.id.list){
  all.genes.names=rownames(expr.df)
  plasmodb.gene.id.names=names(plasmodb.gene.id.list)
  new.names=c()
  for (m in 1:length(all.genes.names)){
    temp.gn=all.genes.names[m]
    if(temp.gn %in% plasmodb.gene.id.names){
      new.names=append(x = new.names,values = temp.gn,after = length(new.names))
      next
    }
    id.parts=unlist(strsplit(x=temp.gn,split = '\\.'))
    id.parts=unique(id.parts[!grepl(pattern = '^\\d',x = id.parts)])
    if(length(id.parts)==1){
      new.names=append(x = new.names,values = as.character(id.parts),after = length(new.names))
      next
    }
    id.part=plasmodb.gene.id.list[id.parts][plasmodb.gene.id.list[id.parts]==max(plasmodb.gene.id.list[id.parts])]
    if(all(!is.na(id.part))){
      new.names=append(x = new.names,values = as.character(names(id.part)),after = length(new.names))
      next
    }
    new.names=append(x = new.names,values = names(plasmodb.gene.id.list[id.parts][!is.na(plasmodb.gene.id.list[id.parts])]),
                     after = length(new.names))
  }
  rownames(expr.df)=new.names
  return(expr.df)
}
temp.function=function(plasmodb.gene.sized.df){
  all.genes.names=rownames(expr.df)
  plasmodb.gene.id.names=names(plasmodb.gene.id.list)
  new.names=c()
  for (m in 1:length(all.genes.names)){
    temp.gn=all.genes.names[m]
    if(temp.gn %in% plasmodb.gene.id.names){
      new.names=append(x = new.names,values = temp.gn,after = length(new.names))
      next
    }
    id.parts=unlist(strsplit(x=temp.gn,split = '\\.'))
    id.parts=unique(id.parts[!grepl(pattern = '^\\d',x = id.parts)])
    if(length(id.parts)==1){
      new.names=append(x = new.names,values = as.character(id.parts),after = length(new.names))
      next
    }
    id.part=plasmodb.gene.id.list[id.parts][plasmodb.gene.id.list[id.parts]==max(plasmodb.gene.id.list[id.parts])]
    if(!is.na(id.part)){
      new.names=append(x = new.names,values = as.character(names(id.part)),after = length(new.names))
      next
    }
    new.names=append(x = new.names,values = names(plasmodb.gene.id.list[id.parts][!is.na(plasmodb.gene.id.list[id.parts])]),after = length(new.names))
  }
  rownames(expr.df)=new.names
  return(expr.df)
}
run.pca=function(expression.df,log.expression=F,pseudo.value=1){
  samples=as.character(colnames(expression.df))
  expression.df=filter.genes.with.zero.variance(df = expression.df)
  if(log.expression){
    log.expression.df=log2(expression.df+pseudo.value)
    expression.df=log.expression.df
  }
  samples.pca=prcomp(t(expression.df),retx=T,center=T,scale.=T)
  return(samples.pca)
}
#Plots PCA object and colors and shapes according the meta columns
plot.pca.object=function(pca.object,meta.df=NULL,plot.pch=19,
                         first.pc='PC1',second.pc='PC2',
                         plot.pairs=F,pairs.cut.off=5,
                         plot.variance.bar=F,show.legend=T,
                         col.factor='project',title.str='PCA plot',
                         cex.factor=1.2,label.sample.names=F){
  pca.summary=summary(pca.object)$importance
  if(plot.variance.bar){
    ylim.vec=c(0,max(100*pca.summary[2,])+10)
    barplot(100*pca.summary[2,][1:20],
            main='Components explained variance',
            ylab='Proportion(%)',las=2,cex.names = .5,ylim = ylim.vec)
    abline(h=100/length(colnames(pca.summary)),
           b=1)
  }
  first.pc.explained.var=
    round(as.numeric(as.character(pca.summary[2,
                                              first.pc]))*100,2)
  second.pc.explained.var=
    round(as.numeric(as.character(pca.summary[2,
                                              second.pc]))*100,2)
  pca.scores.df=pca.object$x
  samples=rownames(pca.scores.df)
  meta.df=meta.df[samples,]
  # kmeans.clusters=kmeans(x =pca.scores.df[,1:6],centers = 3 )
  # 
  # kmeans.clusters=kmeans.clusters$cluster
  # 
  # meta.df$kmeans.cluster=kmeans.clusters[samples]
  pc.scores.meta.df=data.frame(pca.scores.df,meta.df)
  x.points=as.numeric(as.character(pc.scores.meta.df[,first.pc]))
  y.points=as.numeric(as.character(pc.scores.meta.df[,second.pc]))
  first.pc.explained.var=paste('(Explained var: ',
                               paste(first.pc.explained.var,'%)',
                                     sep=''),sep='')
  second.pc.explained.var=paste('(Explained var: ',
                                paste(second.pc.explained.var,'%)',
                                      sep=''),sep='')
  first.pc.lab=paste(first.pc,first.pc.explained.var,'')
  second.pc.lab=paste(second.pc,second.pc.explained.var,'')
  col.factor.vec=as.character(pc.scores.meta.df[,col.factor])
  col.fact.map=get.rainbow.color.vec(in.vec = col.factor.vec)
  pca.col.vec=col.fact.map[col.factor.vec]
  plot(x=x.points,y=y.points,col=pca.col.vec,frame.plot = F,
       xlab=first.pc.lab,ylab=second.pc.lab,
       main=title.str,pch=plot.pch,cex=cex.factor)
  if(label.sample.names){
    text(x=x.points,y=y.points,labels=rownames(pc.scores.meta.df),cex = .3)
  }
  if(show.legend){
    legend('topright',legend =names(col.fact.map) ,cex = 0.5,
           fill=col.fact.map,border = NA,box.lty = 0)
  }
  if(plot.pairs){
    pca.pair.mat=pca.scores.df
    total.no.pc=ncol(pca.scores.df)
    pairs.cut.off=ifelse(total.no.pc>=pairs.cut.off,pairs.cut.off,total.no.pc)
    pca.pair.mat=pca.scores.df[,1:pairs.cut.off]
    pairs(pca.pair.mat,pch=plot.pch,col=pca.col.vec,
          main=title.str)
  }
  return(pc.scores.meta.df)
}
plot.all.cells.average.cell.to.pop.cor.per.timepoint=function(sc.rpkm.df,sc.meta.df,pop.rpkm.df,pop.meta.df,
                                                              in.col.ramp=NULL,out.dir='/Users/mtanga/Downloads/aggregate_sc_vs_bulk_correlation_analysis'){
  sc.rpkm.df=filter.none.expressed.samples(input.data = filter.none.expressed.genes(input.data = sc.rpkm.df))
  pop.rpkm.df=filter.none.expressed.samples(input.data = filter.none.expressed.genes(input.data = pop.rpkm.df))
  sc.meta.df=sc.meta.df[colnames(sc.rpkm.df),]
  pop.meta.df=pop.meta.df[colnames(pop.rpkm.df),]
  sc.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = sc.meta.df)
  pop.meta.df=change.development.stage.to.abbrev.stage.meta.df(meta.df = pop.meta.df)
  sc.timepoints=c(sc.meta.df$development.stage) 
  pop.timepoints=c(pop.meta.df$development.stage)
  sc.meta.list=split(x=sc.meta.df,f=sc.meta.df$development.stage)
  pop.meta.list=split(x=pop.meta.df,f=pop.meta.df$development.stage)
  development.stages=intersect(names(sc.meta.list),names(pop.meta.list))
  development.stages=intersect(std.stages.ordered.str,development.stages)
  stages.col.list=get.col.factor(col.factor =development.stages )
  stages.col.str=stages.col.list$legend.col
  stages.name.str=stages.col.list$legend.str
  len.development.stages=length(development.stages)
  aggregate.sc.rpkm.df=generate.aggregate.per.timepoint(input.data = sc.rpkm.df,meta.data = sc.meta.df)
  intersect.genes.vec=intersect(rownames(pop.rpkm.df),rownames(aggregate.sc.rpkm.df))
  merged.expr.df=cbind(aggregate.sc.rpkm.df[intersect.genes.vec,],pop.rpkm.df[intersect.genes.vec,])
  merged.cor.df=cor(merged.expr.df,method = 'spearman')  
  out.mat=matrix(nrow = len.development.stages,ncol = len.development.stages)
  pvalue.mat=matrix(nrow = len.development.stages,ncol = len.development.stages)
  #dir.create(out.dir,mode = 0777,recursive = T)
  for(n in 1:len.development.stages){
    temp.aggregate=development.stages[n]
    temp.aggregate.sc.rpkm.df=subset(aggregate.sc.rpkm.df,select = temp.aggregate)
    for (m in 1:len.development.stages){
      temp.stage=development.stages[m]
      temp.pop.exp.df=subset(pop.rpkm.df,select = rownames(pop.meta.list[[temp.stage]]))
      temp.pop.exp.df=filter.none.expressed.genes(input.data = temp.pop.exp.df)
      temp.intersect.genes.vec=intersect(rownames(temp.pop.exp.df),rownames(temp.aggregate.sc.rpkm.df))
      temp.pop.exp.df=temp.pop.exp.df[temp.intersect.genes.vec,]
      filt.temp.aggregate.sc.rpkm.df=temp.aggregate.sc.rpkm.df[temp.intersect.genes.vec,]
      temp.merged.expr.df=cbind(temp.pop.exp.df,filt.temp.aggregate.sc.rpkm.df)
      colnames(temp.merged.expr.df)=c(colnames(temp.pop.exp.df),temp.aggregate)
      temp.merged.cor.res=Hmisc::rcorr(x = as.matrix(temp.merged.expr.df),type = 'spearman')
      temp.merged.cor.df=temp.merged.cor.res$r
      temp.merged.pval.df=temp.merged.cor.res$P
      out.mat[n,m]=mean(temp.merged.cor.df[temp.aggregate,colnames(temp.pop.exp.df)]) 
      #out.mat[n,m]=median(temp.merged.cor.df[temp.aggregate,colnames(temp.pop.exp.df)],trim=0.1) 
      name.str=paste(out.dir,paste(c(paste(c('sc_aggregate',temp.aggregate,'bulk',temp.stage,'correlations'),
                                           collapse ='_'),'.tab'),collapse = ''),sep ='/') 
      write.table(x = temp.merged.cor.df,file = name.str,sep='\t',quote = F,col.names = NA)         
      name.str=paste(out.dir,paste(c(paste(c('sc_aggregate',temp.aggregate,'bulk',temp.stage,'pvalues'),
                                           collapse ='_'),'.tab'),collapse = ''),sep ='/') 
      write.table(x = temp.merged.pval.df,file = name.str,sep='\t',quote = F,col.names = NA) 
      pvalue.mat[n,m]=mean(temp.merged.pval.df[temp.aggregate,colnames(temp.pop.exp.df)])
    }
  }
  rownames(out.mat)= development.stages
  colnames(out.mat)= development.stages
  reverse.ordered.development.stages=rev(c("S" ,"ES", "T", "ET" , "LR"  , "R"))
  out.mat=out.mat[,reverse.ordered.development.stages]
  pheatmap::pheatmap(mat =out.mat, cluster_rows = F,cluster_cols = F,
                     main='Average cell vs. bulk',cellwidth = 50,
                     cellheight = 50,fontsize_row = 10,fontsize_col = 10,
                     border_color = 'lightgray')
  pval.fname=paste(out.dir,'aggregate_sc_vs_bulk_mean_pvalues.tab',sep ='/')
  write.table(x = pvalue.mat,file = pval.fname,sep='\t',quote = F,col.names = NA)
  cor.fname=paste(out.dir,'aggregate_sc_vs_bulk_mean_correlation.tab',sep ='/')
  write.table(x = out.mat,file = cor.fname,sep='\t',quote = F,col.names = NA)
}
generate.aggregate.per.timepoint=function(input.data,meta.data){
  sc.meta.list=split(x=meta.data,f=meta.data$development.stage)
  mean.expr=lapply(sc.meta.list,function(item){ 
    subset.samples=intersect(rownames(item),colnames(input.data)) 
    return(as.numeric(apply(input.data[,subset.samples],1,mean)))
  })
  return(as.data.frame(mean.expr,row.names = rownames(input.data)))
}
plot.heatmap = function(in.df,title.str='Test plot',cluster_rows=T,cluster_cols=T,show_rownames = T,
                        show_colnames = T,cellwidth=1,cellheight = 1,fontsize_row = 5,
                        fontsize_col = 5,fontsize = 10,border_color = 'yellow',column.annotation.df=NULL,
                        col.pelette=colorRampPalette(c('navy','gray','firebrick1'))(100)){
  in.matrix=as.matrix(in.df)
  if(is.null(column.annotation.df)){
    pheatmap(mat = in.matrix,main =title.str,color = col.pelette,annotation_col = NULL,
             cellwidth = cellwidth,cellheight = cellheight,
             fontsize_row = fontsize_row,cluster_cols = cluster_cols,
             cluster_rows = cluster_rows,fontsize_col = fontsize_col,
             show_rownames = show_rownames,show_colnames = show_colnames,
             border_color = border_color)
  }
  else{
    pheatmap(mat = in.matrix,main =title.str,color = col.pelette,annotation_col = column.annotation.df,
               cellwidth = cellwidth,cellheight = cellheight,
               fontsize_row = fontsize_row,cluster_cols = cluster_cols,
               cluster_rows = cluster_rows,fontsize_col = fontsize_col,
               show_rownames = show_rownames,show_colnames = show_colnames,
               border_color = border_color) 
    }
}

get.rainbow.color.vec=function(in.vec){
  col.factor=in.vec
  col.categories=sort(unique(col.factor))
  color.names=rainbow(n = length(col.categories))
  col.list=c()
  for (m in  1 : length(col.categories)){
    col.list=append(col.list,color.names[m],after = length(col.list))
  }
  names(col.list)=col.categories
  return(col.list)
}
               
get.col.factor=function(col.factor){
  col.categories=unique(col.factor)
  #all.colors=c('green','blue','red','orange','purple')
  col.list=list()
  color.names=rainbow(length(col.categories))
  #color.names=colorpanel(n=length(col.categories),'green','blue','red')
  #color.names=sample(colours(),length(col.categories))
  for (m in  1 : length(col.categories)){
    col.item=as.character(col.categories[m])
    col.list[[col.item]]=color.names[m] 
  }
  out.col=as.character(unlist(lapply(col.factor,function(col.factor.item){
    col.name=as.character(col.list[[col.factor.item]])
  })))
  legend.text=names(col.list)
  legend.cols=as.character(unlist(col.list))
  out.list=list(col.str=out.col,legend.str=legend.text,legend.cols=legend.cols)
  return(out.list)
}

get.pch.factor=function(pch.factor){
  col.categories=unique(col.factor)
  all.colors=c('green','blue','red','orange','purple')
  col.list=list()
  color.names=rainbow(length(col.categories))
  #color.names=colorpanel(n=length(col.categories),'green','blue','red')
  #color.names=sample(colours(),length(col.categories))
  for (m in  1 : length(col.categories)){
    col.item=as.character(col.categories[m])
    col.list[[col.item]]=color.names[m] 
  }
  out.col=as.character(unlist(lapply(col.factor,function(col.factor.item){
    col.name=as.character(col.list[[col.factor.item]])
  })))
  legend.text=names(col.list)
  legend.cols=as.character(unlist(col.list))
  out.list=list(col.str=out.col,legend.str=legend.text,legend.cols=legend.cols)
  return(out.list)
}
get.color.brewer.list=function(in.vec){
  out.list=list()
  col.factor=in.vec
  col.categories=unique(col.factor)
  color.names=brewer.pal(n = length(col.categories),name = 'Set3')
  #color.names=brewer.pal(n = length(col.categories),name = 'Pastel1')
  col.list=list()
  for (m in  1 : length(col.categories)){
    col.item=as.character(col.categories[m])
    col.list[[col.item]]=color.names[m] 
  }
  out.col=as.character(unlist(lapply(col.factor,function(col.factor.item){
    col.name=as.character(col.list[[col.factor.item]])
  })))
  legend.text=names(col.list)
  legend.cols=as.character(unlist(col.list))
  out.list=list(col.str=out.col,legend.str=legend.text,legend.cols=legend.cols)
  return(out.list)
}
get.color.list.for.pheatmap=function(in.vec){
  out.list=list()
  col.factor=in.vec
  col.categories=sort(unique(col.factor))
  color.names=brewer.pal(n = length(col.categories),name = 'Set3')
  #color.names=rainbow(n =length(col.categories) )
  col.list=c()
  for (m in  1 : length(col.categories)){
    col.list=append(col.list,color.names[m],after = length(col.list))
  }
  names(col.list)=col.categories
  return(col.list)
}
get.color.rainbow.list.for.pheatmap=function(in.vec){
  out.list=list()
  col.factor=in.vec
  col.categories=sort(unique(col.factor))
  color.names=rainbow(n =length(col.categories) )
  col.list=c()
  for (m in  1 : length(col.categories)){
    col.list=append(col.list,color.names[m],after = length(col.list))
  }
  names(col.list)=col.categories
  return(col.list)
}
get.distinct.colour.list.for.pheatmap=function(in.vec){
  out.list=list()
  col.factor=in.vec
  col.categories=sort(unique(col.factor))
  color.names=rainbow(n =length(col.categories) )
  col.list=randomcoloR::randomColor(count = length(col.categories))
  names(col.list)=col.categories
  return(col.list)
}
get.color.pastel.list.for.pheatmap=function(in.vec){
  out.list=list()
  col.factor=in.vec
  col.categories=sort(unique(col.factor))
  #brewer.pal(n = length(col.categories),name = 'Pastel')
  color.names=brewer.pal(n = length(col.categories), name = "Pastel1")
  #color.names=rainbow(n =length(col.categories) )
  col.list=c()
  for (m in  1 : length(col.categories)){
    col.list=append(col.list,color.names[m],after = length(col.list))
  }
  names(col.list)=col.categories
  return(col.list)
}
                             
#Installs a package given the name
install.bioconductor.package=function(package.name){
  source("https://bioconductor.org/biocLite.R")
  biocLite(package.name)
}
test.function=function(in_dir){
  all.files=list.files(in_dir,full.names = T)
  samples.vec=c()
  batch.vec=c()
  for (f in all.files){
    if(grepl(pattern = '.bam$',x = f)){
      f.parts=unlist(strsplit(x = f,split='/'))
      t.sample=sub(pattern = '_unique.bam',replacement = '',x = f.parts[length(f.parts)])
      samples.vec=append(samples.vec,values = t.sample,after = length(samples.vec))
      batch.vec=append(x = batch.vec,values = f.parts[length(f.parts)-1],after =length(batch.vec) )
    }
  }
  out.df=data.frame(sample=samples.vec,batch=batch.vec)
  return(out.df)
}
