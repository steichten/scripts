options(echo=T)
library(fields)
args=commandArgs(trailingOnly=T)
print(args)

#prerecs
library(edmr)
library(methylKit)
library(GenomicRanges)
library(mixtools)
library(data.table)
library(reshape2)


#EDIT AS NEEDED###########################
context=args[1]
min.cov=args[2]
assembly='Bd21Control'
##########################################
#list of all sam files created by bismark (sorted!)
a=dir(pattern="*.sam")

biglist=as.list(a)

#make one gigantic list of all the sam read-ins
for(i in 1:length(a)){
biglist[[i]]=read.bismark(a[i],sample.id=strsplit(a[i],'_')[[1]][1],assembly=assembly,save.context=context,read.context=context,treatment=c(0),mincov=min.cov)
}

aa=combn(a,2) #combinations of filenames
bb=combn(1:length(biglist),2)

#make methylRawList sets
for(i in 1:(length(aa)/2)){
#
myobj=new("methylRawList",list(biglist[[bb[1,i]]],biglist[[bb[2,i]]]),treatment=c(0,1))

#only keep cytosines in which all samples meet the coverage requirement
meth=unite(myobj)

contrast=paste(strsplit(aa[1,i],'_')[[1]][1],strsplit(aa[2,i],'_')[[1]][1],sep='v')
# cluster all samples using correlation distance and return a tree object for plclust
#hc = clusterSamples(meth, dist="correlation", method="ward", plot=T)

# principal component anlaysis of all samples.
#PCASamples(meth)

# calculate differential methylation p-values and q-values
myDiff=calculateDiffMeth(meth)

write.table(getData(myDiff),paste("1_MethylKit_results/",contrast,"_mincov_",min.cov,"_",context,"_",Sys.Date(),".bed",sep=''),sep='\t',row.names=F,quote=F,col.names=F)
# get differentially methylated regions with 25% difference and qvalue<0.01
#myDiff25p=get.methylDiff(myDiff,difference=25,qvalue=0.01)


#DMRs with eDMR

# fitting the bimodal normal distribution to CpGs distribution
myMixmdl=myDiff.to.mixmdl(myDiff,plot=F,main=contrast)

# plot cost function and the determined distance cutoff
#plotCost(myMixmdl, main="cost function")

# calculate all DMRs candidate
mydmr=edmr(myDiff, mode=1, ACF=TRUE)

dmr.output <- data.frame(seqnames=seqnames(mydmr),
  starts=start(mydmr)-1, #coordinate 0-based like bed file
  ends=end(mydmr),
  names=c(rep(contrast,length(mydmr))),
  scores=c(rep(".", length(mydmr))),
  strands=strand(mydmr),
  mean.meth.diff=elementMetadata(mydmr)$mean.meth.diff,
  num.c=elementMetadata(mydmr)$num.CpGs,
  num.dmcs=elementMetadata(mydmr)$num.DMCs,
  dmr.pvalue=elementMetadata(mydmr)$DMR.pvalue,
  dmr.qvalue=elementMetadata(mydmr)$DMR.qvalue)
  
write.table(dmr.output,paste("2_eDMR/",contrast,"_mincov_",min.cov,'_',context,'_',Sys.Date(),".bed",sep=''),sep='\t',row.names=F,quote=F,col.names=F)



# further filtering the DMRs
mysigdmr=filter.dmr(mydmr)

dmr.output <- data.frame(seqnames=seqnames(mysigdmr),
  starts=start(mysigdmr)-1,
  ends=end(mysigdmr),
  names=c(rep(contrast,length(mysigdmr))),
  scores=c(rep(".", length(mysigdmr))),
  strands=strand(mysigdmr),
  mean.meth.diff=elementMetadata(mysigdmr)$mean.meth.diff,
  num.c=elementMetadata(mysigdmr)$num.CpGs,
  num.dmcs=elementMetadata(mysigdmr)$num.DMCs,
  dmr.pvalue=elementMetadata(mysigdmr)$DMR.pvalue,
  dmr.qvalue=elementMetadata(mysigdmr)$DMR.qvalue)
  
write.table(dmr.output,paste("3_eDMR_filtered/",contrast,"_mincov_",min.cov,'_',context,"_filtered_",Sys.Date(),".bed",sep=''),sep='\t',row.names=F,quote=F,col.names=F)

}

