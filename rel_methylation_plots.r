options(echo=T)
library(fields)
args=commandArgs(trailingOnly=T)
print(args)

#read in files#
cpg=read.delim(paste(args[1],'_CpG_',args[2],'.1k.bed',sep=''),head=F)
chg=read.delim(paste(args[1],'_CHG_',args[2],'.1k.bed',sep=''),head=F)
chh=read.delim(paste(args[1],'_CHH_',args[2],'.1k.bed',sep=''),head=F)

#remove scaffolds
cpg.sub=subset(cpg,cpg$V1!='Mt' & cpg$V1!='chrMt' & cpg$V1!='Pt' & cpg$V1!='chrPt')
chg.sub=subset(chg,chg$V1!='Mt' & chg$V1!='chrMt' & chg$V1!='Pt' & chg$V1!='chrPt')
chh.sub=subset(chh,chh$V1!='Mt' & chh$V1!='chrMt' & chh$V1!='Pt' & chh$V1!='chrPt')

f1.end=as.numeric(args[3])
f2.end=as.numeric(args[4])+1
cpg.sub=subset(cpg.sub,cpg.sub[,f1.end + f2.end]!= -1)
chg.sub=subset(chg.sub,chg.sub[,f1.end + f2.end]!= -1)
chh.sub=subset(chh.sub,chh.sub[,f1.end + f2.end]!= -1)

#CpG
real.dist=matrix(ifelse(cpg.sub[,f1.end + 6]=='+',-1*cpg.sub[,f1.end+f2.end],cpg.sub[,f1.end + f2.end]),ncol=1)
cpg.sub=cbind(cpg.sub,real.dist)
rel.dist=matrix(ifelse(cpg.sub$real.dist==0,ifelse(cpg.sub[,f1.end + 6]=="-",((cpg.sub[,f1.end + 3] - (cpg.sub[,2]))/(cpg.sub[,f1.end + 3] - cpg.sub[,f1.end + 2]))*1000,(((cpg.sub[,2]) - cpg.sub[,f1.end + 2])/(cpg.sub[,f1.end + 3] - cpg.sub[, f1.end + 2]))*1000),ifelse(cpg.sub$real.dist>0,cpg.sub$real.dist + 1000,cpg.sub$real.dist)),ncol=1)
cpg.sub=cbind(cpg.sub,rel.dist)
fixy=ifelse(cpg.sub$rel.dist < 0 & cpg.sub$real.dist==0,0,ifelse(cpg.sub$rel.dist >1000 & cpg.sub$real.dist==0,1000,cpg.sub$rel.dist))
cpg.sub$rel.dist=fixy
cpg.bin=stats.bin(cpg.sub$rel.dist,cpg.sub$V4,N=200)
p.cpg.bin=cbind(matrix(cpg.bin$centers,ncol=1),cpg.bin$stats["mean",])

#CHG
real.dist=matrix(ifelse(chg.sub[,f1.end + 6]=='+',-1*chg.sub[,f1.end + f2.end],chg.sub[,f1.end + f2.end]),ncol=1)
chg.sub=cbind(chg.sub,real.dist)
rel.dist=matrix(ifelse(chg.sub$real.dist==0,ifelse(chg.sub[,f1.end + 6]=="-",((chg.sub[,f1.end + 3] - (chg.sub$V2))/(chg.sub[,f1.end + 3] - chg.sub[,f1.end + 2]))*1000,(((chg.sub$V2) - chg.sub[,f1.end + 2])/(chg.sub[,f1.end + 3] - chg.sub[,f1.end + 2]))*1000),ifelse(chg.sub$real.dist>0,chg.sub$real.dist + 1000,chg.sub$real.dist)),ncol=1)
chg.sub=cbind(chg.sub,rel.dist)
fixy=ifelse(chg.sub$rel.dist < 0 & chg.sub$real.dist==0,0,ifelse(chg.sub$rel.dist >1000 & chg.sub$real.dist==0,1000,chg.sub$rel.dist))
chg.sub$rel.dist=fixy
chg.bin=stats.bin(chg.sub$rel.dist,chg.sub$V4,N=200)
p.chg.bin=cbind(matrix(chg.bin$centers,ncol=1),chg.bin$stats["mean",])

#CHH
real.dist=matrix(ifelse(chh.sub[,f1.end + 6]=='+',-1*chh.sub[,f1.end + f2.end],chh.sub[,f1.end + f2.end]),ncol=1)
chh.sub=cbind(chh.sub,real.dist)
rel.dist=matrix(ifelse(chh.sub$real.dist==0,ifelse(chh.sub[,f1.end + 6]=="-",((chh.sub[,f1.end + 3] - (chh.sub$V2))/(chh.sub[,f1.end + 3] - chh.sub[,f1.end + 2]))*1000,(((chh.sub$V2) - chh.sub[,f1.end + 2])/(chh.sub[,f1.end + 3] - chh.sub[,f1.end + 2]))*1000),ifelse(chh.sub$real.dist>0,chh.sub$real.dist + 1000,chh.sub$real.dist)),ncol=1)
chh.sub=cbind(chh.sub,rel.dist)
fixy=ifelse(chh.sub$rel.dist < 0 & chh.sub$real.dist==0,0,ifelse(chh.sub$rel.dist >1000 & chh.sub$real.dist==0,1000,chh.sub$rel.dist))
chh.sub$rel.dist=fixy
chh.bin=stats.bin(chh.sub$rel.dist,chh.sub$V4,N=200)
p.chh.bin=cbind(matrix(chh.bin$centers,ncol=1),chh.bin$stats["mean",])

#create plots
pdf(paste(args[1],'_',args[2],'_methylation.pdf',sep=''),h=10,w=12)
plot(x=NULL,y=NULL,xlim=c(-1000,2000),ylim=c(0,100),xlab='',ylab='% methylation',main=paste(args[1],' methylation over ',args[2],sep=''))
lines(p.cpg.bin,col=1,lwd=2)
lines(p.chg.bin,col=2,lwd=2)
lines(p.chh.bin,col=3,lwd=2)
abline(v=0,lty=2)
abline(v=1000,lty=2)
legend('topright',c(paste(args[1],' - CpG',sep=''),paste(args[1],' - CHG',sep=''),paste(args[1],' - CHH',sep='')),col=c(1,2,3),lwd=2,lty=1)
dev.off()
#####################################################################################

out=cbind(p.cpg.bin,p.chg.bin[,2],p.chh.bin[,2])
colnames(out)=c('pos',paste(args[1],'-CpG',sep=''),paste(args[1],'-CHG',sep=''),paste(args[1],'-CHH',sep=''))
write.table(out,paste(args[1],'_',args[2],'values.txt',sep=''),sep='\t')