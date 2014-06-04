################################
library(fields)

#read in files
gene.cov=read.delim('allmiseq_to_genes.1k.bed',head=F)
gene.bga.cov=read.delim('allmiseq_to_genes_bga.1k.bed',head=F)
gaut.cov=read.delim('gaut1.bd21leaf_bismark_gene.1k.bed',head=F)
gaut.bga.cov=read.delim('gaut1.bd21leaf_bismark_gene_bga.1k.bed',head=F)

#remove scaffolds
gene.cov.chrom=subset(gene.cov,gene.cov$V1=='Bd1' | gene.cov$V1=='Bd2' | gene.cov$V1=='Bd3' | gene.cov$V1=='Bd4' | gene.cov$V1=='Bd5')
gene.bga.cov.chrom=subset(gene.bga.cov,gene.bga.cov$V1=='Bd1' | gene.bga.cov$V1=='Bd2' | gene.bga.cov$V1=='Bd3' | gene.bga.cov$V1=='Bd4' | gene.bga.cov$V1=='Bd5')
gaut.cov.chrom=subset(gaut.cov,gaut.cov$V1=='Bd1' | gaut.cov$V1=='Bd2' | gaut.cov$V1=='Bd3' | gaut.cov$V1=='Bd4' | gaut.cov$V1=='Bd5')
gaut.bga.cov.chrom=subset(gaut.bga.cov,gaut.bga.cov$V1=='Bd1' | gaut.bga.cov$V1=='Bd2' | gaut.bga.cov$V1=='Bd3' | gaut.bga.cov$V1=='Bd4' | gaut.bga.cov$V1=='Bd5')

#get real (strand corrected) distances and relitive distances
real.dist=matrix(ifelse(gene.cov.chrom$V8=='+',-1*gene.cov.chrom$V10,gene.cov.chrom$V10),ncol=1)
gene.cov.chrom=cbind(gene.cov.chrom,real.dist)

real.dist=matrix(ifelse(gene.bga.cov.chrom$V8=='+',-1*gene.bga.cov.chrom$V10,gene.bga.cov.chrom$V10),ncol=1)
gene.bga.cov.chrom=cbind(gene.bga.cov.chrom,real.dist)

real.dist=matrix(ifelse(gaut.cov.chrom$V8=='+',-1*gaut.cov.chrom$V10,gaut.cov.chrom$V10),ncol=1)
gaut.cov.chrom=cbind(gaut.cov.chrom,real.dist)

real.dist=matrix(ifelse(gaut.bga.cov.chrom$V8=='+',-1*gaut.bga.cov.chrom$V10,gaut.bga.cov.chrom$V10),ncol=1)
gaut.bga.cov.chrom=cbind(gaut.bga.cov.chrom,real.dist)

#rel distance
rel.dist=matrix(ifelse(gene.cov.chrom$real.dist==0,ifelse(gene.cov.chrom$V8=="-",((gene.cov.chrom$V7 - (gene.cov.chrom$V2))/(gene.cov.chrom$V7 - gene.cov.chrom$V6))*1000,(((gene.cov.chrom$V2) - gene.cov.chrom$V6)/(gene.cov.chrom$V7 - gene.cov.chrom$V6))*1000),ifelse(gene.cov.chrom$real.dist>0,gene.cov.chrom$real.dist + 1000,gene.cov.chrom$real.dist)),ncol=1)
gene.cov.chrom=cbind(gene.cov.chrom,rel.dist)

rel.dist=matrix(ifelse(gene.bga.cov.chrom$real.dist==0,ifelse(gene.bga.cov.chrom$V8=="-",((gene.bga.cov.chrom$V7 - (gene.bga.cov.chrom$V2))/(gene.bga.cov.chrom$V7 - gene.bga.cov.chrom$V6))*1000,(((gene.bga.cov.chrom$V2) - gene.bga.cov.chrom$V6)/(gene.bga.cov.chrom$V7 - gene.bga.cov.chrom$V6))*1000),ifelse(gene.bga.cov.chrom$real.dist>0,gene.bga.cov.chrom$real.dist + 1000,gene.bga.cov.chrom$real.dist)),ncol=1)
gene.bga.cov.chrom=cbind(gene.bga.cov.chrom,rel.dist)

rel.dist=matrix(ifelse(gaut.cov.chrom$real.dist==0,ifelse(gaut.cov.chrom$V8=="-",((gaut.cov.chrom$V7 - (gaut.cov.chrom$V2))/(gaut.cov.chrom$V7 - gaut.cov.chrom$V6))*1000,(((gaut.cov.chrom$V2) - gaut.cov.chrom$V6)/(gaut.cov.chrom$V7 - gaut.cov.chrom$V6))*1000),ifelse(gaut.cov.chrom$real.dist>0,gaut.cov.chrom$real.dist + 1000,gaut.cov.chrom$real.dist)),ncol=1)
gaut.cov.chrom=cbind(gaut.cov.chrom,rel.dist)

rel.dist=matrix(ifelse(gaut.bga.cov.chrom$real.dist==0,ifelse(gaut.bga.cov.chrom$V8=="-",((gaut.bga.cov.chrom$V7 - (gaut.bga.cov.chrom$V2))/(gaut.bga.cov.chrom$V7 - gaut.bga.cov.chrom$V6))*1000,(((gaut.bga.cov.chrom$V2) - gaut.bga.cov.chrom$V6)/(gaut.bga.cov.chrom$V7 - gaut.bga.cov.chrom$V6))*1000),ifelse(gaut.bga.cov.chrom$real.dist>0,gaut.bga.cov.chrom$real.dist + 1000,gaut.bga.cov.chrom$real.dist)),ncol=1)
gaut.bga.cov.chrom=cbind(gaut.bga.cov.chrom,rel.dist)

#quickfix for bullshit regions that cause incorrect relative positioning
#fixy=ifelse(gene.cov.chrom$rel.dist < 0 & gene.cov.chrom$real.dist==0,0,ifelse(gene.cov.chrom$rel.dist >1000 & gene.cov.chrom$real.dist==0,1000,gene.cov.chrom$rel.dist))
#gene.cov.chrom$rel.dist=fixy

#fixy=ifelse(gene.bga.cov.chrom$rel.dist < 0 & gene.bga.cov.chrom$real.dist==0,0,ifelse(gene.bga.cov.chrom$rel.dist >1000 & gene.bga.cov.chrom$real.dist==0,1000,gene.bga.cov.chrom$rel.dist))
#gene.bga.cov.chrom$rel.dist=fixy

#fixy=ifelse(gaut.cov.chrom$rel.dist < 0 & gaut.cov.chrom$real.dist==0,0,ifelse(gaut.cov.chrom$rel.dist >1000 & gaut.cov.chrom$real.dist==0,1000,gaut.cov.chrom$rel.dist))
#gaut.cov.chrom$rel.dist=fixy

#fixy=ifelse(gaut.bga.cov.chrom$rel.dist < 0 & gaut.bga.cov.chrom$real.dist==0,0,ifelse(gaut.bga.cov.chrom$rel.dist >1000 & gaut.bga.cov.chrom$real.dist==0,1000,gaut.bga.cov.chrom$rel.dist))
#gaut.bga.cov.chrom$rel.dist=fixy

#quickfixv2 - just remove the coverage tracks that cause odd positioning with overlap
fixy=ifelse(gene.cov.chrom$rel.dist < 0 & gene.cov.chrom$real.dist==0,0,ifelse(gene.cov.chrom$rel.dist >1000 & gene.cov.chrom$real.dist==0,0,1))
gene.cov.chrom=cbind(gene.cov.chrom,fixy)
gene.cov.chrom=subset(gene.cov.chrom,gene.cov.chrom$fixy==1)

fixy=ifelse(gene.bga.cov.chrom$rel.dist < 0 & gene.bga.cov.chrom$real.dist==0,0,ifelse(gene.bga.cov.chrom$rel.dist >1000 & gene.bga.cov.chrom$real.dist==0,0,1))
gene.bga.cov.chrom=cbind(gene.bga.cov.chrom,fixy)
gene.bga.cov.chrom=subset(gene.bga.cov.chrom,gene.bga.cov.chrom$fixy==1)

fixy=ifelse(gaut.cov.chrom$rel.dist < 0 & gaut.cov.chrom$real.dist==0,0,ifelse(gaut.cov.chrom$rel.dist >1000 & gaut.cov.chrom$real.dist==0,0,1))
gaut.cov.chrom=cbind(gaut.cov.chrom,fixy)
gaut.cov.chrom=subset(gaut.cov.chrom,gaut.cov.chrom$fixy==1)

fixy=ifelse(gaut.bga.cov.chrom$rel.dist < 0 & gaut.bga.cov.chrom$real.dist==0,0,ifelse(gaut.bga.cov.chrom$rel.dist >1000 & gaut.bga.cov.chrom$real.dist==0,0,1))
gaut.bga.cov.chrom=cbind(gaut.bga.cov.chrom,fixy)
gaut.bga.cov.chrom=subset(gaut.bga.cov.chrom,gaut.bga.cov.chrom$fixy==1)

#perform the stat binning
cov.gene.test=stats.bin(gene.cov.chrom$rel.dist,gene.cov.chrom$V4,N=100)
p.cov.gene.test=cbind(matrix(cov.gene.test$centers,ncol=1),cov.gene.test$stats["mean",])

cov.gene.gene.bga=stats.bin(gene.bga.cov.chrom$rel.dist,gene.bga.cov.chrom$V4,N=100)
p.cov.gene.gene.bga=cbind(matrix(cov.gene.gene.bga$centers,ncol=1),cov.gene.gene.bga$stats["mean",])

cov.gaut.test=stats.bin(gaut.cov.chrom$rel.dist,gaut.cov.chrom$V4,N=100)
p.cov.gaut.test=cbind(matrix(cov.gaut.test$centers,ncol=1),cov.gaut.test$stats["mean",])

cov.gaut.gaut.bga=stats.bin(gaut.bga.cov.chrom$rel.dist,gaut.bga.cov.chrom$V4,N=100)
p.cov.gaut.gaut.bga=cbind(matrix(cov.gaut.gaut.bga$centers,ncol=1),cov.gaut.gaut.bga$stats["mean",])

pdf('coverage_comparison_fix2_13-5-14.pdf',h=8,w=8)
plot(x=NULL,y=NULL,xlim=c(-1000,2000),ylim=c(0,30),xlab='',ylab='readcount',main='Alignment coverage over genes')
lines(p.cov.gene.test,col=1,lwd=2)
lines(p.cov.gene.gene.bga,col=2,lwd=2)
lines(p.cov.gaut.test,col=3,lwd=2)
lines(p.cov.gaut.gaut.bga,col=4,lwd=2)
abline(v=0,lty=2)
abline(v=1000,lty=2)
legend('topright',c('allmiseq','allmiseq_bga','gaut','gaut_bga'),col=c(1,2,3,4),lwd=2,lty=1)
#####################################################################################


#te avg coverage plotting

repeats.cov=read.delim('allmiseq_to_repeats.1k.bed',head=F)
repeats_bga.cov=read.delim('allmiseq_to_repeats_bga.1k.bed',head=F)
gaut.repeats.cov=read.delim('gaut1.bd21leaf_bismark_repeat.1k.bed',head=F)
gaut.repeats_bga.cov=read.delim('gaut1.bd21leaf_bismark_repeat_bga.1k.bed',head=F)

#remove scaffolds
repeats.cov.chrom=subset(repeats.cov,repeats.cov$V1=='Bd1' | repeats.cov$V1=='Bd2' | repeats.cov$V1=='Bd3' | repeats.cov$V1=='Bd4' | repeats.cov$V1=='Bd5')
repeats_bga.cov.chrom=subset(repeats_bga.cov,repeats_bga.cov$V1=='Bd1' | repeats_bga.cov$V1=='Bd2' | repeats_bga.cov$V1=='Bd3' | repeats_bga.cov$V1=='Bd4' | repeats_bga.cov$V1=='Bd5')
gaut.repeats.cov.chrom=subset(gaut.repeats.cov,gaut.repeats.cov$V1=='Bd1' | gaut.repeats.cov$V1=='Bd2' | gaut.repeats.cov$V1=='Bd3' | gaut.repeats.cov$V1=='Bd4' | gaut.repeats.cov$V1=='Bd5')
gaut.repeats_bga.cov.chrom=subset(gaut.repeats_bga.cov,gaut.repeats_bga.cov$V1=='Bd1' | gaut.repeats_bga.cov$V1=='Bd2' | gaut.repeats_bga.cov$V1=='Bd3' | gaut.repeats_bga.cov$V1=='Bd4' | gaut.repeats_bga.cov$V1=='Bd5')

#grab +/- 1kb
repeats.cov.chrom.sub=subset(repeats.cov.chrom,repeats.cov.chrom[,ncol(repeats.cov.chrom)]<=1000 & repeats.cov.chrom[,ncol(repeats.cov.chrom)]>= -1000)
repeats_bga.cov.chrom.sub=subset(repeats_bga.cov.chrom,repeats_bga.cov.chrom[,ncol(repeats_bga.cov.chrom)]<=1000 & repeats_bga.cov.chrom[,ncol(repeats_bga.cov.chrom)]>= -1000)
gaut.repeats.cov.chrom.sub=subset(gaut.repeats.cov.chrom,gaut.repeats.cov.chrom[,ncol(gaut.repeats.cov.chrom)]<=1000 & gaut.repeats.cov.chrom[,ncol(gaut.repeats.cov.chrom)]>= -1000)
gaut.repeats_bga.cov.chrom.sub=subset(gaut.repeats_bga.cov.chrom,gaut.repeats_bga.cov.chrom[,ncol(gaut.repeats_bga.cov.chrom)]<=1000 & gaut.repeats_bga.cov.chrom[,ncol(gaut.repeats_bga.cov.chrom)]>= -1000)
#get real (strand corrected) distances and relitive distances
real.dist=matrix(ifelse(repeats.cov.chrom.sub$V8=='+',-1*repeats.cov.chrom.sub$V10,repeats.cov.chrom.sub$V10),ncol=1)
repeats.cov.chrom.sub=cbind(repeats.cov.chrom.sub,real.dist)

real.dist=matrix(ifelse(repeats_bga.cov.chrom.sub$V8=='+',-1*repeats_bga.cov.chrom.sub$V10,repeats_bga.cov.chrom.sub$V10),ncol=1)
repeats_bga.cov.chrom.sub=cbind(repeats_bga.cov.chrom.sub,real.dist)

real.dist=matrix(ifelse(gaut.repeats.cov.chrom.sub$V8=='+',-1*gaut.repeats.cov.chrom.sub$V10,gaut.repeats.cov.chrom.sub$V10),ncol=1)
gaut.repeats.cov.chrom.sub=cbind(gaut.repeats.cov.chrom.sub,real.dist)

real.dist=matrix(ifelse(gaut.repeats_bga.cov.chrom.sub$V8=='+',-1*gaut.repeats_bga.cov.chrom.sub$V10,gaut.repeats_bga.cov.chrom.sub$V10),ncol=1)
gaut.repeats_bga.cov.chrom.sub=cbind(gaut.repeats_bga.cov.chrom.sub,real.dist)

#rel distance
rel.dist=matrix(ifelse(repeats.cov.chrom.sub$real.dist==0,ifelse(repeats.cov.chrom.sub$V8=="-",((repeats.cov.chrom.sub$V7 - (repeats.cov.chrom.sub$V2))/(repeats.cov.chrom.sub$V7 - repeats.cov.chrom.sub$V6))*1000,(((repeats.cov.chrom.sub$V2) - repeats.cov.chrom.sub$V6)/(repeats.cov.chrom.sub$V7 - repeats.cov.chrom.sub$V6))*1000),ifelse(repeats.cov.chrom.sub$real.dist>0,repeats.cov.chrom.sub$real.dist + 1000,repeats.cov.chrom.sub$real.dist)),ncol=1)
repeats.cov.chrom.sub=cbind(repeats.cov.chrom.sub,rel.dist)

rel.dist=matrix(ifelse(repeats_bga.cov.chrom.sub$real.dist==0,ifelse(repeats_bga.cov.chrom.sub$V8=="-",((repeats_bga.cov.chrom.sub$V7 - (repeats_bga.cov.chrom.sub$V2))/(repeats_bga.cov.chrom.sub$V7 - repeats_bga.cov.chrom.sub$V6))*1000,(((repeats_bga.cov.chrom.sub$V2) - repeats_bga.cov.chrom.sub$V6)/(repeats_bga.cov.chrom.sub$V7 - repeats_bga.cov.chrom.sub$V6))*1000),ifelse(repeats_bga.cov.chrom.sub$real.dist>0,repeats_bga.cov.chrom.sub$real.dist + 1000,repeats_bga.cov.chrom.sub$real.dist)),ncol=1)
repeats_bga.cov.chrom.sub=cbind(repeats_bga.cov.chrom.sub,rel.dist)

rel.dist=matrix(ifelse(gaut.repeats.cov.chrom.sub$real.dist==0,ifelse(gaut.repeats.cov.chrom.sub$V8=="-",((gaut.repeats.cov.chrom.sub$V7 - (gaut.repeats.cov.chrom.sub$V2))/(gaut.repeats.cov.chrom.sub$V7 - gaut.repeats.cov.chrom.sub$V6))*1000,(((gaut.repeats.cov.chrom.sub$V2) - gaut.repeats.cov.chrom.sub$V6)/(gaut.repeats.cov.chrom.sub$V7 - gaut.repeats.cov.chrom.sub$V6))*1000),ifelse(gaut.repeats.cov.chrom.sub$real.dist>0,gaut.repeats.cov.chrom.sub$real.dist + 1000,gaut.repeats.cov.chrom.sub$real.dist)),ncol=1)
gaut.repeats.cov.chrom.sub=cbind(gaut.repeats.cov.chrom.sub,rel.dist)

rel.dist=matrix(ifelse(gaut.repeats_bga.cov.chrom.sub$real.dist==0,ifelse(gaut.repeats_bga.cov.chrom.sub$V8=="-",((gaut.repeats_bga.cov.chrom.sub$V7 - (gaut.repeats_bga.cov.chrom.sub$V2))/(gaut.repeats_bga.cov.chrom.sub$V7 - gaut.repeats_bga.cov.chrom.sub$V6))*1000,(((gaut.repeats_bga.cov.chrom.sub$V2) - gaut.repeats_bga.cov.chrom.sub$V6)/(gaut.repeats_bga.cov.chrom.sub$V7 - gaut.repeats_bga.cov.chrom.sub$V6))*1000),ifelse(gaut.repeats_bga.cov.chrom.sub$real.dist>0,gaut.repeats_bga.cov.chrom.sub$real.dist + 1000,gaut.repeats_bga.cov.chrom.sub$real.dist)),ncol=1)
gaut.repeats_bga.cov.chrom.sub=cbind(gaut.repeats_bga.cov.chrom.sub,rel.dist)
#quickfix for bullshit
fixy=ifelse(repeats.cov.chrom.sub$rel.dist < 0 & repeats.cov.chrom.sub$real.dist==0,0,ifelse(repeats.cov.chrom.sub$rel.dist >1000 & repeats.cov.chrom.sub$real.dist==0,1000,repeats.cov.chrom.sub$rel.dist))
repeats.cov.chrom.sub$rel.dist=fixy

fixy=ifelse(repeats_bga.cov.chrom.sub$rel.dist < 0 & repeats_bga.cov.chrom.sub$real.dist==0,0,ifelse(repeats_bga.cov.chrom.sub$rel.dist >1000 & repeats_bga.cov.chrom.sub$real.dist==0,1000,repeats_bga.cov.chrom.sub$rel.dist))
repeats_bga.cov.chrom.sub$rel.dist=fixy

fixy=ifelse(gaut.repeats.cov.chrom.sub$rel.dist < 0 & gaut.repeats.cov.chrom.sub$real.dist==0,0,ifelse(gaut.repeats.cov.chrom.sub$rel.dist >1000 & gaut.repeats.cov.chrom.sub$real.dist==0,1000,gaut.repeats.cov.chrom.sub$rel.dist))
gaut.repeats.cov.chrom.sub$rel.dist=fixy

fixy=ifelse(gaut.repeats_bga.cov.chrom.sub$rel.dist < 0 & gaut.repeats_bga.cov.chrom.sub$real.dist==0,0,ifelse(gaut.repeats_bga.cov.chrom.sub$rel.dist >1000 & gaut.repeats_bga.cov.chrom.sub$real.dist==0,1000,gaut.repeats_bga.cov.chrom.sub$rel.dist))
gaut.repeats_bga.cov.chrom.sub$rel.dist=fixy

cov.repeats.test=stats.bin(repeats.cov.chrom.sub$rel.dist,repeats.cov.chrom.sub$V4,N=100)
p.cov.repeats.test=cbind(matrix(cov.repeats.test$centers,ncol=1),cov.repeats.test$stats["mean",])

cov.repeats.repeats_bga=stats.bin(repeats_bga.cov.chrom.sub$rel.dist,repeats_bga.cov.chrom.sub$V4,N=100)
p.cov.repeats.repeats_bga=cbind(matrix(cov.repeats.repeats_bga$centers,ncol=1),cov.repeats.repeats_bga$stats["mean",])

cov.gaut.repeats.test=stats.bin(gaut.repeats.cov.chrom.sub$rel.dist,gaut.repeats.cov.chrom.sub$V4,N=100)
p.cov.gaut.repeats.test=cbind(matrix(cov.gaut.repeats.test$centers,ncol=1),cov.gaut.repeats.test$stats["mean",])

cov.gaut.repeats.gaut.repeats_bga=stats.bin(gaut.repeats_bga.cov.chrom.sub$rel.dist,gaut.repeats_bga.cov.chrom.sub$V4,N=100)
p.cov.gaut.repeats.gaut.repeats_bga=cbind(matrix(cov.gaut.repeats.gaut.repeats_bga$centers,ncol=1),cov.gaut.repeats.gaut.repeats_bga$stats["mean",])

plot(x=NULL,y=NULL,xlim=c(-1000,2000),ylim=c(0,30),xlab='',ylab='readcount',main='Alignment coverage over repeats')
lines(p.cov.repeats.test,col=1,lwd=2)
lines(p.cov.repeats.repeats_bga,col=2,lwd=2)
lines(p.cov.gaut.repeats.test,col=3,lwd=2)
lines(p.cov.gaut.repeats.gaut.repeats_bga,col=4,lwd=2)
abline(v=0,lty=2)
abline(v=1000,lty=2)
legend('topright',c('allmiseq','allmiseq_bga','gaut','gaut_bga'),col=c(1,2,3,4),lwd=2,lty=1)
####################################################################
dev.off()
#