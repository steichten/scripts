


library(fields)

data=read.delim('window100_B73130M_to_MTEC.output',head=F)

real.dist=matrix(ifelse(data$V13=='+' ,-1* data$V22,data$V22),ncol=1)

data=cbind(data,real.dist)

rel.dist=matrix(ifelse(data$V22==0,ifelse(data$V13=="-",((data$V12 - (data$V2+50))/(data$V19))*1000,(((data$V2+50) - data$V11)/data$V19)*1000),ifelse(data$real.dist>0,data$real.dist + 1000,data$real.dist)),ncol=1)


data=cbind(data,rel.dist)

###



data.sub=subset(data,data$rel.dist< 2000 & data$rel.dist > -1000)
class.list=names(table(data$V20))

cg.output=matrix(ncol=length(class.list),nrow=60)

for(i in 1:length(class.list)){
data.sub.1=subset(data.sub,data.sub$V20==class.list[i])

stats.test=stats.bin(data.sub.1$rel.dist,data.sub.1$V8,N=60)
cg.output[,i]=stats.test$stats["mean",]*100

}

centers=matrix(stats.test$centers,ncol=1)
cg.output=cbind(centers,cg.output)
colnames(cg.output)=c('window_center',class.list)

pdf('testA.pdf',width=8,height=11)
par(mfrow=c(3,1))
#CG
#plot
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,100),xlab="",ylab='% methylation',main='B73 CG methylation across MTEC classifications')

for(i in 2:ncol(cg.output)){
lines(cbind(cg.output[,1],cg.output[,i]),col=i,lwd=2)
}
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')

legend("topright",c('LINE_RIL','LINE_RIX','LTR_RLC','LTR_RLG','LTR_RLX','TIR_DTA','TIR_DTC','TIR_DTH','TIR_DTM','TIR_DTT'),col=c(2:11),lty=1,lwd=2)
#


#CHG
#plot
chg.output=matrix(ncol=length(class.list),nrow=60)

for(i in 1:length(class.list)){
data.sub.1=subset(data.sub,data.sub$V20==class.list[i])

stats.test=stats.bin(data.sub.1$rel.dist,data.sub.1$V4,N=60)
chg.output[,i]=stats.test$stats["mean",]*100

}

centers=matrix(stats.test$centers,ncol=1)
chg.output=cbind(centers,chg.output)
colnames(chg.output)=c('window_center',class.list)

plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,100),xlab="",ylab='% methylation',main='B73 CHG methylation across MTEC classifications')


for(i in 2:ncol(chg.output)){
lines(cbind(chg.output[,1],chg.output[,i]),col=i,lwd=2)
}
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')

legend("topright",c('LINE_RIL','LINE_RIX','LTR_RLC','LTR_RLG','LTR_RLX','TIR_DTA','TIR_DTC','TIR_DTH','TIR_DTM','TIR_DTT'),col=c(2:11),lty=1,lwd=2)
#

#CHH
#plot
chh.output=matrix(ncol=length(class.list),nrow=60)

for(i in 1:length(class.list)){
data.sub.1=subset(data.sub,data.sub$V20==class.list[i])

stats.test=stats.bin(data.sub.1$rel.dist,data.sub.1$V6,N=60)
chh.output[,i]=stats.test$stats["mean",]*100

}

centers=matrix(stats.test$centers,ncol=1)
chh.output=cbind(centers,chh.output)
colnames(chh.output)=c('window_center',class.list)


plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,40),xlab="",ylab='% methylation',main='B73 CHH methylation across MTEC classifications')


for(i in 2:ncol(chh.output)){
lines(cbind(chh.output[,1],chh.output[,i]),col=i,lwd=2)
}
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')

legend("topright",c('LINE_RIL','LINE_RIX','LTR_RLC','LTR_RLG','LTR_RLX','TIR_DTA','TIR_DTC','TIR_DTH','TIR_DTM','TIR_DTT'),col=c(2:11),lty=1,lwd=2)
#

dev.off()
#

write.table(cg.output,'cg.outputA.txt',sep='\t',row.names=F,quote=F)
write.table(chg.output,'chg.outputA.txt',sep='\t',row.names=F,quote=F)
write.table(chh.output,'chh.outputA.txt',sep='\t',row.names=F,quote=F)