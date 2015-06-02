data=read.delim('chia_IBD_5geno_data.out',head=F)


first.chg=matrix(ifelse(data$V24=='B73',data$V4,ifelse(data$V24=='MO17',data$V10,ifelse(data$V24=='CML322',data$V7,ifelse(data$V24=='OH43',data$V13,ifelse(data$V24=='TX303',data$V16,NA))))),ncol=1)

first.chh=matrix(ifelse(data$V24=='B73',data$V5,ifelse(data$V24=='MO17',data$V11,ifelse(data$V24=='CML322',data$V8,ifelse(data$V24=='OH43',data$V14,ifelse(data$V24=='TX303',data$V17,NA))))),ncol=1)

first.cg=matrix(ifelse(data$V24=='B73',data$V6,ifelse(data$V24=='MO17',data$V12,ifelse(data$V24=='CML322',data$V9,ifelse(data$V24=='OH43',data$V15,ifelse(data$V24=='TX303',data$V18,NA))))),ncol=1)


last.chg=matrix(ifelse(data$V25=='B73',data$V4,ifelse(data$V25=='MO17',data$V10,ifelse(data$V25=='CML322',data$V7,ifelse(data$V25=='OH43',data$V13,ifelse(data$V25=='TX303',data$V16,NA))))),ncol=1)

last.chh=matrix(ifelse(data$V25=='B73',data$V5,ifelse(data$V25=='MO17',data$V11,ifelse(data$V25=='CML322',data$V8,ifelse(data$V25=='OH43',data$V14,ifelse(data$V25=='TX303',data$V17,NA))))),ncol=1)

last.cg=matrix(ifelse(data$V25=='B73',data$V6,ifelse(data$V25=='MO17',data$V12,ifelse(data$V25=='CML322',data$V9,ifelse(data$V25=='OH43',data$V15,ifelse(data$V25=='TX303',data$V18,NA))))),ncol=1)

ibd.genos=cbind(first.chg,first.chh,first.cg,last.chg,last.chh,last.cg)

colnames(ibd.genos)=c('first.chg','first.chh','first.cg','last.chg','last.chh','last.cg')
data=cbind(data,ibd.genos)

diff.chg=matrix(data$last.chg - data$first.chg,ncol=1)
diff.chh=matrix(data$last.chh - data$first.chh,ncol=1)
diff.cg=matrix(data$last.cg - data$first.cg,ncol=1)

plot(density(diff.chg,na.rm=T),main='difference in last-first methylation in IBD')
lines(density(diff.chh,na.rm=T),col='red')
lines(density(diff.cg,na.rm=T),col='blue')
legend('topright',c('chg','chh','cg'),col=c('black','red','blue'),lty=1)


tapply(data$first.chg,data$V24,mean,na.rm=T)
tapply(data$first.chh,data$V24,mean,na.rm=T)
tapply(data$first.cg,data$V24,mean,na.rm=T)
tapply(data$last.chg,data$V24,mean,na.rm=T)
tapply(data$last.chh,data$V24,mean,na.rm=T)
tapply(data$last.cg,data$V24,mean,na.rm=T)


#geno methylation in IBD vs elsewhere
output=matrix(NA,ncol=4,nrow=5)
output[,1]=c('B73','Mo17','Oh43','CML322','Tx303')

b1=subset(data,data$V24=='B73' | data$V25=='B73')
output[1,2]=mean(b1$V4,na.rm=T)
output[1,3]=mean(b1$V5,na.rm=T)
output[1,4]=mean(b1$V6,na.rm=T)

m1=subset(data,data$V24=='MO17' | data$V25=='MO17')
output[2,2]=mean(m1$V10,na.rm=T)
output[2,3]=mean(m1$V11,na.rm=T)
output[2,4]=mean(m1$V12,na.rm=T)
o1=subset(data,data$V24=='OH43' | data$V25=='OH43')
output[3,2]=mean(o1$V13,na.rm=T)
output[3,3]=mean(o1$V14,na.rm=T)
output[3,4]=mean(o1$V15,na.rm=T)
c1=subset(data,data$V24=='CML322' | data$V25=='CML322')
output[4,2]=mean(c1$V7,na.rm=T)
output[4,3]=mean(c1$V8,na.rm=T)
output[4,4]=mean(c1$V9,na.rm=T)
t1=subset(data,data$V24=='TX303' | data$V25=='TX303')
output[5,2]=mean(t1$V16,na.rm=T)
output[5,3]=mean(t1$V17,na.rm=T)
output[5,4]=mean(t1$V18,na.rm=T)

colnames(output)=c('genotype','CHG','CHH','CG')

write.table(output,'IBD_geno_avg_methylation.txt',sep='\t',row.names=F,quote=F)
#


chr10=read.delim('temp_chr10_grab.txt',head=F)

par(mfrow=c(3,5))

#B_CHG
dif1=chr10$V4-chr10$V7
dif2=chr10$V4-chr10$V10
dif3=chr10$V4-chr10$V13
dif4=chr10$V4-chr10$V16
plot(density(dif1,na.rm=T),main='B73_CHG')
lines(density(dif2,na.rm=T))
lines(density(dif3,na.rm=T))
lines(density(dif4,na.rm=T))

a1=length(subset(dif1,is.na(dif1)==F & (dif1>50 | dif1< -50)))
a2=length(subset(dif2,is.na(dif2)==F & (dif2>50 | dif2< -50)))
a3=length(subset(dif3,is.na(dif3)==F & (dif3>50 | dif3< -50)))
a4=length(subset(dif4,is.na(dif4)==F & (dif4>50 | dif4< -50)))

#CML322_CHG
dif5=chr10$V7-chr10$V10
dif6=chr10$V7-chr10$V13
dif7=chr10$V7-chr10$V16
dif8=chr10$V7-chr10$V4
plot(density(dif5,na.rm=T),main='CML322_CHG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a5=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a6=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a7=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a8=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#MO17_CHG
dif5=chr10$V10-chr10$V7
dif6=chr10$V10-chr10$V13
dif7=chr10$V10-chr10$V16
dif8=chr10$V10-chr10$V4
plot(density(dif5,na.rm=T),main='Mo17_CHG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a9=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a10=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a11=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a12=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#OH43_CHG
dif5=chr10$V13-chr10$V7
dif6=chr10$V13-chr10$V10
dif7=chr10$V13-chr10$V16
dif8=chr10$V13-chr10$V4
plot(density(dif5,na.rm=T),main='Oh43_CHG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a13=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a14=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a15=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a16=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#Tx303_CHG
dif5=chr10$V16-chr10$V7
dif6=chr10$V16-chr10$V10
dif7=chr10$V16-chr10$V13
dif8=chr10$V16-chr10$V4
plot(density(dif5,na.rm=T),main='Tx303_CHG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a17=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a18=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a19=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a20=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#
#B_CHh
dif1=chr10$V5-chr10$V8
dif2=chr10$V5-chr10$V11
dif3=chr10$V5-chr10$V14
dif4=chr10$V5-chr10$V17
plot(density(dif1,na.rm=T),main='B73_CHH')
lines(density(dif2,na.rm=T))
lines(density(dif3,na.rm=T))
lines(density(dif4,na.rm=T))

a21=length(subset(dif1,is.na(dif1)==F & (dif1>50 | dif1< -50)))
a22=length(subset(dif2,is.na(dif2)==F & (dif2>50 | dif2< -50)))
a23=length(subset(dif3,is.na(dif3)==F & (dif3>50 | dif3< -50)))
a24=length(subset(dif4,is.na(dif4)==F & (dif4>50 | dif4< -50)))
#CML322_CHh
dif5=chr10$V8-chr10$V11
dif6=chr10$V8-chr10$V14
dif7=chr10$V8-chr10$V17
dif8=chr10$V8-chr10$V5
plot(density(dif5,na.rm=T),main='CML322_CHH')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a25=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a26=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a27=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a28=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#MO17_CHh
dif5=chr10$V11-chr10$V8
dif6=chr10$V11-chr10$V14
dif7=chr10$V11-chr10$V17
dif8=chr10$V11-chr10$V5
plot(density(dif5,na.rm=T),main='Mo17_CHH')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a29=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a30=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a31=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a32=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#OH43_CHh
dif5=chr10$V14-chr10$V8
dif6=chr10$V14-chr10$V11
dif7=chr10$V14-chr10$V17
dif8=chr10$V14-chr10$V5
plot(density(dif5,na.rm=T),main='Oh43_CHH')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a33=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a34=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a35=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a36=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#Tx303_CHh
dif5=chr10$V17-chr10$V8
dif6=chr10$V17-chr10$V11
dif7=chr10$V17-chr10$V14
dif8=chr10$V17-chr10$V5
plot(density(dif5,na.rm=T),main='Tx303_CHH')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a37=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a38=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a39=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a40=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))
#

#B_Cg
dif1=chr10$V6-chr10$V9
dif2=chr10$V6-chr10$V12
dif3=chr10$V6-chr10$V15
dif4=chr10$V6-chr10$V18
plot(density(dif1,na.rm=T),main='B73_CG')
lines(density(dif2,na.rm=T))
lines(density(dif3,na.rm=T))
lines(density(dif4,na.rm=T))

a41=length(subset(dif1,is.na(dif1)==F & (dif1>50 | dif1< -50)))
a42=length(subset(dif2,is.na(dif2)==F & (dif2>50 | dif2< -50)))
a43=length(subset(dif3,is.na(dif3)==F & (dif3>50 | dif3< -50)))
a44=length(subset(dif4,is.na(dif4)==F & (dif4>50 | dif4< -50)))

#CML322_Cg
dif5=chr10$V9-chr10$V12
dif6=chr10$V9-chr10$V15
dif7=chr10$V9-chr10$V18
dif8=chr10$V9-chr10$V6
plot(density(dif5,na.rm=T),main='CML322_CG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a45=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a46=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a47=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a48=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#MO17_Cg
dif5=chr10$V12-chr10$V9
dif6=chr10$V12-chr10$V15
dif7=chr10$V12-chr10$V18
dif8=chr10$V12-chr10$V6
plot(density(dif5,na.rm=T),main='Mo17_CG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a49=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a50=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a51=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a52=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#OH43_Cg
dif5=chr10$V15-chr10$V9
dif6=chr10$V15-chr10$V12
dif7=chr10$V15-chr10$V18
dif8=chr10$V15-chr10$V6
plot(density(dif5,na.rm=T),main='Oh43_CG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a53=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a54=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a55=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a56=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#Tx303_Cg
dif5=chr10$V18-chr10$V9
dif6=chr10$V18-chr10$V12
dif7=chr10$V18-chr10$V15
dif8=chr10$V18-chr10$V6
plot(density(dif5,na.rm=T),main='Tx303_CG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a57=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a58=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a59=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a60=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))
#

out=c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60)

chg=out[1:20]
chh=out[21:40]
cg=out[41:60]

exp.chg.diff=length(subset(diff.chg,is.na(diff.chg)==F & (diff.chg>50 | diff.chg< -50)))
exp.chh.diff=length(subset(diff.chh,is.na(diff.chh)==F & (diff.chh>50 | diff.chh< -50)))
exp.cg.diff=length(subset(diff.cg,is.na(diff.cg)==F & (diff.cg>50 | diff.cg< -50)))

par(mfrow=c(1,3))
boxplot(chg,main='dist of chg difference on chr10')
abline(h=exp.chg.diff,col=2)
boxplot(chh,main='dist of chh difference on chr10')
abline(h=exp.chh.diff,col=2)
boxplot(cg,ylim=c(10000,40000),main='dist of cg difference on chr10')
abline(h=exp.cg.diff,col=2)