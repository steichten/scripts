data=read.delim('chia_IBD_5geno_data.out.txt',head=F)

#grab the methylation data from the first genotype in IBD
first.chg=matrix(ifelse(data$V24=='B73',data$V4,ifelse(data$V24=='MO17',data$V10,ifelse(data$V24=='CML322',data$V7,ifelse(data$V24=='OH43',data$V13,ifelse(data$V24=='TX303',data$V16,NA))))),ncol=1)

first.chh=matrix(ifelse(data$V24=='B73',data$V5,ifelse(data$V24=='MO17',data$V11,ifelse(data$V24=='CML322',data$V8,ifelse(data$V24=='OH43',data$V14,ifelse(data$V24=='TX303',data$V17,NA))))),ncol=1)

first.cg=matrix(ifelse(data$V24=='B73',data$V6,ifelse(data$V24=='MO17',data$V12,ifelse(data$V24=='CML322',data$V9,ifelse(data$V24=='OH43',data$V15,ifelse(data$V24=='TX303',data$V18,NA))))),ncol=1)

#get the last
last.chg=matrix(ifelse(data$V25=='B73',data$V4,ifelse(data$V25=='MO17',data$V10,ifelse(data$V25=='CML322',data$V7,ifelse(data$V25=='OH43',data$V13,ifelse(data$V25=='TX303',data$V16,NA))))),ncol=1)

last.chh=matrix(ifelse(data$V25=='B73',data$V5,ifelse(data$V25=='MO17',data$V11,ifelse(data$V25=='CML322',data$V8,ifelse(data$V25=='OH43',data$V14,ifelse(data$V25=='TX303',data$V17,NA))))),ncol=1)

last.cg=matrix(ifelse(data$V25=='B73',data$V6,ifelse(data$V25=='MO17',data$V12,ifelse(data$V25=='CML322',data$V9,ifelse(data$V25=='OH43',data$V15,ifelse(data$V25=='TX303',data$V18,NA))))),ncol=1)

ibd.genos=cbind(first.chg,first.chh,first.cg,last.chg,last.chh,last.cg)

colnames(ibd.genos)=c('first.chg','first.chh','first.cg','last.chg','last.chh','last.cg')
data=cbind(data,ibd.genos)

#get the difference as last geno - first geno
diff.chg=matrix(data$last.chg - data$first.chg,ncol=1)
diff.chh=matrix(data$last.chh - data$first.chh,ncol=1)
diff.cg=matrix(data$last.cg - data$first.cg,ncol=1)

#plot out the differences
plot(density(diff.chg,na.rm=T),main='difference in last-first methylation in IBD')
lines(density(diff.chh,na.rm=T),col='red')
lines(density(diff.cg,na.rm=T),col='blue')
legend('topright',c('chg','chh','cg'),col=c('black','red','blue'),lty=1)

#develop the proprtion of windows showing differential methylation
exp.chg.diff=(length(subset(diff.chg,is.na(diff.chg)==F & (diff.chg>50 | diff.chg< -50))))/length(subset(diff.chg,is.na(diff.chg)==F))
exp.chh.diff=(length(subset(diff.chh,is.na(diff.chh)==F & (diff.chh>50 | diff.chh< -50))))/length(subset(diff.chh,is.na(diff.chh)==F))
exp.cg.diff=(length(subset(diff.cg,is.na(diff.cg)==F & (diff.cg>50 | diff.cg< -50))))/length(subset(diff.cg,is.na(diff.cg)==F))

#geno methylation in IBD vs chr10


#read in chr10 data
chr10=read.delim('temp_chr10_grab.txt',head=F)

par(mfrow=c(3,5))

#here I am getting the difference in b73 CHG compared to chg in the other 4 genotypes. I made a few temp plots as we go
#B_CHG
dif1=chr10$V4-chr10$V7
dif2=chr10$V4-chr10$V10
dif3=chr10$V4-chr10$V13
dif4=chr10$V4-chr10$V16
plot(density(dif1,na.rm=T),main='B73_CHG')
lines(density(dif2,na.rm=T))
lines(density(dif3,na.rm=T))
lines(density(dif4,na.rm=T))

a1=length(subset(dif1,is.na(dif1)==F))
a2=length(subset(dif2,is.na(dif2)==F))
a3=length(subset(dif3,is.na(dif3)==F))
a4=length(subset(dif4,is.na(dif4)==F))

#no for the other genos...
#CML322_CHG
dif5=chr10$V7-chr10$V10
dif6=chr10$V7-chr10$V13
dif7=chr10$V7-chr10$V16
dif8=chr10$V7-chr10$V4
plot(density(dif5,na.rm=T),main='CML322_CHG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a5=length(subset(dif5,is.na(dif5)==F))
a6=length(subset(dif6,is.na(dif6)==F))
a7=length(subset(dif7,is.na(dif7)==F))
a8=length(subset(dif8,is.na(dif8)==F))

#MO17_CHG
dif5=chr10$V10-chr10$V7
dif6=chr10$V10-chr10$V13
dif7=chr10$V10-chr10$V16
dif8=chr10$V10-chr10$V4
plot(density(dif5,na.rm=T),main='Mo17_CHG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a9=length(subset(dif5,is.na(dif5)==F))
a10=length(subset(dif6,is.na(dif6)==F))
a11=length(subset(dif7,is.na(dif7)==F))
a12=length(subset(dif8,is.na(dif8)==F))

#OH43_CHG
dif5=chr10$V13-chr10$V7
dif6=chr10$V13-chr10$V10
dif7=chr10$V13-chr10$V16
dif8=chr10$V13-chr10$V4
plot(density(dif5,na.rm=T),main='Oh43_CHG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a13=length(subset(dif5,is.na(dif5)==F))
a14=length(subset(dif6,is.na(dif6)==F))
a15=length(subset(dif7,is.na(dif7)==F))
a16=length(subset(dif8,is.na(dif8)==F))

#Tx303_CHG
dif5=chr10$V16-chr10$V7
dif6=chr10$V16-chr10$V10
dif7=chr10$V16-chr10$V13
dif8=chr10$V16-chr10$V4
plot(density(dif5,na.rm=T),main='Tx303_CHG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a17=length(subset(dif5,is.na(dif5)==F))
a18=length(subset(dif6,is.na(dif6)==F))
a19=length(subset(dif7,is.na(dif7)==F))
a20=length(subset(dif8,is.na(dif8)==F))

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

a21=length(subset(dif1,is.na(dif1)==F))
a22=length(subset(dif2,is.na(dif2)==F))
a23=length(subset(dif3,is.na(dif3)==F))
a24=length(subset(dif4,is.na(dif4)==F))
#CML322_CHh
dif5=chr10$V8-chr10$V11
dif6=chr10$V8-chr10$V14
dif7=chr10$V8-chr10$V17
dif8=chr10$V8-chr10$V5
plot(density(dif5,na.rm=T),main='CML322_CHH')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a25=length(subset(dif5,is.na(dif5)==F))
a26=length(subset(dif6,is.na(dif6)==F))
a27=length(subset(dif7,is.na(dif7)==F))
a28=length(subset(dif8,is.na(dif8)==F))

#MO17_CHh
dif5=chr10$V11-chr10$V8
dif6=chr10$V11-chr10$V14
dif7=chr10$V11-chr10$V17
dif8=chr10$V11-chr10$V5
plot(density(dif5,na.rm=T),main='Mo17_CHH')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a29=length(subset(dif5,is.na(dif5)==F))
a30=length(subset(dif6,is.na(dif6)==F))
a31=length(subset(dif7,is.na(dif7)==F))
a32=length(subset(dif8,is.na(dif8)==F))

#OH43_CHh
dif5=chr10$V14-chr10$V8
dif6=chr10$V14-chr10$V11
dif7=chr10$V14-chr10$V17
dif8=chr10$V14-chr10$V5
plot(density(dif5,na.rm=T),main='Oh43_CHH')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a33=length(subset(dif5,is.na(dif5)==F))
a34=length(subset(dif6,is.na(dif6)==F))
a35=length(subset(dif7,is.na(dif7)==F))
a36=length(subset(dif8,is.na(dif8)==F))

#Tx303_CHh
dif5=chr10$V17-chr10$V8
dif6=chr10$V17-chr10$V11
dif7=chr10$V17-chr10$V14
dif8=chr10$V17-chr10$V5
plot(density(dif5,na.rm=T),main='Tx303_CHH')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a37=length(subset(dif5,is.na(dif5)==F))
a38=length(subset(dif6,is.na(dif6)==F))
a39=length(subset(dif7,is.na(dif7)==F))
a40=length(subset(dif8,is.na(dif8)==F))
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

a41=length(subset(dif1,is.na(dif1)==F))
a42=length(subset(dif2,is.na(dif2)==F))
a43=length(subset(dif3,is.na(dif3)==F))
a44=length(subset(dif4,is.na(dif4)==F))

#CML322_Cg
dif5=chr10$V9-chr10$V12
dif6=chr10$V9-chr10$V15
dif7=chr10$V9-chr10$V18
dif8=chr10$V9-chr10$V6
plot(density(dif5,na.rm=T),main='CML322_CG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a45=length(subset(dif5,is.na(dif5)==F))
a46=length(subset(dif6,is.na(dif6)==F))
a47=length(subset(dif7,is.na(dif7)==F))
a48=length(subset(dif8,is.na(dif8)==F))

#MO17_Cg
dif5=chr10$V12-chr10$V9
dif6=chr10$V12-chr10$V15
dif7=chr10$V12-chr10$V18
dif8=chr10$V12-chr10$V6
plot(density(dif5,na.rm=T),main='Mo17_CG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a49=length(subset(dif5,is.na(dif5)==F))
a50=length(subset(dif6,is.na(dif6)==F))
a51=length(subset(dif7,is.na(dif7)==F))
a52=length(subset(dif8,is.na(dif8)==F))

#OH43_Cg
dif5=chr10$V15-chr10$V9
dif6=chr10$V15-chr10$V12
dif7=chr10$V15-chr10$V18
dif8=chr10$V15-chr10$V6
plot(density(dif5,na.rm=T),main='Oh43_CG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a53=length(subset(dif5,is.na(dif5)==F))
a54=length(subset(dif6,is.na(dif6)==F))
a55=length(subset(dif7,is.na(dif7)==F))
a56=length(subset(dif8,is.na(dif8)==F))

#Tx303_Cg
dif5=chr10$V18-chr10$V9
dif6=chr10$V18-chr10$V12
dif7=chr10$V18-chr10$V15
dif8=chr10$V18-chr10$V6
plot(density(dif5,na.rm=T),main='Tx303_CG')
lines(density(dif6,na.rm=T))
lines(density(dif7,na.rm=T))
lines(density(dif8,na.rm=T))

a57=length(subset(dif5,is.na(dif5)==F))
a58=length(subset(dif6,is.na(dif6)==F))
a59=length(subset(dif7,is.na(dif7)==F))
a60=length(subset(dif8,is.na(dif8)==F))
#
#so I made a long list of the number of windows that we have for each geno-methylation type combination. 
total.windows=c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60)

#now to get the number of windows that we would call 'differentialy methylated'

#B_CHG
dif1=chr10$V4-chr10$V7
dif2=chr10$V4-chr10$V10
dif3=chr10$V4-chr10$V13
dif4=chr10$V4-chr10$V16
a1=length(subset(dif1,is.na(dif1)==F & (dif1>50 | dif1< -50)))
a2=length(subset(dif2,is.na(dif2)==F & (dif2>50 | dif2< -50)))
a3=length(subset(dif3,is.na(dif3)==F & (dif3>50 | dif3< -50)))
a4=length(subset(dif4,is.na(dif4)==F & (dif4>50 | dif4< -50)))

#CML322_CHG
dif5=chr10$V7-chr10$V10
dif6=chr10$V7-chr10$V13
dif7=chr10$V7-chr10$V16
dif8=chr10$V7-chr10$V4
a5=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a6=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a7=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a8=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#MO17_CHG
dif5=chr10$V10-chr10$V7
dif6=chr10$V10-chr10$V13
dif7=chr10$V10-chr10$V16
dif8=chr10$V10-chr10$V4
a9=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a10=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a11=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a12=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#OH43_CHG
dif5=chr10$V13-chr10$V7
dif6=chr10$V13-chr10$V10
dif7=chr10$V13-chr10$V16
dif8=chr10$V13-chr10$V4
a13=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a14=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a15=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a16=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#Tx303_CHG
dif5=chr10$V16-chr10$V7
dif6=chr10$V16-chr10$V10
dif7=chr10$V16-chr10$V13
dif8=chr10$V16-chr10$V4
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
a21=length(subset(dif1,is.na(dif1)==F & (dif1>50 | dif1< -50)))
a22=length(subset(dif2,is.na(dif2)==F & (dif2>50 | dif2< -50)))
a23=length(subset(dif3,is.na(dif3)==F & (dif3>50 | dif3< -50)))
a24=length(subset(dif4,is.na(dif4)==F & (dif4>50 | dif4< -50)))
#CML322_CHh
dif5=chr10$V8-chr10$V11
dif6=chr10$V8-chr10$V14
dif7=chr10$V8-chr10$V17
dif8=chr10$V8-chr10$V5
a25=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a26=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a27=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a28=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#MO17_CHh
dif5=chr10$V11-chr10$V8
dif6=chr10$V11-chr10$V14
dif7=chr10$V11-chr10$V17
dif8=chr10$V11-chr10$V5
a29=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a30=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a31=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a32=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#OH43_CHh
dif5=chr10$V14-chr10$V8
dif6=chr10$V14-chr10$V11
dif7=chr10$V14-chr10$V17
dif8=chr10$V14-chr10$V5
a33=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a34=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a35=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a36=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#Tx303_CHh
dif5=chr10$V17-chr10$V8
dif6=chr10$V17-chr10$V11
dif7=chr10$V17-chr10$V14
dif8=chr10$V17-chr10$V5
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
a41=length(subset(dif1,is.na(dif1)==F & (dif1>50 | dif1< -50)))
a42=length(subset(dif2,is.na(dif2)==F & (dif2>50 | dif2< -50)))
a43=length(subset(dif3,is.na(dif3)==F & (dif3>50 | dif3< -50)))
a44=length(subset(dif4,is.na(dif4)==F & (dif4>50 | dif4< -50)))

#CML322_Cg
dif5=chr10$V9-chr10$V12
dif6=chr10$V9-chr10$V15
dif7=chr10$V9-chr10$V18
dif8=chr10$V9-chr10$V6
a45=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a46=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a47=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a48=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#MO17_Cg
dif5=chr10$V12-chr10$V9
dif6=chr10$V12-chr10$V15
dif7=chr10$V12-chr10$V18
dif8=chr10$V12-chr10$V6
a49=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a50=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a51=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a52=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#OH43_Cg
dif5=chr10$V15-chr10$V9
dif6=chr10$V15-chr10$V12
dif7=chr10$V15-chr10$V18
dif8=chr10$V15-chr10$V6
a53=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a54=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a55=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a56=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))

#Tx303_Cg
dif5=chr10$V18-chr10$V9
dif6=chr10$V18-chr10$V12
dif7=chr10$V18-chr10$V15
dif8=chr10$V18-chr10$V6
a57=length(subset(dif5,is.na(dif5)==F & (dif5>50 | dif5< -50)))
a58=length(subset(dif6,is.na(dif6)==F & (dif6>50 | dif6< -50)))
a59=length(subset(dif7,is.na(dif7)==F & (dif7>50 | dif7< -50)))
a60=length(subset(dif8,is.na(dif8)==F & (dif8>50 | dif8< -50)))
#

#so here are all the diff.met windows per chr10 contrast
diff=c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60)

#make the proportion
prop=diff/total.windows

#split by the met type
chg=prop[1:20]
chh=prop[21:40]
cg=prop[41:60]

contrast.names=matrix(c('B-CML','B-Mo17','B-Oh43','B-Tx','CML-Mo17','CML-Oh43','CML-Tx','CML-B','Mo17-CML','Mo17-Oh43','Mo17-Tx','Mo17-B','Oh43-CML','Oh43-Mo17','Oh43-Tx','Oh43-B','Tx-CML','Tx-Mo17','Tx-Oh43','Tx-B'),ncol=1)
#name.list=c(paste(contrast.names,'CHG',sep='_'),paste(contrast.names,'CHH',sep='_'),paste(contrast.names,'CG',sep='_'))

out.table=cbind(contrast.names,matrix(chg,ncol=1),matrix(chh,ncol=1),matrix(cg,ncol=1))
colnames(out.table)=c('contrast','chg','chh','cg')
write.table(out.table,'chia_IBD_chr10_propDiff.txt',sep='\t',row.names=F,quote=F)

#make the plot comparing all of the contrast proportions with our IBD proportion

par(mfrow=c(1,3))
boxplot(chg,ylim=c(0,0.1),main='dist of prop chg difference on chr10')
abline(h=exp.chg.diff,col=2)
boxplot(chh,ylim=c(0,0.01),main='dist of prop chh difference on chr10')
abline(h=exp.chh.diff,col=2)
boxplot(cg,ylim=c(0,0.1),main='dist of prop cg difference on chr10')
abline(h=exp.cg.diff,col=2)



##################
#see if there are any concurrent windows (or at least close to each other)

#subset the dataset based on our filter
data.sub.cg=subset(data,data$diff.cg < -50 | data$diff.cg > 50)
data.sub.chg=subset(data,data$diff.chg < -50 | data$diff.chg > 50)
data.sub.chh=subset(data,data$diff.chh < -50 | data$diff.chh > 50)

write.table(data.sub.cg,'data.sub.cg.txt',sep='\t',row.names=F,quote=F)
write.table(data.sub.chg,'data.sub.chg.txt',sep='\t',row.names=F,quote=F)
write.table(data.sub.chh,'data.sub.chh.txt',sep='\t',row.names=F,quote=F)

#add dist_from_last and uid column in excel. Then reimport

table(table(data.sub.cg$uid))
table(table(data.sub.chg$uid))
table(table(data.sub.chh$uid))

#then something along these lines to divide by the ibd geno pairs.
v1=matrix(paste(data.sub.cg$V24,data.sub.cg$V25,sep='_'),ncol=1)
ddd=subset(data.sub.cg,data.sub.cg$v1=='CML322_B73')
table(table(ddd[,5]))
