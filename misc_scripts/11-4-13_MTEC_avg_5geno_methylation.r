#MTEC familyname average methylation across 5 geno data
#11-4-13
#SRE
#This goes through the file (MTEC_5geno_intersect_all.out) containg all 5geno methylation 100bp windows (window100_5geno_data.bed) that have intersected the MTEC NSclass file (MTEC_ind_spread_NSclass.bed)

#V24 is for family, V25 is for 'name'
#prep
data=read.delim('MTEC_5geno_intersect_all.out',head=F)

classification=matrix(names(table(data$V25)),ncol=1)
colnames(classification)='MTEC_name'

#B73
temp.chg=tapply(data$V4,data$V25,mean,na.rm=T)
temp.chh=tapply(data$V5,data$V25,mean,na.rm=T)
temp.cg=tapply(data$V6,data$V25,mean,na.rm=T)

temp.chg=matrix(temp.chg,ncol=1)
temp.chh=matrix(temp.chh,ncol=1)
temp.cg=matrix(temp.cg,ncol=1)

colnames(temp.chg)='B73_chg'
colnames(temp.chh)='B73_chh'
colnames(temp.cg)='B73_cg'

output=cbind(classification,temp.chg,temp.chh,temp.cg)

#CML322
temp.chg=tapply(data$V7,data$V25,mean,na.rm=T)
temp.chh=tapply(data$V8,data$V25,mean,na.rm=T)
temp.cg=tapply(data$V9,data$V25,mean,na.rm=T)

temp.chg=matrix(temp.chg,ncol=1)
temp.chh=matrix(temp.chh,ncol=1)
temp.cg=matrix(temp.cg,ncol=1)

colnames(temp.chg)='CML322_chg'
colnames(temp.chh)='CML322_chh'
colnames(temp.cg)='CML322_cg'

output=cbind(output,temp.chg,temp.chh,temp.cg)


#Mo17
temp.chg=tapply(data$V10,data$V25,mean,na.rm=T)
temp.chh=tapply(data$V11,data$V25,mean,na.rm=T)
temp.cg=tapply(data$V12,data$V25,mean,na.rm=T)

temp.chg=matrix(temp.chg,ncol=1)
temp.chh=matrix(temp.chh,ncol=1)
temp.cg=matrix(temp.cg,ncol=1)

colnames(temp.chg)='Mo17_chg'
colnames(temp.chh)='Mo17_chh'
colnames(temp.cg)='Mo17_cg'

output=cbind(output,temp.chg,temp.chh,temp.cg)

#Oh43
temp.chg=tapply(data$V13,data$V25,mean,na.rm=T)
temp.chh=tapply(data$V14,data$V25,mean,na.rm=T)
temp.cg=tapply(data$V15,data$V25,mean,na.rm=T)

temp.chg=matrix(temp.chg,ncol=1)
temp.chh=matrix(temp.chh,ncol=1)
temp.cg=matrix(temp.cg,ncol=1)

colnames(temp.chg)='Oh43_chg'
colnames(temp.chh)='Oh43_chh'
colnames(temp.cg)='Oh43_cg'

output=cbind(output,temp.chg,temp.chh,temp.cg)

#Tx303
temp.chg=tapply(data$V16,data$V25,mean,na.rm=T)
temp.chh=tapply(data$V17,data$V25,mean,na.rm=T)
temp.cg=tapply(data$V18,data$V25,mean,na.rm=T)

temp.chg=matrix(temp.chg,ncol=1)
temp.chh=matrix(temp.chh,ncol=1)
temp.cg=matrix(temp.cg,ncol=1)

colnames(temp.chg)='Tx303_chg'
colnames(temp.chh)='Tx303_chh'
colnames(temp.cg)='Tx303_cg'

output=cbind(output,temp.chg,temp.chh,temp.cg)

write.table(output,'5GENO_WINDOW100_AVG_MTEC_NAME_11-5.txt',sep='\t',row.names=F,quote=F)
