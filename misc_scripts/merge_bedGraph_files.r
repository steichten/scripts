#5geno - merge bedgraph files from jawon into single 100bp window file
#10-30-13
#SRE
#make 5 files containg CHG,CHH,and CG methylation for each of the genotypes.

#make a folder of the bedgraph files
a=dir('bedgraph')


#B73
data1=read.delim(paste('bedgraph/',a[1],sep=''),head=F)
colnames(data1)[4]=a[1]


for(i in 2:3){
data2=read.delim(paste('bedgraph/',a[i],sep=''),head=F)

temp1=merge(data1,data2,by=c('V1','V2','V3'),all=T)
colnames(temp1)[i+3]=a[i]
data1=temp1

}

write.table(data1,'bedgraph/window100_5geno_B73.bed',sep='\t',row.names=F,quote=F)

#CML322
data1=read.delim(paste('bedgraph/',a[4],sep=''),head=F)
colnames(data1)[4]=a[4]


for(i in 5:6){
data2=read.delim(paste('bedgraph/',a[i],sep=''),head=F)

temp1=merge(data1,data2,by=c('V1','V2','V3'),all=T)
colnames(temp1)[i]=a[i]
data1=temp1

}

write.table(data1,'bedgraph/window100_5geno_CML322.bed',sep='\t',row.names=F,quote=F)

#Mo17
data1=read.delim(paste('bedgraph/',a[7],sep=''),head=F)
colnames(data1)[4]=a[7]


for(i in 8:9){
data2=read.delim(paste('bedgraph/',a[i],sep=''),head=F)

temp1=merge(data1,data2,by=c('V1','V2','V3'),all=T)
colnames(temp1)[i-3]=a[i]
data1=temp1

}

write.table(data1,'bedgraph/window100_5geno_Mo17.bed',sep='\t',row.names=F,quote=F)




#Oh43
data1=read.delim(paste('bedgraph/',a[10],sep=''),head=F)
colnames(data1)[4]=a[10]


for(i in 11:12){
data2=read.delim(paste('bedgraph/',a[i],sep=''),head=F)

temp1=merge(data1,data2,by=c('V1','V2','V3'),all=T)
colnames(temp1)[i-6]=a[i]
data1=temp1

}

write.table(data1,'bedgraph/window100_5geno_Oh43.bed',sep='\t',row.names=F,quote=F)


#Tx303
data1=read.delim(paste('bedgraph/',a[13],sep=''),head=F)
colnames(data1)[4]=a[13]


for(i in 14:15){
data2=read.delim(paste('bedgraph/',a[i],sep=''),head=F)

temp1=merge(data1,data2,by=c('V1','V2','V3'),all=T)
colnames(temp1)[i-9]=a[i]
data1=temp1

}

write.table(data1,'bedgraph/window100_5geno_Tx303.bed',sep='\t',row.names=F,quote=F)
#