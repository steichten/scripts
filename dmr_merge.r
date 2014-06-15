options(echo=T)
library(reshape2)
args=commandArgs(trailingOnly=T)
print(args)

############
# quick script to grab bedtools results, sum across C's in DMR window, and make final table
############


#grab all the individual bedtools results (*.dmr)

filelist=dir(pattern="*.dmr$")

#read them in and add a row with the samplename. Rbind them all together
tes=read.delim(filelist[1],head=F)
group=rep(strsplit(filelist[1],"_")[[1]][1],nrow(tes))
tes=cbind(tes,group)
for(i in 2:length(filelist)){
	ss=read.delim(filelist[i],head=F)
	group=rep(strsplit(filelist[i],"_")[[1]][1],nrow(ss))
    ss=cbind(ss,group)
	tes=rbind(tes,ss)
	}

#use dcast (reshape2) to get it into a summary table	
t1=dcast(tes, V1 + V2 + V3 ~ group,value.var='V4')
colnames(t1)[4:ncol(t1)]=paste(names(t1[4:ncol(t1)]),"_prop",sep='')
t2=dcast(tes, V1 + V2 + V3 ~ group,value.var='V5')
colnames(t2)[4:ncol(t2)]=paste(names(t2[4:ncol(t2)]),"_met",sep='')
t3=dcast(tes, V1 + V2 + V3 ~ group,value.var='V6')
colnames(t3)[4:ncol(t3)]=paste(names(t3[4:ncol(t3)]),"_unmet",sep='')

#make a table of it all
tout=cbind(t1,t2[,4:ncol(t2)],t3[,4:ncol(t3)])

#write it out, and write out a version with only rows with data for all samples
write.table(tout,'all_dmr_metvalues.txt',sep='\t',row.names=F,quote=F)
toutnarm=na.omit(tout)
write.table(toutnarm,'all_dmr_metvalues_noNA.txt',sep='\t',row.names=F,quote=F)

#
