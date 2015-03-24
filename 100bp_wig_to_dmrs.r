options(echo=T)
library(reshape2)
args=commandArgs(trailingOnly=T)
print(args)

# 1. grab wig files
# 2. take each pairwise sample of a context and ID windows showing differences of at least X%
# 3. collapse all windows that are adjacent to each other. Give IDs

#define argments
context=args[1]
difference=as.numeric(args[2])
coverage=as.numeric(args[3])
sitecounts=as.numeric(args[4])

print(args)
#grab all the wig files, select those for your context
a=dir(pattern="*.wig")
a=subset(a,grepl(context,a)==T)
biglist=as.list(a)

aa=combn(a,2)
bb=combn(1:length(biglist),2) #combinations of elements in the biglist


out=NULL
#loop through all pairwise combinations#########################
for(i in 1:(length(aa)/2)){

file1=read.delim(aa[1,i],head=F,skip=1)
file2=read.delim(aa[2,i],head=F,skip=1)

#take windows where there is coverage for both samples
merged=merge(file1,file2,by=c('V1','V2','V3'))

#ID windows that show the selected difference
test.diff=matrix(ifelse(abs(merged$V4.x - merged$V4.y) >= difference,1,0),ncol=1)

merged=cbind(merged,test.diff)

diff.windows=subset(merged,merged$test.diff==1)

#select windows that meet the coverage threshold
diff.windows.cov=subset(diff.windows,diff.windows$V7.x>=coverage & diff.windows$V7.y>=coverage)

#select windows that also meet the sitecount threshold (looking at at least <s> CG/CHG/CHH sites with coverage in the window)
diff.windows.cov=subset(diff.windows.cov,diff.windows.cov$V8.x >=sitecounts & diff.windows.cov$V8.y >=sitecounts)

if(nrow(diff.windows.cov)==0){ next}

diff.windows.cov=diff.windows.cov[with(diff.windows.cov, order(diff.windows.cov[,1],diff.windows.cov[,2])),]
group.id=c(1,rep(NA,nrow(diff.windows.cov)-1))

for(q in 2:nrow(diff.windows.cov)){
	group.id[q]=ifelse(diff.windows.cov[q,2] - diff.windows.cov[q-1,2]<=100,group.id[q-1],group.id[q-1]+1)
	}
group.id=matrix(group.id,ncol=1)
diff.windows.cov=cbind(diff.windows.cov,group.id)

calling=matrix(rep(paste(aa[1,i],"vs",aa[2,i],sep=''),nrow(diff.windows.cov)),ncol=1)

diff.windows.cov=cbind(diff.windows.cov,calling)

out=rbind(out,diff.windows.cov)

}

colnames(out)=c('chr','start','stop','prop1','met1','unmet1','total1','prop2','met2','unmet2','total2','difference.pass','group.id','contrast')

#sort out
out=out[with(out,order(out[,1],out[,2])),]

dmr.id=c(1,rep(NA,nrow(out)-1))


#UPDATED
for(q in 2:nrow(out)){
	dmr.id[q]=ifelse(out[q,2] - out[q-1,2]<=100 & out[q,2] - out[q-1,2] >= 0 & ((out[q,4] + out[q-1,4] <= 200 - 2*difference) | (out[q,4] + out[q-1,4] >= 2*difference)),dmr.id[q-1],dmr.id[q-1]+1)
	}
dmr.id=matrix(dmr.id,ncol=1)

out=cbind(out,dmr.id)

write.table(out,paste('100bp_',context,'_',difference,'diff','_',coverage,'cov2.txt',sep=''),sep='\t',row.names=F)

#collapse that shiz

collapsed.dmrs=matrix(NA,ncol=5,nrow=max(dmr.id))
for(i in 1:max(dmr.id)){
	subs=subset(out,out$dmr.id==i)
	chr=as.character(subs[1,1])
	starts=min(subs[,2])
	stops=max(subs[,3])
	dmrid=i
	size=stops-starts
	collapsed.dmrs[i,]=c(chr,starts,stops,dmrid,size)
	}
	
write.table(collapsed.dmrs,paste('100bp_',context,'_',difference,'diff','_',coverage,'collapsed2.bed',sep=''),sep='\t',row.names=F,col.names=F,quote=F)