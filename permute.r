options(echo=F)
args=commandArgs(trailingOnly=T)
file=args[1]
input=args[2]
contrast=args[3]


infile=read.delim(file,head=F)

pdf(paste(file,'_',input,'_',contrast,'.pdf',sep=''),width=10,height=10)
in.rand=infile[2:nrow(infile),]
in.real=infile[1,]

plot(density(in.rand[,2]),main=paste(file,'_',input,'_',contrast,'.pdf',sep=''),xlim=c(0,max(infile[,2])+10))
abline(v=in.real[,2],col='red')
dev.off()