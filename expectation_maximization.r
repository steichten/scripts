#####################
#Conducting EM
#
# SRE 1-27-2012
# Modified from Tieming
#
# Look for ###EDIT in code for lines to modify
#####################




#####
#PART 1 - run Mclust to get optimal values for furthur EM
#####

library(mclust)
#dataset (vector, no NAs)
x=c(rnorm(1000,mean=-1),rnorm(500,mean=2))

#the number of clusters (distributions within your dataset)
col=2

#conduct EM at your set number of distributions. "V" sets it as variable variances
data.mclust=Mclust(x,G=col,modelNames="V")

#OUTPUT FILE
#these are the probibilities for all data values  for all distributions created
results=cbind(x,data.mclust$z)

#pull out values of mean, variance, proportion for visualization
#m0: initial values of mean vector mu;
m0=data.mclust$parameters$mean

#v0: initial values of var vector sigma^2;
v0=data.mclust$parameters$variance$sigmasq

#p0: initial values of weights for normal dist's;
p0=data.mclust$parameters$pro

para0=list(m0=m0,v0=v0,p0=p0)
em=list(para0=para0)

#####
#PART 2 - Plot EM distributions
#####

xrng <- range(x)
yrng <- range(dnorm(x,em$para0$m0[2],sqrt(em$para0$v0[2]))*em$para0$p0[2])
h <- max(density(x)$y)
yrng <- c(yrng[1],max(yrng[2],h))
par(cex.main=1.5,cex.axis=1.25,cex.lab=1.7)

#####EDIT
plot(x=NULL, y=NULL,xlim=xrng,ylim=yrng,xlab=paste(colnames(data)[i], "distributions"),ylab="Density",
 main="Raw Data Density, Estimated Mixture Density, \n and Estimated Component Densities")


par(lwd=2)
col1=NULL
nn=length(em$para0$m0)
for (q in 1:nn){ 
	curve(dnorm(x,em$para0$m0[q],sqrt(em$para0$v0[q]))*em$para0$p0[q],
	col=(3*q),add=T)
col1=c(col1,3*q)
}


####EDIT this as needed to obtain the full combined density if you are doing more than 3 subdistributions


curve(dnorm(x,em$para0$m0[1],sqrt(em$para0$v0[1]))*em$para0$p0[1]+
	dnorm(x,em$para0$m0[2],sqrt(em$para0$v0[2]))*em$para0$p0[2],
	col="red",add=T,lty="dashed")

#######################################

lines(density(x)$x, density(x)$y, lty="dashed")

#####EDIT (as necessary)
legend("topright",c("Raw data density","Estimated Mixture Desnity",
		"Estimated Component Density 1","Estimated Component Density 2"),
		col=c('black','red',col1),
		lty=c("dashed","dashed",rep("solid",nn)),cex=0.7)
	

#####END