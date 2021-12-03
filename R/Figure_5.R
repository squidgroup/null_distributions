rm(list=ls())

library(beeswarm)
library(scales)
wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/functions.R"))

hist(boot::inv.logit(rnorm(100000,0,sqrt(0.2))), breaks=50)


load( file=paste0(wd,"Data/Intermediate/GLMM_sim.Rdata"))
alt <- as.data.frame(
	t(sapply(results, function(x) c(x$param,x$actual[1:6])))
)

p_alt <- as.data.frame(
	t(sapply(results, function(x) c(x$param,
		freq = p_func(x$actual["freq"],x$null[,"freq"]),
		mean=p_func(x$actual["mean"],x$null[,"mean"]),
		median=p_func(x$actual["median"],x$null[,"median"]),
		mode= p_func(x$actual["mode"],x$null[,"mode"])
		# ,x$actual["LRT"]
		)))
)

load( file=paste0(wd,"Data/Intermediate/GLMM_sim_null.Rdata"))
null <- as.data.frame(
	t(sapply(results, function(x) c(x$param,x$actual[1:6])))
)
null$sim_var <- 0

p_null <- as.data.frame(
	t(sapply(results, function(x) c(x$param,
		freq = p_func(x$actual["freq"],x$null[,"freq"]),
		mean=p_func(x$actual["mean"],x$null[,"mean"]),
		median=p_func(x$actual["median"],x$null[,"median"]),
		mode= p_func(x$actual["mode"],x$null[,"mode"])
		# ,x$actual["LRT"]
		)))
)
p_null$sim_var <- 0

long_alt <- Stack(alt,c("freq","mean","median","mode"), group.name="type",value.name = "estimate")
long_null <- Stack(null,c("freq","mean","median","mode"), group.name="type",value.name = "estimate")

p_long_alt <- Stack(p_alt,c("freq","mean","median","mode"), group.name="type",value.name = "estimate")#,"LRT"
p_long_null <- Stack(p_null,c("freq","mean","median","mode"), group.name="type",value.name = "estimate")#,"LRT"

long <- rbind(long_alt,long_null)
p_long <- rbind(p_long_alt,p_long_null)


setEPS()
pdf(paste0(wd,"Figures/Fig5.pdf"), height=6, width=8)

{
par(mfrow=c(2,1), mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))
# boxplot(estimate~N_within+type,bp2_long,boxwex=0.5,xlab="Within group sample size", names=rep(c(2,4,8),3), ylab="Estimate", pch=19, cex=0.5)
beeswarm(estimate~sim_var+type,long, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Between group variance", ylab="Estimate", labels=rep(c(0,0.2),4))#, 
mtext(c("ML","Post. Mean","Post. Median","Post. Mode"),side=3,adj=c(0.1,0.35,0.65 ,0.9))
means<-aggregate(estimate~sim_var+type,long,mean)$estimate
CIs<-aggregate(estimate~sim_var+type,long,function(x) sd(x)/sqrt(length(x)))$estimate*qnorm(0.975)
points(means, pch=19, cex=1.5, col="blue")
arrows(1:9,means+CIs,1:9,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)
abline(h=0.2, col="red", lwd=2)
abline(v=c(2.5,4.5,6.5))
# mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))

# beeswarm(estimate~sim_var+type,p_long, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Between group variance", ylab="Estimate")#, labels=c("ML","Posterior Mean","Posterior Median","Posterior Mode"))
boxplot(estimate~sim_var+type,p_long, pch=19, cex=0.5,xlab="Between group variance", ylab=bquote(P[perm]),boxwex=0.5, names=rep(c(0,0.2),4))#, 
abline(v=c(2.5,4.5,6.5))
mtext(c("ML","Post. Mean","Post. Median","Post. Mode"),side=3,adj=c(0.1,0.35,0.65 ,0.9))

}

dev.off()

apply(p_alt[,c("mean","median","mode","freq","LRT")],2, function(x) mean(x<0.05))
apply(p_null[,c("mean","median","mode","freq","LRT")],2, function(x) sum(x<0.05))








	

{
	par(mfrow=c(2,1), mar=c(5,5,2,1), cex.axis=0.75)
	# boxplot(estimate~N_within+type,bp2_long,boxwex=0.5,xlab="Within group sample size", names=rep(c(2,4,8),3), ylab="Estimate", pch=19, cex=0.5)
	beeswarm(estimate~type,long_alt, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Estimate Type", ylab="Estimate", labels=c("ML","Posterior Mean","Posterior Median","Posterior Mode"))
	means<-aggregate(estimate~type,long_alt,mean)$estimate
	CIs<-aggregate(estimate~type,long_alt,function(x) sd(x)/sqrt(length(x)))$estimate*qnorm(0.975)
	points(means, pch=19, cex=2, col="blue")
	arrows(1:9,means+CIs,1:9,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)
	abline(h=0.2, col="red", lwd=2)
	# abline(v=c(3.5,6.5))
	# mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))
	
	beeswarm(estimate~type,long_null, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Estimate Type", ylab="Estimate", labels=c("ML","Posterior Mean","Posterior Median","Posterior Mode"))
	means<-aggregate(estimate~type,long_null,mean)$estimate
	CIs<-aggregate(estimate~type,long_null,function(x) sd(x)/sqrt(length(x)))$estimate*qnorm(0.975)
	points(means, pch=19, cex=2, col="blue")
	arrows(1:9,means+CIs,1:9,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)
	abline(h=0, col="red", lwd=2)
	# abline(v=c(3.5,6.5))
	# mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))

}

par(mfrow=c(3,1))
	boxplot(estimate~N_within+type,long,boxwex=0.5,xlab="Within group sample size", ylab="Estimate", pch=19, cex=0.5)
	# abline(h=0.2,col="red",lwd=2)
	points(colMeans(actual[,c("mean","median","mode","freq")]), col="blue", pch=19)
colMeans(actual[,c("mean","median","mode","freq")])
apply(actual[,c("mean","median","mode","freq")],2,median)

	boxplot(estimate~type,long_null,boxwex=0.5,xlab="Within group sample size", ylab="Estimate", pch=19, cex=0.5)
	abline(h=0.2,col="red",lwd=2)
	points(colMeans(actual[,c("mean","median","mode","freq")]), col="blue", pch=19)
colMeans(actual[,c("mean","median","mode","freq")])
apply(actual[,c("mean","median","mode","freq")],2,median)




p_perm <- as.data.frame(
	t(sapply(results, function(x) c(x$param,
		freq = p_func(x$actual["freq"],x$null[,"freq"]),
		mean=p_func(x$actual["mean"],x$null[,"mean"]),
		median=p_func(x$actual["median"],x$null[,"median"]),
		mode= p_func(x$actual["mode"],x$null[,"mode"]),
		x$actual["LRT"]
		)))
)

p_long<-Stack(p_perm,c("mean","median","mode","freq","LRT"), group.name="type",value.name = "estimate")

	boxplot(estimate~N_within+type,p_long,boxwex=0.5,xlab="Within group sample size", ylab="p-value", pch=19, cex=0.5)
	# abline(h=c(0.25,0.5,0.75))

apply(p_perm[,c("mean","median","mode","freq","LRT")],2, function(x) mean(x<0.05))

boxplot(sapply(results, function(x) c(x$actual["LRT"])))
