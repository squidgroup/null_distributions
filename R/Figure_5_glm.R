rm(list=ls())

library(beeswarm)
library(scales)
wd <- "~/github/bayes_perm/"
source(paste0(wd,"R/00_functions.R"))

hist(boot::inv.logit(rnorm(100000,0,sqrt(0.2))), breaks=50)

# load(paste0(wd,"Data/Intermediate/GLMM_sim_null.Rdata"))
load( file=paste0(wd,"Data/Intermediate/GLMM_sim.Rdata"))

alt <- as.data.frame(
	t(sapply(results0.2, function(x) c(x$param,x$actual)))
)
head(alt)

p_alt <- as.data.frame(
	t(sapply(results0.2, function(x) c(x$param,
		mean_perm=p_func(x$actual["mean"],x$perm[,"mean"]),
		median_perm=p_func(x$actual["median"],x$perm[,"median"]),
		mode1_perm= p_func(x$actual["mode1"],x$perm[,"mode1"]),
		mode0.1_perm= p_func(x$actual["mode0.1"],x$perm[,"mode0.1"]),
		mean_boot=p_func(x$actual["mean"],x$boot[,"mean"]),
		median_boot=p_func(x$actual["median"],x$boot[,"median"]),
		mode1_boot= p_func(x$actual["mode1"],x$boot[,"mode1"]),
		mode0.1_boot= p_func(x$actual["mode0.1"],x$boot[,"mode0.1"])
		# ,x$actual["LRT"]
		)))
)
head(p_alt)
# load( file=paste0(wd,"Data/Intermediate/GLMM_sim_null.Rdata"))
null <- as.data.frame(
	t(sapply(results0, function(x) c(x$param,x$actual)))
)

p_null <- as.data.frame(
	t(sapply(results0, function(x) c(x$param,
		mean_perm=p_func(x$actual["mean"],x$perm[,"mean"]),
		median_perm=p_func(x$actual["median"],x$perm[,"median"]),
		mode1_perm= p_func(x$actual["mode1"],x$perm[,"mode1"]),
		mode0.1_perm= p_func(x$actual["mode0.1"],x$perm[,"mode0.1"]),
		mean_boot=p_func(x$actual["mean"],x$boot[,"mean"]),
		median_boot=p_func(x$actual["median"],x$boot[,"median"]),
		mode1_boot= p_func(x$actual["mode1"],x$boot[,"mode1"]),
		mode0.1_boot= p_func(x$actual["mode0.1"],x$boot[,"mode0.1"])
		# ,x$actual["LRT"]
		)))
)

long_alt <- Stack(alt,c("mean","median","mode1","mode0.1"), group.name="type",value.name = "estimate")
long_null <- Stack(null,c("mean","median","mode1","mode0.1"), group.name="type",value.name = "estimate")

p_long_alt_boot <- Stack(p_alt,c("mean_boot","median_boot","mode1_boot","mode0.1_boot"), group.name="type",value.name = "estimate")#,"LRT"
p_long_null_boot <- Stack(p_null,c("mean_boot","median_boot","mode1_boot","mode0.1_boot"), group.name="type",value.name = "estimate")#,"LRT"
p_long_alt_perm <- Stack(p_alt,c("mean_perm","median_perm","mode1_perm","mode0.1_perm"), group.name="type",value.name = "estimate")#,"LRT"
p_long_null_perm <- Stack(p_null,c("mean_perm","median_perm","mode1_perm","mode0.1_perm"), group.name="type",value.name = "estimate")#,"LRT"


long <- rbind(long_alt,long_null)
p_long_boot <- rbind(p_long_alt_boot,p_long_null_boot)
p_long_perm <- rbind(p_long_alt_perm,p_long_null_perm)



setEPS()
pdf(paste0(wd,"Figures/Fig5_glmm.pdf"), height=6, width=8)

{
par(mfrow=c(2,1), mar=c(1,4,3,1), cex.axis=0.75, mgp=c(2,0.5,0))
# boxplot(estimate~N_within+type,bp2_long,boxwex=0.5,xlab="Within group sample size", names=rep(c(2,4,8),3), ylab="Estimate", pch=19, cex=0.5)
beeswarm(estimate~sim_var+type,long, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="", ylab="Estimate", labels=rep(c(0,0.2),4), xaxt="n")#, 
mtext(c("Post. Mean","Post. Median","Post. Mode 1","Post. Mode 0.1"),side=3,adj=c(0.1,0.35,0.65 ,0.9))
means<-aggregate(estimate~sim_var+type,long,mean)$estimate
CIs<-aggregate(estimate~sim_var+type,long,function(x) sd(x)/sqrt(length(x)))$estimate*qnorm(0.975)
points(means, pch=19, cex=1.5, col="blue")
arrows(1:9,means+CIs,1:9,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)
arrows(c(1,3,5,7)-0.4,0,c(1,3,5,7)+0.4,0,code=0, angle=90, length=0.05, col="red", lwd=2)
arrows(c(2,4,6,8)-0.4,0.2,c(2,4,6,8)+0.4,0.2,code=0, angle=90, length=0.05, col="red", lwd=2)
abline(v=c(2.5,4.5,6.5))
# mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))
mtext("a)",2,padj=-8, las=1, line=2)

par(mar=c(3,4,1,1))
beeswarm(estimate~sim_var+type,p_long_boot, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Among group variance", ylab="p-value", labels=rep(c(0,0.2),4))#, labels=c("ML","Posterior Mean","Posterior Median","Posterior Mode"))
# boxplot(estimate~sim_var+type,p_long, pch=19, cex=0.5,xlab="Between group variance", ylab=bquote(P[perm]),boxwex=0.5, names=rep(c(0,0.2),4))#, 
abline(v=c(2.5,4.5,6.5))
mtext("b)",2,padj=-8, las=1, line=2)

}

dev.off()
