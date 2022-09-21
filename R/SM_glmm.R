rm(list=ls())

library(beeswarm)
library(scales)
wd <- "~/github/bayes_perm/"
source(paste0(wd,"R/00_functions.R"))

# hist(boot::inv.logit(rnorm(100000,0,sqrt(0.2))), breaks=50)

load( file=paste0(wd,"Data/Intermediate/GLMM_sim.Rdata"))

alt <- as.data.frame(
	t(sapply(results0.2, function(x) c(x$param,x$actual)))
)

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
		)))
)


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
		)))
)


p <- rbind(p_null,p_alt)

pch=19
cex=0.5
col=alpha(1,0.7)


setEPS()
pdf(paste0(wd,"Figures/FigSM_bp_glmm.pdf"), height=11, width=11)
{


par(mfrow=c(2,2), cex.lab=1.5)
plot(mean_perm~mean_boot,p, pch=pch, cex=cex, col=col, main="Posterior Mean", xlab="Parametric Bootstrap", ylab="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(p[,"mean_boot"],p[,"mean_perm"]),3)), cex=1.5)

plot(median_perm~median_boot,p, pch=pch, cex=cex, col=col, main="Posterior Median", xlab="Parametric Bootstrap", ylab="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(p[,"median_boot"],p[,"median_perm"]),3)), cex=1.5)

plot(mode1_perm~mode1_boot,p, pch=pch, cex=cex, col=col, main="Posterior Mode (adjust=1)", xlab="Parametric Bootstrap", ylab="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(p[,"mode1_boot"],p[,"mode1_perm"]),3)), cex=1.5)

plot(mode0.1_perm~mode0.1_boot,p, pch=pch, cex=cex, col=col, main="Posterior Mode (adjust=0.1)", xlab="Parametric Bootstrap", ylab="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(p[,"mode0.1_boot"],p[,"mode0.1_perm"]),3)), cex=1.5)
}

dev.off()



setEPS()
pdf(paste0(wd,"Figures/FigSM_metric_glmm.pdf"), height=11, width=11)
{
par(mfrow=c(2,2), cex.lab=1.5)
plot(mean_boot~median_boot,p, pch=pch, cex=cex, col=col, xlab="Median",ylab="Mean")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(p[,"mean_boot"],p[,"median_boot"]),3)), cex=1.5)


plot(mean_boot~mode1_boot,p, pch=pch, cex=cex, col=col, xlab="Mode-1",ylab="Mean")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(p[,"mean_boot"],p[,"mode1_boot"]),3)), cex=1.5)



plot(mean_boot~mode0.1_boot,p, pch=pch, cex=cex, col=col, xlab="Mode-0.1",ylab="Mean")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(p[,"mean_boot"],p[,"mode0.1_boot"]),3)), cex=1.5)



plot(mode1_boot~mode0.1_boot,p, pch=pch, cex=cex, col=col, xlab="Mode-0.1",ylab="Mode-1")#col=alpha(1:4,0.3)[as.factor(p$N_group)])
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(p[,"mode1_boot"],p[,"mode0.1_boot"]),3)), cex=1.5)


}
dev.off()
