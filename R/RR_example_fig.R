
load("~/github/null_distributions/Data/Intermediate/RR_example_results.RData")
library(beeswarm)
library(scales)

pdf("~/github/null_distributions/Figures/Fig7_RRexample.pdf", height=6, width=8)
{
	par(mfrow=c(2,1),mar=c(1,4,3.5,1), cex.axis=1, mgp=c(2,0.5,0))
beeswarm(est~type, res, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="", ylab="Variance estimate",  ylim=c(0,0.15), xaxt="n")
mtext("a)",2,padj=-8, las=1, line=2)
abline(h=obs_est_sim[4], col="red")

par(mar=c(3.5,4,1,1))
beeswarm(est~type, res_GT, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap", ylab="Variance estimate",  ylim=c(0,0.15), xlab="", xaxt="n")#
mtext("b)",2,padj=-8, las=1, line=2)
abline(h=obs_est_real[4], col="red")
 axis(1, c(1,2,3,4,5,6), c("Both","ID", "ID", "x", "x within ID", "y"))
  # axis(1, c(1.5,4.5), c("Simulation","Permutation"),line=1.3, tick=FALSE)
axis(1, c(1,1.5,2), c("","Simulation",""),lwd.ticks=0,line=1.7)
axis(1, c(3,4.5,6), c("","Permutation",""),lwd.ticks=0,line=1.7)

  
}
dev.off()
