rm(list=ls())
wd <- "~/github/null_distributions/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Intermediate/fay_results.Rdata"))

results_sum <- t(sapply(results,colMeans))
 # results_sum <- t(sapply(results,function(x)colMeans(sqrt(x))))
# results_sum <- t(sapply(results,function(x)apply(x,2,median)))

setEPS()
pdf(paste0(wd,"Figures/FigSM_fay.pdf"), height=6, width=6)
{
par(mar=c(4,4,1,1))
plot(results_sum[,"n"],results_sum[,"mean"], ylim=c(0,0.3), type="b", pch=19, col=1, ylab="Mean Estimate ", xlab="Number of Individuals")
points(results_sum[,"n"],results_sum[,"median"], type="b", pch=19, col=2)
points(results_sum[,"n"],results_sum[,"mode0.1"], type="b", pch=19, col=3)
points(results_sum[,"n"],results_sum[,"mode1"], type="b", pch=19, col=4)
abline(h=0.2^2, col="grey")
legend("topright",c("mean","median","mode-0.1","mode-1"), pch=19, col=c(1,2,3,4), bty="n")

}
dev.off()

library(beeswarm)
beeswarm(mode1~n,do.call(rbind,results))
beeswarm(median~n,do.call(rbind,results))

boxplot(ESS~n,do.call(rbind,results))

