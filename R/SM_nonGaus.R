rm(list=ls())

library(scales)
library(beeswarm)

wd <- "~/github/null_distributions/"

source(paste0(wd,"R/00_functions.R"))

convert_func <- function(x){
	param <- if(length(x$param)==4) x$param[2:4] else x$param
	# param <-  x$param[2:4]
	names(param) <-c("ICC", "N_group", "N_within")
	c(param,
		median = as.numeric(x$actual["median"]),
		mean = as.numeric(x$actual["mean"]),
		mode1 = as.numeric(x$actual["mode1"]),
		mode0.1 = as.numeric(x$actual["mode0.1"]),
		median_p = p_func(x$actual["median"],x$boot[,"median"]),
		mean_p = p_func(x$actual["mean"],x$boot[,"mean"]),
		mode1_p = p_func(x$actual["mode1"],x$boot[,"mode1"]),
		mode0.1_p = p_func(x$actual["mode0.1"],x$boot[,"mode0.1"])
)}

# files_p <- list.files(paste0(wd,"Data/Intermediate"))
# results_p <- files_p[grep("gaus_perm",files_p)]
# list_names_p <- gsub("gaus_perm_|.Rdata","",results_p)

# p_perm<-as.data.frame(do.call(rbind,lapply(results_p, function(i){
#   load(paste0(wd,"Data/Intermediate/",i))
#   t(sapply(perm_out, convert_func))
# })))

# head(p_perm)

files <- list.files(paste0(wd,"Data/Intermediate"))

results_bern <- files[grep("bern_boot",files)]
p_bern<-as.data.frame(do.call(rbind,lapply(results_bern, function(i){
  load(paste0(wd,"Data/Intermediate/",i))
  # print(length(boot_out))
  t(sapply(boot_out, convert_func))
})))

results_pois <- files[grep("pois_boot",files)]
p_pois<-as.data.frame(do.call(rbind,lapply(results_pois, function(i){
	## adjust ICC by total latent var
  load(paste0(wd,"Data/Intermediate/",i))
    boot_out<-lapply(boot_out, function(x) {
  	x$param["ICC"] <- x$param["ICC"]*0.2
		x
	})
  t(sapply(boot_out, convert_func))
})))

long_bern <- Stack(p_bern,c("mean","median","mode1","mode0.1"), group.name="type",value.name = "estimate")
long_pois <- Stack(p_pois,c("mean","median","mode1","mode0.1"), group.name="type",value.name = "estimate")


setEPS()
pdf(paste0(wd,"Figures/FigSM_nonGaus_samp_dist.pdf"), height=10, width=10)
{
par(mfrow=c(2,1), mar=c(4,4,4,1), cex.axis=1.5,cex.axis=1, mgp=c(2,0.5,0))
beeswarm(estimate~ICC+type,long_bern, pch=19, cex=0.1, col=alpha(1,0.2),method = "compactswarm",corral="wrap",xlab="", ylab="Estimate", labels=c(0,0.2,0.4,0.8), main = "Bernoulli")#, 
mtext(c("Mean","Median","Mode 1","Mode 0.1"),side=3,adj=c(0.1,0.35,0.65 ,0.9))
means<-aggregate(estimate~ICC+type,long_bern,mean)$estimate
CIs<-aggregate(estimate~ICC+type,long_bern,function(x) sd(x)/sqrt(length(x)))$estimate*qnorm(0.975)
points(means, pch=19, cex=1, col="blue")
arrows(1:16,means+CIs,1:16,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)
arrows(1:16-0.4,c(0,0.2,0.4,0.8),1:16+0.4,c(0,0.2,0.4,0.8),code=0, angle=90, length=0.05, col="red", lwd=2)
abline(v=c(4.5,8.5,12.5))
mtext("a)",2,padj=-14, las=1, line=2.5)

beeswarm(estimate~ICC+type,long_pois, pch=19, cex=0.1, col=alpha(1,0.2),method = "compactswarm",corral="wrap",xlab="", ylab="Estimate", labels=c(0,0.02,0.04,0.08), main="Poisson")#, 
mtext(c("Mean","Median","Mode 1","Mode 0.1"),side=3,adj=c(0.1,0.35,0.65 ,0.9))
means<-aggregate(estimate~ICC+type,long_pois,mean)$estimate
CIs<-aggregate(estimate~ICC+type,long_pois,function(x) sd(x)/sqrt(length(x)))$estimate*qnorm(0.975)
points(means, pch=19, cex=1, col="blue")
arrows(1:16,means+CIs,1:16,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)
arrows(1:16-0.4,c(0,0.02,0.04,0.08),1:16+0.4,c(0,0.02,0.04,0.08),code=0, angle=90, length=0.05, col="red", lwd=2)
abline(v=c(4.5,8.5,12.5))
mtext("b)",2,padj=-14, las=1, line=2.5)

}
dev.off()


setEPS()
pdf(paste0(wd,"Figures/FigSM_nonGaus_p_dist.pdf"), height=8, width=6)
{
par(mfrow=c(2,1), mar=c(4,4,4,1), cex.axis=1.5,cex.axis=1, mgp=c(2,0.5,0))
beeswarm(median_p~ICC,p_bern, pch=19, cex=0.1, col=alpha(1,0.2),method = "compactswarm",corral="wrap",xlab="Simulated variance on latent scale", ylab="P-value", labels=c(0,0.2,0.4,0.8), main = "Bernoulli")
mtext("a)",2,padj=-11, las=1, line=2.5)

beeswarm(median_p~ICC,p_pois, pch=19, cex=0.1, col=alpha(1,0.2),method = "compactswarm",corral="wrap",xlab="Simulated variance on latent scale", ylab="P-value", labels=c(0,0.02,0.04,0.08), main="Poisson")
mtext("b)",2,padj=-11, las=1, line=2.5)
}
dev.off()


##### correlation of p's

p_cor_plot <- function(m1,m2,data,...){
	plot(data[,paste0(m1,"_p")],data[,paste0(m2,"_p")], pch=19, cex=0.2, col=scales::alpha(1,0.3), xlab=m1,ylab=m2,...)
	abline(0,1,col="red")
	text(0.1,0.9,paste0("r=",round(cor(data[,paste0(m1,"_p")],data[,paste0(m2,"_p")]),3)), cex=1.5)
}

setEPS()
pdf(paste0(wd,"Figures/FigSM_nonGaus_metric_comp.pdf"), height=15, width=9)
{
par(mfcol=c(4,2), mar=c(3,4,2,1), cex.lab=1.5, cex.axis=1, mgp=c(2,0.5,0))
p_cor_plot("median","mean",p_bern,main="Bernoulli")
p_cor_plot("mode1","mean",p_bern)
p_cor_plot("mode0.1","mean",p_bern)
p_cor_plot("mode0.1","mode1",p_bern)

p_cor_plot("median","mean",p_pois,main="Poisson")
p_cor_plot("mode1","mean",p_pois)
p_cor_plot("mode0.1","mean",p_pois)
p_cor_plot("mode0.1","mode1",p_pois)
}
dev.off()
