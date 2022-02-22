
rm(list=ls())

library("beeswarm")

library(scales)
wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/functions.R"))

files <- list.files(paste0(wd,"Data/Intermediate"))
results <- files[grep("sims",files)]
list_names <- gsub("sims_|.Rdata","",results)

## get all results into a list
all <- list()
for(i in 1:length(results)){
	load(paste0(wd,"Data/Intermediate/",results[i]))
	all[[list_names[i]]] <- sim_dat
	rm("sim_dat")
}



###
# -- FIGURE 2
###


actual <- as.data.frame(do.call(rbind,lapply(all, function(z){
	t(sapply(z, function(x) c(x$param,x$actual[1:6])))
})))

boxplot(mean~N_within+N_group+ICC,actual,boxwex=0.5, pch=19, cex=0.5)


{
	bias_plot0 <- subset(actual, ICC ==0 & N_group==40)[,1:6]
	bias_plot2 <- subset(actual, ICC ==0.2 & N_group==40)[,1:6]

	bp2_long<-Stack(bias_plot2, col2stack = c("freq","mean","mode"),value.name = "estimate",group.name = "type")  
	bp0_long<-Stack(bias_plot0, col2stack = c("freq","mean","mode"),value.name = "estimate",group.name = "type")

	par(mfrow=c(2,1), mar=c(5,5,2,1), cex.axis=0.75)
	# boxplot(estimate~N_within+type,bp2_long,boxwex=0.5,xlab="Within group sample size", names=rep(c(2,4,8),3), ylab="Estimate", pch=19, cex=0.5)
	beeswarm(estimate~N_within+type,bp2_long, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Within group sample size", labels=rep(c(2,4,8),3), ylab="Estimate")
	means<-aggregate(estimate~N_within+type,bp2_long,mean)$estimate
	CIs<-aggregate(estimate~N_within+type,bp2_long,function(x) sd(x)/sqrt(length(x)))$estimate*qnorm(0.975)
	points(means, pch=19, cex=2, col="blue")
	arrows(1:9,means+CIs,1:9,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)
	abline(h=0.2, col="red", lwd=2)
	abline(v=c(3.5,6.5))
	mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))

	# boxplot(estimate~N_within+type,bp0_long,boxwex=0.5,xlab="Within group sample size", ylab="Estimate", pch=19, cex=0.5, names=rep(c(2,4,8),3), ylim=c(0,0.3))#
	beeswarm(estimate~N_within+type,bp0_long, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Within group sample size", labels=rep(c(2,4,8),3), ylab="Estimate")
	means<-aggregate(estimate~N_within+type,bp0_long,mean)$estimate
	CIs<-aggregate(estimate~N_within+type,bp0_long,function(x) sd(x)/sqrt(length(x)))$estimate*qnorm(0.975)
	points(means, pch=19, cex=2, col="blue")
	arrows(1:9,means+CIs,1:9,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)
	abline(h=0, col="red", lwd=2)
	abline(v=c(3.5,6.5))
	mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))


}

setEPS()
pdf(paste0(wd,"Figures/Fig2.pdf"), height=8, width=8)
{
	bias_plot0 <- subset(actual, ICC ==0 & N_within==2)[,1:7]
	bp0_long<-Stack(bias_plot0, col2stack = c("freq","mean","median","mode"),value.name = "estimate",group.name = "type")

	bias_plot2 <- subset(actual, ICC==0.2 & N_within==2)[,1:7]
	bp2_long<-Stack(bias_plot2, col2stack = c("freq","mean","median","mode"),value.name = "estimate",group.name = "type")

	par(mfrow=c(2,1), mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))


	# boxplot(estimate~N_group+type,bp0_long,boxwex=0.5,xlab="Number of Groups", ylab="Estimate", pch=19, cex=0.5, names=rep(c(20,40,80,160),3),  ylim=c(0,1))#
	beeswarm(estimate~N_group+type,bp0_long, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Number of Groups", labels=rep(c(20,40,80,160),3), ylab="Estimate",  ylim=c(0,1))
	abline(h=0, col="red", lwd=2)
	abline(v=c(4.5, 8.5, 12.5))
	mtext(c("ML","Post. Mean","Post. Median","Post. Mode"),side=3,adj=c(0.12,0.375,0.625 ,0.9))
	mtext("a)",2,padj=-12, las=1, line=2)
	means<-aggregate(estimate~N_group+type,bp0_long,mean)$estimate
	CIs<-aggregate(estimate~N_group+type,bp0_long,function(x) sd(x)/sqrt(length(x)))$estimate*qnorm(0.975)
	points(means, pch=19, cex=1.5, col="blue")
	arrows(1:16,means+CIs,1:16,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)

	# boxplot(estimate~N_group+type,bp2_long,boxwex=0.5,xlab="Number of Groups",  ylab="Estimate", pch=19, cex=0.5, names=rep(c(20,40,80,160),3),  ylim=c(0,1))
	beeswarm(estimate~N_group+type,bp2_long, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Number of Groups", labels=rep(c(20,40,80,160),3), ylab="Estimate",  ylim=c(0,1))
	abline(h=0.2, col="red", lwd=2)
	abline(v=c(4.5, 8.5, 12.5))
	mtext(c("ML","Post. Mean","Post. Median","Post. Mode"),side=3,adj=c(0.12,0.375,0.625 ,0.9))
	mtext("b)",2,padj=-12, las=1, line=2)
	means<-aggregate(estimate~N_group+type,bp2_long,mean)$estimate
	CIs<-aggregate(estimate~N_group+type,bp2_long,function(x) sd(x)/sqrt(length(x)))$estimate*qnorm(0.975)
	points(means, pch=19, cex=1.5, col="blue")
	arrows(1:16,means+CIs,1:16,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)
}
dev.off()



## check ESS
 boxplot(ESS~N_group+N_within+ICC,actual, las=2)
 abline(h=500)

