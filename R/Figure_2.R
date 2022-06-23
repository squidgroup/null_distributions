
rm(list=ls())

library("beeswarm")

library(scales)
wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))

files <- list.files(paste0(wd,"Data/Intermediate"))
results <- files[grep("gaus_sims",files)]
list_names <- gsub("gaus_sims_|.Rdata","",results)

## get all results into a list
# all <- list()
# for(i in 1:length(results)){
# 	load(paste0(wd,"Data/Intermediate/",results[i]))
# 	all[[list_names[i]]] <- out
# 	rm("out")
# }

actual<-as.data.frame(do.call(rbind,lapply(results, function(y){
	load(paste0(wd,"Data/Intermediate/",y))
	do.call(rbind,lapply(out, function(x){
	  z <- c(ICC=x$ICC,N_group=x$N_group, N_within=x$N_within,x$summary)
	  names(z) <- c("ICC","N_group","N_within","freq","mode0.1","mode1","median","mean","LCI","UCI","ESS")	
	  z
	}))
} )))

nrow(actual)

actual$mean_bias<-actual$mean-actual$ICC
actual$median_bias<-actual$median-actual$ICC
actual$mode1_bias<-actual$mode1-actual$ICC
actual$mode0.1_bias<-actual$mode0.1-actual$ICC

bias<-aggregate(cbind(mode0.1_bias,mode1_bias,mean_bias,median_bias)~N_group+N_within+ICC,actual,mean)
abs_bias<-aggregate(cbind(mode0.1_bias,mode1_bias,mean_bias,median_bias)~N_group+N_within+ICC,actual,function(x) mean(abs(x)))
accuracy<-aggregate(cbind(mode0.1_bias,mode1_bias,mean_bias,median_bias)~N_group+N_within+ICC,actual,sd)
precision<-aggregate(cbind(mode0.1,mode1,mean,median)~N_group+N_within+ICC,actual,function(x)1/sd(x))

plot(bias$mean_bias-bias$median_bias);abline(h=0)
plot(accuracy$mean_bias-accuracy$median_bias);abline(h=0)


setEPS()
pdf(paste0(wd,"Figures/Fig_bias.pdf"), height=8, width=10)
{
	par(mfrow=c(4,1), mar=c(2,4,3,1))

plot(bias$mean_bias, pch=19, ylim=c(-0.2,0.2), ylab="Bias", xaxt="n")
abline(h=0)
points(bias$median_bias, pch=19,col="red")
points(bias$mode1_bias, pch=19,col="blue")
points(bias$mode0.1_bias, pch=19,col="green")
	abline(v=(1:5)*4+0.5)
axis(1,1:24,rep(c(20,40,80,160),6))
mtext(c("","0","","0.2","","0.4",""),side=3,adj=c(seq(0,1,length.out=7)),line=1)
mtext(c("","2","4","2","4","2","4"),side=3,adj=c(seq(0,1,length.out=7)-1/12))
legend("bottomleft",c("mean","median","mode-0.1","mode-1"), pch=19, col=c(1,2,3,4))

plot(abs_bias$mean_bias, pch=19, ylim=c(0,0.4), ylab="Absolute Bias", xaxt="n")
abline(h=0)
points(abs_bias$median_bias, pch=19,col="red")
points(abs_bias$mode1_bias, pch=19,col="blue")
points(abs_bias$mode0.1_bias, pch=19,col="green")
	abline(v=(1:5)*4+0.5)
axis(1,1:24,rep(c(20,40,80,160),6))

plot(accuracy$mean_bias, pch=19, ylim=c(0,0.4), ylab="RMSE", xaxt="n")
abline(h=0)
points(accuracy$median_bias, pch=19,col="red")
points(accuracy$mode1_bias, pch=19,col="blue")
points(accuracy$mode0.1_bias, pch=19,col="green")
	abline(v=(1:5)*4+0.5)
axis(1:24,rep(c(20,40,80,160),4))
axis(1,1:24,rep(c(20,40,80,160),6))

plot(precision$mean, pch=19, ylim=c(0,100), ylab="Precision", xaxt="n")
abline(h=0)
points(precision$median, pch=19,col="red")
points(precision$mode1, pch=19,col="blue")
points(precision$mode0.1, pch=19,col="green")
	abline(v=(1:5)*4+0.5)
axis(1,1:24,rep(c(20,40,80,160),6))

}
dev.off()

aggregate(cbind(mean,median)~N_group+N_within+ICC,actual,mean)



###
# -- FIGURE 2
###


actual <- as.data.frame(do.call(rbind,lapply(all, function(z){
	t(sapply(z, function(x) c(x$param,x$actual[1:6])))
})))
str(actual)
boxplot(mean~N_within+N_group+ICC,actual,boxwex=0.5, pch=19, cex=0.5)
head(actual)

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
pdf(paste0(wd,"Figures/Fig2.pdf"), height=8, width=10)
{
	bias_plot0 <- subset(actual, ICC ==0 & N_within==2)[,1:8]
	bp0_long<-Stack(bias_plot0, col2stack = c("freq","mean","median","mode0.1","mode1"),value.name = "estimate",group.name = "type")

	bias_plot2 <- subset(actual, ICC==0.2 & N_within==2)[,1:8]
	bp2_long<-Stack(bias_plot2, col2stack = c("freq","mean","median","mode0.1","mode1"),value.name = "estimate",group.name = "type")

	par(mfrow=c(2,1), mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))


	# boxplot(estimate~N_group+type,bp0_long,boxwex=0.5,xlab="Number of Groups", ylab="Estimate", pch=19, cex=0.5, names=rep(c(20,40,80,160),3),  ylim=c(0,1))#
	beeswarm(estimate~N_group+type,bp0_long, pch=19, cex=0.2, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Number of Groups", labels=rep(c(20,40,80,160),3), ylab="Estimate",  ylim=c(0,1))
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


bp0_long$bias<-abs(bp0_long$estimate - bp0_long$ICC)
bp2_long$bias<-abs(bp2_long$estimate - bp2_long$ICC)

{
	par(mfrow=c(2,1), mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))

	# boxplot(estimate~N_group+type,bp0_long,boxwex=0.5,xlab="Number of Groups", ylab="Estimate", pch=19, cex=0.5, names=rep(c(20,40,80,160),3),  ylim=c(0,1))#
	beeswarm(bias~N_group+type,bp0_long, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Number of Groups", labels=rep(c(20,40,80,160),3), ylab="bias",  ylim=c(0,1))
	abline(h=0, col="red", lwd=2)
	abline(v=c(4.5, 8.5, 12.5))
	mtext(c("ML","Post. Mean","Post. Median","Post. Mode"),side=3,adj=c(0.12,0.375,0.625 ,0.9))
	mtext("a)",2,padj=-12, las=1, line=2)
	means<-aggregate(bias~N_group+type,bp0_long,mean)$bias
	CIs<-aggregate(bias~N_group+type,bp0_long,function(x) sd(x)/sqrt(length(x)))$bias*qnorm(0.975)
	points(means, pch=19, cex=1.5, col="blue")
	arrows(1:16,means+CIs,1:16,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)

	# boxplot(bias~N_group+type,bp2_long,boxwex=0.5,xlab="Number of Groups",  ylab="bias", pch=19, cex=0.5, names=rep(c(20,40,80,160),3),  ylim=c(0,1))
	beeswarm(bias~N_group+type,bp2_long, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Number of Groups", labels=rep(c(20,40,80,160),3), ylab="bias",  ylim=c(0,1))
	abline(h=0, col="red", lwd=2)
	abline(v=c(4.5, 8.5, 12.5))
	mtext(c("ML","Post. Mean","Post. Median","Post. Mode"),side=3,adj=c(0.12,0.375,0.625 ,0.9))
	mtext("b)",2,padj=-12, las=1, line=2)
	means<-aggregate(bias~N_group+type,bp2_long,mean)$bias
	CIs<-aggregate(bias~N_group+type,bp2_long,function(x) sd(x)/sqrt(length(x)))$bias*qnorm(0.975)
	points(means, pch=19, cex=1.5, col="blue")
	arrows(1:16,means+CIs,1:16,means-CIs,code=0, angle=90, length=0.05, col="blue", lwd=2)
}

