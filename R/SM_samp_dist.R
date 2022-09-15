
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
	  z <- c(x$param,x$summary)
	  # ICC=x$ICC,N_group=x$N_group, N_within=x$N_within
	  names(z) <- c("pop","ICC","N_group","N_within","freq","mode0.1","mode1","median","mean","LCI","UCI","ESS")	
	  z
	}))
} )))

nrow(actual)

actual$mean_bias<-actual$mean-actual$ICC
actual$median_bias<-actual$median-actual$ICC
actual$mode1_bias<-actual$mode1-actual$ICC
actual$mode0.1_bias<-actual$mode0.1-actual$ICC

bias<-aggregate(cbind(mode0.1_bias,mode1_bias,mean_bias,median_bias)~ICC+N_group+N_within,actual,mean)
abs_bias<-aggregate(cbind(mode0.1_bias,mode1_bias,mean_bias,median_bias)~ICC+N_group+N_within,actual,function(x) mean(abs(x)))
rmse<-aggregate(cbind(mode0.1_bias,mode1_bias,mean_bias,median_bias)~ICC+N_group+N_within,actual,function(x) sqrt(mean(x^2)) )
# accuracy<-aggregate(cbind(mode0.1_bias,mode1_bias,mean_bias,median_bias)~ICC+N_group+N_within,actual,sd)
precision<-aggregate(cbind(mode0.1,mode1,mean,median)~ICC+N_group+N_within,actual,function(x)1/sd(x))


agg <- aggregate(cbind(mean,median,mode1,mode0.1)~ICC+N_group+N_within,actual,mean)

line_coords <- (1:5)*4+0.5

setEPS()
pdf(paste0(wd,"Figures/FigSM_samp_dist.pdf"), height=8, width=12)
{
par(mfrow=c(4,1), mar=c(0,4,4.5,0))
beeswarm(mean~ICC+N_group+N_within,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",labels=rep(c(0,0.1,0.2,0.4)), ylab="Posterior Mean",  ylim=c(0,1))
# abline(h=0)
points(agg$mean, pch=19, col="blue")
points(agg$ICC, pch=19, col="red")

abline(v=line_coords, lty=c(2,2,1,2,2))
# axis(1,1:24,rep(c(0,0.1,0.2,0.4),6))
axis(3,(1:6) *4 -1.5,rep(c(20,40,80),2), tick=FALSE, line=-0.5)
axis(3,(1:3) *4 -1.5,c("",2,""), lwd.ticks=0, line=2, padj=1)
axis(3,(4:6) *4 -1.5,c("",4,""), lwd.ticks=0, line=2, padj=1)


par( mar=c(1.5,4,3,0))

beeswarm(median~ICC+N_group+N_within,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",labels=rep(c(0,0.1,0.2,0.4)), ylab="Posterior Median",  ylim=c(0,1))
# abline(h=0)
points(agg$median, pch=19, col="blue")
points(agg$ICC, pch=19, col="red")
abline(v=line_coords, lty=c(2,2,1,2,2))

par( mar=c(3,4,1.5,0))

beeswarm(mode1~ICC+N_group+N_within,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",labels=rep(c(0,0.1,0.2,0.4)), ylab="Posterior Mode (adjust=1)",  ylim=c(0,1))
# abline(h=0)
points(agg$mode1, pch=19, col="blue")
points(agg$ICC, pch=19, col="red")
abline(v=line_coords, lty=c(2,2,1,2,2))

par( mar=c(4.5,4,0,0))

beeswarm(mode0.1~ICC+N_group+N_within,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="ICC", ylab="Posterior Mode (adjust=0.1",  ylim=c(0,1),labels=rep(c(0,0.1,0.2,0.4),3))
# abline(h=0)
points(agg$mode0.1, pch=19, col="blue")
points(agg$ICC, pch=19, col="red")
abline(v=line_coords, lty=c(2,2,1,2,2))

}

dev.off()


