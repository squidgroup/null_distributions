
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

agg <- aggregate(cbind(mean,median,mode1,mode0.1)~N_group+N_within+ICC,actual,mean)

line_coords <- (1:7)*3+0.5
line_lty <- c(2,1,2,1,2,1,2)

setEPS()
pdf(paste0(wd,"Figures/FigSM_samp_dist.pdf"), height=8, width=12)
{
par(mfrow=c(4,1), mar=c(0,4,4.5,0))
beeswarm(mean~N_group+N_within+ICC,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",labels=c(20,40,80), ylab="Posterior Mean",  ylim=c(0,1))
# abline(h=0)
points(agg$mean, pch=19, col="orange")
# points(agg$ICC, pch="-", col="red")
arrows(1:24 -0.2,agg$ICC,1:24 +0.2,agg$ICC,code=0,col="red")


abline(v=line_coords, lty=line_lty)
# axis(1,1:24,rep(c(0,0.1,0.2,0.4),6))
# axis(3,(1:6) *4 -1.5,rep(c(20,40,80),2), tick=FALSE, line=-0.5)
# axis(3,(1:3) *4 -1.5,c("",2,""), lwd.ticks=0, line=2, padj=1)
# axis(3,(4:6) *4 -1.5,c("",4,""), lwd.ticks=0, line=2, padj=1)

axis(3,1:8 *3 -1,rep(c(2,4),4), tick=FALSE, line=-0.5, cex.axis=1.25)

axis(3,c(2,3.5,5),c("",0,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)
axis(3,c(2,3.5,5) + 6,c("",0.1,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)
axis(3,c(2,3.5,5) + 12,c("",0.2,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)
axis(3,c(2,3.5,5) + 18,c("",0.4,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)



par( mar=c(1.5,4,3,0))

beeswarm(median~N_group+N_within+ICC,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",labels=c(20,40,80), ylab="Posterior Median",  ylim=c(0,1))
# abline(h=0)
points(agg$median, pch=19, col="orange")
# points(agg$ICC, pch="-", col="red")
arrows(1:24 -0.2,agg$ICC,1:24 +0.2,agg$ICC,code=0,col="red")
abline(v=line_coords, lty=line_lty)

par( mar=c(3,4,1.5,0))

beeswarm(mode1~N_group+N_within+ICC,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",labels=c(20,40,80), ylab="Posterior Mode (adjust=1)",  ylim=c(0,1))
# abline(h=0)
points(agg$mode1, pch=19, col="orange")
# points(agg$ICC, pch="-", col="red")
arrows(1:24 -0.2,agg$ICC,1:24 +0.2,agg$ICC,code=0,col="red")
abline(v=line_coords, lty=line_lty)

par( mar=c(4.5,4,0,0))

beeswarm(mode0.1~N_group+N_within+ICC,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="ICC", ylab="Posterior Mode (adjust=0.1",  ylim=c(0,1),labels=c(20,40,80))
# abline(h=0)
points(agg$mode0.1, pch=19, col="orange")
# points(agg$ICC, pch="-", col="red")
arrows(1:24 -0.2,agg$ICC,1:24 +0.2,agg$ICC,code=0,col="red")

abline(v=line_coords, lty=line_lty)

}

dev.off()


