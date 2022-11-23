
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
par(mfrow=c(4,1), mar=c(0,4,3,0), cex.lab=1.4)
beeswarm(mean~N_group+N_within+ICC,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xaxt="n", ylab="Posterior Mean",  ylim=c(0,1))
# abline(h=0)
points(agg$mean, pch=19, col="orange")
# points(agg$ICC, pch="-", col="red")
arrows(1:24 -0.2,agg$ICC,1:24 +0.2,agg$ICC,code=0,col="red")


abline(v=line_coords, lty=line_lty)

axis(3,c(1,3.5,6),c("",0,""), lwd.ticks=0, line=1, padj=1, cex.axis=1.25)
axis(3,c(1,3.5,6) + 6,c("",0.1,""), lwd.ticks=0, line=1, padj=1, cex.axis=1.25)
axis(3,c(1,3.5,6) + 12,c("",0.2,""), lwd.ticks=0, line=1, padj=1, cex.axis=1.25)
axis(3,c(1,3.5,6) + 18,c("",0.4,""), lwd.ticks=0, line=1, padj=1, cex.axis=1.25)
mtext("ICC", side=3, line=-2, outer=TRUE, adj=0.05)


par( mar=c(1,4,2,0))

beeswarm(median~N_group+N_within+ICC,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xaxt="n", ylab="Posterior Median",  ylim=c(0,1))
# abline(h=0)
points(agg$median, pch=19, col="orange")
# points(agg$ICC, pch="-", col="red")
arrows(1:24 -0.2,agg$ICC,1:24 +0.2,agg$ICC,code=0,col="red")
abline(v=line_coords, lty=line_lty)

par( mar=c(2,4,1,0))

beeswarm(mode1~N_group+N_within+ICC,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xaxt="n", ylab="Posterior Mode - 1",  ylim=c(0,1))
# abline(h=0)
points(agg$mode1, pch=19, col="orange")
# points(agg$ICC, pch="-", col="red")
arrows(1:24 -0.2,agg$ICC,1:24 +0.2,agg$ICC,code=0,col="red")
abline(v=line_coords, lty=line_lty)

par( mar=c(3,4,0,0))

beeswarm(mode0.1~N_group+N_within+ICC,actual, pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="", ylab="Posterior Mode - 0.1",  ylim=c(0,1),xaxt="n")
# abline(h=0)
points(agg$mode0.1, pch=19, col="orange")
# points(agg$ICC, pch="-", col="red")
arrows(1:24 -0.2,agg$ICC,1:24 +0.2,agg$ICC,code=0,col="red")

abline(v=line_coords, lty=line_lty)

axis(1,1:24,rep(c(20,40,80),8), cex.axis=1)
axis(1,1:8 *3 -1,rep(c(2,4),4), tick=FALSE, line=1, cex.axis=1)

mtext("N among", side=1, line=-2, outer=TRUE, adj=0, cex=0.9)
mtext("N within", side=1, line=-1, outer=TRUE, adj=0, cex=0.9)

}
dev.off()


