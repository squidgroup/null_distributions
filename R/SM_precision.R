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

precision<-aggregate(cbind(mode0.1,mode1,mean,median)~N_group+N_within+ICC,actual,function(x)1/sd(x))

line_coords <- (1:7)*3+0.5

setEPS()
pdf(paste0(wd,"Figures/FigSM_precision.pdf"), height=5, width=11)
{

# line_coords <- (1:7)*3+0.5
par(mfrow=c(1,1), mar=c(4,4,3,1))

plot(precision$mean, pch=20, ylim=c(0,70), ylab="Precision", xlab="", xaxt="n")
abline(h=0)
abline(v=line_coords, lty=c(2,1,2,1,2,1,2))

points(precision$median, pch=17,col="red")
points(precision$mode1, pch=15,col="blue")
points(precision$mode0.1, pch=18,col="green")

axis(3,c(1,3.5,6),c("",0,""), lwd.ticks=0, line=1, padj=1)
axis(3,c(1,3.5,6) + 6,c("",0.1,""), lwd.ticks=0, line=1, padj=1)
axis(3,c(1,3.5,6) + 12,c("",0.2,""), lwd.ticks=0, line=1, padj=1)
axis(3,c(1,3.5,6) + 18,c("",0.4,""), lwd.ticks=0, line=1, padj=1)
mtext("ICC", side=3, line=-2, outer=TRUE, adj=0.05)


axis(1,1:24,rep(c(20,40,80),8))
axis(1,1:8 *3 -1,rep(c(2,4),4), tick=FALSE, line=1.5)

mtext("N among", side=1, line=-3, outer=TRUE, adj=0)
mtext("N within", side=1, line=-1.5, outer=TRUE, adj=0)
 legend("topleft",c("mean","median","mode-1","mode-0.1"), pch=c(20,17,15,18), col=c(1,2,4,3), bty="n")

}

dev.off()
