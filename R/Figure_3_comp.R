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

actual$mean_rel_bias<-(actual$mean-actual$ICC)/actual$ICC
actual$median_rel_bias<-(actual$median-actual$ICC)/actual$ICC
actual$mode1_rel_bias<-(actual$mode1-actual$ICC)/actual$ICC
actual$mode0.1_rel_bias<-(actual$mode0.1-actual$ICC)/actual$ICC





# bias<-aggregate(cbind(mode0.1_bias,mode1_bias,mean_bias,median_bias)~ICC+N_group+N_within,actual,mean)
# abs_bias<-aggregate(cbind(mode0.1_bias,mode1_bias,mean_bias,median_bias)~ICC+N_group+N_within,actual,function(x) mean(abs(x)))
# relative_bias<-aggregate(cbind(mode0.1_rel_bias,mode1_rel_bias,mean_rel_bias,median_rel_bias)~ICC+N_group+N_within,actual,mean)

# abs_rel_bias<-aggregate(cbind(mode0.1_rel_bias,mode1_rel_bias,mean_rel_bias,median_rel_bias)~ICC+N_group+N_within,actual,function(x) mean(abs(x)))



# setEPS()
# pdf(paste0(wd,"Figures/Fig3_bias.pdf"), height=8, width=11)
# {

# line_coords <- (1:5)*4+0.5
# par(mfrow=c(2,1), mar=c(1,4,4,1))

# plot(bias$mean_bias, pch=19, ylim=c(-0.2,0.2), ylab="Bias", xaxt="n")
# abline(h=0)
# abline(v=line_coords, lty=c(2,2,1,2,2))
# # axis(1,1:24,rep(c(0,0.1,0.2,0.4),6))
# axis(3,(1:6) *4 -1.5,rep(c(20,40,80),2), tick=FALSE, line=-0.5)
# axis(3,(1:3) *4 -1.5,c("",2,""), lwd.ticks=0, line=2, padj=1)
# axis(3,(4:6) *4 -1.5,c("",4,""), lwd.ticks=0, line=2, padj=1)

# points(bias$median_bias, pch=17,col="red")
# points(bias$mode1_bias, pch=19,col="blue")
# points(bias$mode0.1_bias, pch=19,col="green")


# par(mar=c(4,4,1,1))
# plot(abs_bias$mean_bias, pch=19, ylim=c(0,0.4), ylab="Absolute Bias", xaxt="n", xlab="ICC")
# abline(h=0)
# abline(v=line_coords, lty=c(2,2,1,2,2))
# points(abs_bias$median_bias, pch=17,col="red")
# points(abs_bias$mode1_bias, pch=19,col="blue")
# points(abs_bias$mode0.1_bias, pch=19,col="green")
	
#  axis(1,1:24,rep(c(0,0.1,0.2,0.4),6))
#  legend("topright",c("mean","median","mode-0.1","mode-1"), pch=19, col=c(1,2,3,4))

# }

# dev.off()



# setEPS()
# pdf(paste0(wd,"Figures/Fig3_bias.pdf"), height=8.5, width=9)
# {
# pt.cex=1.5
# line_coords <- (1:5)*4+0.5
# par(mfrow=c(3,1), mar=c(0,5,4,2.5), cex.lab=1.5)

# plot(bias$mean_bias, cex=pt.cex,pch=19, ylim=c(-0.2,0.2), ylab="Bias", xaxt="n")
# abline(h=0)
# abline(v=line_coords, lty=c(2,2,1,2,2))
# points(bias$median_bias, cex=pt.cex,pch=17,col="red")
# points(bias$mode1_bias, cex=pt.cex,pch=15,col="blue")
# points(bias$mode0.1_bias, cex=pt.cex,pch=18,col="green")
# axis(4)
# mtext("a)",2,padj=-9.5, las=1, line=3)

# axis(3,(1:6) *4 -1.5,rep(c(20,40,80),2), tick=FALSE, line=-0.5, cex.axis=1.25)
# axis(3,(1:3) *4 -1.5,c("",2,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)
# axis(3,(4:6) *4 -1.5,c("",4,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)


# par(mar=c(2,5,2,2.5))

# plot(relative_bias$mean_rel_bias*100, cex=pt.cex,pch=19, ylim=c(-150,150), ylab="Relative Bias (%)", xaxt="n", xlab="ICC")
# abline(h=0)
# abline(v=line_coords, lty=c(2,2,1,2,2))
# points(relative_bias$median_rel_bias*100, cex=pt.cex,pch=17,col="red")
# points(relative_bias$mode1_rel_bias*100, cex=pt.cex,pch=15,col="blue")
# points(relative_bias$mode0.1_rel_bias*100, cex=pt.cex,pch=18,col="green")
# axis(4)
# mtext("b)",2,padj=-9.5, las=1, line=3)

# legend("topright",c("mean","median","mode-0.1","mode-1"), cex=1.25,pch=19, col=c(1,2,3,4))



# par(mar=c(4,5,0,2.5))

# plot(abs_rel_bias$mean_rel_bias*100, cex=pt.cex,pch=19, ylim=c(0,150), ylab="Absolute Relative Bias (%)", xaxt="n", xlab="ICC")
# abline(h=0)
# abline(v=line_coords, lty=c(2,2,1,2,2))
# points(abs_rel_bias$median_rel_bias*100, cex=pt.cex,pch=17,col="red")
# points(abs_rel_bias$mode1_rel_bias*100, cex=pt.cex,pch=15,col="blue")
# points(abs_rel_bias$mode0.1_rel_bias*100, cex=pt.cex,pch=18,col="green")
# axis(4)
# mtext("c)",2,padj=-9.5, las=1, line=3)

# axis(1,1:24,rep(c(0,0.1,0.2,0.4),6), cex.axis=1.2)


# }

# dev.off()





bias<-aggregate(cbind(mode0.1_bias,mode1_bias,mean_bias,median_bias)~N_group+N_within+ICC,actual,mean)
relative_bias<-aggregate(cbind(mode0.1_rel_bias,mode1_rel_bias,mean_rel_bias,median_rel_bias)~N_group+N_within+ICC,actual,mean)
abs_rel_bias<-aggregate(cbind(mode0.1_rel_bias,mode1_rel_bias,mean_rel_bias,median_rel_bias)~N_group+N_within+ICC,actual,function(x) mean(abs(x)))


setEPS()
pdf(paste0(wd,"Figures/Fig3_bias.pdf"), height=8.5, width=9)
{


pt.cex=1.5
line_coords <- (1:7)*3+0.5
line_lty <- c(2,1,2,1,2,1,2)

par(mfrow=c(3,1), mar=c(0,5,4,2.5), cex.lab=1.5)

plot(bias$mean_bias, cex=pt.cex,pch=19, ylim=c(-0.2,0.2), ylab="Bias", xaxt="n")
abline(h=0)
abline(v=line_coords, lty=line_lty)
points(bias$median_bias, cex=pt.cex,pch=17,col="red")
points(bias$mode1_bias, cex=pt.cex,pch=15,col="blue")
points(bias$mode0.1_bias, cex=pt.cex,pch=18,col="green")
axis(4)
mtext("a)",2,padj=-9.5, las=1, line=3)

# axis(3,(1:6) *4 -1.5,rep(c(20,40,80),2), tick=FALSE, line=-0.5, cex.axis=1.25)
# axis(3,(1:3) *4 -1.5,c("",2,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)
# axis(3,(4:6) *4 -1.5,c("",4,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)

axis(3,1:8 *3 -1,rep(c(2,4),4), tick=FALSE, line=-0.5, cex.axis=1.25)

axis(3,c(2,3.5,5),c("",0,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)
axis(3,c(2,3.5,5) + 6,c("",0.1,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)
axis(3,c(2,3.5,5) + 12,c("",0.2,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)
axis(3,c(2,3.5,5) + 18,c("",0.4,""), lwd.ticks=0, line=2, padj=1, cex.axis=1.25)

# mtext(c("","0","","0.1","","0.2","","0.4",""),side=3,adj=c(seq(0,1,length.out=9)),line=1)
# mtext(c("","2","4","2","4","2","4","2","4"),side=3,adj=c(seq(0,1,length.out=9)-1/16))


par(mar=c(2,5,2,2.5))

plot(relative_bias$mean_rel_bias*100, cex=pt.cex,pch=19, ylim=c(-150,150), ylab="Relative Bias (%)", xaxt="n", xlab="ICC")
abline(h=0)
abline(v=line_coords, lty=line_lty)
points(relative_bias$median_rel_bias*100, cex=pt.cex,pch=17,col="red")
points(relative_bias$mode1_rel_bias*100, cex=pt.cex,pch=15,col="blue")
points(relative_bias$mode0.1_rel_bias*100, cex=pt.cex,pch=18,col="green")
axis(4)
mtext("b)",2,padj=-9.5, las=1, line=3)

legend("topright",c("mean","median","mode-0.1","mode-1"), cex=1.25,pch=c(19,17,15,18), col=c(1,2,3,4))



par(mar=c(4,5,0,2.5))

plot(abs_rel_bias$mean_rel_bias*100, cex=pt.cex,pch=19, ylim=c(0,150), ylab="Absolute Relative Bias (%)", xaxt="n", xlab="ICC")
abline(h=0)
abline(v=line_coords, lty=line_lty)
points(abs_rel_bias$median_rel_bias*100, cex=pt.cex,pch=17,col="red")
points(abs_rel_bias$mode1_rel_bias*100, cex=pt.cex,pch=15,col="blue")
points(abs_rel_bias$mode0.1_rel_bias*100, cex=pt.cex,pch=18,col="green")
axis(4)
mtext("c)",2,padj=-9.5, las=1, line=3)


dev.off()


# axis(1,1:24,rep(c(0,0.1,0.2,0.4),6), cex.axis=1.2)
axis(1,1:24,rep(c(20,40,80),8), cex.axis=1.2)


}