
rm(list=ls())

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
	boxplot(estimate~N_within+type,bp2_long,boxwex=0.5,xlab="Within group sample size", names=rep(c(2,4,8),3), ylab="Estimate", pch=19, cex=0.5)
	abline(h=0.2, col="red", lwd=2)
	abline(v=c(3.5,6.5))
	mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))

	boxplot(estimate~N_within+type,bp0_long,boxwex=0.5,xlab="Within group sample size", ylab="Estimate", pch=19, cex=0.5, names=rep(c(2,4,8),3), ylim=c(0,0.3))#
	abline(h=0, col="red", lwd=2)
	abline(v=c(3.5,6.5))
	mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))


}

setEPS()
pdf(paste0(wd,"Figures/Fig2.pdf"), height=6, width=8)
{
	bias_plot0 <- subset(actual, ICC ==0 & N_within==2)[,1:6]
	bp0_long<-Stack(bias_plot0, col2stack = c("freq","mean","mode"),value.name = "estimate",group.name = "type")

	bias_plot2 <- subset(actual, ICC==0.2 & N_within==2)[,1:6]
	bp2_long<-Stack(bias_plot2, col2stack = c("freq","mean","mode"),value.name = "estimate",group.name = "type")

	par(mfrow=c(2,1), mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.75,0))
	boxplot(estimate~N_group+type,bp2_long,boxwex=0.5,xlab="Number of Groups",  ylab="Estimate", pch=19, cex=0.5, names=rep(c(20,40,80,160),3),  ylim=c(0,1))
	abline(h=0.2, col="red", lwd=2)
	abline(v=c(4.5,8.5))
	mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))
	mtext("a)",2,padj=-9, las=1, line=2)

	boxplot(estimate~N_group+type,bp0_long,boxwex=0.5,xlab="Number of Groups", ylab="Estimate", pch=19, cex=0.5, names=rep(c(20,40,80,160),3),  ylim=c(0,1))#
	abline(h=0, col="red", lwd=2)
	abline(v=c(4.5,8.5))
	mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))
	mtext("b)",2,padj=-9, las=1, line=2)

}
dev.off()



## check ESS
 boxplot(ESS~N_group+N_within+ICC,actual, las=2)
 abline(h=500)

