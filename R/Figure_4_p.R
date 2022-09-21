
rm(list=ls())

library(scales)
library(beeswarm)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))
	

# x<-perm_out[[1]]
# colnames(x$perm)
#  tnames <- c("freq","mode0.1","mode1","median","mean","LCI","UCI","ESS")

 
# lapply(results, function(i){
#   load(paste0(wd,"Data/Intermediate/",i))
#   for(j in 1:length(PB_out)){
#   	names(PB_out[[j]]$actual)<-c("freq","mode0.1","mode1","median","mean","LCI","UCI","ESS")
#   }
#   save(PB_out, file=paste0(wd,"Data/Intermediate/",i))
# })
# check<-lapply(results, function(i){
#   load(paste0(wd,"Data/Intermediate/",i))
#   t(sapply(PB_out, function(x){
# 		c(a=all(tnames == names(x$actual) ),p=all( tnames == colnames(x$perm) ),b=all( tnames == colnames(x$boot)))
# 	}))
# })
# all(check)

files <- list.files(paste0(wd,"Data/Intermediate"))
results <- files[grep("gaus_perm",files)]
list_names <- gsub("gaus_perm_|.Rdata","",results)

results_b <- files[grep("gaus_boot",files)]
list_names_b <- gsub("gaus_boot_|.Rdata","",results_b)



p_perm<-as.data.frame(do.call(rbind,lapply(results, function(i){
  load(paste0(wd,"Data/Intermediate/",i))
  t(sapply(perm_out, function(x){
		names(x$param) <-c("pop", "ICC", "N_group", "N_within")
		c(x$param,
		freq_perm = p_func(x$actual["freq"],x$perm[,"freq"]),
		mean_perm = p_func(x$actual["mean"],x$perm[,"mean"]),
		median_perm = p_func(x$actual["median"],x$perm[,"median"]),
		mode0.1_perm = p_func(x$actual["mode0.1"],x$perm[,"mode0.1"]),
		mode1_perm = p_func(x$actual["mode1"],x$perm[,"mode1"])
		# median_ratio_boot = x$actual["median"]/median(x$boot[,"median"]),
		# mode1_ratio_boot = x$actual["mode1"]/median(x$boot[,"mode1"])

		)
	}))
})))

p_boot<-as.data.frame(do.call(rbind,lapply(results_b, function(i){
  load(paste0(wd,"Data/Intermediate/",i))
  t(sapply(boot_out, function(x){
		c(
		freq_boot = p_func(x$actual["freq"],x$boot[,"freq"]),
		mean_boot = p_func(x$actual["mean"],x$boot[,"mean"]),
		median_boot = p_func(x$actual["median"],x$boot[,"median"]),
		mode0.1_boot = p_func(x$actual["mode0.1"],x$boot[,"mode0.1"]),
		mode1_boot = p_func(x$actual["mode1"],x$boot[,"mode1"])
		# median_ratio_boot = x$actual["median"]/median(x$boot[,"median"]),
		# mode1_ratio_boot = x$actual["mode1"]/median(x$boot[,"mode1"])

		)
	}))
})))

nrow(p_boot)
nrow(p_perm)


head(p_perm)
head(p_boot)

all<-cbind(p_perm,p_boot)
head(all)

###
# -- FIGURE 3
###


# {
# ICCs <- c(0,0.1,0.2,0.4)
# par(mfrow=c(4,1),mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))
# for(j in seq_along(ICCs)){
# 	plot_dat <- subset(all, ICC==ICCs[j])
# 	beeswarm(mean_boot~N_group + N_within,plot_dat, pch=19, cex=0.3, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Number of Groups", labels=rep(c(20,40,80),3), ylab=bquote(P[perm]),  ylim=c(0,1))
# 	abline(v=3.5)
# # mtext(paste0(letters[j],")"),2,padj=-5, las=1, line=2)

# }
# }




setEPS()
pdf(paste0(wd,"Figures/Fig4_p.pdf"), height=6, width=8)

{
par(mfrow=c(2,1),mar=c(0.5,4,4,1), cex.axis=0.75, mgp=c(2,0.5,0))
	beeswarm(median_boot~N_group + ICC,subset(all, N_within==2), pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="", labels=rep(c(20,40,80),3), ylab="p-value",  ylim=c(0,1), xaxt="n")
	abline(v=c(3.5,6.5,9.5))
	axis(3,c(2,5,8,11),c(0,0.1,0.2,0.4), tick=FALSE)
	mtext("ICC",side=3, line=2)
	mtext("a)",2,padj=-8, las=1, line=2)

	par(mar=c(4,4,0.5,1))
	beeswarm(median_boot~ N_group + ICC,subset(all, N_within==4), pch=19, cex=0.1, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Number of Groups", labels=rep(c(20,40,80),3), ylab="p-value",  ylim=c(0,1))	
	abline(v=c(3.5,6.5,9.5))
	mtext("b)",2,padj=-8, las=1, line=2)
}
dev.off()