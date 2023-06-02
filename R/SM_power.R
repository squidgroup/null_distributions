rm(list=ls())

library(scales)
library(beeswarm)

wd <- "~/github/null_distributions/"

source(paste0(wd,"R/00_functions.R"))

files <- list.files(paste0(wd,"Data/Intermediate"))
results <- files[grep("gaus_sims",files)]
list_names <- gsub("gaus_sims_|.Rdata","",results)

actual<-as.data.frame(do.call(rbind,lapply(results, function(y){
	load(paste0(wd,"Data/Intermediate/",y))
	do.call(rbind,lapply(out, function(x){
	  z <- c(x$param,x$summary)
	  # ICC=x$ICC,N_group=x$N_group, N_within=x$N_within
	  names(z) <- c("pop","ICC","N_group","N_within","freq","mode0.1","mode1","median","mean","LCI","UCI","ESS")	
	  z
	}))
} )))


files_p <- list.files(paste0(wd,"Data/Intermediate"))
results_p <- files_p[grep("gaus_perm",files_p)]
list_names_p <- gsub("gaus_perm_|.Rdata","",results_p)

results_b <- files[grep("gaus_boot",files)]
list_names_b <- gsub("gaus_boot_|.Rdata","",results_b)

p_perm<-as.data.frame(do.call(rbind,lapply(results_p, function(i){
  load(paste0(wd,"Data/Intermediate/",i))
  t(sapply(perm_out, function(x){
		names(x$param) <-c("pop", "ICC", "N_group", "N_within")
		c(x$param,
		freq_perm = p_func(x$actual["freq"],x$perm[,"freq"]),
		mean_perm = p_func(x$actual["mean"],x$perm[,"mean"]),
		median_perm = p_func(x$actual["median"],x$perm[,"median"]),
		mode0.1_perm = p_func(x$actual["mode0.1"],x$perm[,"mode0.1"]),
		mode1_perm = p_func(x$actual["mode1"],x$perm[,"mode1"])
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
		)
	}))
})))


all<-cbind(p_perm,p_boot)

pch=19
cex=0.2
col=alpha(1,0.3)



ICCs <- c(0,0.1,0.2,0.4)
power_dat <- aggregate(cbind(mode1_boot,mode0.1_boot,mean_boot,median_boot,median_perm)~N_group+ N_within+ICC,all,function(x) mean(x<0.05))

setEPS()
pdf(paste0(wd,"Figures/FigSM_power_comp.pdf"), height=5, width=10)

	{
		par(mfrow=c(1,2),mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))
		plot(NA,xlim=c(20,80),xlab="Number of Groups", ylab="Power / FPR",  ylim=c(0,1))
		abline(h=0.05, col="grey")
		# abline(h=c(0.025,0.075), col="grey")
		legend("topleft",c("median","mean","mode1","mode0.1"), pch=c(18,19,17,15), col="grey",pt.bg=c("grey"), bty="n")
mtext("a)",2,padj=-15, las=1, line=2)

text(rep(65,4),c(0.1,0.25,0.55,0.97),c("ICC=0","ICC=0.1","ICC=0.2","ICC=0.4"), col=1:4)


		for(j in seq_along(ICCs)){
			points(mean_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==2, col=j, pch=19, type="b")
			points(median_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==2, col=j, pch=18, type="b")
			points(mode1_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==2, col=j, pch=17, type="b")
			points(mode0.1_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==2, col=j, pch=15, type="b")
		
			# lines(mean_boot~N_group,power_dat, subset=ICC==ICCs[j], col=j, lwd=2)
		}
		
			plot(NA,xlim=c(20,80),xlab="Number of Groups", ylab="Power / FPR",  ylim=c(0,1))
			abline(h=0.05, col="grey")
			# abline(h=c(0.028,0.081), col="grey")
			for(j in seq_along(ICCs)){
				points(mean_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==4, col=j, pch=19, type="b")
				points(median_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==4, col=j, pch=18, type="b")
				points(mode1_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==4, col=j, pch=17, type="b")
				points(mode0.1_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==4, col=j, pch=15, type="b")
			# lines(mean_boot~N_group,power_dat, subset=ICC==ICCs[j], col=j, lwd=2)
			}
			mtext("b)",2,padj=-15, las=1, line=2)

		}
dev.off()