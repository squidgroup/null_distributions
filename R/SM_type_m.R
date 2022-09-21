rm(list=ls())

library(scales)
library(beeswarm)

wd <- "~/github/bayes_perm/"

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
		mean = as.numeric(x$actual["mean"]),
		median = as.numeric(x$actual["median"]),
		mode1 = as.numeric(x$actual["mode1"]),
		mode0.1 = as.numeric(x$actual["mode0.1"]),
		mean_median_dif = as.numeric(x$actual["mean"] - x$actual["median"]),
		mean_mode_dif = as.numeric(x$actual["mean"] - x$actual["mode1"]),
		freq_perm = p_func(x$actual["freq"],x$perm[,"freq"]),
		mean_perm = p_func(x$actual["mean"],x$perm[,"mean"]),
		median_perm = p_func(x$actual["median"],x$perm[,"median"]),
		mode0.1_perm = p_func(x$actual["mode0.1"],x$perm[,"mode0.1"]),
		mode1_perm = p_func(x$actual["mode1"],x$perm[,"mode1"])
		)
	}))
})))

par(mfrow=c(1,2))
plot(mean_median_dif~median_perm,p_perm, pch=19, cex=0.2, col=scales::alpha(1,0.5))
abline(v=c(0.05))
plot(mean_mode_dif~median_perm,p_perm, pch=19, cex=0.2, col=scales::alpha(1,0.5))
abline(h=c(0.007))

p_boot<-as.data.frame(do.call(rbind,lapply(results_b, function(i){
  load(paste0(wd,"Data/Intermediate/",i))
  t(sapply(boot_out, function(x){
		c(
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
power_dat <- aggregate(cbind(median_boot, mean_boot, mode1_boot,mode0.1_boot)~N_group+ N_within+ICC,all,function(x) mean(x<0.05))

type_m <- aggregate(cbind(median,mean,mode1,mode0.1)~N_group+ N_within+ICC,subset(all, mode1_perm<0.05),function(x)mean(abs(x)))

setEPS()
pdf(paste0(wd,"Figures/FigSM_type_m.pdf"), height=9, width=9)

{
par(mfrow=c(2,2), cex.lab=1.5)
plot((type_m$mean)/type_m$ICC~ power_dat$mean_boot, col=c(1:4)[factor(type_m$ICC)], pch=19, ylim=c(1,5.5), xlab="Power",ylab="Type M error", main="Mean")
plot((type_m$median)/type_m$ICC~ power_dat$median_boot, col=c(1:4)[factor(type_m$ICC)], pch=19, ylim=c(1,5.5), xlab="Power",ylab="Type M error", main="Median")
plot((type_m$mode1)/type_m$ICC~ power_dat$mode1_boot, col=c(1:4)[factor(type_m$ICC)], pch=19, ylim=c(1,5.5), xlab="Power",ylab="Type M error", main="Mode-1")
plot((type_m$mode0.1)/type_m$ICC~ power_dat$mode0.1_boot, col=c(1:4)[factor(type_m$ICC)], pch=19, ylim=c(1,5.5), xlab="Power",ylab="Type M error", main="Mode-0.1")

}
dev.off()