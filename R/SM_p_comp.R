
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


setEPS()
pdf(paste0(wd,"Figures/FigSM_sim_perm.pdf"), height=11, width=11)
{


par(mfrow=c(2,2), cex.lab=1.5)
plot(mean_perm~mean_boot,all, pch=pch, cex=cex, col=col, main="Posterior Mean", xlab="Simulation", ylab="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"mean_boot"],all[,"mean_perm"]),3)), cex=1.5)

plot(median_perm~median_boot,all, pch=pch, cex=cex, col=col, main="Posterior Median", xlab="Simulation", ylab="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"median_boot"],all[,"median_perm"]),3)), cex=1.5)

plot(mode1_perm~mode1_boot,all, pch=pch, cex=cex, col=col, main="Posterior Mode (adjust=1)", xlab="Simulation", ylab="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"mode1_boot"],all[,"mode1_perm"]),3)), cex=1.5)

plot(mode0.1_perm~mode0.1_boot,all, pch=pch, cex=cex, col=col, main="Posterior Mode (adjust=0.1)", xlab="Simulation", ylab="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"mode0.1_boot"],all[,"mode0.1_perm"]),3)), cex=1.5)
}

dev.off()


setEPS()
pdf(paste0(wd,"Figures/FigSM_metric_comp.pdf"), height=18, width=9)
{
par(mfrow=c(4,2), cex.lab=1.5)
plot(mean_boot~median_boot,all, pch=pch, cex=cex, col=col, xlab="Median",ylab="Mean", main="Simulation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"mean_boot"],all[,"median_boot"]),3)), cex=1.5)
plot(mean_perm~median_perm,all, pch=pch, cex=cex, col=col, xlab="Median",ylab="Mean", main="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"mean_perm"],all[,"median_perm"]),3)), cex=1.5)

plot(mean_boot~mode1_boot,all, pch=pch, cex=cex, col=col, xlab="Mode-1",ylab="Mean", main="Simulation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"mean_boot"],all[,"mode1_boot"]),3)), cex=1.5)
plot(mean_perm~mode1_perm,all, pch=pch, cex=cex, col=col, xlab="Mode-1",ylab="Mean", main="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"mean_perm"],all[,"mode1_perm"]),3)), cex=1.5)


plot(mean_boot~mode0.1_boot,all, pch=pch, cex=cex, col=col, xlab="Mode-0.1",ylab="Mean", main="Simulation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"mean_boot"],all[,"mode0.1_boot"]),3)), cex=1.5)
plot(mean_perm~mode0.1_perm,all, pch=pch, cex=cex, col=col, xlab="Mode-0.1",ylab="Mean", main="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"mean_perm"],all[,"mode0.1_perm"]),3)), cex=1.5)


plot(mode1_boot~mode0.1_boot,all, pch=pch, cex=cex, col=col, xlab="Mode-0.1",ylab="Mode-1", main="Simulation")#col=alpha(1:4,0.3)[as.factor(all$N_group)])
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"mode1_boot"],all[,"mode0.1_boot"]),3)), cex=1.5)
plot(mode1_perm~mode0.1_perm,all, pch=pch, cex=cex, col=col, xlab="Mode-0.1",ylab="Mode-1", main="Permutation")
abline(0,1, col="red")
text(0.1,0.9,paste0("r=",round(cor(all[,"mode1_perm"],all[,"mode0.1_perm"]),3)), cex=1.5)

}
dev.off()

