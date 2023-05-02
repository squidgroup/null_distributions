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




power_dist <- NULL

within_groups <- c(2,4)
Ns <- c(20,40,80)
ICCs <- c(0.1,0.2,0.4)


for(i in ICCs){
	for(j in Ns){
		for(k in within_groups){
			null <- subset(actual, ICC==0 & N_within==k & N_group==j)$median
			alt <- subset(actual, N_within==k & N_group==j & ICC==i)$median
			power_dist<- rbind(power_dist,c(ICC=i, N_group=j, N_within=k, power=mean(sapply(alt, function(x) p_func(x,null))<0.05)))
		}
	}
}
power_dist <- as.data.frame(power_dist)
## looks weird for n_within=2, ICC=0.1 and 0.2
## plot out histograms

hist(subset(actual, ICC==0.1 & N_within==2 & N_group==40)$median)
hist(subset(actual, ICC==0.1 & N_within==2 & N_group==80)$median)



setEPS()
pdf(paste0(wd,"Figures/Fig6_power.pdf"), height=5, width=10)

{
ICCs <- c(0,0.1,0.2,0.4)
	power_dat <- aggregate(cbind(mode1_boot,mode0.1_boot,mean_boot,median_boot,median_perm)~N_group+ N_within+ICC,all,function(x) mean(x<0.05))
par(mfrow=c(1,2),mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))
plot(NA,xlim=c(20,80),xlab="Number of Groups", ylab="Power",  ylim=c(0,1),)
abline(h=0.05, col="grey")
# legend("topleft",c("power","perm","boot","ICC=0","ICC=0.1","ICC=0.2","ICC=0.4"), pch=c(17,18,19,rep(NA,4)), col=c("grey","grey","grey",1:4),pt.bg=c("grey"), lty=c(0,0,0,1,1,1,1),lwd=2)
 legend("topleft",c("power","perm","sim"), pch=c(17,18,19), col="grey",pt.bg=c("grey"), bty="n")
	# mtext("a)",3,line=1,adj=0)
mtext("a)",2,padj=-15, las=1, line=2)

# abline(h=c(0.025,0.075), col="grey")
for(j in seq_along(ICCs)){
	points(median_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==2, col=j, pch=19, type="b")
	points(median_perm~N_group,power_dat, subset=ICC==ICCs[j]&N_within==2, col=j, pch=18, type="b")
		points(power~N_group,power_dist, subset=ICC==ICCs[j]&N_within==2, col=j, pch=17, type="b")
}

text(rep(65,4),c(0.1,0.25,0.55,0.97),c("ICC=0","ICC=0.1","ICC=0.2","ICC=0.4"), col=1:4)

plot(NA,xlim=c(20,80),xlab="Number of Groups", ylab="Power",  ylim=c(0,1))
abline(h=0.05, col="grey")
# abline(h=c(0.028,0.081), col="grey")
for(j in seq_along(ICCs)){
	points(median_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==4, col=j, pch=19, type="b")
		points(median_perm~N_group,power_dat, subset=ICC==ICCs[j]&N_within==4, col=j, pch=18, type="b")
		points(power~N_group,power_dist, subset=ICC==ICCs[j]&N_within==4, col=j, pch=17, type="b")
}
mtext("b)",2,padj=-15, las=1, line=2)


}
dev.off()




power_dat <- aggregate(cbind(mean_boot,mean_perm)~N_group+ N_within+ICC,all,function(x) mean(x<0.05))

power_dat2 <- aggregate(cbind(mean_boot,mean_perm)~N_group+ N_within+ICC,all,mean)

power_dat

plot(power_dat$mean_boot,power_dat2$mean_boot)
plot(power_dat$mean_perm,power_dat2$mean_perm, ylab="Mean p-value", xlab="Power")

alpha <- 1

mu = seq(0,alpha,0.01)
power = pbeta(0.05,1,(alpha/mu -alpha))#1-(1-0.05)^(1/mu -1)

lines(mu~power)
