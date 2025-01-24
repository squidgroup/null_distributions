rm(list=ls())

library(scales)
library(beeswarm)

wd <- "~/github/null_distributions/"

source(paste0(wd,"R/00_functions.R"))

files <- list.files(paste0(wd,"Data/Intermediate"))



results_gaus <- files[grep("gaus_sim",files)]
se_gaus<-as.data.frame(do.call(rbind,lapply(results_gaus, function(i){
  load(paste0(wd,"Data/Intermediate/",i))

  # boot_out<-lapply(out, function(x){
  # 	x$data
  # 	x$param <- as.numeric(strsplit(gsub("gaus_boot_|.Rdata","",i),"_")[[1]])
		
	# })
	 t(sapply(out, function(x){
		param <-  x$param[2:4]
	 	names(param) <-c("ICC", "N_group", "N_within")
		c(param, x$summary,se=sd(x$post$sigma2_ID))
  }))
})))
head(se_gaus)

results_gaus <- files[grep("gaus_boot",files)]
p_gaus<-as.data.frame(do.call(rbind,lapply(results_gaus, function(i){
  load(paste0(wd,"Data/Intermediate/",i))
  t(sapply(boot_out, function(x){
		c(median2 =as.numeric(x$actual["median"]), 
			median_p = p_func(x$actual["median"],x$boot[,"median"]),
			mean_p = p_func(x$actual["mean"],x$boot[,"mean"])

			)
  }))
})))

head(p_gaus)

results <- cbind(se_gaus,p_gaus)
results$cv <- results$se/results$median
results$z_median <- results$median/results$se
results$z_mean <- results$mean/results$se
results$z_mode1 <- results$mode1/results$se
results$z_mode0.1 <- results$mode0.1/results$se

results$median_zp <- 2*(1-pnorm(results$z_median))
results$mean_zp <- 2*(1-pnorm(results$z_mean))


power_dat <- aggregate(cbind(mean_p,median_p,median_zp,mean_zp)~N_group+ N_within+ICC,results,function(x) mean(x<0.05))

pdf("~/github/null_distributions/Figures/FigSM_z_p.pdf", height=5, width=9)
{
par(mfrow=c(1,2))

plot(results$z_mean,results$mean_p, pch=19, cex=0.05, col=scales::alpha(1,0.3), ylab="P-value", xlab="Z-ratio")
abline(h=0.05, v=qnorm(0.95), col="grey")

plot(mean_p~mean_zp,power_dat, ylab="Power p-value", xlab="Power z-ratio", pch=19)
abline(0,1)
}
dev.off()


{
par(mfrow=c(2,2))
	plot(results$z_median,results$median_p, pch=19, cex=0.05, col=scales::alpha(1,0.5))
	abline(h=0.05, v=qnorm(0.95))
plot(results$z_mean,results$mean_p, pch=19, cex=0.05, col=scales::alpha(1,0.5))
abline(h=0.05, v=qnorm(0.95))
plot(results$z_mode1,results$median_p, pch=19, cex=0.05, col=scales::alpha(1,0.5))
abline(h=0.05, v=qnorm(0.95))
plot(results$z_mode0.1,results$median_p, pch=19, cex=0.05, col=scales::alpha(1,0.5))
abline(h=0.05, v=qnorm(0.95))
}

cor(log(results$z_median),results$median_p)

{
par(mfrow=c(2,2))
	plot((results$median_zp),results$median_p, pch=19, cex=0.05, col=scales::alpha(1,0.5))
	abline(h=0.05,v=0.05)
plot(results$mean_zp,results$median_p, pch=19, cex=0.05, col=scales::alpha(1,0.5))
abline(h=0.05,v=0.05)
plot(results$mean_zp,results$median_zp, pch=19, cex=0.05, col=scales::alpha(1,0.5))
abline(0,1)
# plot(results$z_mode0.1,results$median_p, pch=19, cex=0.05, col=scales::alpha(1,0.5))
# abline(h=0.05,v=0.05)
}



