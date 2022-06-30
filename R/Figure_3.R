
rm(list=ls())

library(scales)
library(beeswarm)

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))

files <- list.files(paste0(wd,"Data/Intermediate"))
results <- files[grep("gaus_perm",files)]
list_names <- gsub("gaus_perm_|.Rdata","",results)


i=15
	

x<-perm_out[[1]]
colnames(x$perm)
 tnames <- c("freq","mode0.1","mode1","median","mean","LCI","UCI","ESS")

 
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

all$total_n <- all$N_group*all$N_within
head(all)
plot(log(median_ratio_boot.median)~median_boot, all, col=c(1:4)[as.factor(all$ICC)])
abline(v=0.05, h=1)
aggregate(log(median_ratio_boot.median)~ICC + N_within +N_group,all,var)
boxplot(median_ratio_boot.median~ICC,all)
boxplot(log(median_ratio_boot.median)~total_n+ICC,all, ylim=c(-5,5))
boxplot(log(mode1_ratio_boot.mode1)~total_n+ICC,all)

boxplot(median_boot~total_n+ICC,all)

plot(all[,"mean_perm"],all[,"mean_boot"])
plot(all[,"median_perm"],all[,"median_boot"])

plot(all[,"mean_perm"],all[,"median_perm"])
plot(all[,"mean_boot"],all[,"median_boot"])
plot(all[,"mean_boot"],all[,"mode0.1_boot"])
plot(all[,"mean_perm"],all[,"mode0.1_perm"])

plot(all[,"mean_boot"],all[,"freq_boot"])
plot(all[,"mean_perm"],all[,"freq_perm"])

plot(all[,"freq_boot"],all[,"freq_perm"])

plot(all[,"mean_boot"],all[,"mode1_boot"])
plot(all[,"median_boot"],all[,"mode1_boot"])


## get all results into a list
all <- list()
for(i in 1:length(results)){
	load(paste0(wd,"Data/Intermediate/",results[i]))
	all[[list_names[i]]] <- sim_dat
	rm("sim_dat")
}

###
# -- FIGURE 3
###
names(all[[1]])

p_perm <- as.data.frame(do.call(rbind,lapply(all, function(z){
	t(sapply(z, function(x) c(x$param,
		freq = p_func(x$actual["freq"],x$null[,"freq"]),
		mean=p_func(x$actual["mean"],x$null[,"mean"]),
		median=p_func(x$actual["median"],x$null[,"median"]),
		mode= p_func(x$actual["mode"],x$null[,"mode"])
		)))
})))
# p_perm<-subset(p_perm,ICC!=0.8 & N_within!=16)
head(p_perm)

# setEPS()
# pdf(paste0(wd,"Figures/Fig3.pdf"), height=6, width=8)
# ICCs <- c(0,0.2,0.4)
# {
# par(mfrow=c(3,1),mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))
# for(j in seq_along(ICCs)){
# 	plot_dat <- subset(p_perm, ICC==ICCs[j])
# 	# boxplot(mean~N_group + N_within,plot_dat,ylim=c(0,1), ylab=bquote(P[perm]), xlab="Number of Groups",boxwex=0.5, names=rep(c(20,40,80,160),3), pch=19, cex=0.5, col=rep(grey(c(1,0.8,0.6,0.4)),3))
# 	beeswarm(mean~N_group + N_within,plot_dat, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Number of Groups", labels=rep(c(20,40,80,160),3), ylab=bquote(P[perm]),  ylim=c(0,1))

# 	abline(v=c(4.5,8.5))
# 	# abline(h=c(0.25,0.75))
# 	mtext(c("2","4","8"),side=3,adj=c(0.15,0.5, 0.85), cex=0.8)
# mtext(paste0(letters[j],")"),2,padj=-5, las=1, line=2)

# }
# }
# dev.off()



p_perm<-subset(p_perm,ICC%in%c(0,0.2) & N_within==2)


setEPS()
pdf(paste0(wd,"Figures/Fig3.pdf"), height=6, width=8)
ICCs <- c(0,0.1,0.2,0.4)
{
par(mfrow=c(4,1),mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))
for(j in seq_along(ICCs)){
	plot_dat <- subset(all, ICC==ICCs[j])
	beeswarm(mean_boot~N_group + N_within,plot_dat, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Number of Groups", labels=rep(c(20,40,80),3), ylab=bquote(P[perm]),  ylim=c(0,1))
	abline(v=3.5)
# mtext(paste0(letters[j],")"),2,padj=-5, las=1, line=2)

}
}
dev.off()

ICCs <- c(0,0.1,0.2,0.4)

{
	power_dat <- aggregate(cbind(mode1_boot,mode0.1_boot,mean_boot,median_boot,median_perm)~N_group+ N_within+ICC,all,function(x) mean(x<=0.05))
par(mfrow=c(2,1),mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))
plot(NA,xlim=c(20,80),xlab="Number of Groups", ylab="Power",  ylim=c(0,1))
abline(h=0.05, col="grey")
abline(h=c(0.025,0.075), col="grey")
for(j in seq_along(ICCs)){
	points(median_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==2, col=j, pch=19, type="b")
		points(median_perm~N_group,power_dat, subset=ICC==ICCs[j]&N_within==2, col=j, pch=18, type="b")

	# points(mode0.1_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==2, col=j, pch=18, type="b")
	# 	points(mode1_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==2, col=j, pch=17, type="b")

	# lines(mean_boot~N_group,power_dat, subset=ICC==ICCs[j], col=j, lwd=2)
}

plot(NA,xlim=c(20,80),xlab="Number of Groups", ylab="Power",  ylim=c(0,1))
abline(h=0.05, col="grey")
abline(h=c(0.028,0.081), col="grey")
for(j in seq_along(ICCs)){
	points(median_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==4, col=j, pch=19, type="b")
		points(median_perm~N_group,power_dat, subset=ICC==ICCs[j]&N_within==4, col=j, pch=18, type="b")

		# points(mode0.1_boot~N_group,power_dat, subset=ICC==ICCs[j]&N_within==4, col=j, pch=18, type="b")

	# lines(mean_boot~N_group,power_dat, subset=ICC==ICCs[j], col=j, lwd=2)
}
}


par(mfrow=c(2,2))
plot(mode1_boot~mode0.1_boot,p_perm)
plot(mean_boot~mode1_boot,p_perm)
plot(mean_boot~mode0.1_boot,p_perm)
plot(mean_boot~median_boot,p_perm)
plot(mean_boot~freq_boot,p_perm)
plot(mode0.1_boot~freq_boot,p_perm)
# cor(p_perm$mean,p_perm$mode)
