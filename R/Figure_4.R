
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



#### power from permutation tests

power <- as.data.frame(t(sapply(all,function(z){
	c(
		z[[1]]$param,LRT=mean(sapply(z, function(x) x$actual["LRT"])< 0.05),
		freq=mean(sapply(z, function(x) p_func(x$actual["freq"],x$null[,"freq"]) < 0.05)),
		mean=mean(sapply(z, function(x) p_func(x$actual["mean"],x$null[,"mean"]) < 0.05)),
		mode=mean(sapply(z, function(x) p_func(x$actual["mode"],x$null[,"mode"]) < 0.05))
	)
})))
power <- power[order(power$ICC,power$N_group,power$N_within),]
power<-subset(power,ICC!=0.8 & N_within!=16)


#### power from null distribution

actual <- as.data.frame(do.call(rbind,lapply(all, function(z){
	t(sapply(z, function(x) c(x$param,x$actual[1:6])))
})))

null_dist <- subset(actual, ICC==0)
power_dist <- NULL

within_groups <- c(2,4,8)
ICCs <- c(0.2,0.4)
Ns <- c(20,40,80,160)

for(i in ICCs){
	for(j in Ns){
		for(k in within_groups){
			null <- subset(null_dist, N_within==k & N_group==j)$mean
			alt <- subset(actual, N_within==k & N_group==j & ICC==i)$mean
			power_dist<- rbind(power_dist,c(ICC=i, N_group=j, N_within=k, power=mean(sapply(alt, function(x) p_func(x,null))<0.05)))
		}
	}
}
power_dist <- as.data.frame(power_dist)


cols <- c("black","red","blue")
setEPS()
pdf(paste0(wd,"Figures/Fig4.pdf"), height=6, width=4)

{
par(mar=c(4,4,1,1), mfrow=c(2,1), cex.axis=0.75, mgp=c(2,0.5,0))
for(j in seq_along(ICCs)){
	plot_dat <- subset(power, ICC==ICCs[j])
	plot_dat2 <- subset(power_dist, ICC==ICCs[j])
	plot(mean~N_group,plot_dat,xlim=c(20,160),ylim=c(0,1), ylab="Power", xlab="Number of Groups")
	points(power~N_group,plot_dat2, pch=17)
	for(i in seq_along(within_groups)){
		lines(mean~N_group,subset(plot_dat, N_within==within_groups[i]), col=cols[i])	
		lines(power~N_group,subset(plot_dat2, N_within==within_groups[i]), col=cols[i])	
	}
	mtext(paste0(letters[j],")"),2,padj=-9, las=1, line=3)
	# mtext("a)",2,padj=-12, las=1, line=2)

}

}
dev.off()

