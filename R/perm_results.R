
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



#### check power

actual <- as.data.frame(do.call(rbind,lapply(all, function(z){
	t(sapply(z, function(x) c(x$param,x$actual[1:6])))
})))


null_dist <- subset(actual, ICC==0)
par(mfrow=c(3,4), mar=c(1,1,1,1))
for(i in c(2,4,8)){
	for(j in c(20,40,80,160)){
		try(hist(subset(null_dist, N_within==i & N_group==j)$mean, xlim=c(0,1), main="", xlab="Posterior Mean"))
	}
}



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

{
# layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0),ncol=2))
within_groups <-c(2,4,8)
par(mfrow=c(3,1),mar=c(5,5,1,1))

plot_dat <- subset(power, ICC==0)
plot(mean~N_group,plot_dat,xlim=c(20,160),ylim=c(0,0.2), ylab="False Positive Rate", xlab="Number of Groups")
abline(h=0.05, lty=2, col="grey")
	for(i in seq_along(within_groups)){
		lines(mean~N_group,subset(plot_dat, N_within==within_groups[i]), col=i)	
	}
polygon(c(0,180,190,0),c(0.01,0.01,0.09,0.09), border=NA, col=alpha("grey",0.2))

## with this sample size (100 sims) we would probably expect anything between 0 and 10
# binom.test(10,100,0.05)
# binom.test(0,100,0.05)

for(j in c(0.2,0.4)){
	
	plot_dat <- subset(power, ICC==j)
	plot(mean~N_group,plot_dat,xlim=c(20,160),ylim=c(0,1), ylab="Power", xlab="Number of Groups")

	for(i in seq_along(within_groups)){
		lines(mean~N_group,subset(plot_dat, N_within==within_groups[i]), col=i)	
	}
}
}


power_dist <- NULL

# i=2;j=20;k=0.2
for(i in c(0.2,0.4,0.8)){
	for(j in c(20,40,80,160)){
		for(k in c(2,4,8,16)){
			null <- subset(null_dist, N_within==k & N_group==j)$mean
			alt <- subset(actual, N_within==k & N_group==j & ICC==i)$mean
			power_dist<- rbind(power_dist,c(ICC=i, N_group=j, N_within=k, power=mean(sapply(alt, function(x) p_func(x,null))<0.05)))
		}
	}
}
power_dist <- as.data.frame(power_dist)

par(mar=c(5,5,1,1), mfrow=c(2,1))
for(j in c(0.2,0.4)){
	plot_dat <- subset(power, ICC==j)
	plot_dat2 <- subset(power_dist, ICC==j)
	plot(mean~N_group,plot_dat,xlim=c(20,160),ylim=c(0,1), ylab="Power", xlab="Number of Groups")
	points(power~N_group,plot_dat2, pch=17)
	for(i in seq_along(within_groups)){
		lines(mean~N_group,subset(plot_dat, N_within==within_groups[i]), col=i)	
		lines(power~N_group,subset(plot_dat2, N_within==within_groups[i]), col=i)	
	}
}



###
# -- FIGURE 4
###
par(mfrow=c(1,1),mar=c(5,5,1,1))
plot(freq~mean,power, xlab="Bayesian Power", ylab="ML Power"); abline(0,1, col="grey")
points(LRT~mean,power, pch=19, col="red")




# for(i in results){
# 	load(paste0(wd,"Data/Intermediate/",i))
# 	sim_dat<-lapply(sim_dat, function(x) {
# 		names(x$actual) <- c("freq","mode","mean","LCI","UCI","ESS")
# 		colnames(x$null) <- c("freq","mode","mean","LCI","UCI","ESS")
# 		x
# 	})

# 	save(sim_dat, file=paste0(wd,"Data/Intermediate/",i))
# 	rm("sim_dat")
# }