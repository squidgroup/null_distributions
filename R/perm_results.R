
rm(list=ls())

library(scales)
wd <- "~/github/bayes_perm/"

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
# -- FIGURE 3
###

p_func <- function(actual,null) mean(actual<null)

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

# plot(mode~mean,power); abline(0,1)

}
}

p_perm <- as.data.frame(do.call(rbind,lapply(all, function(z){
	t(sapply(z, function(x) c(x$param,
		freq = p_func(x$actual["freq"],x$null[,"freq"]),
		mean=p_func(x$actual["mean"],x$null[,"mean"]),
		mode= p_func(x$actual["mode"],x$null[,"mode"])
		)))
})))
p_perm<-subset(p_perm,ICC!=0.8 & N_within!=16)

par(mfrow=c(3,1),mar=c(5,5,1,1))
for(j in c(0,0.2,0.4)){
	plot_dat <- subset(p_perm, ICC==j)
	boxplot(mean~N_group + N_within,plot_dat,ylim=c(0,1), ylab="P_perm", xlab="Number of Groups",boxwex=0.5, names=rep(c(20,40,80,160),3), pch=19, cex=0.5)
	abline(v=c(4.5,8.5))
	# abline(h=c(0.25,0.75))
	mtext(c("2","4","8"),side=3,adj=c(0.15,0.5, 0.85))

}

# p_perm <- as.data.frame(t(sapply(all,function(z){
# 	c(
# 		z[[1]]$param,LRT=mean(sapply(z, function(x) x$actual["LRT"])< 0.05),
# 		freq=mean(sapply(z, function(x) p_func(x$actual["freq"],x$null[,"freq"]) < 0.05)),
# 		mean=mean(sapply(z, function(x) p_func(x$actual["mean"],x$null[,"mean"]) < 0.05)),
# 		mode=mean(sapply(z, function(x) p_func(x$actual["mode"],x$null[,"mode"]) < 0.05))
# 	)
# })))




###
# -- FIGURE 4
###
par(mfrow=c(1,1),mar=c(5,5,1,1))
plot(freq~mean,power, xlab="Bayesian Power", ylab="ML Power"); abline(0,1, col="grey")
points(LRT~mean,power, pch=19, col="red")




###
# -- FIGURE 2
###

joelStack <- function(x, col2stack, value.name="value", group.name="group", levels=NULL){
	if(is.character(col2stack)) {
		col2 <- NA
		for(i in 1:length(col2stack)) col2[i] <- which(colnames(x) %in% col2stack[i])
		col2stack <- col2
		}
	y <- stack(x[,col2stack])
	colnames(y) <- c(value.name,group.name)
	n.var <- length(col2stack)
	n.obs <- nrow(x)
	if(!is.null(levels)) y[,2] <- rep(levels,rep(n.obs,n.var))
	col.names <- colnames(x)[-c(col2stack)]
	z <- as.data.frame(x[, -c(col2stack)])
	colnames(z) <- col.names
	z.new <- z
	for (i in 1:(n.var - 1)) z.new <- rbind(z.new,z)
	return(cbind(z.new,y))
}


actual <- as.data.frame(do.call(rbind,lapply(all, function(z){
	t(sapply(z, function(x) c(x$param,x$actual[1:6])))
})))

bias_plot0 <- subset(actual, ICC ==0 & N_group==80)[,1:6]
bias_plot2 <- subset(actual, ICC ==0.2 & N_group==40)[,1:6]

bp2_long<-joelStack(bias_plot2, col2stack = c("freq","mean","mode"),value.name = "estimate",group.name = "type")  
bp0_long<-joelStack(bias_plot0, col2stack = c("freq","mean","mode"),value.name = "estimate",group.name = "type")    


{
	par(mfrow=c(2,1), mar=c(5,5,2,1), cex.axis=0.75)
	boxplot(estimate~N_within+type,bp2_long,boxwex=0.5,xlab="Within group sample size", names=rep(c(2,4,8,16),3), ylab="Estimate", pch=19, cex=0.5)
	abline(h=0.2, col="red", lwd=2)
	abline(v=c(4.5,8.5))
	mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))

	boxplot(estimate~N_within+type,bp0_long,boxwex=0.5,xlab="Within group sample size", ylab="Estimate", pch=19, cex=0.5, xlim=c(0.5,12.5))#, names=rep(c(2,4,8,16),3)
	abline(h=0, col="red", lwd=2)
	abline(v=c(4.5,8.5))
	mtext(c("ML","Posterior Mean","Posterior Mode"),side=3,adj=c(0.15,0.5, 0.85))


}


## check ESS
 boxplot(ESS~N_group+N_within+ICC,actual, las=2)
 abline(h=500)


#### check power


null_dist <- subset(actual, ICC==0)
par(mfrow=c(4,4), mar=c(1,1,1,1))
for(i in c(2,4,8,16)){
	for(j in c(20,40,80,160)){
		try(hist(subset(null_dist, N_within==i & N_group==j)$mean, xlim=c(0,1), main="", xlab="Posterior Mean"))
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
	points(power~N_group,plot_dat2, pch=19)
	for(i in seq_along(within_groups)){
		lines(mean~N_group,subset(plot_dat, N_within==within_groups[i]), col=i)	
		lines(power~N_group,subset(plot_dat2, N_within==within_groups[i]), col=i)	
	}
}

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