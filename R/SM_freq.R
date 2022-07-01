
rm(list=ls())

library("beeswarm")

library(scales)
wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/00_functions.R"))

files <- list.files(paste0(wd,"Data/Intermediate"))
results <- files[grep("gaus_sims",files)]
list_names <- gsub("gaus_sims_|.Rdata","",results)

## get all results into a list
# all <- list()
# for(i in 1:length(results)){
# 	load(paste0(wd,"Data/Intermediate/",results[i]))
# 	all[[list_names[i]]] <- out
# 	rm("out")
# }

actual<-as.data.frame(do.call(rbind,lapply(results, function(y){
	load(paste0(wd,"Data/Intermediate/",y))
	do.call(rbind,lapply(out, function(x){
	  z <- c(x$param,x$summary)
	  # ICC=x$ICC,N_group=x$N_group, N_within=x$N_within
	  names(z) <- c("pop","ICC","N_group","N_within","freq","mode0.1","mode1","median","mean","LCI","UCI","ESS")	
	  z
	}))
} )))
head(actual)

nrow(actual)

actual$freq_bias<-actual$freq-actual$ICC
actual$mean_bias<-actual$mean-actual$ICC
actual$median_bias<-actual$median-actual$ICC
actual$mode1_bias<-actual$mode1-actual$ICC
actual$mode0.1_bias<-actual$mode0.1-actual$ICC

actual$mean_rmse<-actual$mean-actual$ICC

rmse<-aggregate(cbind(freq_bias,mode0.1_bias,mode1_bias,mean_bias,median_bias)~N_group+N_within+ICC,actual,function(x) sqrt(mean(x^2)) )
bias<-aggregate(cbind(freq_bias,mode0.1_bias,mode1_bias,mean_bias,median_bias)~N_group+N_within+ICC,actual,mean)
abs_bias<-aggregate(cbind(freq_bias,mode0.1_bias,mode1_bias,mean_bias,median_bias)~N_group+N_within+ICC,actual,function(x) mean(abs(x)))

# plot(abs_bias$median_bias,rmse$median_bias);abline(0,1)
# plot(abs_bias$mode1_bias,rmse$mode1_bias);abline(0,1)
# plot(abs_bias$mean_bias,rmse$mean_bias);abline(0,1)
# plot(abs_bias$mode0.1_bias,rmse$mode0.1_bias);abline(0,1)
# plot(abs_bias$freq_bias,rmse$freq_bias);abline(0,1)


line_coords <- (1:7)*3+0.5
	par(mfrow=c(4,1), mar=c(2,4,3,1))

plot(bias$mean_bias, pch=19, ylim=c(-0.2,0.2), ylab="Bias", xaxt="n")
abline(h=0)
points(bias$freq_bias, pch=19,col="purple")

points(bias$median_bias, pch=19,col="red")
points(bias$mode1_bias, pch=19,col="blue")
points(bias$mode0.1_bias, pch=19,col="green")
	abline(v=line_coords)
axis(1,1:24,rep(c(20,40,80,160),6))

mtext(c("","0","","0.1","","0.2","","0.4",""),side=3,adj=c(seq(0,1,length.out=9)),line=1)
mtext(c("","2","4","2","4","2","4","2","4"),side=3,adj=c(seq(0,1,length.out=9)-1/16))
legend("bottomleft",c("mean","median","mode-0.1","mode-1"), pch=19, col=c(1,2,3,4))

plot(abs_bias$mean_bias, pch=19, ylim=c(0,0.4), ylab="Absolute Bias", xaxt="n")
abline(h=0)
points(abs_bias$freq_bias, pch=19,col="purple")
points(abs_bias$median_bias, pch=19,col="red")
points(abs_bias$mode1_bias, pch=19,col="blue")
points(abs_bias$mode0.1_bias, pch=19,col="green")
	abline(v=line_coords)
axis(1,1:24,rep(c(20,40,80,160),6))



plot(rmse$mean_bias, pch=19, ylim=c(0,0.4), ylab="Root mean squared error", xaxt="n")
abline(h=0)
points(rmse$freq_bias, pch=19,col="purple")
points(rmse$median_bias, pch=19,col="red")
points(rmse$mode1_bias, pch=19,col="blue")
points(rmse$mode0.1_bias, pch=19,col="green")
	abline(v=line_coords)
axis(1,1:24,rep(c(20,40,80,160),6))


plot(rmse2$mean_bias, pch=19, ylim=c(0,0.4), ylab="Root mean squared error", xaxt="n")
abline(h=0)
points(rmse2$freq_bias, pch=19,col="purple")
points(rmse2$median_bias, pch=19,col="red")
points(rmse2$mode1_bias, pch=19,col="blue")
points(rmse2$mode0.1_bias, pch=19,col="green")
	abline(v=line_coords)
axis(1,1:24,rep(c(20,40,80,160),6))
