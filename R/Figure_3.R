
rm(list=ls())

library(scales)
library(beeswarm)

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
ICCs <- c(0,0.2)
{
par(mfrow=c(2,1),mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))
for(j in seq_along(ICCs)){
	plot_dat <- subset(p_perm, ICC==ICCs[j])
	beeswarm(mean~N_group + N_within,plot_dat, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap",xlab="Number of Groups", labels=rep(c(20,40,80,160),3), ylab=bquote(P[perm]),  ylim=c(0,1))

# mtext(paste0(letters[j],")"),2,padj=-5, las=1, line=2)

}
}
dev.off()

{
	power_dat <- aggregate(mean~N_group+ICC,p_perm,function(x) sum(x<=0.05)/100)
par(mfrow=c(2,1),mar=c(4,4,2,1), cex.axis=0.75, mgp=c(2,0.5,0))
plot(NA,xlim=c(0,160),xlab="Number of Groups", ylab="Power",  ylim=c(0,1))
abline(h=0.05, col="grey")
abline(h=c(0.025,0.075), col="grey")
for(j in seq_along(ICCs)){
	points(mean~N_group,power_dat, subset=ICC==ICCs[j], col=j, pch=19)
	lines(mean~N_group,power_dat, subset=ICC==ICCs[j], col=j, lwd=2)
}
}


par(mfrow=c(2,2))
plot(mean~mode,p_perm)
plot(mean~median,p_perm)
plot(mean~freq,p_perm)
plot(mode~freq,p_perm)
cor(p_perm$mean,p_perm$mode)
