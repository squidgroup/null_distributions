rm(list=ls())

wd <- "~/github/bayes_perm/"

load(file="~/Dropbox/0_Presentations/20220207_SEED/examples.Rdata")

data<-res[[3]]$post

setEPS()
pdf(paste0(wd,"Figures/Fig2_mode.pdf"), height=6, width=10)

v <- c(1,2,5,10)
{
par(mfcol=c(2,4),mgp=c(2.5,1,0), mar=c(4,4,3,1))
for(i in 1:4){
	dx <- density(data, cut=0, adjust = 1/v[i], bw="SJ")
#kernel="biweight"
	pm <- dx$x[which.max(dx$y)]
	# print(pm)
	plot(dx$y~dx$x,type="l", main=paste0("adjust=",1/v[i]), xlab="Posterior samples",ylab="Density")
	# if(i%in%c(1,10)) box("figure",col="blue",lwd=4)
	abline(v=0, col="grey",lwd=1)
	abline(v=pm, col="red",lwd=1)
	text(0.475,0.75*max(dx$y),round(pm,3))
	mtext(paste0(letters[i],")"),3,line=1,adj=0)
	hist(data,breaks=seq(min(dx$x),max(dx$x),length.out=10*v[i]+1), main=paste0("breaks=",10*v[i]),xlab="Posterior samples")
	mtext(paste0(letters[i+4],")"),3,line=1,adj=0)
}}

dev.off()
