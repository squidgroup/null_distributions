load(file="~/Dropbox/0_Presentations/20220207_SEED/examples.Rdata")
load(file="/Users/joelpick/github/bayes_perm/Data/Intermediate/figure1_data.Rdata")
ls()
figure1_data.Rdata

data<-res[[3]]$post

par(mfcol=c(2,5))
for(i in c(0.5,1,2,5,10)){
	dx <- density(data, cut=0, adjust = 1/i, bw="SJ")
#kernel="biweight"
	pm <- dx$x[which.max(dx$y)]
	# print(pm)
	plot(dx$y~dx$x,type="l", main=paste0("adjust=",1/i), xlab="Posterior samples",ylab="Density")
	if(i%in%c(1,10)) box("figure",col="blue",lwd=4)
	abline(v=0, col="grey",lwd=1)
	abline(v=pm, col="red",lwd=1)
	text(0.475,0.75*max(dx$y),round(pm,3))
	hist(data,breaks=seq(min(dx$x),max(dx$x),length.out=10*i+1), main=paste0("breaks=",10*i),xlab="Posterior samples")
}