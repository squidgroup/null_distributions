Modes2 <- function (x, adjust = 0.1, min.size = 0.1) {
    if (missing(x)) 
        stop("The x argument is required.")
    x <- as.vector(as.numeric(as.character(x)))
    x <- x[is.finite(x)]
    if (is.constant(x)) 
        return(list(modes = NA, mode.dens = NA, size = 1))
    length(density(x)$y)
    dens <- density(x, adjust=adjust)
    dens.y.diff <- dens$y[-1] - dens$y[-length(dens$y)]
    incr <- dens.y.diff
    incr[which(dens.y.diff > 0)] <- 1
    incr[which(dens.y.diff <= 0)] <- 0
    begin <- 1
    count <- 1
    for (i in 2:length(incr)) {
        if (incr[i] != incr[i - 1]) {
            count <- count + 1
            begin <- c(begin, i)
        }
    }
    begin <- c(begin, length(incr))
    size <- modes <- mode.dens <- rep(0, count/2)
    init <- 1
    sumdens <- sum(dens$y)
    if (incr[1] == 0) {
        size[1] <- sum(dens$y[1:begin[2]])/sumdens
        init <- 2
    }
    j <- init
    for (i in init:length(size)) {
        size[i] <- sum(dens$y[begin[j]:begin[j + 2]])/sumdens
        kde <- dens
        kde$x <- kde$x[begin[j]:begin[j + 2]]
        kde$y <- kde$y[begin[j]:begin[j + 2]]
        modes[i] <- kde$x[kde$y == max(kde$y)][1]
        mode.dens[i] <- kde$y[kde$y == max(kde$y)][1]
        j <- j + 2
    }
    size <- size[order(mode.dens, decreasing = TRUE)]
    modes <- modes[order(mode.dens, decreasing = TRUE)]
    mode.dens <- mode.dens[order(mode.dens, decreasing = TRUE)]
    if (any(size < 0.1)) {
        modes <- modes[-which(size < min.size)]
        mode.dens <- mode.dens[-which(size < min.size)]
        size <- size[-which(size < min.size)]
    }
    if (sum(size) > 1) 
        size <- size/sum(size)
    return(list(modes = modes, mode.dens = mode.dens, size = size))
}




load(file="~/Dropbox/0_Presentations/20220207_SEED/examples.Rdata")

data<-res[[2]]$post
hist(data, breaks=150, freq=FALSE)
library(multimode)
library(MCMCglmm)
modetest(data)

locmodes(data,mod0=1,display=TRUE,lowsup=0.000001)
mean(data)

library(LaplacesDemon)

lapply(res, function(x){
	data<-x$post
	c(unimodal=is.unimodal(data),modes=Modes(data)$modes, postmode=posterior.mode(as.mcmc(matrix(data))), mean=mean(data))

})

lapply(res, function(x){
	data<-x$post
	c(unimodal=is.unimodal(data),modes=Modes2(data, adjust=0.1)$modes, postmode=posterior.mode(as.mcmc(matrix(data)), adjust=0.1), mean=mean(data))

})

lapply(res, function(x){
	data<-x$post
	Modes2(data, adjust=1)$modes
})
sapply(res, function(x){
	mean(x$post)
})
sapply(res, function(x){
	posterior.mode(as.mcmc(matrix(x$post)))
})
sapply(res, function(x){
	posterior.mode(as.mcmc(matrix(x$post)),adjust=1)
})



data<-res[[3]]$post
dx <- density(data, adjust = 0.2, bw="SJ")
#kernel="biweight"
dx$x[which.max(dx$y)]
plot(dx$y,type="l")
hist(data,breaks=50)

par(mfcol=c(2,4))
for(i in c(1,2,5,10)){
	dx <- density(data, adjust = 1/i, bw="SJ")
#kernel="biweight"
dx$x[which.max(dx$y)]
plot(dx$y,type="l")
hist(data,breaks=10*i)

}


data<-res[[3]]$post

par(mfcol=c(2,4))
for(i in c(1,2,5,100)){
	dx <- density(data, cut=0, adjust = 1/i, bw="SJ")
#kernel="biweight"
	pm <- dx$x[which.max(dx$y)]
	print(pm)
	plot(dx$y~dx$x,type="l", main=paste0("adjust=",1/i), xlab="Posterior samples",ylab="Density")
	abline(v=pm, col="red",lwd=1)
	abline(v=0, col="grey",lwd=1)
	hist(data,breaks=10*i, main=paste0("breaks=",10*i),xlab="Posterior samples")

}

length(data)




density(data, adjust = 10)$y


abline(v=posterior.mode(as.mcmc(matrix(data)), adjust = 1), col="blue")
abline(v=posterior.mode(as.mcmc(matrix(data))), col="red")

plot(density(data, adjust = 0.1)$y)