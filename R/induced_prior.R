
## prior on sd to induced prior on variance
n <- 100000
x<-seq(0,5,0.1)

par(mfrow=c(1,2))
## compare density and histogram
sd_prior <- runif(n,0,5)
sd_priorD <- dunif(x,0,5)
hist(sd_prior, freq=FALSE)
lines(sd_priorD~x,col="red")

##compare density and histogram of transformed values
var_prior <- sd_prior^2
var_priorD <- sd_priorD/(2*x) ## multiplying by the Jacobian of the transform
hist(var_prior, freq=FALSE, breaks=100)
lines(var_priorD~I(x^2),col="red")


## prior on variance to induced prior on sd

x<-seq(0,25,0.1)

par(mfrow=c(1,2))
## compare density and histogram
var_prior <- runif(n,0,25)
var_priorD <- dunif(x,0,25)
hist(var_prior, freq=FALSE)
lines(var_priorD~x,col="red")

##compare density and histogram of transformed values
sd_prior <- sqrt(var_prior)
sd_priorD <- var_priorD*(2*sqrt(x)) ## multiplying by the Jacobian of the transform (or divide by the derivative??)
hist(sd_prior, freq=FALSE, breaks=100)
lines(sd_priorD~I(sqrt(x)),col="red")




x<-seq(0,5,0.1)

par(mfrow=c(1,2))
## compare density and histogram
var_prior <- rnorm(n,0,2)
var_priorD <- dnorm(x,0,2)
hist(var_prior, freq=FALSE)
lines(var_priorD~x,col="red")

##compare density and histogram of transformed values
sd_prior <- sqrt(var_prior)
sd_priorD <- var_priorD*(2*sqrt(x))*2 ## multiplying by the Jacobian of the transform - x2 because of the folder distribution
hist(sd_prior, freq=FALSE, breaks=100)
lines(sd_priorD~I(sqrt(x)),col="red")




x<-seq(0,20,0.1)

par(mfrow=c(1,2))
## compare density and histogram
var_prior <- rcauchy(n,0,5)
var_priorD <- dcauchy(x,0,5)
hist(var_prior, freq=FALSE, xlim=c(0,20), breaks=1000000)
lines(var_priorD~x,col="red")

##compare density and histogram of transformed values
sd_prior <- sqrt(var_prior)
sd_priorD <- var_priorD*(2*sqrt(x))*2 ## multiplying by the Jacobian of the transform - x2 because of the folder distribution
hist(sd_prior, freq=FALSE, xlim=c(0,10), breaks=10000)
lines(sd_priorD~I(sqrt(x)),col="red")

## --------


## --------
par(mfrow=c(1,2))
x<-seq(0,5,0.01)
sd_prior <- dunif(x,0,5)
plot(sd_prior~x,type="l")

var_prior <- sd_prior/(2*x)
plot(var_prior~I(x^2),type="l", xlim=c(0,1))


sd_prior <- dnorm(x,0,1)
plot(sd_prior~x,type="l")

var_prior <- sd_prior/(2*x)
plot(var_prior~I(x^2),type="l", xlim=c(0,1))

## --------


par(mfrow=c(1,1))

x<-seq(0.1,5,0.01)


sd_prior_N <- dnorm(x,0,1)
var_prior_N <- sd_prior_N/(2*x)
plot(var_prior_N~I(x^2),type="l", xlim=c(0,1))

sd_prior_U <- dunif(x,0,5)
var_prior_U <- sd_prior_U/(2*x)
lines(var_prior_U~I(x^2),type="l", col="blue")

sd_prior_C <- dcauchy(x,0,2)
var_prior_C <- sd_prior_C/(2*x)
lines(var_prior_C~I(x^2),type="l", col="red")




# hist(plogis(rnorm(10000,0,0.01)))
# hist(rnorm(10000,0.5,0.01))
