
rm(list=ls())

devtools::load_all("~/github/squidSim/R")

options(mc.cores = parallel::detectCores())

library(lme4)
library(parallel)
library(rstan)
rstan_options("auto_write" = TRUE)
library("beeswarm")

wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/functions.R"))


p_mode <- function(x, adjust, ...) {
  dx <- density(x, adjust = adjust, cut=0, bw="SJ")
  dx$x[which.max(dx$y)]
}

post_summary <- function(post){
  c(
    mode0.1 = p_mode(post,adjust=0.1),
    mode1 = p_mode(post,adjust=1),
    median = median(post),
    mean = mean(post) 
  )
}


prior_sim <- function(N_pop, ICC, N_group, N_within, mc.cores=8){
#N_pop=100; ICC=0.2; N_group=80; N_within=2; mc.cores=8
  squid_dat<-simulate_population(
    data_structure= make_structure(paste0("ID(",N_group,")"),repeat_obs=N_within),
    parameters= list( ID=if(ICC==0){list(vcov=0.01, beta=0)}else{list(vcov=ICC)}, residual=list(vcov=1-ICC)),
    n_pop=N_pop
    )
  
  dat <- get_population_data(squid_dat, list=TRUE)
  dist_M <- mclapply(1:N_pop,function(i){
    #i=1
    
    LMM_stan <- stan_model(file = paste0(wd,"stan/simple_LMM.stan"))

    stan_dat1 <- list(
        N = nrow(dat[[i]]),
        N_ID = length(unique(dat[[i]][,ID])),
        y = dat[[i]][,y],
        ID = dat[[i]][,ID],
        cauchy_scale=2)

    stan_mod1 <- sampling(LMM_stan, data=stan_dat1, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID"), refresh=0)
    cat(i, " ")
    extract(stan_mod1)["sigma2_ID"][[1]]
    
  },mc.cores=mc.cores)
  return(dist_M)
}

out<-prior_sim(N_pop=1000, ICC=0.2, N_group=80, N_within=2, mc.cores=8)

apply(sapply(out,post_summary),1,function(x) sum(x>0.015))

pairs(t(sapply(out,post_summary)))

out2<-do.call(rbind,lapply(out,function(x)data.frame(estimate=post_summary(x),type=c("mode0.1","mode1","median","mean")))
)

  beeswarm(estimate~type,out2, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap", ylab="Estimate")
points(with(out2,tapply(estimate,type,mean)), pch=19,col="red")
abline(h=0.2,col="blue")

with(out2,tapply(estimate,type,function(x)mean(x-0.2)))
with(out2,tapply(estimate,type,function(x)sd(x-0.2)))




out0<-prior_sim(N_pop=1000, ICC=0, N_group=80, N_within=2, mc.cores=8)

apply(sapply(out0,post_summary),1,function(x) sum(x>0.01))

pairs(t(sapply(out0,post_summary)))

out0_2<-do.call(rbind,lapply(out0,function(x)data.frame(estimate=post_summary(x),type=c("mode0.1","mode1","median","mean")))
)

  beeswarm(estimate~type,out0_2, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap", ylab="Estimate")
points(with(out0_2,tapply(estimate,type,median)), pch=19,col="red")
abline(h=0,col="blue")

with(out0_2,tapply(estimate,type,mean))
with(out0_2,tapply(estimate,type,function(x)sd(x-0.2)))




out0.1<-prior_sim(N_pop=1000, ICC=0.1, N_group=80, N_within=2, mc.cores=8)

apply(sapply(out0.1,post_summary),1,function(x) sum(x>0.01))

pairs(t(sapply(out0.1,post_summary)))

out0.1_2<-do.call(rbind,lapply(out0.1,function(x)data.frame(estimate=post_summary(x),type=c("mode0.1","mode1","median","mean")))
)

  beeswarm(estimate~type,out0.1_2, pch=19, cex=0.5, col=alpha(1,0.3),method = "compactswarm",corral="wrap", ylab="Estimate")
points(with(out0.1_2,tapply(estimate,type,median)), pch=19,col="red")
abline(h=0.1,col="blue")

with(out0.1_2,tapply(estimate,type,function(x)mean(x-0.1)))
with(out0.1_2,tapply(estimate,type,function(x)sd(x-0.1)))
with(out0.1_2,tapply(estimate,type,sd))

save(out,out0.1,out0,file=paste0(wd,"Data/Intermediate/modes.Rdata"))
