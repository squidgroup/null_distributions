#Functions
ran_y<-function(x, n_perm){
  dran<-list()
  for(i in 1:n_perm){
    x$y<-sample(x$y,nrow(x))
    dran[[i]]<-x
  }
  dran
}

ran_ID<-function(x, n_perm){
  dran<-list()
  for(i in 1:n_perm){
    x$individual<-sample(x$individual,nrow(x))
    dran[[i]]<-x
  }
  dran
}

ran_x<-function(x, n_perm){
  dran<-list()
  for(i in 1:n_perm){
    x$environment<-sample(x$environment,nrow(x))
    dran[[i]]<-x
  }
  dran
}

ran_x_ID<-function(x, n_perm){
  dran<-list()
  x_ind<-list()
  n.ind<-max(x$individual)
  for(i in 1:n_perm){
    for(j in 1:n.ind){
      dI<-x[x$individual==j,]
      x_ind[[j]]<-sample(dI$environment,nrow(dI))
    }
    x$environment<-unlist(x_ind)
    dran[[i]]<-x
  }
  dran
}

##Function to analyze the simulations
mod<-function(x){
  
  stan_data = list(n = nrow(x),
                   n_ind=max(x$individual),
                   individual=x$individual,
                   x=x$environment,
                   y=x$y)
  
  
  params<-c("beta_0", "beta_1", "sigma2_I2")
  ni <- 1500 ##Number of iterations
  nc <- 1  ##Number of chains
  nt <- 1 ## Thinning interval
  nb <- 500 ## Number of iterations to discard
  
  stan_data = list(n = nrow(x),
                   n_ind=max(x$individual),
                   individual=x$individual,
                   x=x$environment,
                   y=x$y)
  
  md <- stan("Stan/RR_example.stan", data = stan_data, pars = params, chains = nc, iter = ni,  warmup = nb, thin = nt, cores=1)
  round(summary(md)$summary[3,1],3)
}
