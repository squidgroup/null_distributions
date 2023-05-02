
## from adapted from MCMCglmm::posterior.mode
post_mode <- function(x, adjust=1, cut=0, bw="SJ", ...) {
  dx <- density(x, adjust = adjust, cut=cut, bw=bw)
  dx$x[which.max(dx$y)]
}

stan_out <- function(model){
	out <- c(
		post_mode(extract(model)$sigma2_ID,adjust=0.1),
		post_mode(extract(model)$sigma2_ID,adjust=1), 
		median(extract(model)$sigma2_ID),
		summary(model)$summary["sigma2_ID",c(1,4,8,9)])
	names(out) <- c("mode0.1","mode1","median","mean","LCI","UCI","ESS")
	out	
}

gaussian_mods <- function(dat){
		data <- dat[,c("y","ID")]
		data$ID <- as.numeric(as.factor(data$ID))
		
		modF <- suppressMessages(lmer(y~1+(1|ID),data))

		stan_dat <- list(
			N = nrow(data),
			N_ID = length(unique(data[,"ID"])),
			y = data[,"y"],
			ID = as.numeric(as.factor(data[,"ID"])),
			cauchy_scale=2)

		stan_mod <- sampling(LMM_stan, data=stan_dat, chains=1,iter=5000, warmup=2000, pars=c("sigma2_ID","sigma2_E"), refresh=0)
		
		list(summary = c(freq=as.numeric(summary(modF)$varcor),stan_out(stan_mod)), data=data, post = as.data.frame(extract(stan_mod, permuted=FALSE)[,,1:2]))
}

bern_mods <- function(dat){
		data <- dat[,c("y","ID")]
		data$ID <- as.numeric(as.factor(data$ID))
		
		stan_dat <- list(
			N = nrow(data),
			N_ID = length(unique(data[,"ID"])),
			y = data[,"y"],
			ID = as.numeric(as.factor(data[,"ID"])))

		stan_mod <- sampling(bern_stan, data=stan_dat, chains=1,iter=5000, warmup=2000, pars=c("beta_0","sigma2_ID"), refresh=0)
		
		list(summary = stan_out(stan_mod), data=data, post = as.data.frame(extract(stan_mod, permuted=FALSE)[,,1:2]))
}

pois_mods <- function(dat){
		data <- dat[,c("y","ID")]
		data$ID <- as.numeric(as.factor(data$ID))
		
		stan_dat <- list(
			N = nrow(data),
			N_ID = length(unique(data[,"ID"])),
			y = data[,"y"],
			ID = as.numeric(as.factor(data[,"ID"])))

		stan_mod <- sampling(pois_stan, data=stan_dat, chains=1,iter=5000, warmup=2000, pars=c("beta_0","sigma2_ID","sigma2_E"), refresh=0)
		
		list(summary = stan_out(stan_mod), data=data, post = as.data.frame(extract(stan_mod, permuted=FALSE)[,,1:3]))
}

Stack <- function(x, col2stack, value.name="value", group.name="group", levels=NULL){
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

p_func <- function(actual,null) mean(actual<null)
p_func_LRT <- function(actual_p,null_p) mean(
	qchisq(actual_p*2, df=1, lower.tail=FALSE)<
	qchisq(null_p*2, df=1, lower.tail=FALSE))

stan_out_RR <- function(model){
	out <- cbind(MCMCglmm::posterior.mode(coda::as.mcmc(as.data.frame(extract(model)))[,c("sigma2_int","sigma2_slope")]), 
		apply(as.data.frame(extract(model))[,c("sigma2_int","sigma2_slope")],2,median),
		summary(model)$summary[c("sigma2_int","sigma2_slope"),c(1,4,8,9)])
	colnames(out) <- c("mode","median","mean","LCI","UCI","ESS")
	out	
}


print_func <- function(x) {
	x<-round(x,3)
	paste0(x["mean"]," (",x["mode1"],") \n [",x["LCI"],",",x["UCI"],"]")
}


## convert from expected to latent scale
exp2lat<- function(mean,cov){
	mean <- as.matrix(mean)
	cov <- as.matrix(cov)
	
	mean_out <- rep(NA,nrow(cov))
	for(i in 1:nrow(cov)) mean_out[i] <- log(mean[i]^2/sqrt(mean[i]^2+cov[i,i]))
	
	cov_out <- matrix(NA,nrow(cov),ncol(cov))
	for(i in 1:nrow(cov)){
		for(j in 1:ncol(cov)) cov_out[i,j] <- log(1 + cov[i,j]/(mean[i]*mean[j]))
	}
	return(list(mean=mean_out,cov=cov_out))
}

lat2exp<- function(mean,cov){
	mean <- as.matrix(mean)
	cov <- as.matrix(cov)
	
	mean_out <- rep(NA,nrow(cov))
	for(i in 1:nrow(cov)) mean_out[i] <- exp(mean[i]+ cov[i,i]/2)

	cov_out <- matrix(NA,nrow(cov),ncol(cov))
	for(i in 1:nrow(cov)){
		for(j in 1:ncol(cov)) cov_out[i,j] <- exp(mean[i]+mean[j] + (cov[i,i]+cov[j,j])/2)*(exp(cov[i,j])-1)
	}
	return(list(mean=mean_out,cov=cov_out))
}