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



stan_out <- function(model){
	out <- c(MCMCglmm::posterior.mode(coda::as.mcmc(as.data.frame(extract(model)))[,"sigma2_ID"]), 
		median(extract(model)$sigma2_ID),
		summary(model)$summary["sigma2_ID",c(1,4,8,9)])
	names(out) <- c("mode","median","mean","LCI","UCI","ESS")
	out	
}
stan_out_RR <- function(model){
	out <- cbind(MCMCglmm::posterior.mode(coda::as.mcmc(as.data.frame(extract(model)))[,c("sigma2_int","sigma2_slope")]), 
		apply(as.data.frame(extract(model))[,c("sigma2_int","sigma2_slope")],2,median),
		summary(model)$summary[c("sigma2_int","sigma2_slope"),c(1,4,8,9)])
	colnames(out) <- c("mode","median","mean","LCI","UCI","ESS")
	out	
}


print_func <- function(x) {
	x<-round(x,3)
	paste0(x["mean"]," (",x["mode"],") \n [",x["LCI"],",",x["UCI"],"]")
}