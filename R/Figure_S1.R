rm(list=ls())

library(scales)
wd <- "~/github/bayes_perm/"

source(paste0(wd,"R/functions.R"))

files <- list.files(paste0(wd,"Data/Intermediate"))
results <- files[grep("sims",files)]
list_names <- gsub("sims_|.Rdata","",results)


all <- list()
for(i in 1:length(results)){
	load(paste0(wd,"Data/Intermediate/",results[i]))
	all[[list_names[i]]] <- sim_dat
	rm("sim_dat")
}



p_perm <- as.data.frame(do.call(rbind,lapply(all, function(z){
	t(sapply(z, function(x) c(x$param,
		freq = p_func(x$actual["freq"],x$null[,"freq"]),
		LRT=p_func_LRT(x$actual["LRT"],x$null[,"LRT"])
		)))
})))
p_perm<-subset(p_perm,ICC!=0.8 & N_within!=16)
head(p_perm)

plot(freq~LRT,p_perm, pch=19, cex=0.5, col=alpha(1,0.2))
cor(p_perm$freq,p_perm$LRT)
hist(p_perm$LRT)


# devtools::install("~/github/squidSim")
library(squidSim)
library(lme4)
squid_dat <- simulate_population(
	data_structure = make_structure("ID(100)",repeat_obs=4),
	parameters = list(
		ID=list(vcov=0.2),
		residual = list(vcov=0.3)
		)
	)

dat <- get_population_data(squid_dat)

		modF <- suppressMessages(lmer(y~1+(1|ID),dat, REML=TRUE))
		as.numeric(summary(modF)$varcor)
		modF <- suppressMessages(lmer(y~1+(1|ID),dat, REML=FALSE))
		as.numeric(summary(modF)$varcor)
		modF_null <- lm(y~1,dat)
		LRT <- anova(modF,modF_null)$P[2]/2

		anova(modF,modF_null)$Chi[2]


		qchisq(LRT*2, df=1, lower.tail=FALSE)