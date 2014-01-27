	# (1=gill, 2=heart, 3=liver, 4=spleen, 5=kidney)	
	
params<- c("psi_z","psi_organ","p","sifAB","sifBC","sifAC","sifAD","sifBD","sifCD")


inits<- function(t){
	list(theta=runif(15),eta=runif(20), beta=runif(20)
	)}

	
	# BUNDLE DATA UP
	combos<-as.matrix(expand.grid(c(0,1),c(0,1),c(0,1),c(0,1)))
	dat<- list(yyy=analysis_dat,nfish=26,combos=combos)
	
	out <- bugs(data=dat,inits=inits,parameters=params,model = filename,
		n.chains = 3,n.iter = 500000,n.burnin = 100000, debug=TRUE,n.thin=20,	
		bugs.directory="C:/Documents and Settings/mcolvin/My Documents/WinBUGS14 - 1",	
		working.directory=getwd())	
		
	out <- bugs(data=dat,inits=inits,parameters=params,model = filename,
		n.chains = 3,n.iter = 5000,n.burnin = 1000, debug=TRUE,	bugs.directory="C:/Users/colvinm/Documents/WinBUGS14 - 1",
		working.directory=getwd())	
	out$mean

	save(file="./output/ms_mcmc.Rdata", out)

	
	
	
	
	
	
	
	
	### DIAGNOSTICS AND GOF
		out_gof <- as.mcmc.list(out)		
		plot(out_gof, ask=T)
		summary(out_gof) 		# Notice that rss and rss.sim have similar posteriors
								# indicating good fit
		##### Calculate a Bayesian P-value (aka Posterior predictive check)
		gofMat <- as.matrix(out_gof)
		PrGreater   <- sum(gofMat[,"rss"] > gofMat[,"rss_sim"]) / nrow(gofMat)
		# 0.48
		Pr2tail <- 2 * min(PrGreater, 1 - PrGreater)	# No evidence of lack of fit
		# 0.96	
		#### Generic convergence diagnostics
		gelman.diag(out_gof)
		gelman.plot(out_gof)
		geweke.diag(out_gof)
		geweke.plot(out_gof)
		crosscorr(out_gof)
		crosscorr.plot(out_gof)
		autocorr.diag(out_gof)
		autocorr.plot(out_gof)
	### END 
	