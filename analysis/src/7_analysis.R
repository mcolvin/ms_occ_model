	# (1=gill, 2=heart, 3=liver, 4=spleen, 5=kidney)
	######################################################################
	# PM [done] (1=gill, 2=heart, 3=liver, 4=spleen, 5=kidney)
	######################################################################
	param <- c('pred_organ_psi','pred_organ_p','pred.fish','prev','prev_org','rss','rss_sim',"y.sim")
	dat<- list(y=pm_dat,nfish=nfish,no.orgs=5)
	lat.org<- matrix(1,nrow=nfish, ncol=5)
	lat.org[,c(1:3)]<- 0
	inits_base <- function(t) 
		{
		list(lat.org=lat.org,alpha0=0 ,beta=runif(5), gamma=runif(5),lat.fish = rep(1, nfish))	
		list(lat.org=lat.org,alpha0=0 ,beta=runif(5), gamma=runif(5),lat.fish = rep(1, nfish))	
		list(lat.org=lat.org,alpha0=0 ,beta=runif(5), gamma=runif(5),lat.fish = rep(1, nfish))	
		}	
	out <- bugs(data=dat,inits=inits_base,parameters=param,model = filename,
		n.chains = 3,n.iter = 130000,n.burnin = 30000,	debug=TRUE,	bugs.directory=bugsdir,n.thin=20)
	write.csv(out$summary,"./dat/pm_out.csv")	
	ests<- list(pm= out$summary)
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
	