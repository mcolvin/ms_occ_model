model<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		## j indexes number of organs
		for (j in 1:no.orgs) 
			{
			## conditional organ occupancy
			logit(organ.psi0[i,j]) <- beta[j] 
			organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
			## latent organ occupancy variable
			lat.org[i,j] ~ dbern(organ.psi[i,j])
			lat.org2[i,j] <-  lat.org[i,j] + 1

			for (h in 1:3) 			## h indexes number of repeated samples (visits) per organ
				{     
				## detection probability 
				logit(theta[i,j,h,2]) <- gamma[j] 
				theta2[i,j,h,1] <- 0
				theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
				y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
				}
			}## end no.spc
	}  ## end n.obs

	
# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy

for(o in 1:5)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict occupancy
	logit(pred_organ_p[o]) <- gamma[o] # predict detection
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
beta0 ~ dnorm(0, 0.37)
gamma0 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	#beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	beta[z] ~ dnorm(beta0, beta.tau)
	gamma[z] ~ dnorm(gamma0, gamma.tau)
	#gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}

## Priors for the random effects from the precision (tau)
beta.tau <- pow(beta.ss,-2)
beta.ss ~ dunif(0,6)

gamma.tau <- pow(gamma.ss,-2)
gamma.ss ~ dunif(0,6)

}
filename <- paste(getwd(),"base.bug", sep="/")
write.model(model, filename)

model_sim<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		## j indexes number of organs
		for (j in 1:no.orgs) 
			{
			## conditional organ occupancy
			logit(organ.psi0[i,j]) <- beta[j] 
			organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
			## latent organ occupancy variable
			lat.org[i,j] ~ dbern(organ.psi[i,j])
			lat.org2[i,j] <-  lat.org[i,j] + 1

			for (h in 1:reps) 			## h indexes number of repeated samples (visits) per organ
				{     
				## detection probability 
				logit(theta[i,j,h,2]) <- gamma[j] 
				theta2[i,j,h,1] <- 0
				theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
				y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
				}
			}## end no.spc
	}  ## end n.obs

	
# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy

for(o in 1:5)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict occupancy
	logit(pred_organ_p[o]) <- gamma[o] # predict detection
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
beta0 ~ dnorm(0, 0.37)
gamma0 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:5) 
	{
	#beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	beta[z] ~ dnorm(beta0, beta.tau)
	gamma[z] ~ dnorm(gamma0, gamma.tau)
	#gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}

## Priors for the random effects from the precision (tau)
beta.tau <- pow(beta.ss,-2)
beta.ss ~ dunif(0,6)

gamma.tau <- pow(gamma.ss,-2)
gamma.ss ~ dunif(0,6)

}
filename_sim <- paste(getwd(),"base_sim.bug", sep="/")
write.model(model_sim, filename_sim)

model_NS_M01<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		## j indexes number of organs
		for (j in 1:no.orgs) 
			{
			gill_occ[i,j]<- lat.org[i,j]*X[i,j,1]	# ARE THE GILLS OCCUPIED MATRIX		
			gill_occ_mat[i,j]<- min(sum(gill_occ[i,]),1)# FILLS ROWS WITH 1 IF GILLS ARE OCCUPIED, THEN X[i,j,2] DOES THE RIGHT MULTIPLICATION
			## conditional organ occupancy
			logit(organ.psi0[i,j]) <- beta[j]  + beta2*gill_occ_mat[i,j]*X[i,j,2]
			### Ignore this it helps keeps the program in track, note multiplied by 
			#     fish level latent variable
			organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
			## latent organ occupancy variable
			lat.org[i,j] ~ dbern(organ.psi[i,j])
			lat.org2[i,j] <-  lat.org[i,j] + 1
			## h indexes number of repeated samples (visits) per organ
			for (h in 1:3) 
				{     
				## detection probability 
				logit(theta[i,j,h,2]) <- gamma[j] 
				theta2[i,j,h,1] <- 0
				theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
				y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
				}
			}## end no.spc
	}  ## end n.obs

	
# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy
for(o in 1:5)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict occupancy
	logit(pred_organ_p[o]) <- gamma[o] # predict detection
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
#beta0 ~ dnorm(0, 0.37)
#gamma0 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	#beta[z] ~ dnorm(beta0, beta.tau)
	#gamma[z] ~ dnorm(gamma0, gamma.tau)
	gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}
beta2 ~  dnorm(0, 0.37)
## Priors for the random effects from the precision (tau)
#beta.tau <- pow(beta.ss,-2)
#beta.ss ~ dunif(0,6)

#gamma.tau <- pow(gamma.ss,-2)
#gamma.ss ~ dunif(0,6)

}
filename_NS_M01 <- paste(getwd(),"base_NS_01.bug", sep="/")
write.model(model_NS_M01, filename_NS_M01)


model_NS_M01_cut<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		## j indexes number of organs
		for (j in 1:no.orgs) 
			{
			#gill_occ[i,j]<- lat.org[i,j]*X[i,j,1]	# ARE THE GILLS OCCUPIED MATRIX		
			gill_occ[i,j]<- lat.org.cut[i,j]*X[i,j,1]	# ARE THE GILLS OCCUPIED MATRIX		
			gill_occ_mat[i,j]<- min(sum(gill_occ[i,]),1)# FILLS ROWS WITH 1 IF GILLS ARE OCCUPIED, THEN X[i,j,2] DOES THE RIGHT MULTIPLICATION
			## conditional organ occupancy
			logit(organ.psi0[i,j]) <- beta[j]  + beta2*gill_occ_mat[i,j]*X[i,j,2]
			### Ignore this it helps keeps the program in track, note multiplied by 
			#     fish level latent variable
			organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
			## latent organ occupancy variable
			lat.org[i,j] ~ dbern(organ.psi[i,j])
			lat.org2[i,j] <-  lat.org[i,j] + 1
			lat.org.cut[i,j]<- cut(lat.org[i,j])
			## h indexes number of repeated samples (visits) per organ
			for (h in 1:3) 
				{     
				## detection probability 
				logit(theta[i,j,h,2]) <- gamma[j] 
				theta2[i,j,h,1] <- 0
				theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
				y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
				}
			}## end no.spc
	}  ## end n.obs

	
# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy
for(o in 1:5)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict occupancy
	logit(pred_organ_p[o]) <- gamma[o] # predict detection
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
beta0 ~ dnorm(0, 0.37)
gamma0 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	#beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	beta[z] ~ dnorm(beta0, beta.tau)
	gamma[z] ~ dnorm(gamma0, gamma.tau)
	#gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}
beta2 ~  dnorm(0, 0.37)
## Priors for the random effects from the precision (tau)
beta.tau <- pow(beta.ss,-2)
beta.ss ~ dunif(0,6)

gamma.tau <- pow(gamma.ss,-2)
gamma.ss ~ dunif(0,6)

}
filename_NS_M01_cut <- paste(getwd(),"base_NS_01_cut.bug", sep="/")
write.model(model_NS_M01_cut, filename_NS_M01_cut)

model_M01_sex<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0 + alpha1*X[i,1]
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		## j indexes number of organs
		for (j in 1:no.orgs) 
			{
			## conditional organ occupancy
			logit(organ.psi0[i,j]) <- beta[j]  
			### Ignore this it helps keeps the program in track, note multiplied by 
			#     fish level latent variable
			organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
			## latent organ occupancy variable
			lat.org[i,j] ~ dbern(organ.psi[i,j])
			lat.org2[i,j] <-  lat.org[i,j] + 1
			lat.org.cut[i,j]<- cut(lat.org[i,j])
			## h indexes number of repeated samples (visits) per organ
			for (h in 1:3) 
				{     
				## detection probability 
				logit(theta[i,j,h,2]) <- gamma[j] 
				theta2[i,j,h,1] <- 0
				theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
				y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
				}
			}## end no.spc
	}  ## end n.obs

	
# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy
for(o in 1:5)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict occupancy
	logit(pred_organ_p[o]) <- gamma[o] # predict detection
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
alpha1 ~ dnorm(0, 0.37)
beta0 ~ dnorm(0, 0.37)
gamma0 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	#beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	beta[z] ~ dnorm(beta0, beta.tau)
	gamma[z] ~ dnorm(gamma0, gamma.tau)
	#gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}
beta2 ~  dnorm(0, 0.37)
## Priors for the random effects from the precision (tau)
beta.tau <- pow(beta.ss,-2)
beta.ss ~ dunif(0,6)

gamma.tau <- pow(gamma.ss,-2)
gamma.ss ~ dunif(0,6)

}
filename_M01_sex <- paste(getwd(),"model_M01_sex.bug", sep="/")
write.model(model_M01_sex, filename_M01_sex)

model_M01_len<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0 + alpha1*X[i,2]
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		## j indexes number of organs
		for (j in 1:no.orgs) 
			{
			## conditional organ occupancy
			logit(organ.psi0[i,j]) <- beta[j]  
			### Ignore this it helps keeps the program in track, note multiplied by 
			#     fish level latent variable
			organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
			## latent organ occupancy variable
			lat.org[i,j] ~ dbern(organ.psi[i,j])
			lat.org2[i,j] <-  lat.org[i,j] + 1
			lat.org.cut[i,j]<- cut(lat.org[i,j])
			## h indexes number of repeated samples (visits) per organ
			for (h in 1:3) 
				{     
				## detection probability 
				logit(theta[i,j,h,2]) <- gamma[j] 
				theta2[i,j,h,1] <- 0
				theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
				y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
				}
			}## end no.spc
	}  ## end n.obs

	
# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy
for(o in 1:5)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict occupancy
	logit(pred_organ_p[o]) <- gamma[o] # predict detection
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
alpha1 ~ dnorm(0, 0.37)
beta0 ~ dnorm(0, 0.37)
gamma0 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	#beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	beta[z] ~ dnorm(beta0, beta.tau)
	gamma[z] ~ dnorm(gamma0, gamma.tau)
	#gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}
beta2 ~  dnorm(0, 0.37)
## Priors for the random effects from the precision (tau)
beta.tau <- pow(beta.ss,-2)
beta.ss ~ dunif(0,6)

gamma.tau <- pow(gamma.ss,-2)
gamma.ss ~ dunif(0,6)

}
filename_M01_len <- paste(getwd(),"model_M01_len.bug", sep="/")
write.model(model_M01_len, filename_M01_len)

model_NS_M02<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		## j indexes number of organs
		for (j in 1:no.orgs) 
			{
			kid_occ[i,j]<- lat.org[i,j]*X[i,j,1]	# KIDNEY OCCUPANCY MATRIX		
			kid_occ_mat[i,j]<- min(sum(kid_occ[i,]),1)  # FILLS ROWS WITH 1 IF KIDNEYS ARE OCCUPIED, THEN X[i,j,2] DOES THE RIGHT MULTIPLICATION
			## conditional organ occupancy
			logit(organ.psi0[i,j]) <- beta[j]  + beta2[j]*kid_occ_mat[i,j]*X[i,j,2]
			### Ignore this it helps keeps the program in track, note multiplied by 
			#     fish level latent variable
			organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
			## latent organ occupancy variable
			lat.org[i,j] ~ dbern(organ.psi[i,j])
			lat.org2[i,j] <-  lat.org[i,j] + 1
			## h indexes number of repeated samples (visits) per organ
			for (h in 1:3) 
				{     
				## detection probability 
				logit(theta[i,j,h,2]) <- gamma[j] 
				theta2[i,j,h,1] <- 0
				theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
				y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
				}
			}## end no.spc
	}  ## end n.obs

	
# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy
for(o in 1:5)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict occupancy
	logit(pred_organ_p[o]) <- gamma[o] # predict detection
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
#beta0 ~ dnorm(0, 0.37)
#gamma0 ~ dnorm(0, 0.37)
beta2[z] ~  dnorm(0, 0.37)
### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	beta[z] ~  dnorm(0, 0.37) # organ level occupancy

	#beta[z] ~ dnorm(beta0, beta.tau)
	#gamma[z] ~ dnorm(gamma0, gamma.tau)
	gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}

## Priors for the random effects from the precision (tau)
#beta.tau <- pow(beta.ss,-2)
#beta.ss ~ dunif(0,6)

#gamma.tau <- pow(gamma.ss,-2)
#gamma.ss ~ dunif(0,6)

}
filename_NS_M02 <- paste(getwd(),"base_NS_02.bug", sep="/")
write.model(model_NS_M02, filename_NS_M02)

model_NS_M02a<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		## j indexes number of organs
		for (j in 1:no.orgs) 
			{
			kid_occ[i,j]<- lat.org[i,j]*X[i,j,1]	# KIDNEY OCCUPANCY MATRIX		
			kid_occ_mat[i,j]<- min(sum(kid_occ[i,]),1)  # FILLS ROWS WITH 1 IF KIDNEYS ARE OCCUPIED, THEN X[i,j,2] DOES THE RIGHT MULTIPLICATION
			## conditional organ occupancy
			logit(organ.psi0[i,j]) <- beta[j]+ beta2[j]*kid_occ_mat[i,j]
			### Ignore this it helps keeps the program in track, note multiplied by 
			# fish level latent variable
			organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
			## latent organ occupancy variable
			lat.org[i,j] ~ dbern(organ.psi[i,j])
			lat.org2[i,j] <-  lat.org[i,j] + 1
			## h indexes number of repeated samples (visits) per organ
			for (h in 1:3) 
				{     
				## detection probability 
				logit(theta[i,j,h,2]) <- gamma[j] + gamma1*a[i,j]
				theta2[i,j,h,1] <- 0
				theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
				y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
				}
			}## end no.spc
	}  ## end n.obs

	
# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy
for(o in 1:5)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict occupancy
	logit(pred_organ_p[o]) <- gamma[o] # predict detection
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
gamma1 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	beta2[z] ~  dnorm(0, 0.37)
	gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}
}
filename_NS_M02a <- paste(getwd(),"base_NS_02a.bug", sep="/")
write.model(model_NS_M02a, filename_NS_M02a)

# ns counts to predict det probability
model_NS_M02b<- function ()
	{ 
	for (i in 1:nfish) 
		{
		logit(fish.psi0[i]) <- b*XXX[i,5] + alpha0	
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		lat.fish[i] ~ dbern(fish.psi[i])
		for (j in 1:no.orgs) 
			{
			logit(organ.psi0[i,j]) <- beta[j] 
			organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
				lat.org[i,j] ~ dbern(organ.psi[i,j])
			lat.org2[i,j] <-  lat.org[i,j] + 1
			for (h in 1:3) 
				{     
				logit(theta[i,j,h,2]) <- gamma[j]
				theta2[i,j,h,1] <- 0
				theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
				y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
				}
			}
		}  
	alpha0 ~ dnorm(0, 0.37)
	b ~ dnorm(0, 0.37)
	for (z in 1:no.orgs) 
		{
		beta[z] ~  dnorm(0, 0.37) # organ level occupancy
		gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
		}
	}
filename_NS_M02b <- paste(getwd(),"base_NS_02b.bug", sep="/")
write.model(model_NS_M02b, filename_NS_M02b)

model_NS_M05<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		## j indexes number of organs
		for (j in 1:no.orgs) 
			{
			kid_occ[i,j]<- lat.org[i,j]*X[i,j,1]	# KIDNEY OCCUPANCY MATRIX		
			kid_occ_mat[i,j]<- step(sum(kid_occ[i,])-1)  # FILLS ROWS WITH 1 IF KIDNEYS ARE OCCUPIED, THEN X[i,j,2] DOES THE RIGHT MULTIPLICATION
			## conditional organ occupancy
			logit(organ.psi0[i,j]) <- beta[j]  + beta2*kid_occ_mat[i,j]*X[i,j,2]
			### Ignore this it helps keeps the program in track, note multiplied by 
			#     fish level latent variable
			organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
			## latent organ occupancy variable
			lat.org[i,j] ~ dbern(organ.psi[i,j])
			lat.org2[i,j] <-  lat.org[i,j] + 1
			## h indexes number of repeated samples (visits) per organ
			for (h in 1:3) 
				{     
				## detection probability 
				logit(theta[i,j,h,2]) <- gamma[j] 
				theta2[i,j,h,1] <- 0
				theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
				y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
				}
			}## end no.spc
	}  ## end n.obs

	
# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy
for(o in 1:5)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict occupancy
	logit(pred_organ_p[o]) <- gamma[o] # predict detection
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
#beta0 ~ dnorm(0, 0.37)
#gamma0 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	#beta[z] ~ dnorm(beta0, beta.tau)
	#gamma[z] ~ dnorm(gamma0, gamma.tau)
	gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}
beta2 ~  dnorm(0, 0.37)
## Priors for the random effects from the precision (tau)
#beta.tau <- pow(beta.ss,-2)
#beta.ss ~ dunif(0,6)
#gamma.tau <- pow(gamma.ss,-2)
#gamma.ss ~ dunif(0,6)
}
filename_NS_M05 <- paste(getwd(),"base_NS_05.bug", sep="/")
write.model(model_NS_M05, filename_NS_M05)

#P. minibicornis #########################################################
model_pm<- function ()
		{
		# i indexes total number of fish 
		for (i in 1:nfish) 
			{
		   ###  Fish-level model
			logit(fish.psi0[i]) <- alpha0
			## Ignore this it helps keeps the program in track min-max trick
			fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
			# latent fish level variable
			lat.fish[i] ~ dbern(fish.psi[i])
			## j indexes number of organs, non-kidney
			for (j in 1:5) 
				{
				## conditional organ occupancy
				logit(organ.psi0[i,j]) <- beta[j]
				### Ignore this it helps keeps the program in track, note multiplied by 
				#     fish level latent variable
				organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
				## latent organ occupancy variable
				lat.org[i,j] ~ dbern(organ.psi[i,j])
				lat.org2[i,j] <-  lat.org[i,j] + 1				
				}
			for(j in 6:7)
				{# glom and tubules
				## conditional on kidney occupancy
				logit(organ.psi0[i,j]) <- beta[j]
				### Ignore this it helps keeps the program in track, note multiplied by 
				#     kidney level latent variable
				organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.org[i,5]
				lat.org[i,j] ~ dbern(organ.psi[i,j])
				lat.org2[i,j] <-  lat.org[i,j] + 1			
				} # end j
			} # end i
		
			# detection model			
			for(i in 1:nfish)
				{
				for(j in 1:4)
					{
					for (h in 1:3) 
						{ # detection model  h indexes number of repeated samples (visits) per organ
						logit(theta[i,j,h,2]) <- gamma[j] 
						theta2[i,j,h,1] <- 0
						theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
						y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
						}	## end h
					} ## end j
				for(j in 6:7)
					{
					for (h in 1:3) 
						{ # detection model  h indexes number of repeated samples (visits) per organ
						logit(theta[i,j,h,2]) <- gamma[j] 
						theta2[i,j,h,1] <- 0
						theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
						y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
						}	## end h
					}## end j
			}## end i

# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy
for(o in 1:no.orgs)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict organ occupancy probability
	logit(pred_organ_p[o]) <- gamma[o] # predict orgen detection probability
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
#beta0 ~ dnorm(0, 0.37)
#gamma0 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	#beta[z] ~ dnorm(beta0, beta.tau)
	#gamma[z] ~ dnorm(gamma0, gamma.tau)
	gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}

## Priors for the random effects from the precision (tau)
#beta.tau <- pow(beta.ss,-2)
#beta.ss ~ dunif(0,6)

#gamma.tau <- pow(gamma.ss,-2)
#gamma.ss ~ dunif(0,6)
}

filename_pm <- paste(getwd(),"pm.bug", sep="/")
write.model(model_pm, filename_pm)
#P. minibicornis MODEL 1: TUBULE  OCCUPANCY AS A FUNCTION OF GLOMERULUS OCCUPANCY
	model_pm_01<- function ()
		{
		# i indexes total number of fish 
		for (i in 1:nfish) 
			{
		   ###  Fish-level model
			logit(fish.psi0[i]) <- alpha0
			## Ignore this it helps keeps the program in track min-max trick
			fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
			# latent fish level variable
			lat.fish[i] ~ dbern(fish.psi[i])
			## j indexes number of organs, non-kidney
			for (j in 1:5) 
				{
				## conditional organ occupancy
				logit(organ.psi0[i,j]) <- beta[j]
				### Ignore this it helps keeps the program in track, note multiplied by 
				#     fish level latent variable
				organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
				## latent organ occupancy variable
				lat.org[i,j] ~ dbern(organ.psi[i,j])
				lat.org2[i,j] <-  lat.org[i,j] + 1				
				}
	
				## conditional on kidney occupancy,tubule  occupancy depends on glomerulus occpupancy
				logit(organ.psi0[i,6]) <- beta[6]+ beta2*lat.org[i,7]
				### Ignore this it helps keeps the program in track, note multiplied by 
				#     kidney level latent variable
				organ.psi[i,6]<-max(0.00001,min(0.99999, organ.psi0[i,6]))*lat.org[i,5]
				lat.org[i,6] ~ dbern(organ.psi[i,6])
				lat.org2[i,6] <-  lat.org[i,6] + 1			

				## conditional on kidney occupancy
				logit(organ.psi0[i,7]) <- beta[7] 
				### Ignore this it helps keeps the program in track, note multiplied by 
				#     kidney level latent variable
				organ.psi[i,7]<-max(0.00001,min(0.99999, organ.psi0[i,7]))*lat.org[i,5]
				lat.org[i,7] ~ dbern(organ.psi[i,7])
				lat.org2[i,7] <-  lat.org[i,7] + 1				
			}	
		
			# detection model			
			for(i in 1:nfish)
				{
				for(j in 1:4)
					{
					for (h in 1:3) 
						{ # detection model  h indexes number of repeated samples (visits) per organ
						logit(theta[i,j,h,2]) <- gamma[j] 
						theta2[i,j,h,1] <- 0
						theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
						y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
						}	## end h
					} ## end j
				for(j in 6:7)
					{
					for (h in 1:3) 
						{ # detection model  h indexes number of repeated samples (visits) per organ
						logit(theta[i,j,h,2]) <- gamma[j] 
						theta2[i,j,h,1] <- 0
						theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
						y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
						}	## end h
					}## end j
			}## end i

# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy
for(o in 1:no.orgs)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict organ occupancy probability
	logit(pred_organ_p[o]) <- gamma[o] # predict orgen detection probability
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
#beta0 ~ dnorm(0, 0.37)
#gamma0 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	#beta[z] ~ dnorm(beta0, beta.tau)
	#gamma[z] ~ dnorm(gamma0, gamma.tau)
	gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}
beta2 ~  dnorm(0, 0.37)
## Priors for the random effects from the precision (tau)
#beta.tau <- pow(beta.ss,-2)
#beta.ss ~ dunif(0,6)

#gamma.tau <- pow(gamma.ss,-2)
#gamma.ss ~ dunif(0,6)
}
filename_pm_01 <- paste(getwd(),"pm_M01.bug", sep="/")
write.model(model_pm_01, filename_pm_01)	

# P. minibicornis MODEL 2: EFFECT OF SEX
model_pm_02<- function ()
		{
		# i indexes total number of fish 
		for (i in 1:nfish) 
			{
		   ###  Fish-level model
			logit(fish.psi0[i]) <- alpha0 + alpha1*X[i,1]
			## Ignore this it helps keeps the program in track min-max trick
			fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
			# latent fish level variable
			lat.fish[i] ~ dbern(fish.psi[i])
			## j indexes number of organs, non-kidney
			for (j in 1:5) 
				{
				## conditional organ occupancy
				logit(organ.psi0[i,j]) <- beta[j]
				### Ignore this it helps keeps the program in track, note multiplied by 
				#     fish level latent variable
				organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
				## latent organ occupancy variable
				lat.org[i,j] ~ dbern(organ.psi[i,j])
				lat.org2[i,j] <-  lat.org[i,j] + 1				
				}
	
				## conditional on kidney occupancy
				logit(organ.psi0[i,6]) <- beta[6]
				### Ignore this it helps keeps the program in track, note multiplied by 
				#     kidney level latent variable
				organ.psi[i,6]<-max(0.00001,min(0.99999, organ.psi0[i,6]))*lat.org[i,5]
				lat.org[i,6] ~ dbern(organ.psi[i,6])
				lat.org2[i,6] <-  lat.org[i,6] + 1			

				## conditional on kidney occupancy
				logit(organ.psi0[i,7]) <- beta[7] 
				### Ignore this it helps keeps the program in track, note multiplied by 
				#     kidney level latent variable
				organ.psi[i,7]<-max(0.00001,min(0.99999, organ.psi0[i,7]))*lat.org[i,5]
				lat.org[i,7] ~ dbern(organ.psi[i,7])
				lat.org2[i,7] <-  lat.org[i,7] + 1				
				
			}	
		
			# detection model			
			for(i in 1:nfish)
				{
				for(j in 1:4)
					{
					for (h in 1:3) 
						{ # detection model  h indexes number of repeated samples (visits) per organ
						logit(theta[i,j,h,2]) <- gamma[j] 
						theta2[i,j,h,1] <- 0
						theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
						y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
						}	## end h
					} ## end j
				for(j in 6:7)
					{
					for (h in 1:3) 
						{ # detection model  h indexes number of repeated samples (visits) per organ
						logit(theta[i,j,h,2]) <- gamma[j] 
						theta2[i,j,h,1] <- 0
						theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
						y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
						}	## end h
					}## end j
			}## end i

# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy
for(o in 1:no.orgs)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict organ occupancy probability
	logit(pred_organ_p[o]) <- gamma[o] # predict orgen detection probability
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
alpha1 ~ dnorm(0, 0.37)
beta0 ~ dnorm(0, 0.37)
gamma0 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	#beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	beta[z] ~ dnorm(beta0, beta.tau)
	gamma[z] ~ dnorm(gamma0, gamma.tau)
	#gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}
## Priors for the random effects from the precision (tau)
beta.tau <- pow(beta.ss,-2)
beta.ss ~ dunif(0,6)
gamma.tau <- pow(gamma.ss,-2)
gamma.ss ~ dunif(0,6)
}
filename_pm_02 <- paste(getwd(),"pm_M02.bug", sep="/")
write.model(model_pm_02, filename_pm_02)	

# P. minibicornis MODEL 2: EFFECT OF LENGTH
model_pm_03<- function ()
		{
		# i indexes total number of fish 
		for (i in 1:nfish) 
			{
		   ###  Fish-level model
			logit(fish.psi0[i]) <- alpha0 + alpha1*X[i,2]
			## Ignore this it helps keeps the program in track min-max trick
			fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
			# latent fish level variable
			lat.fish[i] ~ dbern(fish.psi[i])
			## j indexes number of organs, non-kidney
			for (j in 1:5) 
				{
				## conditional organ occupancy
				logit(organ.psi0[i,j]) <- beta[j]
				### Ignore this it helps keeps the program in track, note multiplied by 
				#     fish level latent variable
				organ.psi[i,j]<-max(0.00001,min(0.99999, organ.psi0[i,j]))*lat.fish[i]
				## latent organ occupancy variable
				lat.org[i,j] ~ dbern(organ.psi[i,j])
				lat.org2[i,j] <-  lat.org[i,j] + 1				
				}
	
				## conditional on kidney occupancy
				logit(organ.psi0[i,6]) <- beta[6]
				### Ignore this it helps keeps the program in track, note multiplied by 
				#     kidney level latent variable
				organ.psi[i,6]<-max(0.00001,min(0.99999, organ.psi0[i,6]))*lat.org[i,5]
				lat.org[i,6] ~ dbern(organ.psi[i,6])
				lat.org2[i,6] <-  lat.org[i,6] + 1			

				## conditional on kidney occupancy
				logit(organ.psi0[i,7]) <- beta[7] 
				### Ignore this it helps keeps the program in track, note multiplied by 
				#     kidney level latent variable
				organ.psi[i,7]<-max(0.00001,min(0.99999, organ.psi0[i,7]))*lat.org[i,5]
				lat.org[i,7] ~ dbern(organ.psi[i,7])
				lat.org2[i,7] <-  lat.org[i,7] + 1				
				
			}	
		
			# detection model			
			for(i in 1:nfish)
				{
				for(j in 1:4)
					{
					for (h in 1:3) 
						{ # detection model  h indexes number of repeated samples (visits) per organ
						logit(theta[i,j,h,2]) <- gamma[j] 
						theta2[i,j,h,1] <- 0
						theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
						y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
						}	## end h
					} ## end j
				for(j in 6:7)
					{
					for (h in 1:3) 
						{ # detection model  h indexes number of repeated samples (visits) per organ
						logit(theta[i,j,h,2]) <- gamma[j] 
						theta2[i,j,h,1] <- 0
						theta2[i,j,h,2] <- max(0.00001,min(0.99999, theta[i,j,h,2]))
						y[i,j,h] ~ dbin(theta2[i,j,h,lat.org2[i,j]], 1)
						}	## end h
					}## end j
			}## end i

# derived parameters
prev<- sum(lat.fish[])/nfish
logit(pred.fish) <- alpha0	## predict fish occupancy
for(o in 1:no.orgs)
	{
	logit(pred_organ_psi[o]) <- beta[o]# predict organ occupancy probability
	logit(pred_organ_p[o]) <- gamma[o] # predict orgen detection probability
	}

######   Fixed effect priors
alpha0 ~ dnorm(0, 0.37)
alpha1 ~ dnorm(0, 0.37)
beta0 ~ dnorm(0, 0.37)
gamma0 ~ dnorm(0, 0.37)

### model as random effects from logit normal, more efficient
for (z in 1:no.orgs) 
	{
	#beta[z] ~  dnorm(0, 0.37) # organ level occupancy
	beta[z] ~ dnorm(beta0, beta.tau)
	gamma[z] ~ dnorm(gamma0, gamma.tau)
	#gamma[z] ~ dnorm(0, 0.37)# organ level detection probability
	}
## Priors for the random effects from the precision (tau)
beta.tau <- pow(beta.ss,-2)
beta.ss ~ dunif(0,6)
gamma.tau <- pow(gamma.ss,-2)
gamma.ss ~ dunif(0,6)
}
filename_pm_03 <- paste(getwd(),"pm_M03.bug", sep="/")
write.model(model_pm_03, filename_pm_03)		

### AE MODEL 
model_AE_M01<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		lat.fish2[i] <-  lat.fish[i]  + 1	
		} 
		
	# detection model			
	for(i in 1:nfish)
		{
		for (h in 1:3) 
			{ # detection model  h indexes number of repeated samples (visits) per organ
			logit(theta[i,h,2]) <- gamma
			theta2[i,h,1] <- 0
			theta2[i,h,2] <- max(0.00001,min(0.99999, theta[i,h,2]))
			y[i,h] ~ dbin(theta2[i,h,lat.fish2[i]],1)
			}	## end h
		} ## end I
		
	# derived parameters
	prev<- sum(lat.fish[])/nfish
	logit(pred.fish) <- alpha0	## predict fish occupancy
	logit(p) <- gamma # predict detection
	######   Fixed effect priors
	alpha0 ~ dnorm(0, 0.37)
	gamma ~ dnorm(0, 0.37)# organ level detection probability
	}
filename_AE_M01 <- paste(getwd(),"base_AE_01.bug", sep="/")
write.model(model_AE_M01, filename_AE_M01)


### AE MODEL SEX
model_AE_M01_sex<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0 + alpha1*X[i,1]
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		lat.fish2[i] <-  lat.fish[i]  + 1	
		} 
		
	# detection model			
	for(i in 1:nfish)
		{
		for (h in 1:3) 
			{ # detection model  h indexes number of repeated samples (visits) per organ
			logit(theta[i,h,2]) <- gamma
			theta2[i,h,1] <- 0
			theta2[i,h,2] <- max(0.00001,min(0.99999, theta[i,h,2]))
			y[i,h] ~ dbin(theta2[i,h,lat.fish2[i]],1)
			}	## end h
		} ## end I
		
	# derived parameters
	prev<- sum(lat.fish[])/nfish
	logit(pred.fish) <- alpha0	## predict fish occupancy
	logit(p) <- gamma # predict detection
	######   Fixed effect priors
	alpha0 ~ dnorm(0, 0.37)
	alpha1 ~ dnorm(0, 0.37)
	gamma ~ dnorm(0, 0.37)# organ level detection probability
	}
filename_AE_M01_sex <- paste(getwd(),"base_AE_01_sex.bug", sep="/")
write.model(model_AE_M01_sex, filename_AE_M01_sex)

### AE MODEL LENGTH
model_AE_M01_len<- function ()
	{ 
	# i indexes total number of fish 
	for (i in 1:nfish) 
		{
	   ###  Fish-level model
		logit(fish.psi0[i]) <- alpha0 + alpha1*X[i,2]
		## Ignore this it helps keeps the program in track min-max trick
		fish.psi[i]<-max(0.00001,min(0.99999, fish.psi0[i]))
		# latent fish level variable
		lat.fish[i] ~ dbern(fish.psi[i])
		lat.fish2[i] <-  lat.fish[i]  + 1	
		} 
		
	# detection model			
	for(i in 1:nfish)
		{
		for (h in 1:3) 
			{ # detection model  h indexes number of repeated samples (visits) per organ
			logit(theta[i,h,2]) <- gamma
			theta2[i,h,1] <- 0
			theta2[i,h,2] <- max(0.00001,min(0.99999, theta[i,h,2]))
			y[i,h] ~ dbin(theta2[i,h,lat.fish2[i]],1)
			}	
		} 
	prev<- sum(lat.fish[])/nfish
	logit(p) <- gamma # predict detection
	######   Fixed effect priors
	alpha0 ~ dnorm(0, 0.37)
	alpha1 ~ dnorm(0, 0.37)
	gamma ~ dnorm(0, 0.37)# organ level detection probability
	}
filename_AE_M01_len <- paste(getwd(),"base_AE_01_len.bug", sep="/")
write.model(model_AE_M01_len, filename_AE_M01_len)

