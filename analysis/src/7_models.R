
model<- function ()
	{
	# FISH LEVEL PSIs
	Psi_A<- max(0.0001,min(0.9999, psi_z[1]))
	Psi_aB<- max(0.0001,min(0.9999, psi_z[2]))
	Psi_AB<- max(0.0001,min(0.9999, psi_z[3]))	
	Psi_abC<- max(0.0001,min(0.9999, psi_z[4]))	
	Psi_AbC<- max(0.0001,min(0.9999, psi_z[5]))	
	Psi_aBC<- max(0.0001,min(0.9999, psi_z[6]))	
	Psi_ABC<- max(0.0001,min(0.9999, psi_z[7]))	
	Psi_abcD<- max(0.0001,min(0.9999, psi_z[8]))	
	Psi_AbcD<- max(0.0001,min(0.9999, psi_z[9]))	
	Psi_aBcD<- max(0.0001,min(0.9999, psi_z[10]))	
	Psi_ABcD<- max(0.0001,min(0.9999, psi_z[11]))	
	Psi_abCD<- max(0.0001,min(0.9999, psi_z[12]))	
	Psi_AbCD<- max(0.0001,min(0.9999, psi_z[13]))	
	Psi_aBCD<- max(0.0001,min(0.9999, psi_z[14]))	
	Psi_ABCD<- max(0.0001,min(0.9999, psi_z[15]))	
	for(fish in 1:nfish)
		{				
		psi[fish,1]<-(1-Psi_A)*(1-Psi_aB)*(1-Psi_abC)*(1-Psi_abcD)
		psi[fish,2]<-Psi_A*(1-Psi_AB)*(1-Psi_AbC)*(1-Psi_AbcD)
		psi[fish,3]<-(1-Psi_A)*Psi_aB*(1-Psi_aBC)*(1-Psi_aBcD)   
		psi[fish,4]<-Psi_A*Psi_AB*(1-Psi_ABC)*(1-Psi_ABcD)       
		psi[fish,5]<-(1-Psi_A)*(1-Psi_aB)*Psi_abC*(1-Psi_abCD)    
		psi[fish,6]<-Psi_A*(1-Psi_AB)*Psi_AbC*(1-Psi_AbCD)        
		psi[fish,7]<-(1-Psi_A)*Psi_aB*Psi_aBC*(1-Psi_aBCD)        
		psi[fish,8]<-Psi_A*Psi_AB*Psi_ABC*(1-Psi_ABCD)           
		psi[fish,9]<-(1-Psi_A)*(1-Psi_aB)*(1-Psi_abC)*Psi_abcD   
		psi[fish,10]<-Psi_A*(1-Psi_AB)*(1-Psi_AbC)*Psi_AbcD      
		psi[fish,11]<-(1-Psi_A)*Psi_aB*(1-Psi_aBC)*Psi_aBcD      
		psi[fish,12]<-Psi_A*Psi_AB*(1-Psi_ABC)*Psi_ABcD        
		psi[fish,13]<-(1-Psi_A)*(1-Psi_aB)*Psi_abC*Psi_abCD       
		psi[fish,14]<-Psi_A*(1-Psi_AB)*Psi_AbC*Psi_AbCD         
		psi[fish,15]<-(1-Psi_A)*Psi_aB*Psi_aBC*Psi_aBCD           
		psi[fish,16]<-Psi_A*Psi_AB*Psi_ABC*Psi_ABCD
		yy[fish] ~ dcat(psi[fish,]) # categorical 1 draw from a multinomial (1:4 possible states)	
		for(j in 1:4)
			{
			true_fish[fish,j]<- combos[yy[fish],j]
			true_fish_ones[fish,j]<- combos[yy[fish],j]+1
			}
		}	
	### ORGAN LEVEL	
	for(organ in 1:5)
	{
		for(path in 1:4)
		{
		psi_organ_z[organ,path,1]<- 0
		psi_organ_z[organ,path,2]<- max(0.0001,min(0.9999, psi_organ[organ,path]))		
			for(fish in 1:nfish)
			{
			organ_infect_true[fish,path,organ]~dbern(psi_organ_z[organ,path,true_fish_ones[fish,path]])
			organ_infect_true_ones[fish,path,organ]<-organ_infect_true[fish,path,organ]+1
			}
		}
	}
		
	## DETECTION MODEL
	for(path in 1:4)
		{
		for(organ in 1:5)
			{		
			p_z[organ,path,1]<- 0
			p_z[organ,path,2]<- max(0.0001,min(0.9999, p[organ,path]))	
				for(fish in 1:nfish)
				{
				for(reps in 1:3)
					{					
					yyy[fish,path,organ,reps]~dbern(p_z[organ,path,organ_infect_true_ones[fish,path,organ]])
					}
				}
			}
		}
	## DERIVED PARAMETERS
	psiA_est<-psi[1,2]+psi[1,4]+psi[1,6]+psi[1,8]+psi[1,10]+psi[1,12]+psi[1,14]+psi[1,16]
	psiB_est<-psi[1,3]+psi[1,4]+psi[1,7]+psi[1,8]+psi[1,11]+psi[1,12]+psi[1,15]+psi[1,16]
	psiC_est<-psi[1,5]+psi[1,6]+psi[1,7]+psi[1,8]+psi[1,13]+psi[1,14]+psi[1,15]+psi[1,16]
	psiD_est<-psi[1,9]+psi[1,10]+psi[1,11]+psi[1,12]+psi[1,13]+psi[1,14]+psi[1,15]+psi[1,16]
		
	psiAB_est<-psi[1,4]+psi[1,8]+psi[1,12]+psi[1,16]		
	psiBC_est<-psi[1,7]+psi[1,8]+psi[1,15]+psi[1,16]	
	psiAC_est<-psi[1,6]+psi[1,8]+psi[1,14]+psi[1,16]	
	psiAD_est<-psi[1,10]+psi[1,12]+psi[1,14]+psi[1,16]		
	psiBD_est<-psi[1,11]+psi[1,12]+psi[1,15]+psi[1,16]		
	psiCD_est<-psi[1,13]+psi[1,14]+psi[1,15]+psi[1,16]		
	
	# SPECIES INTERACTION FACTORS
	sifAB<-	psiAB_est/(psiA_est*psiB_est)
	sifBC<- psiBC_est/(psiB_est*psiC_est)
	sifAC<- psiAC_est/(psiA_est*psiC_est)
	sifAD<- psiAD_est/(psiA_est*psiD_est)
	sifBD<- psiBD_est/(psiB_est*psiD_est)
	sifCD<-	psiCD_est/(psiC_est*psiD_est)
	
	
	## PRIORS ##	
		# INFECTION PROBABILITIES (conditional
		for(i in 1:15)
			{
			theta[i]~dnorm(0,0.37)
			logit(psi_z[i])<- theta[i]
			}		
		for(path in 1:4)
			{
			for(organ in 1:5)
				{
				eta[organ,path]~dnorm(0,0.37)
				logit(psi_organ[organ,path])<- eta[organ,path]
				beta[organ,path]~dnorm(0,0.37)				
				logit(p[organ,path])<- beta[organ,path]
				}
			}
			
	}
filename <- paste(getwd(),"ms_occ.bug", sep="/")
write.model(model, filename)	
