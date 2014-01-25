	
	
	
	# [1] FORMAT DATA FROM WIDE TO LONG 
	tmp<- reshape(dat,
		varying = names(dat) [4:ncol(dat)],
		v.names = "severity",
		timevar= "tissue_parasite",
		times = names(dat) [4:ncol(dat)],
		direction = "long")
	tmp$organ<- matrix(unlist(strsplit(tmp$tissue_parasite,"_")),ncol=2, byrow=TRUE)[,2]
	tmp$parasite<- matrix(unlist(strsplit(tmp$tissue_parasite,"_")),ncol=2, byrow=TRUE)[,1]
	tmp$severity<- as.character(tmp$severity)
	
	yy<-unlist(strsplit(tmp$severity,","))
	yy[yy=="."]<- -99
	yy<-as.data.frame(matrix(as.numeric(yy),ncol=3,byrow=TRUE))
	names(yy)<-c("x1","x2","x3")	
	tmp<- cbind(tmp,yy)
	
	tmp2<- tmp[,-c(9,10)]
	tmp2$rep<-1
	names(tmp2)[8]<-"sev"
	
	tmp_app<- tmp[,-c(8,10)]
	tmp_app$rep<-2
	names(tmp_app)[8]<-"sev"
	tmp2<- rbind(tmp2,tmp_app)
	
	tmp_app<- tmp[,-c(8,9)]
	tmp_app$rep<-3
	names(tmp_app)[8]<-"sev"	
	tmp2<- rbind(tmp2,tmp_app)
	
	tmp2<- tmp2[,-c(1,3,5)]
	dat<- tmp2
	dat[dat==-99]<-NA

	
	# (1=gill, 2=heart, 3=liver, 4=spleen, 5=kidney)
	dat$organ_id<- 1
	dat[dat$organ=="glo",]$organ_id<- 5
	dat[dat$organ=="kid",]$organ_id<- 5
	dat[dat$organ=="hrt",]$organ_id<- 2
	dat[dat$organ=="liv",]$organ_id<- 3
	dat[dat$organ=="spl",]$organ_id<- 4
	
	dcast(dat, id+organ_id~parasite,value.var="sev",mean, subset=.(rep==1 & organ_id==1))
	analysis_dat<- array(0, dim=c(26,4,5,3)) 
	
	
	
	
	dat$withinOrgan<- dat$organ
	dat[!(dat$organ %in% c("glo", "tub")),]$withinOrgan<- NA
	dat[dat$organ %in% c("glo", "tub"),]$organ<- 'kid'

	dat$pres<- 0
	dat[dat$sev %in% c(1,2,3,4),]$pres<- 1
	dat[dat==-99]<-NA
	dat[is.na(dat$sev),]$pres<- NA
	
	xx<- dcast(dat, id+rep~parasite+organ+withinOrgan, value.var="sev",mean)		
	
	nfish = 26
	
	# [2] MAKE DATASETS FOR EACH PATHOGEN
	# N. SALMONICOLA DATA (1=gill, 2=heart, 3=liver, 4=spleen, 5=kidney, 6=tub, 7=glom)
	ns_dat<- array(0, dim=c(nfish,5,3)) 
	nano<- xx[,c(1,2,12,13,14)]
	names(nano)[-c(1:2)]<- c("xx1",'xx2','xx5')
	nano$xx3<- nano$xx4<-0
	nano<- nano[,c(1,2,3,4,7,6,5)]
	ns_dat[,,1]<- as.matrix(nano[nano$rep==1,-c(1,2)])
	ns_dat[,,2]<- as.matrix(nano[nano$rep==2,-c(1,2)])
	ns_dat[,,3]<- as.matrix(nano[nano$rep==3,-c(1,2)])
	ns_dat[ns_dat>0]<- 1
	
	# P. MINIBCORNIS (1=gill, 2=heart, 3=liver, 4=spleen, 5=kidney, 6=tub, 7=glom)
	pm_dat<- array(0, dim=c(nfish,5,3)) 
	pm<- xx[,c(1,2,15,16)]
	names(pm)[-c(1,2)]<- c("xx4","xx5")
	pm$xx3<-pm$xx2<-pm$xx1<- 0
	pm<- pm[,c(1,2, 5:7, 3,4)]
	pm_dat[,,1]<- as.matrix(pm[pm$rep==1,-c(1,2)])
	pm_dat[,,2]<- as.matrix(pm[pm$rep==2,-c(1,2)])
	pm_dat[,,3]<- as.matrix(pm[pm$rep==3,-c(1,2)])
	pm_dat[pm_dat>0]<- 1
		
	# R. SALMONINARUM (1=gill, 2=heart, 3=liver, 4=spleen, 5=kidney)	
	bkd_dat<- array(0, dim=c(nfish,5,3)) 
	bkd<- xx[,c(1,2,4:6)]
	names(bkd)[-c(1:2)]<- c("xx5",'xx3','xx4')
	bkd$xx1<- bkd$xx2<-0
	bkd<- bkd[,c(1,2,7,6,4,5,3)]
	bkd_dat[,,1]<- as.matrix(bkd[bkd$rep==1,-c(1,2)])
	bkd_dat[,,2]<- as.matrix(bkd[bkd$rep==2,-c(1,2)])
	bkd_dat[,,3]<- as.matrix(bkd[bkd$rep==3,-c(1,2)])
	bkd_dat[bkd_dat>0]<- 1

	# AE (1=gill, 2=heart, 3=liver, 4=spleen, 5=kidney)
	ae_dat<- array(0, dim=c(nfish,5,3)) 
	ae<- xx[,c(1,2,3)]
	names(ae)[-c(1:2)]<- c("xx1")
	ae$xx5<- ae$xx4<-ae$xx3<-ae$xx2<-0
	ae_dat[,,1]<- as.matrix(ae[ae$rep==1,-c(1,2)])
	ae_dat[,,2]<- as.matrix(ae[ae$rep==2,-c(1,2)])
	ae_dat[,,3]<- as.matrix(ae[ae$rep==3,-c(1,2)])
	ae_dat[ae_dat>0]<- 1
	
	ae_dat<- dcast(xx, id~rep, value.var="AE_gill_NA", mean)[,c(2,3,4)]
	ae_dat[ae_dat>0]<- 1
	ae_dat<- as.matrix(ae_dat)
	
	# CLEAN UP ESTIMATES
	#estimates$organ<- "Fish"
	#for(i in 1:nrow(estimates))
	#	{
	#	estimates$organ[i]<- unlist(strsplit(as.character(estimates$X)[i],"\\["))[2]
	#	estimates$parameter[i]<- unlist(strsplit(as.character(estimates$X)[i],"\\["))[1]
	#	}
	#	estimates[is.na(estimates$organ),]$organ<- "Fish"
	#estimates[estimates$organ=="1]",]$organ<- "Gills"
	#estimates[estimates$organ=="2]",]$organ<- "Heart"
	#estimates[estimates$organ=="3]",]$organ<- "Liver"
	#estimates[estimates$organ=="4]",]$organ<- "Spleen"
	#estimates[estimates$organ=="5]",]$organ<- "Kidney"
	#estimates[estimates$organ=="6]",]$organ<- "Tubules"
	#estimates[estimates$organ=="7]",]$organ<- "Glomerulus"
	
	
	# FISH LEVEL COVARIATES
	ns_counts<- subset (ns_counts, organ%in% c("kidney","Kidney"))
	ns_counts<- subset (ns_counts, necropsyId%in% fishdat$necropsyId)
	ns_counts<- ns_counts[-nrow(ns_counts),] # DROP SECOND 466
	ns_counts<-ns_counts[,c(5,10)]
	
	fishdat<- merge(fishdat, ns_counts, by="necropsyId",all.x=TRUE)
	#fishdat[fishdat$length< 0,]$length<- NA
	fishdat$sex_id<- ifelse(fishdat$sex=="M",0,1)
	fishdat<- fishdat[order(fishdat$rId),]
	XX<- as.matrix(cbind(fishdat$sex_id, fishdat$length))
	
	XXX<- matrix(0, nrow=26, ncol=5)
	XXX[,5]<- fishdat[,7]
	
	####
	
	
	
	
	
	