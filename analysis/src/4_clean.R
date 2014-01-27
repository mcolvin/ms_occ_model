	
	
	
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
	dat$organ_id<- as.numeric(as.character(dat$organ_id))
	dat$parasite<- as.factor(dat$parasite)
	analysis_dat<- array(0, dim=c(26,4,5,3)) 

	for(organ in 1:5)
		{
		for(reps in 1:3)
			{
			fillDat<- dat[dat$organ_id==organ & dat$rep==reps,]
			analysis_dat[,,organ,reps]<- as.matrix(dcast(fillDat, id~parasite,value.var="sev",mean,
				fill=0,drop=FALSE)[,-c(1:4)])
			}
		}
	
analysis_dat[analysis_dat>=1]<-1
nfish = 26
