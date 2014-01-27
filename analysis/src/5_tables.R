	tables<- function(n)
		{
		if(n==1)
			{
			# (1=gill/tubules, 2=heart, 3=liver, 4=spleen, 5=kidney)
			# A = AE, B = BKD, C = Nano, D = Pm
			
			
			
			}
		
		if(n==1)
			{
ests<- out$summary
vals<- as.data.frame(rbind( ests["sifAB",c(1,3,7)],
   ests["sifBC",c(1,3,7)],
	ests["sifAC",c(1,3,7)],
	ests["sifAD",c(1,3,7)],
	ests["sifBD",c(1,3,7)],
	ests["sifCD",c(1,3,7)]))
         
    vals$cols<- c("sifAB", "sifBC", "sifAC","sifAD","sifBD", "sifCD")  
  }
           
	}  