

	# [1] load data
	com3<- odbcConnectExcel2007("./dat/triplicate.xlsx")
	dat<- sqlFetch(com3, "data")
	
	fishdat<- sqlFetch(com3, "fishdata")  
	ns_counts<- sqlFetch(com, "Fish data: parasite counts")
	
	odbcClose(com3)
	
	
	load("./output/ms_mcmc.Rdata")
	

	