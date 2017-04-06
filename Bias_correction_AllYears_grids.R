#--------------------------------------------------------------------------------------------------------------------#
# Bias correction for the ESSENCE precipitation model output (17 members) for each grid cell (8 in total Rhine basin)
# A similar correction was applied in a selected set of GCM from CMIP5.
# two methods of bc: one with circulation patterns (CP) and one without (CP0)
# BC method: advance delta method: P^*=aP^b, with excess correction for quantiles>90th
# For CP method: correction is applied for each circulation patterns 
# this correction is done on daily time steps. the published correction was done on 5 day sums because it was used in hydrological simulations and the purpose was to correct hydrological extremes as well.
# C.S.Photiadou

rm(list=ls())
setwd("/***/SVD_ESS/Bias_output_grids/Winter/AllYears/")
DPRECIPIN <- "/***/SVD_ESS/Results_TXT_grids/Winter/"
DINCP <- "/***/SVD_ESS/Results_TXT/Winter/ESS/"
DCORECTOUT <- "/***/SVD_ESS/Bias_output_grids/Winter/AllYears/"
DCoeffOUT <- "/***/SVD_ESS/Bias_output_grids/Winter/Coeff/"

ncp <- 7 # number of circulaion patterns
nmodW<-8019 # number of days in winter
ngrid <- 8 # number of grids for the Rhine basin
### The correction can be applied on control and validation period for validation purposes
### Here the correction was done on all the years 1961-2004
##Control Period
#Cloperiod <- read.table(paste(DPRECIPIN,"Control_period_W.txt",sep=""))$V1
##Validation Period
#Vloperiod <- read.table(paste(DPRECIPIN,"Validation_period_W.txt",sep=""))$V1
#############################
### Read CHR & CP (observations & circulation patterns)
#############################
	chrinfile <-read.table(paste(DPRECIPIN,"CHR_grids_W.txt",sep=""))
	chrdata_all <- chrinfile[,1:8]
	chrcp <- chrinfile[,9]
	
#############################
### Read ESS & CP (ESSNECE ensemble (17memb) & circulation patterns)
#############################
# to save the a & b coeff for P=aP^b for correction with circulation patterns (CP) & correction without (CP0)
Acoeff_CP <- array(NA,dim=c(7,8))
Bcoeff_CP <- array(NA,dim=c(7,8))

Acoeff_CP0 <- array(NA,dim=c(8))
Bcoeff_CP0 <- array(NA,dim=c(8))

i<-0
files <- list.files( "/***/SVD_ESS/Results_TXT_grids/Winter/ESS/")

for (iess in 1:17){
	i<-i+1
	CP_Final<-array(NA,dim=c(nmodW,8))
	CP0_Final<-array(NA,dim=c(nmodW,8))
	
  # read precipitation of model per member
	modfilein <- read.table(paste(DPRECIPIN,files[iess],sep=""))
	# read cp of model per member
	modcp_order <-  modcpin[order(modcpin[,1]),]
	modcp_all <- modcp_order[[2]]
	modcp <- modfilein[,9]
	moddata_all <- modfilein[,1:8]
	
# Correction without CPs (CP0)
	# select low and high quantiles
	lowquant <- 0.6
	highquant <- 0.95
	k<- 0
	for (igrid in 1:ngrid){
	k<-k+1
	chrdata <- chrdata_all[,igrid]
	moddata <- moddata_all[,igrid]
		
	q60obs <- quantile(chrdata,lowquant)
	q95obs <- quantile(chrdata,highquant)
	q60mod <- quantile(moddata,lowquant)
	q95mod <- quantile(moddata,highquant)

  # calculate a & b coeff
	b <- log(/q60obs) / log(q95mod/q60mod)
	a <- q60obs/(q60mod^b)
			
			Acoeff_CP0[k] <- signif(a,digits=3)

			Bcoeff_CP0[k] <- signif(b,digits=3)

			write.table(Acoeff_CP0,paste(DCoeffOUT,"Acoef_NoCP_",iess,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)	
			write.table(Bcoeff_CP0,paste(DCoeffOUT,"Bcoef_NoCP_",iess,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)	

  # apply excess correction >90th quantile
	excessobs <- mean(chrdata[chrdata>])-
	excessmod <- mean(moddata[moddata>q95mod])-q95mod
	lolow <- moddata <= q95mod
	lohigh <- !lolow
	cormoddata1 <- a*moddata ^ b
	cormoddata2 <- (excessmod/excessobs) * (moddata-q95mod) + a*q95mod^b

	cormoddata <- cormoddata1
	cormoddata[lohigh] <- cormoddata2[lohigh]
	n <- length(cormoddata)

# Correction with CPs (CP) (each CP has its own a & b coeff and excess correction)
	cormoddatacp <- array(NA,dim=n)
        j<-0
      	for (icp in 1:9){
	              j<-j+1

	 	             chrdatacp <- chrdata[chrcp == icp]
	 	             moddatacp <- moddata[modcp == icp]

	 	             q60obs <- quantile(chrdatacp,lowquant)
	 	             q95obs <- quantile(chrdatacp,highquant)
	 	             q60mod <- quantile(moddatacp,lowquant)
	 	             q95mod <- quantile(moddatacp,highquant)
	 	             excessobs <- mean(chrdatacp[chrdatacp>])-
	 	             excessmod <- mean(moddatacp[moddatacp>q95mod])-q95mod

	 	             bcp <- log(/q60obs) / log(q95mod/q60mod)
	 	             acp <- q60obs/(q60mod^b)

	 	             Acoeff_CP[j,k] <- signif(acp,digits=3)
	 	             Bcoeff_CP[j,k] <- signif(bcp,digits=3)

	 	             write.table(Acoeff_CP,paste(DCoeffOUT,"Acoef_CP_",iess,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)	
	 	             write.table(Bcoeff_CP,paste(DCoeffOUT,"Bcoef_CP_",iess,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)	
		 	            
	 	            # apply excess correction for each CP.
	 	            # analysis was made on this. one can apply a mean excess correction over all CP but in this cases the excess correction worked better like this 
	 	             lolow <- moddatacp <= q95mod
	 	             lohigh <- !lolow
	 	             cormoddata1 <- acp*moddatacp ^ bcp
	 	             cormoddata2 <- (excessmod/excessobs) * (moddatacp-q95mod) + acp*q95mod^bcp

	 	             cormoddata_tmp <- cormoddata1
	 	             cormoddata_tmp[lohigh] <- cormoddata2[lohigh]
  
	 	             cormoddatacp[modcp == icp] <- cormoddata_tmp
          }
    	CP0_Final[,k]<-cormoddata
	    CP_Final[,k]<-cormoddatacp
  	}#igrid
	# Final bias corrected model precipitation per member for CP & CP0 method
	write.table(CP_Final,paste(DCORECTOUT,"CP_grid_ESS",iess,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(CP0_Final,paste(DCORECTOUT,"CP0_grid_ESS",iess,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
}#iess
