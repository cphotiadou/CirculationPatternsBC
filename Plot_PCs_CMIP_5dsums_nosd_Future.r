#######################################################################################
######Graphical output
######
#X11(display="")
#old.par <- par(no.readonly = TRUE)

.libPaths(c("/***/PROJECTS/CMIP5_Runs/",.libPaths()))
setwd("/Volumes/TimeMachineBackups/ethz/cmip5/rcp45/CMIP/Winter/CP/")
library(RNetCDF)
library(maps)
library(mapdata)
library(lattice)
library(maptools)
library(fields)
library(abind)
# library(zoo)
rm(list = ls(all=TRUE))

options(show.error.messages = TRUE)

#nRhine<-readShapeSpatial("/Users/photiadou/PROJECTS/CMIP5_Runs/CMIP5/Try_Winter/shapes/rhine_hbv_basins_x134")
# nmod<-8151
nlonERA <-20
nlatERA <-20

DGCMIN <- "/Volumes/TimeMachineBackups/ethz/cmip5/rcp45/day/FUTURE/PSL/"
DGCMOUT <- "/Volumes/TimeMachineBackups/ethz/cmip5/rcp45/CMIP/Winter/CP/"
DTXT   <- "/Users/photiadou/PROJECTS/SVD_ESS/INPUT/latlon/"
DIndex   <- "/Users/photiadou/PROJECTS/SVD_ESS/INPUT/Index_season/"

##Read txt fiels
latERA     <- read.table(paste(DTXT,"latERA_slp.txt",sep=""),header=FALSE)$V1
lonERA     <- read.table(paste(DTXT,"lonERA_slp.txt",sep=""),header=FALSE)$V1
#w_61_f<- read.table(paste(DIndex,"winter_index.txt",sep=""))$V1
#s_61_f<- read.table(paste(DIndex,"s_61_f.txt",sep=""))
###############################################################
                    #### Read s$u
###############################################################
UCP_w <- as.matrix(read.table("/Users/photiadou/PROJECTS/SVD_ESS/OUTPUT/sW$u_W.txt"))
################################################################
                    # #### Read CMIP5 & interpolate into ERA grid
# ###############################################################
files <- list.files("/Volumes/TimeMachineBackups/ethz/cmip5/rcp45/day/FUTURE/PSL/")
names_gcm<-c("access1-0","bcc-csm1","canesm2","ccsm4","cnrm-cm5","csiro-mk3","ec-earth",
			"gfdl-cm3","gfdl-esm2g","gfdl-esm2m","giss-e2-r","inmcm4",
			"ipsl-cm5a-mr","miroc-esm","miroc-esm-chem","miroc5","mpi-esm-lr","mri-cgcm3","noresm1-m")

names_cap <- c("ACCESS1-0","bcc-csm1-1","CanESM2","CCSM4","CMRM-CM5","CSIRO-Mk3-6-0","EC-EARTH",
"GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-R","inmcm4",
"IPSL-CM5A-MR","MIROC-ESM","MIROC-ESM-CHEM","MIROC5","MPI-ESM-LR","MRI-CGCM3","NorESM1-M")

#write.table(names_gcm,paste(DGCMOUT,"names_gcm.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
#######################################
### Define ERA-I grid
	grid.dest <- list(x=lonERA,y=latERA)
	nlonERA <- length(lonERA)
	nlatERA <- length(latERA)
####################################### 
ngcm<-19
################ Loop all GCMs
for (igcm in 1:ngcm){
    rm(list=(ls()[ls()=="df_PC"]))
    rm(list=(ls()[ls()=="dlist_p"]))
    rm(list=objects(pattern="^PC_dates"))
    
		filename <- paste(DGCMIN,files[igcm],sep="")
		print(filename)
		varfileOne <- open.nc(filename,write=FALSE)
		varfileOne_psl <- var.get.nc(varfileOne,"psl")
		lat_pls <- var.get.nc(varfileOne,"lat")
		lon_pls <- var.get.nc(varfileOne,"lon")
		inmode <- var.get.nc(varfileOne,"time")
		timeunits <- att.get.nc(varfileOne, "time", "units")
		dtg <- utcal.nc(timeunits,inmode,"n")
	close.nc(varfileOne)
	###Fix years
	index_yyyy <- which(dtg[,1]>=2081)
	dtg_yyyy <- data.frame(subset(dtg,dtg[,1]>=2081))
	varfileOne_psl_81 <- varfileOne_psl[,,index_yyyy]
	
	### Select Winter
    ### Select winter & winter
    w_81_1 <- which(dtg_yyyy[,2]>=10)
    w_81_2 <- which(dtg_yyyy[,2]<=3)
    w_81_3 <- which(dtg_yyyy[,2]==9 & dtg_yyyy[,3]>=28)
    w_81_f <- sort(c(w_81_1,w_81_2,w_81_3))
    #write.table(w_61_f,paste(DIROUT,"winter_index.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
    #write.table(w_61_f,"/Users/photiadou/PROJECTS/SVD_ESS/INPUT/Index_season/winter_index_all.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
    
    
    winter_index<-array(NA,dim=c(length(w_81_f),3))
    i<-0
    for (ni in 1:3){
        i<-i+1
        winter_index[,ni] <- dtg_yyyy[,ni][w_81_f]
    }

	varfileOne_psl_w <- varfileOne_psl_81[,,w_81_f]
	nmodW <- dim(varfileOne_psl_w)[3]
	newvarSLPintW <- array(0,dim=c(nlonERA,nlatERA,nmodW))
	
	### Interpolate
	for (itim in 1:nmodW){
		pat_obj <- list(x=lon_pls, y=lat_pls, z=varfileOne_psl_w[,,itim])
		pat_intp_obj <- interp.surface.grid(pat_obj,grid.dest)
		newvarSLPintW[,,itim] <- pat_intp_obj$z
	}
	
	#### 2D CMIP5 SLP
	lovalgcm<-array(FALSE,dim=c(nlonERA,nlatERA))
	for (ilon in 1:nlonERA){
		for (ilat in 1:nlatERA){
			if (sum(is.na(newvarSLPintW[ilon,ilat,]))==0){
				lovalgcm[ilon,ilat]<-TRUE
			}
		}
	}
	gcm_athroisma <- sum(lovalgcm)
	GCMPres2dW <- array(0,dim=c(gcm_athroisma,nmodW))
	
	i <-0
	for (ilon in 1:nlonERA){
		for (ilat in 1:nlatERA){
			if (lovalgcm[ilon,ilat]){
				i <-i+1
				GCMPres2dW[i,1:nmodW]= newvarSLPintW[ilon,ilat,1:nmodW]
			}
		}#ilat
	}#ilon
	
	
	PCGCM_W <- PressureCorrectGCM_W <- GCMPres2dW
####WINTER  
##Make correlations for positive U  & negative U
## Take the mean and sd
	MeanPressGCM_W <- apply(PressureCorrectGCM_W,1,mean)
	SDPressESS_W <- apply(PressureCorrectGCM_W,1,sd)
	for (ni in 1:nmodW){
		PressureCorrectGCM_W[,ni]<-(PCGCM_W[,ni]-MeanPressGCM_W)/SDPressESS_W
	}

	write.table(GCMPres2dW,paste(DGCMOUT,names_gcm[igcm],"_2dSLP.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
	
	###Make 5dsums SLP
    ### Running mean & sd
	newGCM<-PressureCorrectGCM_W
	transPressureGCM<-t(newGCM)
	PressureGCM_all <- cbind(winter_index[,1:3],transPressureGCM[,1:400])
	slp_athroisma<-400
    
    lastyear <- max(winter_index[,1])

    
set_out_2100 <- subset(PressureGCM_all,PressureGCM_all[,1]!=lastyear | PressureGCM_all[,2]!=9 & PressureGCM_all[,2]!=10 & PressureGCM_all[,2]!=11 & PressureGCM_all[,2]!=12)
diplo_eksw <- set_out_2100
neksw <- dim(diplo_eksw)[1]
final_index <- set_out_2100[91:neksw,]

	slidetime <-c(1,2,3,4,5,10,180)
		sltime<-slidetime[5]
		 ndata <-dim(final_index)[1]
		 nyrs <- floor(ndata/185)
		 nslot_first<- floor(ndata/sltime)
		 nslot_f <- floor(nslot_first/nyrs)*nyrs
        yyyy<-seq(2081,lastyear,1)
	
	

	enagrid_mean<-array(NA,dim=c(ndata))
	enagrid_std<-array(NA,dim=c(ndata,1))
	olamazi_mean<-array(NA,dim=c(ndata,slp_athroisma))
	olamazi_std<-array(NA,dim=c(ndata,slp_athroisma))
	timemean_f <-array(NA,dim=185)
	timestd_f <-array(NA,dim=185)
	j<-0
		for (mi in 1:slp_athroisma){
			j<-j+1
			i<-0
			for (ni in 1:(length(yyyy)-1)){
				i<-i+1
				plapla<- which(PressureGCM_all[,1]==yyyy[ni] & PressureGCM_all[,2]>=9 | PressureGCM_all[,1]==yyyy[ni+1] & PressureGCM_all[,2]<=3) 
			plo<-PressureGCM_all[,mi+3][plapla]
			n<-length(plo)
			sample_slp<-plo[1:n]
			nslot <- floor(n/sltime)
			Preavg <- plo[1:(nslot*sltime)]
			dim(Preavg) <- c(sltime,nslot)
			Preavgt_all <- apply(Preavg,2,sum)
			nx<-length(Preavgt_all)

			enagrid_mean[(nx*(i-1)+1):(nx*i)]<-Preavgt_all[1:nx]
			}#ni
			olamazi_mean[,j]<-enagrid_mean
		}#mi
	PressureGCM_5dsum<-t(olamazi_mean[1:nslot_f,])


	###############################################################
                    #### Create PC_slp of GCM
	###############################################################
	PC_msl <- t(PressureGCM_5dsum) %*% UCP_w

	###################################################################
                ### Find occurrence of CPs using PCs.
	###################################################################	
	
	ind_Winter<-seq(1,nslot_f,1)
	PC_slp<-cbind(PC_msl,ind_Winter)


	# PC_complete<-array(NA,dim=c(1591,134))
	# k<-0
	# for (isub in 1:134){
		# k<-k+1
		# sdPC<-sd(PC_slp[,isub])*0.5
		# ind_sel<-which(abs(PC_slp[,isub])>=abs(sdPC))
		# PC_sel<-cbind(PC_slp[,isub][ind_sel],PC_slp[,135][ind_sel])
		# DT<-data.frame(MSL=PC_sel[,1],Index=PC_sel[,2])
		# z<-zoo(DT$MSL,DT$Index)
	    # z <- merge(zoo(,as.numeric(1:1591)), z) 
		# zdf<-data.frame(z)
		# PC_complete[,k]<-zdf$z
	# }#isub
	

	ind_phases_maxPC_sd<-array(NA,dim=nslot_f)
		maxPC_sd<-array(NA,dim=nslot_f)

	for(mi in 1:nslot_f){
    	maxPCSLP_sd<-max(abs(PC_msl[mi,]),na.rm=T)
    	ind_phases_maxPC_sd[mi]<-which(abs(PC_msl[mi,])==maxPCSLP_sd)
    	maxPC_sd[mi]<-maxPCSLP_sd
		}#mi

	MaxPC_GCM_sd <- tabulate(ind_phases_maxPC_sd)
	df_indices_sd<-data.frame(days=MaxPC_GCM_sd,Patterns=seq(1,length(MaxPC_GCM_sd),1))
	write.table(df_indices_sd,paste(DGCMOUT,"df_indices_sd_",names_gcm[igcm],".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

	################################################
		### Scaled PCs and sd_PC -> Define Threshold
	################################################

#	PCrest_index <-as.matrix( which(ind_phases_maxPC_sd>=4))
#	PCrest_dates <- data.frame(Index= PCrest_index,CP=rep(c(9),length(PCrest_index)))
#	rest_PC<- maxPC_sd[PCrest_index]

    
	################################################
	### Find Positive & Negative stages of PC1-3
	################################################
	 for (ipc in 1:max(ind_phases_maxPC_sd)){
	 	PC_index <- which(ind_phases_maxPC_sd ==ipc)
	 	PC_msl_ipc<-cbind(PC_msl[,ipc][PC_index],PC_index)

	 	ind_PC_nega <- as.matrix(subset(PC_msl_ipc, PC_msl_ipc[,1]<0))
		PC_nega_dates<-data.frame(Pressure=ind_PC_nega[,1],index=ind_PC_nega[,2],CP=rep(c(paste("-",ipc,sep="")),length(ind_PC_nega[,1])))
	
		ind_PC_posit <- as.matrix(subset(PC_msl_ipc, PC_msl_ipc[,1]>=0))
		PC_posit_dates<-data.frame(Pressure=ind_PC_posit[,1],index=ind_PC_posit[,2],CP=rep(c(paste(ipc,sep="")),length(ind_PC_posit[,1])))
	
		assign(paste("PC_dates",ipc,sep=""),rbind(PC_nega_dates,PC_posit_dates))
     }#ipc
    
    
    dlist_p <- lapply(ls(patt='^PC_dates'),get) 
    
    df_PC<- abind(dlist_p,along=1)
    
    write.table(df_PC,paste(DGCMOUT,"CP_dates_",names_gcm[igcm],".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

}#igcm
