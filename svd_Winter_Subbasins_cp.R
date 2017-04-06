## calculates a SVD between normalized ERA-Interim_mslp fields and normalized CHR08_prec records for the period of 1979-2007.
## All input data are assumed to be normalized (z-variable). A single svd
## is computed by correlating a mslp field with the combined prec arrays
##C.S.Photiadou
##Rhine
rm(list = ls(all=TRUE))
.libPaths(c("/***/SVD_NewDomain/Try_Winter/Daily_svd/",.libPaths()))
library(RNetCDF)
library(maps)
library(mapdata)
library(lattice)
library(fields)

options(show.error.messages = TRUE)

#Set the domain boundaries (this was set after a testing over possible domains); 134 are the number of subcatchments in the Rhine
# nmod=timesteps

nlon_msl <- 20
nlat_msl <- 20
precip_sum <- 134
nmod=5285

##########################
# read model data and do normalization
#########################
		## Read MSLP daily
			varfileMSLPdaily <- open.nc("/***/CMIP5/psl/Europe_MSLP_Winter.nc",write=FALSE)
             varMSLPdaily <-  var.get.nc(varfileMSLPdaily, "msl")
             scfactor.MSLP<- att.get.nc(varfileMSLPdaily,"msl","scale_factor")
			 addoffset.MSLP<- att.get.nc(varfileMSLPdaily,"msl","add_offset")
			 
			 lonE <- var.get.nc(varfileMSLPdaily, "longitude")
			 lonE[lonE>180] <- lonE[lonE>180]-360
			 l <- sort(lonE,index.return=T)
			 lonERA <- l$x

			 latE <- var.get.nc(varfileMSLPdaily, "latitude")
			 la <- sort(latE,index.return=T)
			 latERA <- la$x
			close.nc(varfileMSLPdaily)

            newvarERA <- varMSLPdaily[l$i,la$i,]
#################################################################################
##Correct physical units for pressure: multiply by scale factor and add the offset.
MSLP.daily <- array(NA,dim=c(nlon_msl,nlat_msl,nmod))

for (ilon in 1:nlon_msl){
	for (ilat in 1:nlat_msl){
		for (ni in 1:nmod){
	MSLP.daily[ilon,ilat,ni]<-newvarERA[ilon,ilat,ni]*scfactor.MSLP+addoffset.MSLP
		}
	}
}

## Standarize Pressure 
varMSLP_mean <- apply(MSLP.daily,c(1,2),mean)
varMSLP_std<- apply(MSLP.daily,c(1,2),sd)

var_msl <- array(NA,dim=c(nlon_msl,nlat_msl,nmod))
for (ilon in 1:nlon_msl){
	for (ilat in 1:nlat_msl){
		for (ni in 1:nmod){
			var_msl[ilon,ilat,ni] <- (MSLP.daily[ilon,ilat,ni]-varMSLP_mean[ilon,ilat])/varMSLP_std[ilon,ilat]
		}
	}
}


## Read Precipitation Daily
varfilePRECIPdaily <- open.nc("Winter_CHR08.nc",write=FALSE)
varPRECIPdaily <-  var.get.nc(varfilePRECIPdaily, "precipitation")
#lat_arr <- var.get.nc(varfilePRECIP,"y")
#lon_arr <- var.get.nc(varfilePRECIP,"x") 
close.nc(varfilePRECIPdaily)

## Standarize Precipitation 
PRECIP.mean <- apply(varPRECIPdaily,1,mean)
PRECIP.std <- apply(varPRECIPdaily,1,sd)

var_precip <- array(NA,dim=c(precip_sum,nmod))
for (ni in 1:nmod){
	var_precip[,ni] <- (varPRECIPdaily[,ni]-PRECIP.mean)/PRECIP.std
}
#########################################################################################
###### Make 2D matrices
#### 2D standarized MSLP 
 lovalmsl<-array(FALSE,dim=c(nlon_msl,nlat_msl))
 for (ilon in 1:nlon_msl){
 for (ilat in 1:nlat_msl){
    if (sum(is.na(var_msl[ilon,ilat,]))==0){
        lovalmsl[ilon,ilat]<-TRUE
       }
     }
     }
msl_sum <- sum(lovalmsl)
msl2d <- array(0,dim=c(msl_sum,nmod))

 i <-0
 for (ilon in 1:nlon_msl){
 for (ilat in 1:nlat_msl){
    if (lovalmsl[ilon,ilat]){
       i <-i+1
       msl2d[i,1:nmod]= var_msl[ilon,ilat,1:nmod]
       }
     }
     }


#### 2D MSLP daily original
msl_sum <- sum(lovalmsl)
PressureDaily2d <- array(0,dim=c(msl_sum,nmod))

 i <-0
 for (ilon in 1:nlon_msl){
 for (ilat in 1:nlat_msl){
    if (lovalmsl[ilon,ilat]){
       i <-i+1
       PressureDaily2d[i,1:nmod]= MSLP.daily[ilon,ilat,1:nmod]
       }
     }
     }


#############
######Calculate covariance matrix for standarized (square matrix).
#############

M <- msl2d %*% t(var_precip)

############
######And then do the SVD
############

nu <- precip_sum
nv <- precip_sum


s <- svd(M,nu=nu,nv=nv)

######The covariance matrix M is reconstructed by s$u s$d t(s$v)

######
######Patterns of standarized msl for large scale fields. This is what we plot to get the circulation patterns (u3d_msl)
######
i <- 0
u3d_msl <- array(0,dim=c(nlon_msl,nlat_msl,precip_sum))
for (ilon in 1:nlon_msl){
for (ilat in 1:nlat_msl){
   if (lovalmsl[ilon,ilat]){
      i <- i+1
      u3d_msl[ilon,ilat,1:precip_sum] = s$u[i,1:precip_sum]
    } else {
     u3d_msl[ilon,ilat,1:precip_sum] = NA
    }
  }
  }

#######################################################
### Some checks
###fraction of explained coveriance per pattern pair
####This is also stored in s$d=> covariance matrix
###Find first total variance_totvar: each eigenvalue stored in(s$d) gives a measure of the fraction of the total variance in M expalined by the mode 
###Find the variance fraction_varfrac: divide each eigenvalues by the sum of all eigenvalues (trace of s$d)

totvar <- diag(t(s$u) %*% M %*% s$v)
varfrac <- totvar/sum(totvar)
write.table(totvar,"totvar.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(varfrac,"varfrac.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

####Time series of projection on patterns of pressure and precip
####Time evolution for standarized, average, daily

ts_msl <- t(msl2d) %*% s$u
ts_pressureDaily <- t(PressureDaily2d) %*% s$u
Patt_strengthdaily <- ts_pressureDaily/ norm(ts_pressureDaily)
Scaled_strengthdaily <- scale(ts_pressureDaily)

################################
Cmatrix<-s$u%*%diag(s$d)%*%t(s$v)
write.table(Cmatrix,"Cmatrix.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(M,"M.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(s$u,"s$u.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(s$v,"s$v.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(PressureDaily2d,"PressureDaily2d.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)