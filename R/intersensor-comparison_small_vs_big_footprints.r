### intersensor comparison script to compare single soundings
###
### sat1 is the satellite with the bigger footprint, sat2 the one with the smaller footprint
library(raster)
library(data.table)
library(ncdf4)
library(zoo)
library(solaR)
library(colorspace)
source("intersensor-comparison_helpers.r")


## range of dates to analyze:
date.range <- c(as.Date("2020-04-01"),as.Date("2020-04-30"))
###################

### 1) to collect only the relevant data from NC files:
## First satellite (bigger footprints): here TROPOMI
sat1.name <- "TROPOMI"
sat1.dict    <- list(sif="sif",sif.err="sif_err",
                     clat="lat",clon="lon",
                     sza="sza",vza="vza",
                     cont="NIR", dc="daily_correction_factor",
                     cf="cloud_fraction", time="TIME", pa="phase_angle")
sat1.dict.bnds <- list(lat_bnds="lat_bnds",lon_bnds="lon_bnds")
## needs to be subset of dict:
sat1.wishlist <- c("sif","sif.err","cont","clat","clon","vza","sza","pa","dc","time",paste0("lat",1:4),paste0("lon",1:4))
sat1.filter   <- "sza < 75."
sat1NCfiles   <- list.files("/net/fluo/data2/data/TROPOMI_SIF740nm/original",recursive=F,full.names=TRUE)
sat1NCdates   <- as.Date(basename(sat1NCfiles), format="TROPO_SIF_%Y-%m-%d") ##
sat1.time.origin <-  "1970-01-01 00:00:00 UTC"


###### Second Satellite (smaller footprints), here OCO-2 or OCO-3 
sat2.name <- "OCO-3"
#sat2NCfiles <- list.files("/net/fluo/data2/groupMembers/cfranken/data/kurosu/test_data/daily/oco2/2018", recursive=TRUE, full.names=TRUE, pattern ="*.nc4")
#sat2NCdates <- as.Date(basename(oco2NCfiles), format="oco2_LtSIF_%y%m%d") ## fixed
sat2NCfiles <- list.files("/net/fluo/data2/groupMembers/cfranken/data/kurosu/test_data/daily/oco3/2020", recursive=TRUE, full.names=TRUE, pattern ="*.nc4")
sat2NCdates <- as.Date(basename(sat2NCfiles), format="oco3_LtSIF_%y%m%d") ## fixed

##to check the content and define what to read in:
#nc <- nc_open(oco2NCfiles[1])
#nc_close(nc)
# SIF_740nm comes w/o uncertainty!
## OCO-2
# sat2.dict     <- list(sif="SIF_740nm",
#                       clat="Geolocation/latitude",clon="Geolocation/longitude",
#                       sza="Geolocation/solar_zenith_angle",saa="Geolocation/solar_azimuth_angle",
#                       vza="Geolocation/sensor_zenith_angle", vaa="Geolocation/sensor_azimuth_angle",
#                       qual="Science/sounding_qual_flag",mode="Metadata/MeasurementMode",
#                       cloud="Cloud/cloud_flag_abp",igbp="Science/IGBP_index", time="Geolocation/time_tai93")
### OCO-3
sat2.dict     <- list(sif="SIF_740nm",
                      clat="Geolocation/latitude",clon="Geolocation/longitude",
                      sza="Geolocation/solar_zenith_angle",saa="Geolocation/solar_azimuth_angle",
                      vza="Geolocation/sensor_zenith_angle", vaa="Geolocation/sensor_azimuth_angle",
                      qual="Science/sounding_qual_flag", qual2="Quality_Flag",mode="Metadata/MeasurementMode",
                      cloud="Cloud/cloud_flag_abp",igbp="Science/IGBP_index", time="Geolocation/time_tai93")

## TAI93 Time is seconds since 1993-01-01 00:00:00 UTC...
sat2.time.origin <-  "1993-01-01 00:00:00 UTC"
sat2.dict.bnds   <- list(lat_bnds="Latitude_Corners",lon_bnds="Longitude_Corners")

## needs to be subset of nc.dict:
sat2.wishlist <- c("sif","time","clat","clon","igbp","vza","sza","saa","vaa",paste0("lat",1:4),paste0("lon",1:4))
#sat2.filter   <- "cloud==0 & mode==0 & qual==0"
sat2.filter   <- "cloud==0 & mode==3 & qual==0 & qual2==0" ## !3 for area map!



############################# bundle stuff:
sat1.bundle <- list(dict=sat1.dict, dict.bnds=sat1.dict.bnds, files=sat1NCfiles, wishlist=sat1.wishlist, filter=sat1.filter, NCdates=sat1NCdates, sat1.time.origin=sat1.time.origin)

sat2.bundle <- list(dict=sat2.dict, dict.bnds=sat2.dict.bnds, files=sat2NCfiles, wishlist=sat2.wishlist, filter=sat2.filter, NCdates=sat2NCdates, sat2.time.origin=sat2.time.origin)

##### Actual co-location happens here (be patient, this will take a while)
jdata <- colocate_sat1_sat2(date.range,sat1.bundle,sat2.bundle)
##################


coloc.sat1 <- jdata$sat1
coloc.sat2 <- jdata$sat2
### check dimensions for consistency (should be: number of soundings sat1 < sat2)
if (dim(coloc.sat1)[1] >= dim(coloc.sat2)[1]) print("Oopsie, something went wrong, don't go ahead!") else if (dim(coloc.sat1)[1] < dim(coloc.sat2)[1]) print("So far, so good!")

### Because OCO-2/3 data comes w/o phase angle
coloc.sat2$pa <- compute_phase_angle(coloc.sat2)



################## Once the first part is done 
## compute difference in time of the measurements to narrow down co-located soundings wrt to difference in time and geometry:
coloc.sat1$hm  <- hour(coloc.sat1$utc)+minute(coloc.sat1$utc)/60
coloc.sat2$hm  <- hour(coloc.sat2$utc)+minute(coloc.sat2$utc)/60

################ SOME DECISIONS TO BE MADE HERE:!!!!!!!
## order of narrowing down may be changed arbitrarily (some computations are a bit slow)
## save copy to go back if criteria need to be re-defined?
#coloc.sat1.orig <- coloc.sat1
#coloc.sat2.orig <- coloc.sat2
## start here with
#coloc.sat1 <- coloc.sat1.orig
#coloc.sat2 <- coloc.sat2.orig

### start with defining acceptable difference in aquisition time:
coloc.sat2$dhm  <- unlist(lapply(1:dim(coloc.sat2)[1], function(i) coloc.sat2$hm[i]  - coloc.sat1$hm[which(coloc.sat1$sat1ID==coloc.sat2$sat1ID[i])]))
dev.new()
hist(coloc.sat2$dhm, breaks=100)
dhm.max <- 2. ## define acceptable difference in time [h]
sub.sat2  <- subset(coloc.sat2, abs(dhm) < dhm.max)

#### Filter for similar phase angles
sub.sat2$dpa <- unlist(lapply(1:dim(sub.sat2)[1], function(i) abs(sub.sat2$pa[i]) - abs(coloc.sat1$pa[which(coloc.sat1$sat1ID==sub.sat2$sat1ID[i])])))
dev.new()
hist(sub.sat2$dpa, breaks=100)
dpa.max <- 10. ## define acceptable difference in (absolute) phase angles [°]
sub.sat2  <- subset(sub.sat2, abs(dpa) < dpa.max)

#### Filter for similar VZAs
sub.sat2$dvza <- unlist(lapply(1:dim(sub.sat2)[1], function(i) sub.sat2$vza[i] - coloc.sat1$vza[which(coloc.sat1$sat1ID==sub.sat2$sat1ID[i])]))
dev.new()
hist(sub.sat2$dvza, breaks=100)
dvza.max <- 20. ## define acceptable difference in VZA [°]
sub.sat2  <- subset(sub.sat2, abs(dvza) < dvza.max)

#### Filter for similar VZAs
sub.sat2$dsza <- unlist(lapply(1:dim(sub.sat2)[1], function(i) sub.sat2$sza[i] - coloc.sat1$sza[which(coloc.sat1$sat1ID==sub.sat2$sat1ID[i])]))
dev.new()
hist(sub.sat2$dsza, breaks=100)
dsza.max <- 5. ## define acceptable difference in SZA [°]
sub.sat2  <- subset(sub.sat2, abs(dsza) < dsza.max)


####### add more/other filters here
####### For example: compute coverage of the bigger footprint
## potential strategy: create a polygons from smaller footprints and estimate the fractional coverage of the bigger one
## TBD
######################


### for the remaining soundings, look up the number of (smaller) soundings per (bigger) footprint 
u.sat1IDs          <- unique(sub.sat2$sat1ID)
sub.sat1           <- subset(coloc.sat1, sat1ID %in% u.sat1IDs)
sub.sat1$sat2SIF   <- unlist(lapply(u.sat1IDs, function(ID) mean(sub.sat2$sif[which(sub.sat2$sat1ID==ID)])))
sub.sat1$n.sat2    <- unlist(lapply(u.sat1IDs, function(ID) length(which(sub.sat2$sat1ID==ID))))

#sort(unique(sub.sat1$n.sat2))
#sort(unique(sub.sat1$n.sat2))

#### set minimum number of smaller footprints within the bigger one:
dev.new()
hist(sub.sat1$n.sat2)
n.min <- 1
final <- subset(sub.sat1, n.sat2>n.min) #

#####################################################################################
### quick check what it means to have several small footprints within one big footprint:
unique.n.footprints <- sort(unique(final$n.sat2))

idx.show      <- sample(which(final$n.sat2==max(unique.n.footprints)),1)
sat1ID.show   <- final$sat1ID[idx.show]
sps.sat1.show <- get_sps(subset(final, sat1ID==sat1ID.show))
sps.sat2.show <- get_sps(subset(sub.sat2, sub.sat2$sat1ID==sat1ID.show))

plot(sps.sat1.show, col="blue")
plot(sps.sat2.show, add=TRUE)
##########


### compute mean differences (if it is of interest)
#### it should make sense to look at mean differences of absolute values as -6° and 6° VZA difference would cancel out:
final$mdsza     <- unlist(lapply(final$sat1ID, function(ID) mean(abs(sub.sat2$dsza[which(sub.sat2$sat1ID==ID)]))))
final$mdvza     <- unlist(lapply(final$sat1ID, function(ID) mean(abs(sub.sat2$dvza[which(sub.sat2$sat1ID==ID)]))))
final$mdhm      <- unlist(lapply(final$sat1ID, function(ID) mean(abs(sub.sat2$dhm[which(sub.sat2$sat1ID==ID)]))))
final$mdpa      <- unlist(lapply(final$sat1ID, function(ID) mean(abs(sub.sat2$dpa[which(sub.sat2$sat1ID==ID)]))))
############


###########################################################
## order soundings latitude-wise (to avoid close to zero averaging):
idx.lat.ordered <- order(final$clat)

### aggregation levels of the bigger footprint for plot:
agg.facs <- c(2,5,10)
#names(final)
agg.dict <- list(sat1="sif", sat2="sat2SIF")
agg.sif  <- lapply(agg.facs, function(x,...) get_agg_mean(x, final, idx.lat.ordered, agg.dict))


########### Plot: !!!not all labels are not automized!!!
agg.cols <- c(terrain_hcl(length(agg.facs), h = c(265,80), c = c(60,10), l = c(25,95),power=c(1,1.2))[c(1,2)],rev(terrain_hcl(32, h = c(-100,100), c = c(60,100), l = c(15,95),power=c(2.,0.9)))[1])

dem.def  <- with(final,round(unlist(deming(sat2SIF,sif)),digits=2))
dem.a1   <- with(agg.sif[[1]],round(unlist(deming(sat2,sat1)),digits=2))
dem.a2   <- with(agg.sif[[2]],round(unlist(deming(sat2,sat1)),digits=2))
dem.a3   <- with(agg.sif[[3]],round(unlist(deming(sat2,sat1)),digits=2))

##usual linear regression (considering satellite with smaller footprints as 'truth' here)

p.def <- paste("y=",dem.def[1],"+",dem.def[2],"x",sep="")
p.a1  <- paste("y=",dem.a1[1],"+",dem.a1[2],"x",sep="")
p.a2  <- paste("y=",dem.a2[1],"+",dem.a2[2],"x",sep="")
p.a3  <- paste("y=",dem.a3[1],"+",dem.a3[2],"x",sep="")


leg.def <- substitute(paste(fit.eq," (", R^2,"=",r2,")",sep=""), list(fit.eq=p.def,r2=round(with(final,cor(sat2SIF,sif)^2),digits=2)))
leg.a1  <- substitute(paste(fit.eq," (", R^2,"=",r2,")",sep=""), list(fit.eq=p.a1,r2=round(with(agg.sif[[1]],cor(sat2,sat1)^2),digits=2)))
leg.a2  <- substitute(paste(fit.eq," (", R^2,"=",r2,")",sep=""), list(fit.eq=p.a2,r2=round(with(agg.sif[[2]],cor(sat2,sat1)^2),digits=2)))
leg.a3  <- substitute(paste(fit.eq," (", R^2,"=",r2,")",sep=""), list(fit.eq=p.a3,r2=round(with(agg.sif[[3]],cor(sat2,sat1)^2),digits=2)))


lim        <- c(-1.5,4.)
#title_tmp <- expression("April 2020, "*Delta*paste0("time <",dhm.max,"h, ")*Delta*"phase_angle < 10°")
title_tmp <- expression("April 2020, "*Delta*"time < 1h, "*Delta*"phase_angle < 10°")


dev.new()
#pdf(file="TROPO_vs_OCO-2_agg_Jun2018.pdf",width=6.5,height=6.)
par(oma=c(1,1,1,1),mar=c(3.,3.,3.5,1.5))
with(final, plot(sat2SIF,sif, col=makeTransparent("grey",alpha=0.5),xlim=lim,ylim=lim,xlab="",ylab="",pch=16,cex=3., main=title_tmp))#"#00000033"
abline(0,1,col="black",lwd=2)

with(agg.sif[[1]],  points(sat2, sat1, col=makeTransparent(agg.cols[1],alpha=0.8),xlim=lim,ylim=lim,xlab="",ylab="",pch=15,cex=2.))#"#00000033"
with(agg.sif[[2]],  points(sat2, sat1, col=makeTransparent(agg.cols[2],alpha=0.8),xlim=lim,ylim=lim,xlab="",ylab="",pch=17,cex=1.8))#"#00000033"
with(agg.sif[[3]],  points(sat2, sat1, col=makeTransparent(agg.cols[3],alpha=0.8),xlim=lim,ylim=lim,xlab="",ylab="",pch=18,cex=1.6))#"#00000033"

abline(dem.def[1:2],lty=2,col="darkgrey",lwd=2.5)
abline(dem.a1[1:2],lty=2,col=makeTransparent(agg.cols[1],alpha=1),lwd=2.5)
abline(dem.a2[1:2],lty=2,col=makeTransparent(agg.cols[2],alpha=1),lwd=2.5)
abline(dem.a3[1:2],lty=2,col=makeTransparent(agg.cols[3],alpha=1),lwd=2.5)



mtext(substitute("TROPOMI SIF@740 nm [mW/" ~ m^2 * "/sr/nm]"),2,2.5)
mtext(substitute("OCO-3 SIF@740 nm [mW/" ~ m^2 * "/sr/nm]"),1,2.5)

legend("bottomright",c("1:1 line",as.expression(leg.def),as.expression(leg.a1),as.expression(leg.a2),as.expression(leg.a3)),bty="n",lty=c(1,2,2,2),col=c("black","darkgrey",agg.cols),lwd=rep(2,4))
legend("topleft",c(paste0("Single ", sat1.name," soundings (with >= ",n.min," ",sat2.name," soundings)"),"Aggregation level: 2","Aggregation level: 5","Aggregation level: 10"),bty="n",pch=c(16,15,17,18),col=c("grey",agg.cols),pt.cex=rep(2,4))

dev.off()
###
