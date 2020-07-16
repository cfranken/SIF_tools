################################
### create TROPOMI time series for kml/kmz files and compute the fractional coverage of single TROPOMI footprints
################################
### standalone script, creates also a Quicklook
library(data.table)
library(ncdf4)
library(solaR)
library(rgdal)
library(raster)

sif.data.directory <- "/dir/to/TROPOMI_NCs_ungridded"
kmz.data.dir       <- "/dir/to/kmz_files"
out.dir            <- "/dir/to/output"


inp.files <- grep("_ungridded.nc",list.files(sif.data.directory,pattern="TROPO_SIF_20",recursive=TRUE, full.names=TRUE),value=TRUE)


## list of locations to process:
kmz.files   <- list.files(kmz.data.dir, pattern="kmz", full.names=TRUE)
locationIDs <- substr(basename(kmz.files),1,nchar(basename(kmz.files))-4)
####
## quicklook will show single soundings and averages for: 
agg.times <- list("4 days", "1 week", "10 days", "2 weeks", "20 days", "1 month")
########


####################### A few necessary functions
read.TROPO.nc <- function(inp.file,box=box,cf.min=0,cf.max=0.8,vza.min=0,vza.max=60,pa.min=0,pa.max=180){
  ## read data partially for filtering, read bnds later
  nc <- try(nc_open(inp.file),silent=TRUE)
  if (class(nc)!="try-error"){
    clon    <- ncvar_get(nc,"lon")
    clat    <- ncvar_get(nc,"lat")
    fs      <- ncvar_get(nc,"sif")
    vza     <- ncvar_get(nc,"vza")
    pa      <- ncvar_get(nc,"phase_angle")
    cf      <- ncvar_get(nc,"cloud_fraction")
  nc_close(nc)

  idx.filter <- which(clon >= box[1] & clon <= box[2] & clat >= box[3] & clat <= box[4] & cf >= cf.min & cf <= cf.max & vza >= vza.min & vza <= vza.max & abs(pa) >= pa.min & abs(pa) <= pa.max)#
  if (length(idx.filter) == 0) return(NULL)
  if (length(idx.filter) == 1){
    nc <- nc_open(inp.file)
      lonb    <- ncvar_get(nc,"lon_bnds")[idx.filter,]
      latb    <- ncvar_get(nc,"lat_bnds")[idx.filter,]
      time    <- ncvar_get(nc,"TIME")[idx.filter]
      sza     <- ncvar_get(nc,"sza")[idx.filter]
      corr    <- ncvar_get(nc,"daily_correction_factor")[idx.filter]
      sif_err <- ncvar_get(nc,"sif_err")[idx.filter]
    nc_close(nc)
    ret <- data.table(time=time,clon=clon[idx.filter],clat=clat[idx.filter],sif=fs[idx.filter],sif.err=sif_err,sza=sza,vza=vza[idx.filter],pa=pa[idx.filter],cf=cf[idx.filter],
		      lon1=lonb[1],lon2=lonb[2],lon3=lonb[3],lon4=lonb[4],lat1=latb[1],lat2=latb[2],lat3=latb[3],lat4=latb[4],corr=corr)
    return(ret)
    }
  if (length(idx.filter) >= 2){
    nc <- nc_open(inp.file)
      lonb    <- ncvar_get(nc,"lon_bnds")[idx.filter,]
      latb    <- ncvar_get(nc,"lat_bnds")[idx.filter,]
      time    <- ncvar_get(nc,"TIME")[idx.filter]
      sza     <- ncvar_get(nc,"sza")[idx.filter]
      corr    <- ncvar_get(nc,"daily_correction_factor")[idx.filter]
      sif_err <- ncvar_get(nc,"sif_err")[idx.filter]
    nc_close(nc)
    ret <- data.table(time=time,clon=clon[idx.filter],clat=clat[idx.filter],sif=fs[idx.filter],sif.err=sif_err,sza=sza,vza=vza[idx.filter],pa=pa[idx.filter],cf=cf[idx.filter],
		      lon1=lonb[,1],lon2=lonb[,2],lon3=lonb[,3],lon4=lonb[,4],lat1=latb[,1],lat2=latb[,2],lat3=latb[,3],lat4=latb[,4],corr=corr)
    return(ret)
    }
  }else{
  return(NULL)}
}


divLine <- function(lat1,lon1,lat2,lon2,n){
    dLat     <- (lat2-lat1)/(n-1)
    dLon     <- (lon2-lon1)/(n-1)
	  lats     <- c(lat1, lat1 + seq(1,n-1)*dLat)
    lons     <- c(lon1, lon1 + seq(1,n-1)*dLon)
    return(list(lats=lats,lons=lons))
}

# Divide each polygon into multiple points
getPoints <- function(vert_lat, vert_lon, n){
	# Get reference points for two lines at the extremes:
	  line_1 <- divLine(vert_lat[1],vert_lon[1],vert_lat[2],vert_lon[2],n)
    line_2 <- divLine(vert_lat[4],vert_lon[4],vert_lat[3],vert_lon[3],n)
    lats <- logical(0)
    lons <- logical(0)
    for (i in 1:n){
      tmp <- divLine(line_1$lats[i], line_1$lons[i] ,line_2$lats[i], line_2$lons[i],n)
      lats <- c(lats,tmp$lats)
      lons <- c(lons,tmp$lons)
    }
    #lapply(1:n, function(i,...) divLine(line_1$lats[i], line_1$lons[i] ,line_2$lats[i], line_2$lons[i],n))
    return(list(lats=lats,lons=lons))
}

### to generate quicklook:
get.means <- function(ag.level,site.data,t.vec,...){
t.range <- range(t.vec)
t.seq   <- seq(as.Date(t.range[1]),as.Date(t.range[2]),by=ag.level)

t.start                    <- t.seq[1:(length(t.seq)-1)]
t.start[2:length(t.start)] <- t.start[2:length(t.start)] + 1 #adjust times + 1 sec
t.end                      <- t.seq[2:length(t.seq)]

time <- NA
sif  <- NA
n    <- NA
for (i.t in 1:length(t.start)){
#i.t <- 1
time[i.t] <- mean(c(t.start[i.t],t.end[i.t]))
idx.sif   <- which(t.vec >= t.start[i.t] & t.vec <= t.end[i.t])
##inst. values:
# sif[i.t] <- mean(site.data$sif[idx.sif])
##daily corrected:
sif[i.t] <- mean(site.data$sifDC[idx.sif])
n[i.t]   <- length(idx.sif)
}
ret <- list(time=time, sif=sif, n=n)
return(ret)
}
#############################



dlon <- 0.2
dlat <- 0.2

for (i.site in seq_along(locationIDs)){
#i.site <- 1
locationID <- locationIDs[i.site]

print(paste("Processing ",locationID," (",i.site,"/",length(locationIDs),")",sep=""))
print("********************************************************************")


##def N-E-S-W pixel:
tpoly <- readOGR(kmz.files[i.site])
box   <- c(slot(extent(tpoly), "xmin")-dlon, slot(extent(tpoly), "xmax")+dlon, slot(extent(tpoly), "ymin")-dlat, slot(extent(tpoly), "ymax")+dlat)

siteData <- data.table()
# looping through nc files...

pb <- txtProgressBar(min = 0, max = length(inp.files), style = 3) ##intialize progress bar
for (i in 1:length(inp.files)){
setTxtProgressBar(pb, i) ##update progress bar
tmp <- read.TROPO.nc(inp.files[i],box=box)
if (is.null(tmp) == FALSE){
      l         <- list(siteData,tmp)
	    siteData  <- rbindlist(l)
}
}
close(pb) ##close progress bar

tmpUTC  <- as.POSIXct(siteData$time,origin = '1970-01-01 00:00:00',tz="UTC")
tmpMST  <- local2Solar(tmpUTC,lon=siteData$clon)

siteData$date <- as.character(as.Date(tmpUTC, "%m/%d/%Y"))
siteData$lst  <- hour(tmpMST)+minute(tmpMST)/60


## strategy is to divide every measurement into 100 and count the points inside the shape file!
print("Calc fractional coverage of single TROPOMI footprints hitting the spatial extent...")
n.poly <- dim(siteData)[1]
#set n=10 to get percentage right away:
n <- 10

###### get coordinates of shapefile:
lonlat <- lapply(slot(tpoly, "polygons"), function(x) lapply(slot(x,"Polygons"), function(y) slot(y, "coords")))
pol.x <- lonlat[[1]][[1]][,1]
pol.y <- lonlat[[1]][[1]][,2]


siteData$per <- rep(NA,n.poly)
for (i.poly in 1:n.poly){
vert_lat <- c(siteData$lat1[i.poly],siteData$lat2[i.poly],siteData$lat3[i.poly],siteData$lat4[i.poly])
vert_lon <- c(siteData$lon1[i.poly],siteData$lon2[i.poly],siteData$lon3[i.poly],siteData$lon4[i.poly])
p        <- getPoints(vert_lat, vert_lon, n)
siteData$per[i.poly] <- length(which(unlist(lapply(1:n^2, function(i,...) point.in.polygon(p$lons[i], p$lats[i], pol.x, pol.y)))==1))
}

sub.siteData <- subset(siteData, per!=0)
#dev.new()
#hist(sub.siteData$per, breaks=100)
######################################


file.out <- paste0(out.dir,"/",locationID,".txt")
cat("DATE LST CLON CLAT SIF SIFerr PER DC SZA VZA PA CF LON1 LON2 LON3 LON4 LAT1 LAT2 LAT3 LAT4\n", file=file.out)
write.table(format(subset(sub.siteData,select=c("date","lst","clon","clat","sif","sif.err","per","corr","sza","vza","pa","cf","lon1","lon2","lon3","lon4","lat1","lat2","lat3","lat4")),scientific=FALSE),file=file.out,append=T, row.names = F,col.names = F, quote = F)

sub.siteData$sifDC <- sub.siteData$sif*sub.siteData$corr
t.vec              <- as.Date(sub.siteData$date)
t.range            <- range(t.vec)
names(agg.times)   <- letters[1:length(agg.times)]
agg                <- lapply(agg.times,FUN=function(x) get.means(ag.level=x,site.data=sub.siteData,t.vec=t.vec))

### Plotting
cols <- rainbow(length(agg))
#dev.new(width=15,height=6)
pdf(file=paste(out.dir,"/",locationID,"_SIFdc_TROPOMI_quicklook.pdf",sep=""),width=15,height=6)
plot(t.vec,sub.siteData$sifDC,type="b",pch=4,lwd=0.3,ylim=range(sub.siteData$sifDC),yaxt="n",xaxt="n",xlab="",ylab="",cex=1.)
abline(h=0,col="grey")
for (i in 1:length(agg)){
ag.lev <- names(agg)[i]
lines(as.POSIXlt(eval(parse(text=paste("agg$",ag.lev,"$time",sep=""))), origin="1970-01-01"),eval(parse(text=paste("agg$",ag.lev,"$sif",sep=""))),col=cols[i],type="l",lwd=2.,cex=1.)
}
legend("top",c("single meas.",unlist(agg.times)),fill=c("black",cols),bty="n",ncol=length(agg)+1,cex=1.2)
mtext(substitute("daily-corrected SIF [mW/" ~ m^2 * "/sr/nm]"),2,2.2,cex=1.5)
mtext("Date [Month/Year]",1,3.2,cex=1.5)
axis(1, at=pretty(t.vec), labels=format(pretty(t.vec), "%m/%Y"),cex.axis=1.5)
axis(2, cex.axis=1.5)
title(paste("Aggregation levels of daily-corrected SIF at ",locationID,sep=""))
dev.off()
##############################
}
