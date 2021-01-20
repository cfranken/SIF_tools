### script to create a map by temporal oversampling for specified days
library(raster)
library(data.table)
library(ncdf4)
library(zoo)
library(solaR)
library(colorspace)
library(viridis)
library(fields)

setwd("/home/pkoehler/Documents/Caltech/git_stuff/SIF_tools/R")
source("intersensor-comparison_helpers.r")

#Task: the “Tanana West” zone: 64.6 N – 64.8 N and -148.5 W – -147.7 W 
#all sampling was done from July 6-8 2017, so perhaps a period of time roughly centered on that would be good for examining TROPOMI data from 2018, ’19, and ’20?
box <- c(-148.5,-147.7,64.6,64.8)
# choose 1 month around 7 July for 2018/29/20:
init.date <- "-07-04"
years     <- c("2018","2019","2020")

return.date.seq <- function(init.date,delta.day) seq(init.date-delta.day,init.date+delta.day,by="1 day")

all.days <- unlist(lapply(years, function(a,...) return.date.seq(as.Date(paste0(a,init.date)),4)))
# as.Date(all.days)


## specify stuff to collect data:
sat1.name <- "TROPOMI"
sat1.dict    <- list(sif="sif",sif.err="sif_err",
                     clat="lat",clon="lon",
                     sza="sza",vza="vza",
                     cont="NIR", dc="daily_correction_factor",
                     cf="cloud_fraction", time="TIME", pa="phase_angle")
sat1.dict.bnds <- list(lat_bnds="lat_bnds",lon_bnds="lon_bnds")
## needs to be subset of dict:
sat1.wishlist <- c("sif","sif.err","cont","clat","clon","vza","sza","pa","dc","time",paste0("lat",1:4),paste0("lon",1:4))
## “Tanana West” zone: 64.6 N – 64.8 N and -148.5 W – -147.7 W , but use +/- 0.2 deg to collect all measurements which might have overlap:
sat1.filter   <- "sza < 75. & clat > 64.4 & clat < 65. & clon > -148.7 & clon < -147.5"

sat1NCfiles   <- list.files("/net/fluo/data2/data/TROPOMI_SIF740nm/original",recursive=F,full.names=TRUE)
sat1NCdates   <- as.Date(basename(sat1NCfiles), format="TROPO_SIF_%Y-%m-%d") ##
## subset relevant dates:
idx.sub <- which(sat1NCdates %in% all.days)
sat1NCfiles   <- sat1NCfiles[idx.sub]
sat1NCdates   <- sat1NCdates[idx.sub]

sat1.time.origin <-  "1970-01-01 00:00:00 UTC"

sat1.bundle <- list(dict=sat1.dict, dict.bnds=sat1.dict.bnds, files=sat1NCfiles, wishlist=sat1.wishlist, filter=sat1.filter, NCdates=sat1NCdates, sat1.time.origin=sat1.time.origin)

dd <- collect_dt_NC(sat1.bundle)

######## now mapping by oversampling:
### using some old code here, could be done more efficient and works only for small regions:
r      <- raster(xmn=box[1],xmx=box[2],ymn=box[3],ymx=box[4])
res(r) <- 0.005

print("********")
print("Polygon-Rastering...")
#print(names(dd))
idx     <- 1:length(dd$lon1)
counter <- 0
out <- NA
for (i.len in 1:length(dd$lon1)){
counter <- counter+1
poly.coord <- rbind(c(dd$lon1[i.len],dd$lat1[i.len]),c(dd$lon2[i.len],dd$lat2[i.len]),c(dd$lon3[i.len],dd$lat3[i.len]),c(dd$lon4[i.len],dd$lat4[i.len]))
poly.coord <- rbind(poly.coord,poly.coord[1,])
## polygon:
p.name <- paste("p",i.len, sep="")
assign(p.name,(Polygon(poly.coord,hole=F)))
##polygons:
ps.name <- paste("ps",i.len, sep="")
assign(ps.name,(Polygons(list(eval(parse(text=p.name))),ID=paste("s",i.len,sep=""))))
rm(list=p.name)
}
#out       <- out[is.na(out)==FALSE]
#if (length(out) != 0) idx.good  <- idx[-out] else idx.good  <- idx
idx.good  <- idx

sps   <- SpatialPolygons(mget(paste("ps",idx.good,sep=""),environment()), proj4string=CRS("+proj=longlat +ellps=WGS84"))
#aufräumen:
ps.names <- paste("ps",idx,sep="")
rm(list=c(ps.names,"ps.names","ps.name","p.name"))
tropo.spdf <- SpatialPolygonsDataFrame(sps, data.frame(fs=(dd$sif[idx.good]*dd$dc[idx.good]), row.names=paste("s",idx.good,sep="")))
tropo.sldf <- as(tropo.spdf,"SpatialLinesDataFrame") # Create a lines object for the borders of the polygons
## and rasterizing...
polygon.raster <- rasterize(tropo.sldf,r,'fs',fun=function(x,...) mean(na.omit(x)))
polygon.raster <- rasterize(tropo.spdf,polygon.raster,update = T,'fs',fun=function(x,...) mean(na.omit(x)))


dev.new(height=6, width=8)
lim <- c(0,cellStats(polygon.raster,max))
image.plot(polygon.raster,zlim=lim, col=viridis(64), main=" 0.005 deg resolution")
#map("worldHires",interior=F, add=T,myborder = 0.0000001 ,col="grey39", lwd=0.2)#"grey80"

writeRaster(polygon.raster,"/net/fluo/data2/groupMembers/pkoehler/collected_datasets/ABoVE_Tanana_West_0005deg_sifdc_AndyM.tif", format="GTiff", overwrite=TRUE) 



############################################################################################# 
#### Faster gridding:
### !!!wayyyy faster than the rasterize function!!!
### idea: create polygon from measurement and check for all grid cells which are inside with point.in.polygon
gridres <- 0.005
lons <- seq(box[1]+0.5*gridres,box[2]-0.5*gridres,by=gridres)
lats <- seq(box[3]+0.5*gridres,box[4]-0.5*gridres,by=gridres)
grid <- matrix(NA,nrow=length(lons),ncol=length(lats))

## create gridcell information:
pGrid <- list()
counter <- 0
for (iLon in 1:length(lons)){
    for (iLat in 1:length(lats)){
        counter <- counter + 1
        pGrid$lons[counter] <- lons[iLon] 
        pGrid$lats[counter] <- lats[iLat]
        pGrid$iLat[counter] <- iLat   
        pGrid$iLon[counter] <- iLon 
    }
}

# #### working example:
# i.poly <- 99
# pol.x <- c(dd$lon1[i.poly],dd$lon2[i.poly],dd$lon3[i.poly],dd$lon4[i.poly])
# pol.y <- c(dd$lat1[i.poly],dd$lat2[i.poly],dd$lat3[i.poly],dd$lat4[i.poly])

# ## which grid cells might be covered by the polygon?
# idxGrid <- which(pGrid$lons >= range(pol.x)[1] & pGrid$lons <= range(pol.x)[2] &
#                  pGrid$lats >= range(pol.y)[1] & pGrid$lats <= range(pol.y)[2])
# idxHit <- which(unlist(lapply(idxGrid, function(i,...) point.in.polygon(pGrid$lons[i], pGrid$lats[i], pol.x, pol.y)))==1)

# emptyGrid <- grid
# emptyGrid[] <- 0
# test.poly <- get_sps(dd[i.poly,])
# image(lons,lats,emptyGrid)
# abline(v=lons,col="grey")
# abline(h=lats,col="grey")
# points(pGrid$lons[idxGrid],pGrid$lats[idxGrid])
# points(pGrid$lons[idxGrid[idxHit]],pGrid$lats[idxGrid[idxHit]],col="red",pch="x")
# plot(test.poly,add=T)
# #test <- retXY(pGrid,pol.x,pol.y)
# #points(lons[test$x],lats[test$y],col="red",pch="x")
# #########################################################

retXY <- function(pGrid,pol.x,pol.y){
    idxGrid <- which(pGrid$lons >= range(pol.x)[1] & pGrid$lons <= range(pol.x)[2] &
                     pGrid$lats >= range(pol.y)[1] & pGrid$lats <= range(pol.y)[2])
    idxHit <- which(unlist(lapply(idxGrid, function(i,...) point.in.polygon(pGrid$lons[i], pGrid$lats[i], pol.x, pol.y)))==1)
    return(list(x=pGrid$iLon[idxGrid[idxHit]],y=pGrid$iLat[idxGrid[idxHit]]))
}

sgrid   <- ngrid <- grid
ngrid[] <- 0
sgrid[] <- 0 
n.poly <- dim(dd)[1]
for (i.poly in 1:n.poly){
    pol.x <- c(dd$lon1[i.poly],dd$lon2[i.poly],dd$lon3[i.poly],dd$lon4[i.poly])
    pol.y <- c(dd$lat1[i.poly],dd$lat2[i.poly],dd$lat3[i.poly],dd$lat4[i.poly])

    idxXY <- retXY(pGrid,pol.x,pol.y)
    nIDX  <- length(idxXY$x)
    if (nIDX > 0){
        for (iIDX in 1:nIDX){
               ngrid[idxXY$x[iIDX],idxXY$y[iIDX]] <- ngrid[idxXY$x[iIDX],idxXY$y[iIDX]] + 1#nIDX
               sgrid[idxXY$x[iIDX],idxXY$y[iIDX]] <- sgrid[idxXY$x[iIDX],idxXY$y[iIDX]] + dd$sif[i.poly]*dd$dc[i.poly]
        }
     }
}

ngrid[ngrid==0] <- NA
sgrid <- sgrid/ngrid
lim <- c(0,0.8)
dev.new()
image.plot(lons,lats,sgrid, zlim=lim,col=viridis(64),
           main=paste0("+/-4 days (2018) , n=",dim(dd)[1])) #

t.sgrid <- t(sgrid[,dim(sgrid)[2]:1])
test.raster <- raster(t.sgrid,xmn=box[1],xmx=box[2],ymn=box[3],ymx=box[4])
dev.new()
image.plot(test.raster, zlim=lim,col=viridis(64)) #



####### just to compare end result to slow rasterize function:
# dev.new()
# image.plot(polygon.raster, zlim=lim,col=viridis(64)) #

# diff.raster <- test.raster-polygon.raster
# lim.val <- cellStats(abs(diff.raster),max)

# dev.new()
# image.plot(diff.raster,zlim=c(-lim.val,lim.val),col=diverge_hcl(11)) #
#########