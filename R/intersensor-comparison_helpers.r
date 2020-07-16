### necessary functions to read data and perform satellite intersensor comparisons:

cosd<-function(degrees){
  radians<-cos(degrees*pi/180)
  return(radians)
}

sind<-function(degrees){
  radians<-sin(degrees*pi/180)
  return(radians)
}

collect_dt_NC <- function(bundle,...){

  ncfiles      <- bundle$files
  nc.dict      <- bundle$dict
  nc.dict.bnds <- bundle$dict.bnds
  sel.list     <- bundle$wishlist
  filter       <- bundle$filter

  dt.ts <- list()

  pb <- txtProgressBar(min = 0, max = length(ncfiles), style = 3) ##intialize progress bar
  for (ncfile in ncfiles){
      setTxtProgressBar(pb, which(ncfiles==ncfile)) ##update progress bar
      tmp <- list()
      nc <- try(nc_open(ncfile),silent=TRUE)
        if (class(nc)!="try-error"){
          for (i.var in 1:length(nc.dict)) eval(parse(text=paste0("tmp$",names(nc.dict)[i.var] , " <- ncvar_get(nc, nc.dict$",names(nc.dict)[i.var],")")))
          ## read boundaries
          for (i.bnd in 1:length(nc.dict.bnds)) eval(parse(text=paste0(names(nc.dict.bnds)[i.bnd], " <- ncvar_get(nc, nc.dict.bnds$",names(nc.dict.bnds)[i.bnd],")")))
          nc_close(nc)
          ##attach boundaries to tmp:
          idx.lat <- grep("lat",nc.dict.bnds, ignore.case=TRUE)
          idx.lon <- grep("lon",nc.dict.bnds, ignore.case=TRUE)
          ## dimension for boundaries should be 4, but the dimension can change:
          bnd.dim <- dim(eval(parse(text=names(nc.dict.bnds)[1])))
          if (bnd.dim[1] == 4){
             bnd.dim.struc <- function(i.bnd) return(paste0("[",i.bnd,",]"))
           }else if (bnd.dim[2] == 4){
             bnd.dim.struc <- function(i.bnd) return(paste0("[,",i.bnd,"]"))
           }
          for (i.bnd in 1:4){
            eval(parse(text=paste0("tmp$lat", i.bnd, " <- ", names(nc.dict.bnds)[idx.lat],bnd.dim.struc(i.bnd))))
            eval(parse(text=paste0("tmp$lon", i.bnd, " <- ", names(nc.dict.bnds)[idx.lon],bnd.dim.struc(i.bnd))))
          }
          tmp  <- as.data.table(tmp)
          ### abstracted to switch between filter criteria:
          tmp  <- subset(tmp, eval(parse(text=filter)), select=sel.list) ##
          ##may be empty:
          if (dim(tmp)[1]>1){
          ## 2) Filter for special surfaces?
          #tmp$training.code <- extract(ref.area.ras,cbind(tmp$clon,tmp$clat))
          #tmp  <- subset(tmp, is.na(training.code)==FALSE & training.code %in% c(2,1,0,-2))
          #if (dim(tmp)[1]>1){
              l           <- list(dt.ts,tmp)
              dt.ts      <- rbindlist(l)
          #}
          }
          }
  }
  close(pb) ##close progress bar
  return(dt.ts)
}

compute_phase_angle <- function(sat){
  ## Input is satellite data in data.table format with: sza, vza, saa, vaa
  ## Phase angles are set to negative values if the observational azimuth angle is bigger than the solar azimuth angle
  ## (negative phase angle means the sun is to the right of the satellite)
  necessary.names <- c("saa","sza","vaa","vza")

  if (all(necessary.names %in% names(sat))){

  vaa <- sat$vaa ## viewing azimuth angle
  vza <- sat$vza ## viewing zenith angle
  sza <- sat$sza ## solar zenith angle
  saa <- sat$saa ## solar azimuth angle

  pa  <- phase <- raa <- rep(NA, length(sza))

  phase[vaa > saa] <- -1.
  phase[vaa < saa] <-  1.
  raa              <- vaa - saa

  idx <- which(raa < -180.)
  if (length(idx) > 0) raa[idx] <- raa[idx]+360.
  idx <- which(raa > 180.)
  if (length(idx) > 0) raa[idx] <- raa[idx]-360.
  raa <- abs(raa)

  pa  <- acos(cosd(sza)*cosd(vza)+sind(vza)*sind(sza)*cosd(raa))*180./pi
  pa  <- pa * phase
  return(pa)

} else {
  print("!!! Necessary input is missing, function returns NULL !!!")
  return(NULL)
}
}
##########################

return_sign <- function(x){
  if (x > 0.) return(1)
  if (x < 0.) return(-1)
}

get_agg_mean <- function(agg.fac, dd.in, idx.vec, agg.dict){
  print(paste("Aggregation level:",agg.fac))
  n.agg       <- floor(length(idx.vec)/agg.fac)
  dd.out      <- data.frame(sat1=rep(NA,n.agg),sat2=rep(NA,n.agg))
  idx2av      <- idx.vec[1:agg.fac]

  dd.out$sat1[1]  <- mean(eval(parse(text=paste0("dd.in[idx2av, ]$",agg.dict$sat1))))
  dd.out$sat2[1]  <- mean(eval(parse(text=paste0("dd.in[idx2av, ]$",agg.dict$sat2))))

  for (i.step in 1:(n.agg-1)){
  idx2av                <- idx.vec[(i.step*agg.fac+1):(i.step*agg.fac+1+(agg.fac-1))]
  dd.out$sat1[i.step+1] <- mean(eval(parse(text=paste0("dd.in[idx2av, ]$",agg.dict$sat1))))
  dd.out$sat2[i.step+1] <- mean(eval(parse(text=paste0("dd.in[idx2av, ]$",agg.dict$sat2))))
  }
  return(dd.out)
}

deming <- function(x,y){
##Deming regression:
  #delta <- var(sigma_y)/var(sigma_x)
  delta <- 1.
  s_xx  <- 1./(length(x)-1.)*sum((x-mean(x))^2.)
  s_yy  <- 1./(length(y)-1.)*sum((y-mean(y))^2.)
  s_xy  <- 1./(length(y)-1.)*sum((x-mean(x))*(y-mean(y)))
  slope <- (s_yy-delta*s_xx+sqrt((s_yy-delta*s_xx)^2.+4.*delta*s_xy^2.))/(2.*s_xy)
  intercept <- mean(y)-slope*mean(x)
  ##R_sq:
  y_fit  <- slope*x+intercept
  ss_tot <- sum((y-mean(y))^2.)	#total sum of squares
  ss_err <- sum((y-y_fit)^2.)	#residual sum of squares
  r_sq 	 <- 1-ss_err/ss_tot
  demreg <-list(intercept=intercept,slope=slope,r.sq=r_sq)
  return(demreg)
}

makeTransparent = function(..., alpha=0.5){
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")

  alpha = floor(255*alpha)
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)

  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }

  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)

  return(newColor)

}

colocate_sat1_sat2 <- function(date.range, sat1.bundle, sat2.bundle){
## loop through given date range and look for close soundings
## goal is to get two data sets to be merged in a second step, later on: filter for close time and geometry...
ref.ras      <- raster()
res(ref.ras) <- 0.5

### sat1=bigger footprints, sat2=smaller footprints, 
### important because IDs are attached to the satellite with the bigger footprint
coloc.sat1 <- data.table()
coloc.sat2 <- data.table()
sat1ID     <- 0

date.seq <- seq(date.range[1],date.range[2],by="1 day")

## loop through single days:
for (daydate in date.seq){
#  daydate <- date.seq[1]
  print(paste0("Reading: ", as.Date(daydate),", found so far:", sat1ID))
  #### Sat1 (bigger footprints)
  idx.sat1.ncfiles      <- which(sat1.bundle$NCdates == daydate)
  sat1.bundle.tmp       <- sat1.bundle
  sat1.bundle.tmp$files <- sat1.bundle$files[idx.sat1.ncfiles]
  sat1                  <- collect_dt_NC(sat1.bundle.tmp)

  #### Sat2 (smaller footprints)
  idx.sat2.ncfiles      <- which(sat2.bundle$NCdates == daydate)
  sat2.bundle.tmp       <- sat2.bundle
  sat2.bundle.tmp$files <- sat2.bundle$files[idx.sat2.ncfiles]
  sat2                  <- collect_dt_NC(sat2.bundle.tmp)
  if (length(dim(sat2)[1]) != 0){

  print(paste0("Checking for overlapping regions..."))
  spodf  <- SpatialPointsDataFrame(cbind(sat1$clon,sat1$clat),data.frame(dummy=rep(1,length(sat1$clat))),proj4string=CRS("+proj=longlat +ellps=WGS84"),match.ID=F)
  raster.sat1 <- rasterize(spodf, ref.ras ,'dummy',fun=function(x,...) mean(na.omit(x)))

  idx.rough.hit <- which(is.na(extract(raster.sat1,cbind(sat2$clon,sat2$clat)))==FALSE)
  sat2          <- sat2[idx.rough.hit,]

  ## the same for the other satellite:
  spodf <- SpatialPointsDataFrame(cbind(sat2$clon,sat2$clat),data.frame(dummy=rep(1,length(sat2$clat))),proj4string=CRS("+proj=longlat +ellps=WGS84"),match.ID=F)
  raster.sat2 <- rasterize(spodf, ref.ras ,'dummy',fun=function(x,...) mean(na.omit(x)))
  raster.comb <- raster.sat1*raster.sat2

  idx.rough.hit <- which(is.na(extract(raster.comb,cbind(sat1$clon,sat1$clat)))==FALSE)
  sat1          <- sat1[idx.rough.hit,]

  ## filter soundings crossing the dateline (causing potential issues):
  ## aka: make sure that the sign of lon is the same for all lons
  idx.drop <- NULL
  for (i.sdng in 1:dim(sat1)[1]){
  sign.tmp <-  unlist(lapply(1:4, function(i.bnd) eval(parse(text=paste0("return_sign(sat1$lon",i.bnd,"[i.sdng])")))))
  if ((all(sign.tmp==1) | all(sign.tmp==-1))==FALSE) idx.drop <- c(idx.drop,i.sdng)
  }
  if (length(idx.drop)>=1) sat1 <- sat1[-idx.drop,]

  idx.drop <- NULL
  for (i.sdng in 1:dim(sat2)[1]){
  sign.tmp <-  unlist(lapply(1:4, function(i.bnd) eval(parse(text=paste0("return_sign(sat2$lon",i.bnd,"[i.sdng])")))))
  if ((all(sign.tmp==1) | all(sign.tmp==-1))==FALSE) idx.drop <- c(idx.drop,i.sdng)
  }
  if (length(idx.drop)>=1) sat2 <- sat2[-idx.drop,]

  ## attach time in UTC:
  sat1$utc <- as.POSIXct(sat1$time, origin=sat1.time.origin,tz="UTC")
  sat2$utc <- as.POSIXct(sat2$time, origin=sat2.time.origin ,tz="UTC")

  n.poly  <- dim(sat1)[1]
  print(paste0("Co-location..."))
  pb <- txtProgressBar(min = 0, max = n.poly, style = 3) ##intialize progress bar
  for (i.poly in 1:n.poly){
    setTxtProgressBar(pb, i.poly) ##update progress bar
    #i.poly <- 1
    pol.x          <- c(sat1$lon1[i.poly],sat1$lon2[i.poly],sat1$lon3[i.poly],sat1$lon4[i.poly])
    pol.y          <- c(sat1$lat1[i.poly],sat1$lat2[i.poly],sat1$lat3[i.poly],sat1$lat4[i.poly])

    tmp.hit <- point.in.polygon(sat2$clon,sat2$clat,pol.x,pol.y)
    hit.idx <- which(tmp.hit==1)
    if (length(hit.idx) >= 1){
      ## create joint data set here if there are indeed soundings within the bigger footprint for that specific day...
      sat1ID <- sat1ID+1

      coloc.sat1.tmp        <- sat1[i.poly,]
      coloc.sat1.tmp$sat1ID <- sat1ID
      coloc.sat2.tmp        <- sat2[hit.idx,]
      coloc.sat2.tmp$sat1ID <- rep(sat1ID,length(hit.idx))

      l           <- list(coloc.sat1,coloc.sat1.tmp)
      coloc.sat1 <- rbindlist(l)

      l           <- list(coloc.sat2,coloc.sat2.tmp)
      coloc.sat2  <- rbindlist(l)
    }
  }
  close(pb) ##close progress bar
  ## store intermediate results?:
  #jdata <- list(tropo=coloc.tropo,gome2=coloc.gome2)
  #save.name <- paste0(out.dir,"/co-located_tropo-gome2_",ac.date,".RData")
  #print(paste0("Saving joint data set at: ", save.name))
  #save(jdata, file=save.name)
}
}
return(list(sat1=coloc.sat1,sat2=coloc.sat2))
}

## to check coverage of several soundings
get_sps <- function(df){
######## Create spatial polygons from data frame###############
counter <- 0
idx     <- 1:length(df$clat)
for (i.len in idx){
counter <- counter+1
poly.coord <- rbind(c(df$lon1[i.len],df$lat1[i.len]),c(df$lon2[i.len],df$lat2[i.len]),c(df$lon3[i.len],df$lat3[i.len]),c(df$lon4[i.len],df$lat4[i.len]))
poly.coord <- rbind(poly.coord,poly.coord[1,])
## polygon:
p.name <- paste("p",i.len, sep="")
assign(p.name,(Polygon(poly.coord,hole=F)))
##polygons:
ps.name <- paste("ps",i.len, sep="")
assign(ps.name,(Polygons(list(eval(parse(text=p.name))),ID=paste("s",i.len,sep=""))))
}
sps   <- SpatialPolygons(mget(paste("ps",idx,sep=""),environment()), proj4string=CRS("+proj=longlat +ellps=WGS84"))
#aufräumen, clean-up:
ps.names <- paste("ps",idx,sep="")
rm(list=c(ps.names,"ps.names","ps.name","p.name"))
return(sps)
}


