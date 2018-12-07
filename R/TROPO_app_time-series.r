## TROPOMI-SIF app to extract custom time series and test how filtering affects the data availability and seasonality
##
## Copyright 2018 California Institute of Technology
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## Author: Philipp Koehler, pkoehler@caltech.edu
library(shiny)
library(data.table)
library(ncdf4)
library(solaR)

### directory where ungridded TROPOMI data is stored:
sif.data.directory <- "/net/fluo/data2/projects/TROPOMI/nc_ungridded"
##


# Define UI for app ----
ui <- fluidPage(
  # App title ----
  titlePanel("Extract custom TROPOMI SIF time series from lat/lon input"),
  fluidRow(
    column(3, wellPanel(
      ## Define location
      # Old Black Spruce coordinates serve as default location "where-ever-forest"
      textInput("locationID", label="Location ID (optional --> filename)",value = "where-ever-forest"),
      numericInput("clat", label = "Center latitude (-90° - 90°)",value=53.98),
      numericInput("clon", label = "Center longitude (-180° - 180°)",value=-105.12),
      numericInput("dlat", label = "Delta latitude",value=0.5),
      numericInput("dlon", label = "Delta longitude",value=0.5),
      ## Define date range
      dateRangeInput('dateRange',label = 'Date range (yyyy-mm-dd):',start = as.Date("2018-07-01"), end = as.Date("2018-07-15")),
      #submitButton("Kick it!"),
      actionButton("kick_it","Kick it!"),
      downloadButton("downloadData", "Save Table")
    )),
    mainPanel(
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Request review", 
			    verbatimTextOutput("locText"),
			    verbatimTextOutput("dateRangeText"),
			    verbatimTextOutput("cfRange"),
			    verbatimTextOutput("szaRange"),
			    verbatimTextOutput("vzaRange"),
			    verbatimTextOutput("paRange"),
			    verbatimTextOutput("lstRange"),
			    verbatimTextOutput("dc")),
                  tabPanel("Plot", 
			    plotOutput("plot"),
			    ### Sliders for filtering:
			    h4("Additional Filtering"),
			    column(4,
			    sliderInput("cf_slider", label="Cloud fraction", min = 0,max = 0.8, value = c(0, 0.8)),
			    sliderInput("sza_slider", label="Solar Zenith Angle", min = 0,max = 70, value = c(0, 70))
			    ),
			    column(4, offset=0.1,
			    sliderInput("vza_slider", label="Viewing Zenith Angle", min = 0,max = 60, value = c(0, 60)),
			    sliderInput("pa_slider", label="Phase Angle", min = 0,max = 180, value = c(0, 120))
			    ),
			    column(4, offset=0.5,
			    sliderInput("lst_slider", label="Local Solar Time", min = 0,max = 24, value = c(0., 24.)),
			    checkboxInput("dc_check", "Apply daily correction factor", value = FALSE)
			    )			    
			    ),
                  tabPanel("Table", tableOutput("table"))
      )
    )
  )
)


server <- function(input, output, session) {
  output$locText  <- renderText({
      paste("Creating time series for:", 
	paste(input$clat,"° +/-",input$dlat,"° lat, ",input$clon,"° +/-",input$dlon,"° lon",sep="")
      )
  })
  output$dateRangeText  <- renderText({
      paste("Looking for data between", 
	paste(as.character(input$dateRange), collapse = " and ")
      )
  })
  output$cfRange  <- renderText({
      paste("Accepting cloud fractions from", 
	paste(as.character(input$cf_slider), collapse = " to ")
      )
  })
  output$szaRange  <- renderText({
      paste("Accepting solar zenith angles from ",input$sza_slider[1],"° to ",input$sza_slider[2],"°",sep="")
  })
  output$vzaRange  <- renderText({
      paste("Accepting viewing zenith angles from ",input$vza_slider[1],"° to ",input$vza_slider[2],"°",sep="")
  })
  output$paRange  <- renderText({
      paste("Accepting phase angles from ",input$pa_slider[1],"° to ",input$pa_slider[2],"°",sep="")
  })
  output$lstRange  <- renderText({
      paste("Accepting local solar times from ",input$lst_slider[1],"h to ",input$lst_slider[2],"h",sep="")
  })
 # Checkbox
 updateCheckboxInput(session, "dc_check")
 output$dc <- renderText({
	paste("Apply daily correction factor:",input$dc_check)
  })
  ##########################
  ## Collecting, filtering & plotting the data ----
  get.siteData <- reactive({
    input$kick_it # Re-run when button is clicked
    inp.dates <- isolate(seq(input$dateRange[1],input$dateRange[2],by="1 day"))
    ## all available files:
    inp.files <- grep("_ungridded.nc",list.files(sif.data.directory,pattern="TROPO_SIF_20",recursive=TRUE, full.names=TRUE),value=TRUE)
    ## process only desired time range:
    inp.files <- inp.files[which(as.Date(basename(inp.files),format="TROPO_SIF_%Y-%m-%d") %in% inp.dates)]
    n <- length(inp.files)
    if (n > 0){
      siteData <- data.table()
      withProgress(message = 'Generating plot', value = 0, {
      # looping through time steps
      for (i in 1:n) {
        incProgress(1/n, detail = paste("Looping through file #", i))
	tmp <- isolate(read.TROPO.nc(inp.files[i],box=c(input$clon-input$dlon,input$clon+input$dlon,input$clat-input$dlat,input$clat+input$dlat)))
	if (is.null(tmp) == FALSE){
	    l         <- list(siteData,tmp)
	    siteData  <- rbindlist(l)
	  }
	}
      })
      if (dim(siteData)[1] !=0){ 
      tmpUTC  <- as.POSIXct(siteData$time,origin = '1970-01-01 00:00:00',tz="UTC")
      tmpMST  <- local2Solar(tmpUTC,lon=siteData$clon) 

      siteData$date <- as.character(as.Date(tmpUTC, "%m/%d/%Y"))
      siteData$lst  <- hour(tmpMST)+minute(tmpMST)/60

    return(siteData)
  }else return(NULL)
  }else return(NULL)
  })


  output$plot <- renderPlot({
    input$kick_it # Re-run when button is clicked
    siteData <<- get.siteData()

    if (length(siteData) != 0){ 
    sub.siteData <- subset(siteData, cf >= input$cf_slider[1] & cf<=input$cf_slider[2] & vza>=input$vza_slider[1] & vza<=input$vza_slider[2] &
			   abs(pa)>=input$pa_slider[1] & abs(pa)<=input$pa_slider[2] & sza>=input$sza_slider[1] & sza<=input$sza_slider[2] & 
			   lst>=input$lst_slider[1] & lst<=input$lst_slider[2])
    }else sub.siteData <<- data.table()

    if (dim(sub.siteData)[1] > 0){ 
    
    if (input$dc_check==TRUE) sub.siteData$sif <- sub.siteData$sif*sub.siteData$corr

    ##Averaging:
    day.av <- unlist(lapply(unique(sub.siteData$date),function(day,...) mean(sub.siteData$sif[sub.siteData$date==day])))
    
    #define ylim range:
    ylim <- range(sub.siteData$sif)
    if (ylim[1] < -1) ylim[1] <- -1.
    if (ylim[2] > 4)  ylim[2] <- 4.

    plot(as.Date(sub.siteData$date),sub.siteData$sif,pch=4,lwd=0.3,yaxt="n",xlab="",ylab="",cex=1.,cex.axis=1.5,xlim=isolate(input$dateRange),ylim=ylim,
	 main=paste("n=",dim(sub.siteData)[1],sep=""))#,ylim=range(siteData$sif)xaxt="n",
    abline(h=0,col="grey")
    points(as.Date(unique(sub.siteData$date)),day.av,col="red",pch=17,lwd=2.,cex=1.5)
    legend("top",c("single measurement","daily mean"),fill=c("black","red"),bty="n",ncol=2,cex=1.2)
    if (input$dc_check==TRUE)  ylab <- substitute(SIF[dc]~"@740nm [mW/" ~ m^2 * "/sr/nm]") else ylab <- substitute("SIF@740nm [mW/" ~ m^2 * "/sr/nm]")
    mtext(ylab,2,2.2,cex=1.5)
    #axis(1, cex.axis=1.5)
    axis(2, cex.axis=1.5)
    }else{
    ##empty plot (no data)
    plot(0,0, xaxt = 'n', yaxt = 'n', pch = '', ylab = '', xlab = '')
    text(0,0,"No data...",cex=4)
    }
  })

  ###################################### 
  ## Output in table format to "download"----
  output$table <- renderTable({
    subset(siteData,cf >= input$cf_slider[1] & cf<=input$cf_slider[2] & vza>=input$vza_slider[1] & vza<=input$vza_slider[2] &
		    abs(pa)>=input$pa_slider[1] & abs(pa)<=input$pa_slider[2] & sza>=input$sza_slider[1] & sza<=input$sza_slider[2] & 
		    lst>=input$lst_slider[1] & lst<=input$lst_slider[2],
	  select=c("date","lst","clon","clat","sif","corr","sza","vza","pa","cf","lon1","lon2","lon3","lon4","lat1","lat2","lat3","lat4"))
  })

  ### Download selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$locationID, ".txt", sep = "")
    },
    content = function(file) {
    cat("DATE LST CLON CLAT SIF DC SZA VZA PA CF LON1 LON2 LON3 LON4 LAT1 LAT2 LAT3 LAT4\n", file=file)
    write.table(format(subset(siteData, cf >= input$cf_slider[1] & cf<=input$cf_slider[2] & vza>=input$vza_slider[1] & vza<=input$vza_slider[2] &
			      abs(pa)>=input$pa_slider[1] & abs(pa)<=input$pa_slider[2] & sza>=input$sza_slider[1] & sza<=input$sza_slider[2],
		select=c("date","lst","clon","clat","sif","corr","sza","vza","pa","cf","lon1","lon2","lon3","lon4","lat1","lat2","lat3","lat4")),scientific=FALSE),
		file=file,append=T, row.names = F,col.names = F, quote = F)
    }
  )
}


read.TROPO.nc <- function(inp.file,box=box,cf.min=0,cf.max=0.8,vza.min=0,vza.max=60,pa.min=0,pa.max=180){
  ## read data partially for filtering, read bnds later  
  nc <- try(nc_open(inp.file),silent=TRUE)
  if (class(nc)!="try-error"){
    clon  <- ncvar_get(nc,"lon")
    clat  <- ncvar_get(nc,"lat")
    fs    <- ncvar_get(nc,"sif")
    vza   <- ncvar_get(nc,"vza")
    pa    <- ncvar_get(nc,"phase_angle")
    cf    <- ncvar_get(nc,"cloud_fraction")
  nc_close(nc)

  idx.filter <- which(clon >= box[1] & clon <= box[2] & clat >= box[3] & clat <= box[4] & cf >= cf.min & cf <= cf.max & vza >= vza.min & vza <= vza.max & abs(pa) >= pa.min & abs(pa) <= pa.max)#
  if (length(idx.filter) == 0) return(NULL)
  if (length(idx.filter) == 1){
    nc <- nc_open(inp.file)
      lonb <- ncvar_get(nc,"lon_bnds")[idx.filter,]
      latb <- ncvar_get(nc,"lat_bnds")[idx.filter,]
      time <- ncvar_get(nc,"TIME")[idx.filter] 
      sza  <- ncvar_get(nc,"sza")[idx.filter] 
      corr <- ncvar_get(nc,"daily_correction_factor")[idx.filter]
    nc_close(nc)
    ret <- data.table(time=time,clon=clon[idx.filter],clat=clat[idx.filter],sif=fs[idx.filter],sza=sza,vza=vza[idx.filter],pa=pa[idx.filter],cf=cf[idx.filter],
		      lon1=lonb[1],lon2=lonb[2],lon3=lonb[3],lon4=lonb[4],lat1=latb[1],lat2=latb[2],lat3=latb[3],lat4=latb[4],corr=corr)
    return(ret)
    }
  if (length(idx.filter) >= 2){
    nc <- nc_open(inp.file)
      lonb <- ncvar_get(nc,"lon_bnds")[idx.filter,]
      latb <- ncvar_get(nc,"lat_bnds")[idx.filter,]
      time <- ncvar_get(nc,"TIME")[idx.filter] 
      sza  <- ncvar_get(nc,"sza")[idx.filter] 
      corr <- ncvar_get(nc,"daily_correction_factor")[idx.filter]
    nc_close(nc)
    ret <- data.table(time=time,clon=clon[idx.filter],clat=clat[idx.filter],sif=fs[idx.filter],sza=sza,vza=vza[idx.filter],pa=pa[idx.filter],cf=cf[idx.filter],
		      lon1=lonb[,1],lon2=lonb[,2],lon3=lonb[,3],lon4=lonb[,4],lat1=latb[,1],lat2=latb[,2],lat3=latb[,3],lat4=latb[,4],corr=corr)
    return(ret)
    }
  }else{
  return(NULL)}
}


# Create Shiny app ----
# source("TROPO_app_time-series.r")

app <- shinyApp(ui = ui, server = server)
runApp(app)

