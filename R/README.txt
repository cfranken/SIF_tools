########################################
## How to run TROPO_app_time-series.r ##
########################################

The TROPOMI-SIF extraction tool can be used to subset time series for specific locations and test how filtering/normalizing affects the data availability and seasonality.
The program should be very straight forward to use, you just need R (https://www.r-project.org/) and a few installed packages:
1) "shiny", 
2) "data.table",
3) "ncdf4", 
4) "solaR"
Installing packages in R works by executing following command: install.packages("shiny")). 
You will further need to replace the directory in line 25 with the directory where you stored the ungridded TROPOMI data (available to download at ftp://fluo.gps.caltech.edu/data/tropomi/ungridded/). 
Your browser should open once you execute source("TROPO_app_time-series.r") within R. 