# *------------------------------------------------------------------
# | PROGRAM NAME: 03_ap_model_fit
# | FILE NAME: 03_ap_model_fit.R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  This is a code wrapper to fit the Annual Percentile (AP) model.
# | It fits cumulative probability distributions for annual and monthly flows.
# |
# |
# *------------------------------------------------------------------
# | COMMENTS:               
# |
# |  1:  
# |  2: 
# |  3: 
# |*------------------------------------------------------------------
# | DATA USED:               
# | USGS gauge flow data
# | Annual reconstructions from:
# | Allen, E.B., Rittenour, T.M., DeRose, R.J., Bekker, M.F., Kjelgren, R., Buckley, B.M., 2013. A tree-ring based reconstruction of Logan River streamflow, northern Utah. Water Resources Research 49, 8579–8588. doi:10.1002/2013WR014273.
# |
# | DeRose, R.J., Bekker, M.F., Wang, S.Y., Buckley, B.M., Kjelgren, R.K., Bardsley, T., Rittenour, T.M., Allen, E.B., 2015. A millennium-length reconstruction of Bear River stream flow, Utah. Journal of Hydrology 529, Part 2, 524–534. doi:10.1016/j.jhydrol.2015.01.014.
# |
# |*------------------------------------------------------------------
# | CONTENTS:               
# |
# |  PART 1:  
# |  PART 2: 
# |  PART 3: 
# *-----------------------------------------------------------------
# | UPDATES:               
# |
# |
# *------------------------------------------------------------------

### Clear any existing data or functions.
rm(list=ls())

###########################################################################
## Set the Paths
###########################################################################
### Path for Data and Output	
data_path <- "../../data"
output_path <- "../../output"
global_path <- "../global_func"
function_path <- "./functions"

### Set output location
output_name <- "clim_change"

weber_output_path <- file.path(output_path,"paleo_weber")

write_output_base_path <- file.path(weber_output_path, output_name)

dir.create(write_output_base_path)

###########################################################################
###  Load functions
###########################################################################
### Load these functions for all code
require(colorout)
require(assertthat)

### Load these functions for this unique project
require(ggplot2)
require(tidyverse)
require(staggefuncs)
require(paleoAPR)
require(data.table)

### Load project specific functions
file.sources = list.files(function_path, pattern="*.R", recursive=TRUE)
sapply(file.path(function_path, file.sources),source)

### Load global functions
file.sources = list.files(global_path, pattern="*.R", recursive=TRUE)
sapply(file.path(global_path, file.sources),source)

###########################################################################
## Set Initial Values
###########################################################################
### Set site data
site_id_list <- c("10128500")
site_name_list <- c("Weber River")
recons_file_name_list <- c("weber2014flow.txt")

first_month_wy <- 10 ### Water Year starts on Oct 1
param_cd <- "00060"

monthly_distr <- "gamma"
annual_distr <- "logis"

ref_period <- c(1900,2005)

thresh_level <- 0.5

###########################################################################
###  Set up a loop to run through all site_ides and process the Null Model
###########################################################################
n <- 1


site_id <- site_id_list[n]
site_name <- site_name_list[n]
recons_file_name <- recons_file_name_list[n]

### Create output folders


###########################################################################
###  Read in Data
###########################################################################

### Read in observed flow and fix data type
obs_file_name <- paste0(site_id,"_",param_cd,"_mon_wy.csv")
flow_obs <- read.csv(file.path(weber_output_path,paste0("observed_utah_flow/",obs_file_name)))
flow_obs$date <- as.Date(flow_obs$date)  
#head(flow_obs) # Review data frame

################################################
### Read in monthly and annual parameters to Transform Percentiles
#################################################
percentile_path <- file.path(file.path(weber_output_path, "paleo_reconst"), "norm_fit")

monthly_param <- readRDS(file.path(percentile_path, paste0(site_id,"/",site_id,"_monthly_obs_norm.rds")))


################################################
### Calculate Percentile
#################################################
flow_obs_ts <- data.frame(year=flow_obs$year, water_year=flow_obs$water_year, month=flow_obs$month, flow=flow_obs$monthly_mean)
flow_obs_ts <- paleo.ts(ts=flow_obs_ts, time_scale="monthly")

flow_obs$month_perc <- flow_to_norm(flow_series= flow_obs_ts, dist_object=monthly_param)

###########################################################################
## Generate 50th percentile threshold
###########################################################################
thresh_obs <- flow_obs_ts
thresh_obs$ts$norm <- qnorm(thresh_level)

thresh_obs <- norm_to_flow(norm_series=thresh_obs, dist_object=monthly_param)
flow_obs$threshold <- thresh_obs

###########################################################################
## Save results
###########################################################################
flow_obs$data <- "observed"

### Reorganize
cols_to_retain <- c("site_id", "site_name","data", "year", "month", "date")
flow_obs_df <- subset(flow_obs, select=cols_to_retain)
flow_obs_df$flow_m3s <- flow_obs$monthly_mean
flow_obs_df$thresh_m3s <- flow_obs$threshold
flow_obs_df$norm <- flow_obs$month_perc

### Save to csv
write_location <- file.path(write_output_base_path, paste0(site_id,"_obs_perc_ts.csv"))
write.csv(flow_obs_df, file = write_location, row.names=FALSE)

plot(flow_obs_df$date, flow_obs_df$norm, type="l")

plot(flow_obs_df$thresh_m3s[1:20], type="l")

plot(flow_obs_df$norm[1:50], type="l")

plot(flow_obs_df$flow_m3s[1:50], type="l")
lines(flow_obs_df$thresh_m3s[1:50], col="red")




###########################################################################
###  Read in Weber Climate Change
###########################################################################
weber_cc_path <- file.path(data_path,"weber_climate_change")
weber_cc_flow <- read.csv(file.path(weber_cc_path, "weber_cc_flow.csv"))

### Convert to m3s
for (i in seq(5, 10)){
	weber_cc_flow[,i] <- weber_cc_flow[,i]/35.3146666
}

### Convert dates
weber_cc_flow$date <- as.Date(weber_cc_flow$date, "%m/%d/%Y")

################################################
### Convert to monthly
#################################################
flow_daily <- data.frame(month=weber_cc_flow$mo, year=weber_cc_flow$year, weber_cc_flow[,4:10])
	
### Create a datatable with keys to allow for Monthly and Annual mean calculations
flow_daily <- data.table(flow_daily)
setkey(flow_daily, month, year)

### Calculate mean monthly flow and then re-sort by year and month
flow_monthly <- flow_daily[, lapply(.SD, mean, na.rm=TRUE) , by=c("year", "month")]
### Resort and write over weber_cc_flow
weber_cc_flow <- as.data.frame(flow_monthly[order(date)])


################################################
### Calculate Percentile
#################################################
weber_cc_norm <- weber_cc_flow

for (i in seq(4, 9)){
flow_i_ts <- data.frame(year=weber_cc_norm$year, month=weber_cc_norm$mo, flow=weber_cc_norm[,i])
flow_i_ts <- paleo.ts(ts=flow_i_ts, time_scale="monthly")

weber_cc_norm[,i] <- flow_to_norm(flow_series= flow_i_ts, dist_object=monthly_param)
}

###########################################################################
## Generate 50th percentile threshold
###########################################################################
thresh_cc <- flow_i_ts
thresh_cc$ts$norm <- qnorm(thresh_level)

thresh_cc <- norm_to_flow(norm_series=thresh_cc, dist_object=monthly_param)
weber_cc_flow$threshold <- thresh_cc

###########################################################################
###  Reorganize dataframe
###########################################################################
	for (i in seq(4, 9)){

weber_cc_i_ts <- data.frame(site_id=site_id, site_name=site_name, data=names(weber_cc_norm)[i], year=weber_cc_flow$year, month=weber_cc_flow$mo, date=as.Date(paste0(weber_cc_flow$year,"-", weber_cc_flow$mo,"-15")))

	weber_cc_i_ts$flow_m3s <- weber_cc_flow[,i]
	weber_cc_i_ts$thresh_m3s <- weber_cc_flow$threshold
	weber_cc_i_ts$norm <- weber_cc_norm[,i]

	if(i == 4) {
		weber_cc_ts <- weber_cc_i_ts
	} else {
		weber_cc_ts <- rbind(weber_cc_ts, weber_cc_i_ts)
	}}

###########################################################################
## Save results
###########################################################################
### Save to csv
write_location <- file.path(write_output_base_path, paste0(site_id,"_climchange_perc_ts.csv"))
write.csv(weber_cc_ts, file = write_location, row.names=FALSE)



###########################################################################
###  Read in Weber Paleo
###########################################################################
weber_paleo_path <- file.path(data_path,"weber_climate_change")
weber_paleo_flow <- read.csv(file.path(weber_cc_path, "weber_reconst_postproc_ts.csv"))

###########################################################################
## Generate 50th percentile threshold
###########################################################################
thresh_paleo <- paleo.ts(ts=weber_paleo_flow, time_scale="monthly")
thresh_paleo$ts$norm <- qnorm(thresh_level)

thresh_paleo <- norm_to_flow(norm_series=thresh_paleo, dist_object=monthly_param)
weber_paleo_flow$threshold <- thresh_paleo

###########################################################################
###  Re-organize
###########################################################################
weber_paleo_ts <- data.frame(site_id=site_id, site_name=site_name, data="paleo", year=weber_paleo_flow$year, month=weber_paleo_flow$month, date=as.Date(paste0(weber_paleo_flow$year,"-", weber_paleo_flow$mo,"-15")))

weber_paleo_ts$flow_m3s <- weber_paleo_flow$flow_est
weber_paleo_ts$thresh_m3s <- weber_paleo_flow$threshold
weber_paleo_ts$norm <- weber_paleo_flow$norm_est


###########################################################################
## Save results
###########################################################################
### Save to csv
write_location <- file.path(write_output_base_path, paste0(site_id,"_paleo_perc_ts.csv"))
write.csv(weber_paleo_ts, file = write_location, row.names=FALSE)


