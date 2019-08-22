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
output_name <- "ap_model"

weber_output_path <- file.path(output_path,"paleo_weber")

write_output_base_path <- file.path(file.path(weber_output_path, "paleo_reconst"), output_name)
write_figures_base_path <- file.path(file.path(weber_output_path,"figures"),output_name)

dir.create(write_output_base_path)
dir.create(write_figures_base_path)

###########################################################################
###  Load functions
###########################################################################
### Load these functions for all code
require(colorout)
require(assertthat)

### Load these functions for this unique project
require(ggplot2)
require(monthlypaleo)
require(staggefuncs)

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

###########################################################################
###  Set up a loop to run through all site_ides and process the Null Model
###########################################################################
for (n in seq(1,length(site_id_list))) {

site_id <- site_id_list[n]
site_name <- site_name_list[n]
recons_file_name <- recons_file_name_list[n]

### Create output folders
dir.create(file.path(write_output_base_path, site_id))
dir.create(file.path(write_figures_base_path, site_id))
dir.create(file.path(write_figures_base_path, paste0(site_id, "/pdf")))
dir.create(file.path(write_figures_base_path, paste0(site_id, "/svg")))
dir.create(file.path(write_figures_base_path, paste0(site_id, "/png")))

write_figures_path <- file.path(write_figures_base_path, site_id)
write_output_path <- file.path(write_output_base_path, site_id)

###########################################################################
###  Read in Data
###########################################################################

### Read in observed flow and fix data type
obs_file_name <- paste0(site_id,"_",param_cd,"_mon_wy.csv")
flow_obs <- read.csv(file.path(weber_output_path,paste0("observed_flow/",obs_file_name)))
flow_obs$date <- as.Date(flow_obs$date)  
#head(flow_obs) # Review data frame

### Read in reconst flows (use fread because of large header)
flow_recon <- read_table_wheaders(file.path(data_path,paste0("paleo_flow_annual/",recons_file_name)), sep="\t", na.string="-9999")

flow_recon <- merge(flow_recon, data.frame(age_AD=flow_obs$water_year, flow.obs.m3s=flow_obs$annual_mean), by="age_AD", all.x=TRUE)

################################################
### Calculate Monthly Percentile Fit
#################################################
for (j in 1:12) {
	### Create a test for the month and extract flows from observed for this month
	month_test <- flow_obs$month == j & flow_obs$water_year >= ref_period[1] & flow_obs$water_year <= ref_period[2]
	month_flows <- flow_obs$monthly_mean[month_test]
	
	### Run the distribution fit
	plot_name <- paste0("month_",j)
	month_fit <- perc_fit(flows=month_flows, distr=monthly_distr, plot_name=plot_name, write_folder=write_figures_path)
	
	if (j ==1) {
		fit_param <- month_fit
	} else {
		fit_param <- rbind(fit_param, month_fit)	
	}
}

### Write fit results to CSV file
write_location <- file.path(write_output_path, paste0(site_id,"_param_month_",monthly_distr,".csv"))
write.csv(fit_param, file = write_location,row.names=FALSE)


################################################
### Calculate Annual Percentile Fit
#################################################
### Create test for years in the reference period
year_test <- flow_recon$age_AD >= ref_period[1] & flow_recon$age_AD <= ref_period[2]

### Process observed annual flows
annual_flows <- flow_recon$flow.obs.m3s[year_test]
annual_fit <- perc_fit(flows=annual_flows, distr=annual_distr, plot_name="observ_annual", write_folder=write_figures_path)
	

### Process reconstructed flows
### For other sites do reconstructed analysis
### create column name and row name for analyis and results
col_name <- "flow.rec.local.m3s"
plot_name <- paste0("annual_",gsub("\\.", "_", col_name))
### Create a vector with annual flows from column col_name within reference years
annual_flows <- flow_recon[,c(col_name)][year_test]
### Perform fit on annual flows
annual_fit <- rbind(annual_fit, perc_fit(flows=annual_flows, distr=annual_distr, plot_name=plot_name, write_folder=write_figures_path))



### Write fit results to CSV file
write_location <- file.path(write_output_path, paste0(site_id,"_param_annual_",annual_distr,".csv"))
write.csv(annual_fit, file = write_location,row.names=FALSE)

}


