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
n <- 1


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

################################################
### Read in monthly and annual parameters to Transform Percentiles
#################################################
percentile_path <- file.path(file.path(weber_output_path, "paleo_reconst"), "ap_model")

monthly_param <- list(param=read.csv(file.path(percentile_path, paste0(site_id,"/",site_id,"_param_month_",monthly_distr,".csv"))), distr=monthly_distr)
monthly_param$param <- monthly_param$param[,c(2,3)]

################################################
### Calculate Percentile
#################################################
flow_obs$month_perc <- NA

for (j in 1:12) {
	### Create a test for the month and extract flows from observed for this month
	month_test <- flow_obs$month == j & flow_obs$water_year >= ref_period[1] & flow_obs$water_year <= ref_period[2]
	month_flows <- flow_obs$monthly_mean[month_test]

	### Read parameters
	month_perc <- pgamma(month_flows, shape= monthly_param$param[j,1], rate=monthly_param$param[j,2])
	
	flow_obs$month_perc[month_test] <- month_perc	
}
	



################################################
### Calculate Percentile
#################################################
flow_obs$month_perc <- NA

for (j in 1:12) {
	### Create a test for the month and extract flows from observed for this month
	month_test <- flow_obs$month == j
	month_flows <- flow_obs$monthly_mean[month_test]

	### Read parameters
	month_perc <- pgamma(month_flows, shape= monthly_param$param[j,1], rate=monthly_param$param[j,2])
	
	flow_obs$month_perc[month_test] <- month_perc	
}




perc_threshold <- 0.5
gap_months <- 2
min_duration <- 3

###########################################################################
## Read in Reconstructed Annual Flows
###########################################################################
recon_df <- flow_obs
### Calculate percentile from standard normal
recon_df$percentile <- recon_df$month_perc
recon_df$flow_rec_m3s <- recon_df$monthly_mean

distr_df <- monthly_param$param
###########################################################################
## Calculate time series of threshold
###########################################################################
	### Create column to hold the threshold percentile converted to flow rate	
	recon_df$thresh_flow <- NA
	
	for (j in seq(1,12)) {
		### Calculate threshold flow from percentile and distribution
		thresh_flow <- qgamma(perc_threshold, shape=distr_df$shape[j], rate=distr_df$rate[j])
		### Find correct places to insert and insert
		month_test <- recon_df$month == j
		recon_df$thresh_flow[month_test] <- thresh_flow
	}


	### Test plot
	plot(as.Date(recon_df$date), recon_df$flow_rec_m3s, type="l", xlim=c(as.Date("1958-01-01"), as.Date("1965-01-01")))
	lines(as.Date(recon_df$date), recon_df$thresh_flow, col="red")

###########################################################################
## Extract droughts based on percentile threshold
###########################################################################

	row_index <- seq(1,dim(recon_df)[1])
	
	### Extract droughts
	drought_occur <- recon_df$percentile < perc_threshold
	drought_occur_diff <- c(NA,diff(drought_occur))
	
	drought_begin_index <- row_index[drought_occur_diff == 1 & !is.na(drought_occur_diff)]
	drought_end_index <- row_index[drought_occur_diff == -1 & !is.na(drought_occur_diff)] - 1
	
	### Catch if the first or last nonNA is in drought
	firstNonNAindex <- min(which(!is.na(drought_occur)))
	lastNonNAindex <- max(which(!is.na(drought_occur)))
	
	if (drought_occur[firstNonNAindex] == TRUE) {
		drought_begin_index <- c(firstNonNAindex, drought_begin_index)
	}
	if (drought_occur[lastNonNAindex] == TRUE) {
		drought_end_index <- c(drought_end_index, lastNonNAindex)
	}	
	
	### Fill gaps
	gap_length <- c(drought_begin_index, NA) - c(NA, drought_end_index) - 1
	### Test where gap is less than the gap
	gap_test <- gap_length <= gap_months
	gap_test_end <- c(gap_test[seq(2,length(gap_test))], NA)
	### Extract the end and start in middle of gap
	drought_begin_index <- drought_begin_index[!gap_test | is.na(gap_test)]
	drought_end_index <- drought_end_index[!gap_test_end | is.na(gap_test_end)]
	### Remove NAs
	na_test <- is.na(drought_begin_index) & is.na(drought_end_index)
	drought_begin_index <- drought_begin_index[!na_test]
	drought_end_index <- drought_end_index[!na_test]
	
	
	### Calculate minimum duration droughts
	drought_dur <- drought_end_index - drought_begin_index + 1
	### Test against min duration
	min_dur_test <- drought_dur >= min_duration
	### Redo beginning and end indices
	drought_begin_index <- drought_begin_index[min_dur_test]
	drought_end_index <- drought_end_index[min_dur_test]
	
	
###########################################################################
## Calculate summary statistics
###########################################################################
	drought_event_df <- data.frame(begin=as.Date(recon_df$date[drought_begin_index]), end=as.Date(recon_df$date[drought_end_index]))
	
	### Calculate drought durations
	drought_event_df$dura_months <- drought_end_index - drought_begin_index + 1
	
	### Create empty columns for results
	drought_event_df$begin_month <- NA
	drought_event_df$end_month <- NA
	drought_event_df$min_flow <- NA
	drought_event_df$min_perc <- NA
	drought_event_df$max_deficit <- NA
	drought_event_df$cum_deficit <- NA
	drought_event_df$deficit_center_mass <- NA
	drought_event_df$perc_center_mass <- NA
	drought_event_df$prior_year_flow <- NA
	drought_event_df$dec_rate_to_max_deficit <- NA
	drought_event_df$inc_rate_from_max_deficit <- NA	
	drought_event_df$dec_rate_to_min_perc <- NA
	drought_event_df$inc_rate_from_min_perc <- NA	
	
		
	### Loop through each event and calculate information
	for (j in seq(1, length(drought_begin_index))) {
		event_subset_full <- recon_df[seq(drought_begin_index[j]-1, drought_end_index[j]+1),]
		event_subset <- event_subset_full[seq(2,dim(event_subset_full)[1]-1),]
		
		### Start and End
		drought_event_df$begin_month[j] <- event_subset$month[1]
		drought_event_df$end_month[j] <- event_subset$month[dim(event_subset)[1]]
		
		### Min flow and percentile
		drought_event_df$min_flow[j] <- min(event_subset$flow_rec_m3s, na.rm=TRUE)
		drought_event_df$min_perc[j] <- min(event_subset$percentile, na.rm=TRUE)
		
		### Calculate Deficits
		event_subset$deficit <- event_subset$thresh_flow - event_subset$flow_rec_m3s
		event_subset$deficit[event_subset$deficit<0] <- 0
		
		drought_event_df$max_deficit[j] <- max(event_subset$deficit, na.rm=TRUE)
		drought_event_df$cum_deficit[j] <- sum(event_subset$deficit, na.rm=TRUE)

		### Calculate Centers of Mass
		event_subset$time_fromstart <- seq(1,dim(event_subset)[1])-1
		event_subset$perc_deficit <- perc_threshold - event_subset$percentile
		event_subset$perc_deficit[event_subset$perc_deficit<0] <- 0
		
		drought_event_df$deficit_center_mass[j] <- center_mass(weight=event_subset$deficit, distance=event_subset$time_fromstart)/max(event_subset$time_fromstart)
		drought_event_df$perc_center_mass[j] <- center_mass(weight=event_subset$perc_deficit, distance=event_subset$time_fromstart)/max(event_subset$time_fromstart)

		### Determine rates of change to and from min flow
		event_subset_full$deficit <- event_subset_full$thresh_flow - event_subset_full$flow_rec_m3s
		event_subset_full$time_fromstart <- seq(1,dim(event_subset_full)[1])-1
		event_subset_full$perc_deficit <- perc_threshold - event_subset_full$percentile
	
		min_date <- which.max(event_subset_full$deficit)
		drought_event_df$dec_rate_to_max_deficit[j] <- -(event_subset_full$deficit[min_date] - event_subset_full$deficit[1] )/(min_date-1)

		length_dates <- length(event_subset_full$flow_rec_m3s)
		drought_event_df$inc_rate_from_max_deficit[j] <- -( event_subset_full$deficit[length_dates] - event_subset_full$deficit[min_date])/(length_dates - min_date)


		### Determine rates of change to and from min percentile	
		min_date <- which.max(event_subset_full$perc_deficit)
		drought_event_df$dec_rate_to_min_perc[j] <- -(event_subset_full$perc_deficit[min_date] - event_subset_full$perc_deficit[1] )/(min_date-1)

		length_dates <- length(event_subset_full$flow_rec_m3s)
		drought_event_df$inc_rate_from_min_perc[j] <- -( event_subset_full$perc_deficit[length_dates] - event_subset_full$perc_deficit[min_date])/(length_dates - min_date)

	}
	
	### Prior year
	prior_wy <- recon_df$water_year[drought_begin_index] - 1
	
	### Loop through each event and calculate information
	for (j in seq(1, length(drought_begin_index))) {
		prior_subset <- recon_df[recon_df$water_year == prior_wy[j],]
		drought_event_df$prior_year_flow[j] <- sum(prior_subset$flow_rec_m3s)
	}
	



###########################################################################
## Plot summary statistics
###########################################################################
	
	drought_event_df$century <- factor(paste0(substr(drought_event_df$begin,1,2),"00s"))
	drought_event_df$timeperiod <- "Prior Period"
	drought_event_df$timeperiod[as.numeric(substr(drought_event_df$begin,1,4)) >= 1900] <- "Obs Period"

	drought_event_df$century <- factor(paste0(substr(drought_event_df$begin,1,2),"00s"))
	drought_event_df$timeperiod <- "Prior"
	drought_event_df$timeperiod[as.numeric(substr(drought_event_df$begin,1,4)) >= 1904] <- "Observed"

require(ggrepel)
require(RColorBrewer)
require(ggsci)
require(lubridate)
	
### Plot time series	
	
char_years_list <- c("1932-08", "1999-11", "1976-06", "2009-07",  "1987-05", "1957-10", "2012-05", "1993-11")

#manual_pal <- c(pal_d3(palette="category20")(4), brewer.pal(7, "Purples")[seq(3,7)])
#manual_pal <- manual_pal[c(1, 5, 6, 2, 7, 4, 8, 9, 3)]
#manual_pal <- c(manual_pal, manual_pal[1], manual_pal[4])

cluster_event_df <- drought_event_df
rownames(cluster_event_df) <- substr(cluster_event_df$begin,1,7)

for (n in seq(1,length(char_years_list))){

row_name <- char_years_list[n]

row_n <- rownames(cluster_event_df) %in% row_name

yup <- cluster_event_df[row_n,]

begin_test <- which(recon_df$date %in% yup$begin) -1
end_test <- which(recon_df$date %in% yup$end) + 1

subset_ts <- recon_df[seq(begin_test,end_test),]

base_date <- as.Date(paste0(subset_ts$water_year[1]-1, "-10-01"))

yup <- as.Date(subset_ts$date) - base_date
subset_ts$relative_date <-  yup + as.Date("0000-10-01")
subset_ts$date <- as.Date(subset_ts$date)

first_date <- min(subset_ts$date)
last_date <- first_date + years(14)

breaks_qtr <- seq(as.Date("1400-1-1"), by = "6 months", length.out = 6000)
labels_year = format(seq(from = min(breaks_qtr), to = max(breaks_qtr), by = "1 years"), "%Y")
labs = c(sapply(labels_year, function(x) {
    c(x, rep("", 1))
    }))    
#labs <- labs[seq(1,30)]
    
theme_ts <- theme_classic_correct()+ theme(axis.text.x = element_text(angle = 30, hjust = 1))

p <- ggplot(subset_ts, aes(x=date, y=percentile))
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="grey40")
p <- p + geom_line(colour="blue")#colour=manual_pal[n])
p <- p  + theme_classic()
#p <- p + scale_x_date(name="Date", date_labels = "%Y", breaks = seq(as.Date("1400-1-1"), by = "24 months", length.out = 1000),  date_minor_breaks = "12 months")
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent)
p <- p + coord_cartesian(ylim=c(0,0.62),  xlim=c(first_date, last_date))
p
#
save_folder <- file.path(write_output_path, "heir_agglom/ts")
dir.create(save_folder, showWarnings=FALSE)

### Save results
ggsave(paste0(row_name, "_obs_ts.png"),  p, width=6, height=3, dpi=300)
ggsave(paste0(row_name, "_obs_ts.pdf"),  p, width=6, height=3)
#ggsave(paste0(row_name, "_obs_ts.svg"),  p, width=6, height=3)

}



	plot(as.Date(recon_df$date), recon_df$flow_rec_m3s, type="l", xlim=c(as.Date("1958-01-01"), as.Date("1965-01-01")))
	lines(as.Date(recon_df$date), recon_df$thresh_flow, col="red")

plot_df <- data.frame(date=recon_df$date, flow=recon_df$monthly_mean, label="Observed")
thresh_df <- data.frame(date=recon_df$date, flow=recon_df$thresh_flow, label="50th Percentile")
yup <- rbind(plot_df, thresh_df)

p <- ggplot(yup, aes(x=date, y=flow*35.3146667, colour=label, linetype=label))
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow (cfs)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,2000),  xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("2012_drought_flow_blank.png",  p, width=12, height=4, dpi=600)


p <- ggplot(yup, aes(x=date, y=flow*35.3146667, colour=label, linetype=label))
p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow (cfs)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,2000),  xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("2012_drought_flow_2013.png",  p, width=12, height=4, dpi=600)


p <- ggplot(yup, aes(x=date, y=flow*35.3146667, colour=label, linetype=label))
p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow (cfs)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,2000),  xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("2012_drought_flow_2014.png",  p, width=12, height=4, dpi=600)


p <- ggplot(yup, aes(x=date, y=flow*35.3146667, colour=label, linetype=label))
p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2015-04-01"), xmax=as.Date("2015-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2015-10-01"), size=0.6, colour="grey50")
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow (cfs)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,2000),  xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("2012_drought_flow_2015.png",  p, width=12, height=4, dpi=600)



p <- ggplot(yup, aes(x=date, y=flow*35.3146667, colour=label, linetype=label))
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow (cfs)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,2000),  xlim=c(as.Date("1975-01-01"), as.Date("1981-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("1976_drought_flow_blank.png",  p, width=12, height=4, dpi=600)




p <- ggplot(yup, aes(x=date, y=flow*35.3146667, colour=label, linetype=label))
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow (cfs)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,2000),  xlim=c(as.Date("1930-01-01"), as.Date("1936-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("1932_drought_flow_blank.png",  p, width=12, height=4, dpi=600)








p <- ggplot(recon_df, aes(x=date, y=percentile))
p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=-Inf, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2015-04-01"), xmax=as.Date("2015-10-01"), ymin=-Inf, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2015-10-01"), size=0.6, colour="grey50")
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="#d95f02")
p <- p + geom_line(colour="black")#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p
#

ggsave("2012_drought_percentile_2015.png",  p, width=12, height=4, dpi=600)



p <- ggplot(recon_df, aes(x=date, y=percentile))
p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=-Inf, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="#d95f02")
p <- p + geom_line(colour="black")#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p
#

ggsave("2012_drought_percentile_2014.png",  p, width=12, height=4, dpi=600)


p <- ggplot(recon_df, aes(x=date, y=percentile))
p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=-Inf, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="#d95f02")
p <- p + geom_line(colour="black")#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p
#

ggsave("2012_drought_percentile_2013.png",  p, width=12, height=4, dpi=600)


p <- ggplot(recon_df, aes(x=date, y=percentile))
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="#d95f02")
p <- p + geom_line(colour="black")#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p
#

ggsave("2012_drought_percentile_blank.png",  p, width=12, height=4, dpi=600)




p <- ggplot(recon_df, aes(x=date, y=percentile))
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="#d95f02")
p <- p + geom_line(colour="black")#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1), xlim=c(as.Date("1975-01-01"), as.Date("1981-01-01")))
p
#

ggsave("1976_drought_percentile_blank.png",  p, width=12, height=4, dpi=600)





p <- ggplot(recon_df, aes(x=date, y=percentile))
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="#d95f02")
p <- p + geom_line(colour="black")#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("1930-01-01"), as.Date("1936-01-01")))
p
#

ggsave("1932_drought_percentile_blank.png",  p, width=12, height=4, dpi=600)






### Read in Weber Storage

weber_stor <- read.csv("/run/media/jhstagge/Data/Documents/work/data/weber_storage/weber_storage.csv")
### Fix date
weber_stor$Date <- as.Date(strptime(weber_stor$Date, "%m/%d/%Y"))

plot(weber_stor$Date, weber_stor$TotStor, type="l")

### Melt to include all reservoirs
weber_melt <- melt(weber_stor, id.vars="Date")

### Add to make upper and lower
weber_summary <- weber_stor[,c(1,2)]
weber_summary$site <- "Total"
colnames(weber_summary) <- c("Date", "Storage", "Site")
upper <- weber_stor$SMH + weber_stor$Wanship + weber_stor$Echo + weber_stor$LstCrk + weber_stor$ECanyon
weber_upper <- weber_summary
weber_lower <- weber_summary

weber_upper$Storage <- upper
weber_upper$Site <- "Upper"
weber_summary <- rbind(weber_summary, weber_upper)
lower <- weber_stor$Causey + weber_stor$Pineview + weber_stor$Willard
weber_lower$Storage <- lower
weber_lower$Site <- "Lower"
weber_summary <- rbind(weber_summary, weber_lower)

p <- ggplot(weber_stor, aes(x=Date, y=TotStorage/1000))
p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2015-04-01"), xmax=as.Date("2015-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2015-10-01"), size=0.6, colour="grey50")
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Storage (1,000 Acre-Feet)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,550), xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("2012_drought_storage_2015.png",  p, width=12, height=4, dpi=600)


p <- ggplot(weber_stor, aes(x=Date, y=TotStorage/1000))
p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Storage (1,000 Acre-Feet)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,550), xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("2012_drought_storage_2014.png",  p, width=12, height=4, dpi=600)



p <- ggplot(weber_stor, aes(x=Date, y=TotStorage/1000))
p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Storage (1,000 Acre-Feet)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,550), xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("2012_drought_storage_2013.png",  p, width=12, height=4, dpi=600)



p <- ggplot(weber_stor, aes(x=Date, y=TotStorage/1000))
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Storage (1,000 Acre-Feet)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,550), xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("2012_drought_storage_blank.png",  p, width=12, height=4, dpi=600)


p <- ggplot(weber_stor, aes(x=Date, y=TotStorage/1000))
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Storage (1,000 Acre-Feet)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,550), xlim=c(as.Date("1975-01-01"), as.Date("1981-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("1976_drought_storage_blank.png",  p, width=12, height=4, dpi=600)


p <- ggplot(weber_stor, aes(x=Date, y=TotStorage/1000))
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#d95f02"), name="")
p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Storage (1,000 Acre-Feet)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,550), xlim=c(as.Date("1931-01-01"), as.Date("1937-01-01")))
p <- p + theme(legend.position = c(0.90, 0.9))
p
#

ggsave("1932_drought_storage_blank.png",  p, width=12, height=4, dpi=600)








p <- ggplot(weber_melt, aes(x=Date, y=value/1000, colour=variable))
p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2015-04-01"), xmax=as.Date("2015-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2015-10-01"), size=0.6, colour="grey50")
p <- p + geom_line()
p <- p + scale_colour_manual(values=cb_pal(9, pal="d3"), name="")
#p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Storage (1,000 Acre-Feet)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,550), xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p <- p + theme(legend.position = c(0.90, 0.8))
p




p <- ggplot(weber_summary, aes(x=Date, y=Storage/1000, colour=Site))
p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_rect(aes(xmin=as.Date("2015-04-01"), xmax=as.Date("2015-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
p <- p + geom_vline(xintercept=as.Date("2015-10-01"), size=0.6, colour="grey50")
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("black", "#1F77B4", "#FF7F0E"), name="")
#p <- p + scale_linetype_manual(values=c("solid", "longdash"), name="")
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Storage (1,000 Acre-Feet)", expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,550), xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p <- p + theme(legend.position = c(0.90, 0.8))
p


#ggsave("2012_drought_storage_2015.png",  p, width=12, height=4, dpi=600)



weber_recon_df <- read.csv("/run/media/jhstagge/Data/Documents/work/output/paleo_weber/apr_model/10128500_apr_model_enso_pca_impute_std_concur_reconst_ts_rec_local.csv")
weber_recon_df$percentile <- pnorm(weber_recon_df$monthly_norm)
recon_orig_date <- as.Date(weber_recon_df$date) 

weber_recon_df$date <- recon_orig_date + years(2012-1706)
weber_recon_df_plot <- subset(weber_recon_df, date < as.Date("2030-01-01"))
#recon_adjust_df <- recon_df 
#recon_adjust_df$date <- recon_adjust_df$date - years(306)

p <- ggplot(weber_recon_df_plot, aes(x=date, y=percentile))
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="black")
p <- p + geom_line(data=recon_df, colour="grey20", alpha=0.5)
p <- p + geom_line(colour="#1F77B4FF", size=0.7)#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2012-01-01"), as.Date("2025-01-01")))
p
#

ggsave("1706_drought_vs_2012.png",  p, width=8, height=4, dpi=600)


weber_recon_df$date <- recon_orig_date + years(2012-1930)
weber_recon_df_plot <- subset(weber_recon_df, date < as.Date("2020-01-01"))
#recon_adjust_df <- recon_df 
#recon_adjust_df$date <- recon_adjust_df$date - years(306)

p <- ggplot(weber_recon_df_plot, aes(x=date, y=percentile))
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="black")
p <- p + geom_line(data=recon_df, colour="grey20", alpha=0.5)
p <- p + geom_line(colour="#FF7F0EFF", size=0.7)#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2012-01-01"), as.Date("2025-01-01")))
p
#

ggsave("1930_drought_vs_2012.png",  p, width=8, height=4, dpi=600)



weber_recon_df$date <- recon_orig_date + years(2012-1590)
weber_recon_df_plot <- subset(weber_recon_df, date < as.Date("2017-01-01"))
#recon_adjust_df <- recon_df 
#recon_adjust_df$date <- recon_adjust_df$date - years(306)

p <- ggplot(weber_recon_df_plot, aes(x=date, y=percentile))

p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="black")
p <- p + geom_line(data=recon_df, colour="grey20", alpha=0.5)
p <- p + geom_line(colour="#1F77B4FF", size=0.7)#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2012-01-01"), as.Date("2025-01-01")))
p
#

ggsave("1590_drought_vs_2012.png",  p, width=8, height=4, dpi=600)


weber_recon_df$date <- recon_orig_date + years(2012-1886)
weber_recon_df_plot <- subset(weber_recon_df, date < as.Date("2017-01-01"))
#recon_adjust_df <- recon_df 
#recon_adjust_df$date <- recon_adjust_df$date - years(306)

p <- ggplot(weber_recon_df_plot, aes(x=date, y=percentile))
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="black")
p <- p + geom_line(data=recon_df, colour="grey20", alpha=0.5)
p <- p + geom_line(colour="#FF7F0EFF", size=0.7)#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2012-01-01"), as.Date("2025-01-01")))
p
#

ggsave("1886_drought_vs_2012.png",  p, width=8, height=4, dpi=600)




weber_recon_df$date <- recon_orig_date + years(2012-1750)
weber_recon_df_plot <- subset(weber_recon_df, date < as.Date("2016-01-01"))
#recon_adjust_df <- recon_df 
#recon_adjust_df$date <- recon_adjust_df$date - years(306)

p <- ggplot(weber_recon_df_plot, aes(x=date, y=percentile))
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="black")
p <- p + geom_line(data=recon_df, colour="grey20", alpha=0.5)
p <- p + geom_line(colour="#2CA02CFF", size=0.7)#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2012-01-01"), as.Date("2025-01-01")))
p
#

ggsave("1750_drought_vs_2012.png",  p, width=8, height=4, dpi=600)




weber_recon_df$date <- recon_orig_date + years(2012-1861)
weber_recon_df_plot <- subset(weber_recon_df, date < as.Date("2014-01-01"))
#recon_adjust_df <- recon_df 
#recon_adjust_df$date <- recon_adjust_df$date - years(306)

p <- ggplot(weber_recon_df_plot, aes(x=date, y=percentile))
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="black")
p <- p + geom_line(data=recon_df, colour="grey20", alpha=0.5)
p <- p + geom_line(colour="#D62728FF", size=0.7)#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2012-01-01"), as.Date("2025-01-01")))
p
#

ggsave("1861_drought_vs_2012.png",  p, width=8, height=4, dpi=600)






	### Some example plots
	p <- ggplot(drought_event_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	p <- p + geom_point(aes(colour=timeperiod))
	#p <- p + geom_text(hjust = 0, nudge_x = 0.5, nudge_y = 0.001, check_overlap = TRUE)
	#p <- p + geom_text_repel( size=2.7, box.padding = unit(0.5, "lines"),
    #point.padding = unit(0.2, "lines"),
    #segment.color = 'grey', segment.alpha=0.4,nudge_x = 0.5,  nudge_y=0.001)
	p <- p + geom_text_repel( size=2.7,segment.color = 'grey', segment.alpha=0.4)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	#p <- p + scale_color_manual(name="Century", values=pal_cb_5_more)
	#p <- p + scale_color_manual(name="Century", values=rev(plasma(7)))
	p <- p + scale_color_manual(name="Time Period", values=c("red", "black"))
	#p <- p + scale_color_viridis(option="magma")
#	p <- p + theme(legend.position = c(0.8, 0.2))
	#p <- p + theme_classic_correct()
	p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))
	#p <- p + theme(axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))
	#	p <- p + theme(axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))

	p
ggsave(file.path(write_output_path, paste0(save_name,"_drought_perc_vs_duration.png")),  p, width=6.5, height=5, dpi=300)



paleo_drought_event_df <- read.csv("/run/media/jhstagge/Data/Documents/work/output/paleo_clustering/drought_details/weber_river_droughts.csv")
drought_event_df$datasource <- "Observed"
paleo_drought_event_df$datasource <- NA
paleo_drought_event_df$datasource[paleo_drought_event_df$timeperiod=="Prior"] <- "Reconstructed"
paleo_drought_event_df$datasource[paleo_drought_event_df$timeperiod=="Observed"] <- "Reconstructed (1900-2000)"
paleo_drought_event_df$datasource <- factor(paleo_drought_event_df$datasource, levels=c("Reconstructed (1900-2000)", "Reconstructed"))
combined_drought_event_df <- rbind(drought_event_df, paleo_drought_event_df)

	### Some example plots
	p <- ggplot(paleo_drought_event_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
	p <- p + geom_text(data=drought_event_df, aes(colour=datasource), size=2.7)
	p <- p + geom_text(aes(colour=datasource), size=2.7)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	#p <- p + scale_color_manual(name="Century", values=pal_cb_5_more)
	#p <- p + scale_color_manual(name="Century", values=rev(plasma(7)))
	p <- p + scale_color_manual(name="Data Source", values=c("black","cadetblue",  "red"))
	#p <- p + scale_color_viridis(option="magma")
#	p <- p + theme(legend.position = c(0.8, 0.2))
	#p <- p + theme_classic_correct()
	p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))
	#p <- p + theme(axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))
	#	p <- p + theme(axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))

	p
	
ggsave(paste0("observed_recon","_drought_perc_vs_duration.png"),  p, width=6.5, height=5, dpi=600)
ggsave(paste0("observed_recon","_drought_perc_vs_duration.pdf"),  p, width=6.5, height=5)
	
	
	
	
	### Some example plots
	p <- ggplot(subset(paleo_drought_event_df, timeperiod=="Observed"), aes(colour=datasource, x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
	#p <- p + geom_text(data=drought_event_df, aes(colour=datasource), size=2.7)
	p <- p + geom_text( size=2.7)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2), expand=c(0,0))
	p <- p + scale_y_continuous(name="Min Flow Percentile", expand=c(0,0))
	#p <- p + scale_color_manual(name="Century", values=pal_cb_5_more)
	#p <- p + scale_color_manual(name="Century", values=rev(plasma(7)))
	p <- p + scale_color_manual(name="Data Source", values=c("red","black"))
	#p <- p + scale_color_viridis(option="magma")
#	p <- p + theme(legend.position = c(0.8, 0.2))
	#p <- p + theme_classic_correct()
	p <- p + theme_classic_new(14)
	p <- p + coord_cartesian(ylim=c(0,50), xlim=c(0,15))
	p <- p + theme(legend.position = c(0.85, 0.85))
	#p <- p + theme(axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))
	#	p <- p + theme(axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))

	p
	
ggsave(paste0("recon","_drought_perc_vs_duration_only_recent.png"),  p, width=6.5, height=5, dpi=600)
#ggsave(paste0("recon","_drought_perc_vs_duration_only_recent.pdf"),  p, width=6.5, height=5)


	### Some example plots
	p <- ggplot(paleo_drought_event_df, aes(colour=datasource, x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
	#p <- p + geom_text(data=drought_event_df, aes(colour=datasource), size=2.7)
	p <- p + geom_text( size=2.7)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2), expand=c(0,0))
	p <- p + scale_y_continuous(name="Min Flow Percentile", expand=c(0,0))
	#p <- p + scale_color_manual(name="Century", values=pal_cb_5_more)
	#p <- p + scale_color_manual(name="Century", values=rev(plasma(7)))
	p <- p + scale_color_manual(name="Data Source", values=c("red","black"))
	#p <- p + scale_color_viridis(option="magma")
#	p <- p + theme(legend.position = c(0.8, 0.2))
	#p <- p + theme_classic_correct()
	p <- p + theme_classic_new(14)
	p <- p + coord_cartesian(ylim=c(0,50), xlim=c(0,15))
	p <- p + theme(legend.position = c(0.85, 0.85))
	#p <- p + theme(axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))
	#	p <- p + theme(axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))

	p
	
ggsave(paste0("recon","_drought_perc_vs_duration_all.png"),  p, width=6.5, height=5, dpi=600)
	
	
	
	
	
	
	
	
	
	
	
	


###########################################################################
###  Read in Weber Climate Change
###########################################################################
weber_cc_path <- file.path(data_path,"weber_climate_change")
weber_cc_flow <- read.csv(file.path(weber_cc_path, "weber_cc_flow.csv"))

weber_cc_perc <- weber_cc_flow

###########################################################################
###  
###########################################################################
require(data.table)

for(i in seq(5,10)){
	weber_cc_subset <- weber_cc_perc
	
	### Convert to m3/s
	weber_cc_perc[,i] <- weber_cc_perc[,i] * ft3_to_m3 
	
	### Create a datatable with keys to allow for Monthly and Annual mean calculations
	weber_cc_subset <- data.table(weber_cc_perc[,c(1,2,3,4,i)])
	col_name <- names(weber_cc_subset)[5]
	names(weber_cc_subset)[5] <- "flow"
	setkey(weber_cc_subset, mo, year)

	### Calculate mean monthly flow and then re-sort by year and month
	### Assigns a date as the mean day of each month (only to be used for plotting)
	weber_cc_subset_monthly <- weber_cc_subset[,list(date=mean(dy), monthly_sum=sum(flow), monthly_mean=mean(flow)), by=c("year", "mo")]
	weber_cc_subset_monthly <- weber_cc_subset_monthly[order(rank(year), mo)]

	names(weber_cc_subset_monthly)[5] <- col_name
		
	### Write results
	if (i ==5){
		weber_cc_monthly <- weber_cc_subset_monthly[,c(1,2,3,5)]
	} else {
		weber_cc_monthly <- cbind(weber_cc_monthly, weber_cc_subset_monthly[,5])
	}

}

weber_cc_melt <- melt(weber_cc_monthly, id.vars=c("year", "mo", "date"))
weber_cc_melt$month_perc <- NA

for (j in 1:12) {
	### Create a test for the month and extract flows from observed for this month
	month_test <- weber_cc_melt$mo == j 
	month_flows <- weber_cc_melt$value[month_test]

	### Read parameters
	month_perc <- pgamma(month_flows, shape= monthly_param$param[j,1], rate=monthly_param$param[j,2])
	
	weber_cc_melt$month_perc[month_test] <- month_perc	
}


weber_cc_melt$date <- as.Date(paste0(weber_cc_melt$year,"-", weber_cc_melt$mo,"-15"))
weber_cc_melt$norm <- qnorm(weber_cc_melt$month_perc)


p <- ggplot(weber_cc_melt, aes(x=date, y=month_perc, colour=variable))
p <- p + geom_line()
p <- p + theme_classic(9)
p



p <- ggplot(weber_cc_melt, aes(x=date, y=norm, colour=variable))
p <- p + geom_line()
p <- p + theme_classic(9)
p










perc_threshold <- 0.5
gap_months <- 2
min_duration <- 3


cc_levels <- levels(weber_cc_melt$variable)

###########################################################################
## Read in Reconstructed Annual Flows
###########################################################################
for (k in seq(1, length(cc_levels))){

scenario <- cc_levels[[k]]
recon_df <- weber_cc_melt[weber_cc_melt$variable == scenario,]

### Calculate percentile from standard normal
recon_df$percentile <- recon_df$month_perc
recon_df$flow_rec_m3s <- recon_df$value
recon_df$month <- recon_df$mo 

distr_df <- monthly_param$param
###########################################################################
## Calculate time series of threshold
###########################################################################
	### Create column to hold the threshold percentile converted to flow rate	
	recon_df$thresh_flow <- NA
	
	for (j in seq(1,12)) {
		### Calculate threshold flow from percentile and distribution
		thresh_flow <- qgamma(perc_threshold, shape=distr_df$shape[j], rate=distr_df$rate[j])
		### Find correct places to insert and insert
		month_test <- recon_df$month == j
		recon_df$thresh_flow[month_test] <- thresh_flow
	}


	### Test plot
	plot(as.Date(recon_df$date), recon_df$flow_rec_m3s, type="l", xlim=c(as.Date("1980-01-01"), as.Date("2000-01-01")))
	lines(as.Date(recon_df$date), recon_df$thresh_flow, col="red")

###########################################################################
## Extract droughts based on percentile threshold
###########################################################################

	row_index <- seq(1,dim(recon_df)[1])
	
	### Extract droughts
	drought_occur <- recon_df$percentile < perc_threshold
	drought_occur_diff <- c(NA,diff(drought_occur))
	
	drought_begin_index <- row_index[drought_occur_diff == 1 & !is.na(drought_occur_diff)]
	drought_end_index <- row_index[drought_occur_diff == -1 & !is.na(drought_occur_diff)] - 1
	
	### Catch if the first or last nonNA is in drought
	firstNonNAindex <- min(which(!is.na(drought_occur)))
	lastNonNAindex <- max(which(!is.na(drought_occur)))
	
	if (drought_occur[firstNonNAindex] == TRUE) {
		drought_begin_index <- c(firstNonNAindex, drought_begin_index)
	}
	if (drought_occur[lastNonNAindex] == TRUE) {
		drought_end_index <- c(drought_end_index, lastNonNAindex)
	}	
	
	### Fill gaps
	gap_length <- c(drought_begin_index, NA) - c(NA, drought_end_index) - 1
	### Test where gap is less than the gap
	gap_test <- gap_length <= gap_months
	gap_test_end <- c(gap_test[seq(2,length(gap_test))], NA)
	### Extract the end and start in middle of gap
	drought_begin_index <- drought_begin_index[!gap_test | is.na(gap_test)]
	drought_end_index <- drought_end_index[!gap_test_end | is.na(gap_test_end)]
	### Remove NAs
	na_test <- is.na(drought_begin_index) & is.na(drought_end_index)
	drought_begin_index <- drought_begin_index[!na_test]
	drought_end_index <- drought_end_index[!na_test]
	
	
	### Calculate minimum duration droughts
	drought_dur <- drought_end_index - drought_begin_index + 1
	### Test against min duration
	min_dur_test <- drought_dur >= min_duration
	### Redo beginning and end indices
	drought_begin_index <- drought_begin_index[min_dur_test]
	drought_end_index <- drought_end_index[min_dur_test]
	
	
###########################################################################
## Calculate summary statistics
###########################################################################
	drought_event_df <- data.frame(begin=as.Date(recon_df$date[drought_begin_index]), end=as.Date(recon_df$date[drought_end_index]))
	
	### Calculate drought durations
	drought_event_df$dura_months <- drought_end_index - drought_begin_index + 1
	
	### Create empty columns for results
	drought_event_df$begin_month <- NA
	drought_event_df$end_month <- NA
	drought_event_df$min_flow <- NA
	drought_event_df$min_perc <- NA
	drought_event_df$max_deficit <- NA
	drought_event_df$cum_deficit <- NA
	drought_event_df$deficit_center_mass <- NA
	drought_event_df$perc_center_mass <- NA
	drought_event_df$prior_year_flow <- NA
	drought_event_df$dec_rate_to_max_deficit <- NA
	drought_event_df$inc_rate_from_max_deficit <- NA	
	drought_event_df$dec_rate_to_min_perc <- NA
	drought_event_df$inc_rate_from_min_perc <- NA	
	
		
	### Loop through each event and calculate information
	for (j in seq(1, length(drought_begin_index))) {
		event_subset_full <- recon_df[seq(drought_begin_index[j]-1, drought_end_index[j]+1),]
		event_subset <- event_subset_full[seq(2,dim(event_subset_full)[1]-1),]
		
		### Start and End
		drought_event_df$begin_month[j] <- event_subset$month[1]
		drought_event_df$end_month[j] <- event_subset$month[dim(event_subset)[1]]
		
		### Min flow and percentile
		drought_event_df$min_flow[j] <- min(event_subset$flow_rec_m3s, na.rm=TRUE)
		drought_event_df$min_perc[j] <- min(event_subset$percentile, na.rm=TRUE)
		
		### Calculate Deficits
		event_subset$deficit <- event_subset$thresh_flow - event_subset$flow_rec_m3s
		event_subset$deficit[event_subset$deficit<0] <- 0
		
		drought_event_df$max_deficit[j] <- max(event_subset$deficit, na.rm=TRUE)
		drought_event_df$cum_deficit[j] <- sum(event_subset$deficit, na.rm=TRUE)

		### Calculate Centers of Mass
		event_subset$time_fromstart <- seq(1,dim(event_subset)[1])-1
		event_subset$perc_deficit <- perc_threshold - event_subset$percentile
		event_subset$perc_deficit[event_subset$perc_deficit<0] <- 0
		
		drought_event_df$deficit_center_mass[j] <- center_mass(weight=event_subset$deficit, distance=event_subset$time_fromstart)/max(event_subset$time_fromstart)
		drought_event_df$perc_center_mass[j] <- center_mass(weight=event_subset$perc_deficit, distance=event_subset$time_fromstart)/max(event_subset$time_fromstart)

		### Determine rates of change to and from min flow
		event_subset_full$deficit <- event_subset_full$thresh_flow - event_subset_full$flow_rec_m3s
		event_subset_full$time_fromstart <- seq(1,dim(event_subset_full)[1])-1
		event_subset_full$perc_deficit <- perc_threshold - event_subset_full$percentile
	
		min_date <- which.max(event_subset_full$deficit)
		drought_event_df$dec_rate_to_max_deficit[j] <- -(event_subset_full$deficit[min_date] - event_subset_full$deficit[1] )/(min_date-1)

		length_dates <- length(event_subset_full$flow_rec_m3s)
		drought_event_df$inc_rate_from_max_deficit[j] <- -( event_subset_full$deficit[length_dates] - event_subset_full$deficit[min_date])/(length_dates - min_date)


		### Determine rates of change to and from min percentile	
		min_date <- which.max(event_subset_full$perc_deficit)
		drought_event_df$dec_rate_to_min_perc[j] <- -(event_subset_full$perc_deficit[min_date] - event_subset_full$perc_deficit[1] )/(min_date-1)

		length_dates <- length(event_subset_full$flow_rec_m3s)
		drought_event_df$inc_rate_from_min_perc[j] <- -( event_subset_full$perc_deficit[length_dates] - event_subset_full$perc_deficit[min_date])/(length_dates - min_date)

	}
	
	### Prior year
	prior_wy <- recon_df$water_year[drought_begin_index] - 1
	
	### Loop through each event and calculate information
	for (j in seq(1, length(drought_begin_index))) {
		prior_subset <- recon_df[recon_df$water_year == prior_wy[j],]
		drought_event_df$prior_year_flow[j] <- sum(prior_subset$flow_rec_m3s)
	}
	
	drought_event_df$scenario <- scenario
	
	if (k ==1) {
		cc_drought_event_df <- drought_event_df
	} else {
		cc_drought_event_df <- rbind(cc_drought_event_df, drought_event_df)
	}
}
	


###########################################################################
## Plot summary statistics
###########################################################################



	### Some example plots
	p <- ggplot(paleo_drought_event_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
	p <- p + geom_text(colour="grey50", alpha=0.5, size=2.7)
	p <- p + geom_text(data=cc_drought_event_df, aes(colour=scenario), size=2.7)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	#p <- p + scale_colour_brewer(name="Scenario", palette="Set1")
	p <- p + scale_color_manual(name="Scenario", values=pal_d3(palette = "category20")(6))
	#p <- p + scale_color_manual(name="Century", values=rev(plasma(7)))
	#p <- p + scale_color_manual(name="Data Source", values=c("black","cadetblue",  "red"))
	#p <- p + scale_color_viridis(option="magma")
#	p <- p + theme(legend.position = c(0.8, 0.2))
	#p <- p + theme_classic_correct()
	p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))
	#p <- p + theme(axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))
	#	p <- p + theme(axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))

	p
	
		
ggsave(paste0("cc_recon","_drought_perc_vs_duration.png"),  p, width=6.5, height=5, dpi=600)
ggsave(paste0("cc_recon","_drought_perc_vs_duration.pdf"),  p, width=6.5, height=5)
	
	p <- p + coord_cartesian(ylim=c(0,16))
	
		
ggsave(paste0("cc_recon","_drought_perc_vs_duration_zoom.png"),  p, width=6.5, height=5, dpi=600)
ggsave(paste0("cc_recon","_drought_perc_vs_duration_zoom.pdf"),  p, width=6.5, height=5)
	
	
cc_drought_event_df$begin_yearonly <- substr(cc_drought_event_df$begin,1,4)
	
		### Some example plots
	p <- ggplot(paleo_drought_event_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
	p <- p + geom_text(colour="grey50", alpha=0.5, size=2.7)
	p <- p + geom_text(data=subset(cc_drought_event_df, begin_yearonly =="1987"), aes(colour=scenario), size=2.7)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	#p <- p + scale_colour_brewer(name="Scenario", palette="Set1")
	p <- p + scale_color_manual(name="Scenario", values=pal_d3(palette = "category20")(6))
	#p <- p + scale_color_manual(name="Century", values=rev(plasma(7)))
	#p <- p + scale_color_manual(name="Data Source", values=c("black","cadetblue",  "red"))
	#p <- p + scale_color_viridis(option="magma")
#	p <- p + theme(legend.position = c(0.8, 0.2))
	#p <- p + theme_classic_correct()
	p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))
	#p <- p + theme(axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))
	#	p <- p + theme(axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))

	p
	
		
ggsave(paste0("cc_recon","_drought_perc_vs_duration_1987.png"),  p, width=6.5, height=5, dpi=600)
ggsave(paste0("cc_recon","_drought_perc_vs_duration_1987.pdf"),  p, width=6.5, height=5)



	
		### Some example plots
	p <- ggplot(paleo_drought_event_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
	p <- p + geom_text(colour="grey50", alpha=0.5, size=2.7)
	p <- p + geom_text(data=subset(cc_drought_event_df, begin_yearonly >= 1999 & begin_yearonly <= 2001), aes(colour=scenario), size=2.7)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	#p <- p + scale_colour_brewer(name="Scenario", palette="Set1")
	p <- p + scale_color_manual(name="Scenario", values=pal_d3(palette = "category20")(6))
	#p <- p + scale_color_manual(name="Century", values=rev(plasma(7)))
	#p <- p + scale_color_manual(name="Data Source", values=c("black","cadetblue",  "red"))
	#p <- p + scale_color_viridis(option="magma")
#	p <- p + theme(legend.position = c(0.8, 0.2))
	#p <- p + theme_classic_correct()
	p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))
	#p <- p + theme(axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))
	#	p <- p + theme(axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))

	p
	

	
ggsave(paste0("cc_recon","_drought_perc_vs_duration_2000.png"),  p, width=6.5, height=5, dpi=600)
ggsave(paste0("cc_recon","_drought_perc_vs_duration_2000.pdf"),  p, width=6.5, height=5)





	drought_event_df$century <- factor(paste0(substr(drought_event_df$begin,1,2),"00s"))
	drought_event_df$timeperiod <- "Prior Period"
	drought_event_df$timeperiod[as.numeric(substr(drought_event_df$begin,1,4)) >= 1900] <- "Obs Period"

	drought_event_df$century <- factor(paste0(substr(drought_event_df$begin,1,2),"00s"))
	drought_event_df$timeperiod <- "Prior"
	drought_event_df$timeperiod[as.numeric(substr(drought_event_df$begin,1,4)) >= 1904] <- "Observed"
	
	### Some example plots
	p <- ggplot(drought_event_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	p <- p + geom_point(aes(colour=timeperiod))
	#p <- p + geom_text(hjust = 0, nudge_x = 0.5, nudge_y = 0.001, check_overlap = TRUE)
	#p <- p + geom_text_repel( size=2.7, box.padding = unit(0.5, "lines"),
    #point.padding = unit(0.2, "lines"),
    #segment.color = 'grey', segment.alpha=0.4,nudge_x = 0.5,  nudge_y=0.001)
	p <- p + geom_text_repel( size=2.7,segment.color = 'grey', segment.alpha=0.4)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	#p <- p + scale_color_manual(name="Century", values=pal_cb_5_more)
	#p <- p + scale_color_manual(name="Century", values=rev(plasma(7)))
	p <- p + scale_color_manual(name="Time Period", values=c("red", "black"))
	#p <- p + scale_color_viridis(option="magma")
#	p <- p + theme(legend.position = c(0.8, 0.2))
	#p <- p + theme_classic_correct()
	p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))
	#p <- p + theme(axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))
	#	p <- p + theme(axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))

	p
ggsave(file.path(write_output_path, paste0(save_name,"_drought_perc_vs_duration.png")),  p, width=6.5, height=5, dpi=300)



	### Some example plots
	p <- ggplot(drought_event_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(aes(colour=timeperiod))
	p <- p + geom_text(aes(colour=timeperiod), size=2.7)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	#p <- p + scale_color_manual(name="Century", values=pal_cb_5_more)
	#p <- p + scale_color_manual(name="Century", values=rev(plasma(7)))
	p <- p + scale_color_manual(name="Time Period", values=c("red", "black"))
	#p <- p + scale_color_viridis(option="magma")
#	p <- p + theme(legend.position = c(0.8, 0.2))
	#p <- p + theme_classic_correct()
	p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))
	#p <- p + theme(axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))
	#	p <- p + theme(axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))

	p
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	###