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
require(staggefuncs)
require(paleoAPR)
#require(data.table)
require(ggrepel)

### Load project specific functions
file.sources = list.files(function_path, pattern="*.R", recursive=TRUE)
sapply(file.path(function_path, file.sources),source)

### Load global functions
file.sources = list.files(global_path, pattern="*.R", recursive=TRUE)
sapply(file.path(global_path, file.sources),source)


		center_mass <- function (weight, distance) {
			moment <- weight * distance

			moment_sum <- sum(moment)
			weight_sum <- sum(weight)
			return(moment_sum/weight_sum)
		}
		

###########################################################################
## Set Initial Values
###########################################################################
### Set site data
site_id_list <- c("10128500")
site_name_list <- c("Weber River")

thresh_level <- 0.5
gap_months <- 2
min_duration <- 3

###########################################################################
###  
###########################################################################
n <- 1


site_id <- site_id_list[n]
site_name <- site_name_list[n]

### Create output folders


###########################################################################
###  Read in Data
###########################################################################
### Read in observed flow
read_location <- file.path(write_output_base_path, paste0(site_id,"_obs_perc_ts.csv"))
flow_obs <- read.csv(file = read_location)

### Read in climate change
read_location <- file.path(write_output_base_path, paste0(site_id,"_climchange_perc_ts.csv"))
flow_cc <- read.csv(file = read_location)

### Read in paleo flow
read_location <- file.path(write_output_base_path, paste0(site_id,"_paleo_perc_ts.csv"))
flow_paleo <- read.csv(file = read_location)

### Combine all data
flow_all <- rbind(flow_obs, flow_cc, flow_paleo)
flow_all$data <- as.factor(flow_all$data)
flow_all$date <- as.Date(flow_all$date)

data_levels <- levels(flow_all$data)



###########################################################################
## Extract droughts based on percentile threshold
###########################################################################
for (i in seq(1, length(data_levels))){

	### Extract data
	data_name <- data_levels[i]
	flow_i <- subset(flow_all, data==data_name)

	row_index <- seq(1,dim(flow_i)[1])
	
	### Calculate water year
	flow_i$wy <- usgs_wateryear(year=flow_i$year, month=flow_i$month)

	### Cut to full water years
	avail_flow <- tapply(flow_i$flow_m3s, flow_i$wy, function(x) sum(!is.na(x)))
	avail_flow <- as.numeric(names(avail_flow)[avail_flow == 12])
	#flow_i <- flow_i[flow_i$wy >= min(avail_flow) & flow_i$wy <= max(avail_flow),]
	flow_i$flow_m3s[flow_i$wy < min(avail_flow) ] <- NA
	flow_i$flow_m3s[flow_i$wy > max(avail_flow) ] <- NA
			
	### Extract droughts
	drought_occur <- flow_i$flow_m3s < flow_i$thresh_m3s
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
	drought_event_df <- data.frame(begin=as.Date(flow_i$date[drought_begin_index]), end=as.Date(flow_i$date[drought_end_index]))
	
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
		event_subset_full <- flow_i[seq(drought_begin_index[j]-1, drought_end_index[j]+1),]
		event_subset <- event_subset_full[seq(2,dim(event_subset_full)[1]-1),]
		
		### Start and End
		drought_event_df$begin_month[j] <- event_subset$month[1]
		drought_event_df$end_month[j] <- event_subset$month[dim(event_subset)[1]]
		
		### Min flow and percentile
		event_subset$percentile <- pnorm(event_subset$norm)
		event_subset_full$percentile <- pnorm(event_subset_full$norm)
	
		drought_event_df$min_flow[j] <- min(event_subset$flow_m3s, na.rm=TRUE)
		drought_event_df$min_perc[j] <- min(event_subset$percentile, na.rm=TRUE)
		
		### Calculate Deficits
		event_subset$deficit <- event_subset$thresh_m3s - event_subset$flow_m3s
		event_subset$deficit[event_subset$deficit<0] <- 0
		
		drought_event_df$max_deficit[j] <- max(event_subset$deficit, na.rm=TRUE)
		drought_event_df$cum_deficit[j] <- sum(event_subset$deficit, na.rm=TRUE)

		### Calculate Centers of Mass
		event_subset$time_fromstart <- seq(1,dim(event_subset)[1])-1
		event_subset$perc_deficit <- thresh_level - event_subset$percentile
		event_subset$perc_deficit[event_subset$perc_deficit<0] <- 0
		
		drought_event_df$deficit_center_mass[j] <- center_mass(weight=event_subset$deficit, distance=event_subset$time_fromstart)/max(event_subset$time_fromstart)
		drought_event_df$perc_center_mass[j] <- center_mass(weight=event_subset$perc_deficit, distance=event_subset$time_fromstart)/max(event_subset$time_fromstart)

		### Determine rates of change to and from min flow
		event_subset_full$deficit <- event_subset_full$thresh_m3s - event_subset_full$flow_m3s
		event_subset_full$time_fromstart <- seq(1,dim(event_subset_full)[1])-1
		event_subset_full$perc_deficit <- thresh_level - event_subset_full$percentile
	
		min_date <- which.max(event_subset_full$deficit)
		drought_event_df$dec_rate_to_max_deficit[j] <- -(event_subset_full$deficit[min_date] - event_subset_full$deficit[1] )/(min_date-1)

		length_dates <- length(event_subset_full$flow_m3s)
		drought_event_df$inc_rate_from_max_deficit[j] <- -( event_subset_full$deficit[length_dates] - event_subset_full$deficit[min_date])/(length_dates - min_date)


		### Determine rates of change to and from min percentile	
		min_date <- which.max(event_subset_full$perc_deficit)
		drought_event_df$dec_rate_to_min_perc[j] <- -(event_subset_full$perc_deficit[min_date] - event_subset_full$perc_deficit[1] )/(min_date-1)

		length_dates <- length(event_subset_full$flow_m3s)
		drought_event_df$inc_rate_from_min_perc[j] <- -( event_subset_full$perc_deficit[length_dates] - event_subset_full$perc_deficit[min_date])/(length_dates - min_date)

	}
	
	### Prior year
	prior_wy <- flow_i$wy[drought_begin_index] - 1
	
	### Loop through each event and calculate information
	for (j in seq(1, length(drought_begin_index))) {
		prior_subset <- flow_i[flow_i$wy == prior_wy[j],]
		drought_event_df$prior_year_flow[j] <- sum(prior_subset$flow_m3s)
	}
	


drought_event_df$data <- data_name

if (i == 1) {
	drought_event_summary <- drought_event_df
} else {
	drought_event_summary <- rbind(drought_event_summary, drought_event_df)
}
}



###########################################################################
## Save Drought Statistics
###########################################################################
save_file <- file.path(write_output_base_path, paste0(site_id,"_drought_details.csv"))

write.csv(drought_event_summary, save_file, row.names = FALSE)








































###########################################################################
## Plot summary statistics
###########################################################################
	
	drought_event_df$century <- factor(paste0(substr(drought_event_df$begin,1,2),"00s"))
	drought_event_df$timeperiod <- "Prior Period"
	drought_event_df$timeperiod[as.numeric(substr(drought_event_df$begin,1,4)) >= 1900] <- "Obs Period"

	drought_event_df$century <- factor(paste0(substr(drought_event_df$begin,1,2),"00s"))
	drought_event_df$timeperiod <- "Prior"
	drought_event_df$timeperiod[as.numeric(substr(drought_event_df$begin,1,4)) >= 1904] <- "Observed"






### Some example plots
	p <- ggplot(subset(drought_event_summary, data=="observed" | data=="paleo" | data=="HDN5"), aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	p <- p + geom_point(aes(colour=data))
	#p <- p + geom_text_repel( size=2.7,segment.color = 'grey', segment.alpha=0.4)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
		p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))





### Some example plots
	p <- ggplot(subset(drought_event_summary, data!="observed" & data!="paleo"& data!="CTN5"), aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	p <- p + geom_point(aes(colour=data))
	p <- p + geom_text_repel( size=2.7,segment.color = 'grey', segment.alpha=0.4)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
		p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))

p


### Some example plots
observed_background <- subset(drought_event_summary, data=="observed")
cc_event_2000 <- subset(test_subset, year(test_subset$begin) >1999 & year(test_subset$begin) < 2001 & data!="observed" & data!="paleo"& data!="CTN5")

	p <- ggplot(observed_background, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	p <- p + geom_point(colour="grey50")
	p <- p + geom_point(data=cc_event_2000, aes(colour=data))
	p <- p + geom_text_repel(data=cc_event_2000, size=2.7,segment.color = 'grey', segment.alpha=0.4)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
		p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))

p



















###########################################################################
## Calculate cumulative loss
###########################################################################

two_year_limit <- 2*sum(flow_all$thresh_m3s[1:12])

### Set up cumulative loss column
flow_all$cum_loss <- NA
	
for (i in data_levels){
	### Extract data
	i_test <- flow_all$data == i
	flow_i <- flow_all[i_test, ]

	### Loop through each month
	for (j in 2:dim(flow_i)[1]) {

		loss_j <- flow_i$flow_m3s[j] - flow_i$thresh_m3s[j]

		if (!is.na(loss_j)) {
			### If previously is NA and there is loss, start a new drought
			if (is.na(flow_i$cum_loss[j-1])){
				if (loss_j < 0) {flow_i$cum_loss[j] <- loss_j}
			### If previously in drought, does it end or continue
			} else {
			### Cumulative loss is previous plus current loss
				cum_loss_j <- flow_i$cum_loss[j-1] + loss_j

				## If current cumulative loss is les than zero, add to column
				if (cum_loss_j < 0) {
					flow_i$cum_loss[j] <- cum_loss_j
				} else if (flow_i$cum_loss[j-1] == 0){
					flow_i$cum_loss[j] <- NA
				} else {
					flow_i$cum_loss[j] <- 0
				}
		}}
		
		if(flow_i$cum_loss[j] < -two_year_limit && is.na(flow_i$cum_loss[j])!=TRUE) {flow_i$cum_loss[j] <- -two_year_limit}
	}
	### Re-insert cumulative loss into the main dataframe
	flow_all$cum_loss[i_test] <- flow_i$cum_loss
}


###########################################################################
## Test Plot
###########################################################################
plot_test <- flow_all$data %in% c("observed", "HDN5", "HWN5", "WDN5", "WWN5", "paleo")

p <- ggplot(flow_all[plot_test,], aes(x=date, y=cum_loss, colour=data)) + geom_line() + theme_classic_new() + scale_colour_manual(values=c("black", "orange", "green",  "yellow", "blue","purple"))
p

ggsave("cum_shortfall_full.png", p, width=10, height=5, dpi=600)

p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01"))) + scale_x_date(breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

ggsave("cum_shortfall_1900s.png", p, width=7, height=5, dpi=600)



###########################################################################
## Extract droughts based on percentile threshold
###########################################################################
for (i in seq(1, length(data_levels))){

data_name <- data_levels[i]


flow_i <- subset(flow_all, data==data_name)


row_index <- seq(1,dim(flow_i)[1])
	
### Extract droughts
#drought_occur <- flow_i$flow_m3s < flow_i$thresh_m3s
flow_i$cum_loss[is.na(flow_i$cum_loss)] <- 0
drought_occur <- flow_i$cum_loss < 0
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
	drought_event_df <- data.frame(begin=as.Date(flow_i$date[drought_begin_index]), end=as.Date(flow_i$date[drought_end_index]))
	
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
		event_subset_full <- flow_i[seq(drought_begin_index[j]-1, drought_end_index[j]+1),]
		event_subset <- event_subset_full[seq(2,dim(event_subset_full)[1]-1),]
		
		### Start and End
		drought_event_df$begin_month[j] <- event_subset$month[1]
		drought_event_df$end_month[j] <- event_subset$month[dim(event_subset)[1]]
		
		### Min flow and percentile
		drought_event_df$min_flow[j] <- min(event_subset$flow_m3s, na.rm=TRUE)
		drought_event_df$min_perc[j] <- min(pnorm(event_subset$norm), na.rm=TRUE)
		
		### Calculate Deficits
		event_subset$deficit <- event_subset$thresh_m3s - event_subset$flow_m3s
		event_subset$deficit[event_subset$deficit<0] <- 0
		
		drought_event_df$max_deficit[j] <- max(event_subset$deficit, na.rm=TRUE)
		drought_event_df$cum_deficit[j] <- max(-event_subset$cum_loss, na.rm=TRUE)
}


drought_event_df$data <- data_name

if (i == 1) {
	drought_event_summary <- drought_event_df
} else {
	drought_event_summary <- rbind(drought_event_summary, drought_event_df)
}
}






p <- ggplot(subset(drought_event_summary, data=="observed" | data=="paleo"), aes(x=dura_months/12, y=cum_deficit, colour=data,  label=substr(begin,1,4)))
	p <- p + geom_point()
	p
	

	p <- ggplot(drought_event_summary, aes(x=dura_months/12, y=cum_deficit,  label=substr(begin,1,4)))
	p <- p + geom_point()





### Some example plots
	p <- ggplot(subset(drought_event_summary, data=="observed" | data=="paleo"), aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	p <- p + geom_point(aes(colour=data))
	#p <- p + geom_text(hjust = 0, nudge_x = 0.5, nudge_y = 0.001, check_overlap = TRUE)
	#p <- p + geom_text_repel( size=2.7, box.padding = unit(0.5, "lines"),
    #point.padding = unit(0.2, "lines"),
    #segment.color = 'grey', segment.alpha=0.4,nudge_x = 0.5,  nudge_y=0.001)
	p <- p + geom_text_repel( size=2.7,segment.color = 'grey', segment.alpha=0.4)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	#p <- p + scale_color_manual(name="Century", values=pal_cb_5_more)
	#p <- p + scale_color_manual(name="Century", values=rev(plasma(7)))
	#p <- p + scale_color_manual(name="Time Period", values=c("red", "black"))
	#p <- p + scale_color_viridis(option="magma")
#	p <- p + theme(legend.position = c(0.8, 0.2))
	#p <- p + theme_classic_correct()
	p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))







### Some example plots
	p <- ggplot(subset(drought_event_summary, data=="observed" | data=="paleo"), aes(x=dura_months/12, y=cum_deficit, label=substr(begin,1,4)))
	p <- p + geom_point(aes(colour=data))
	#p <- p + geom_text(hjust = 0, nudge_x = 0.5, nudge_y = 0.001, check_overlap = TRUE)
	#p <- p + geom_text_repel( size=2.7, box.padding = unit(0.5, "lines"),
    #point.padding = unit(0.2, "lines"),
    #segment.color = 'grey', segment.alpha=0.4,nudge_x = 0.5,  nudge_y=0.001)
	p <- p + geom_text_repel(data=subset(drought_event_summary, data=="observed"),  size=2.7,segment.color = 'grey', segment.alpha=0.4)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	#p <- p + scale_color_manual(name="Century", values=pal_cb_5_more)
	#p <- p + scale_color_manual(name="Century", values=rev(plasma(7)))
	#p <- p + scale_color_manual(name="Time Period", values=c("red", "black"))
	#p <- p + scale_color_viridis(option="magma")
#	p <- p + theme(legend.position = c(0.8, 0.2))
	#p <- p + theme_classic_correct()
	p <- p + theme_classic_new(14)
	p <- p + theme(legend.position = c(0.85, 0.85))

p


		### Calculate Centers of Mass
		event_subset$time_fromstart <- seq(1,dim(event_subset)[1])-1
		event_subset$perc_deficit <- thresh_level - pnorm(event_subset$norm)
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
	
	


flow_x$cum_loss <- NA

for (j in 2:dim(flow_x)[1]) {

	loss_j <- flow_x$flow_m3s[j] - flow_x$thresh_m3s[j]

	if (!is.na(loss_j)) {
	### If previously is NA and there is loss, start a new drought
	if (is.na(flow_x$cum_loss[j-1])){
		if (loss_j < 0) {flow_x$cum_loss[j] <- loss_j}
	### If previously in drought, does it end or continue
	} else {
		### Cumulative loss is previous plus current loss
		cum_loss_j <- flow_x$cum_loss[j-1] + loss_j

		### If current cumulative loss is les than zero, add to column
		if (cum_loss_j < 0) {
			flow_x$cum_loss[j] <- cum_loss_j
		} else if (flow_x$cum_loss[j-1] == 0){
			flow_x$cum_loss[j] <- NA
		} else {
			flow_x$cum_loss[j] <- 0
		}
	}}
}





flow_i$loss <- flow_i$flow_m3s - flow_i$thresh_m3s
flow_i$excess <- flow_i$loss
flow_i$loss[flow_i$loss > 0] <- 0
flow_i$excess[flow_i$excess < 0] <- 0

plot(flow_i$excess, type="l")
lines(-flow_i$loss, col="red")

plot(cumsum(flow_i$excess), type="l")
lines(-cumsum(flow_i$loss), col="red")



