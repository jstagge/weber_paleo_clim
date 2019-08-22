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
output_name <- "percent_ts"

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
require(monthlypaleo)
require(staggefuncs)
require(ggjoy)

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
site_id <- c("10128500")
site_name <- c("Weber River")


###########################################################################
###  Read in data
###########################################################################

### Read to csv
read_location <- file.path(write_output_base_path, paste0(site_id,"_obs_perc_ts.csv"))
weber_obs <- read.csv(read_location)
weber_obs$date <- as.Date(weber_obs$date)

### Read to csv
read_location <- file.path(write_output_base_path, paste0(site_id,"_climchange_perc_ts.csv"))
weber_climchange <- read.csv(read_location)
weber_climchange$date <- as.Date(weber_climchange$date)



acast(aqm, day ~ month ~ variable)
acast(aqm, month ~ variable, mean)
acast(aqm, month ~ variable, mean, margins = TRUE)
dcast(aqm, month ~ variable, mean, margins = c("month", "variable"))

dcast(weber_climchange, year + month ~ month_perc)

 dcast(weber_climchange, year ~ data, month_perc)
 
 
diet + chick ~ time, length

require(reshape2)

aql <- melt(weber_climchange, id.vars = c("month", "day"))
head(aql)


yup <- melt(weber_climchange, id.vars = c("month", "year", "data"), measure.vars = "month_perc")
yup2 <- acast(yup, year + month ~ data, mean) 

head(yup2)

 

yup <- melt(weber_climchange, id.vars = c("month", "year", "data"), measure.vars = "monthly_mean")
yup2 <- acast(yup, year + month ~ data, mean) 
yup2 <- data.frame(yup2)

head(yup2)

obs_ts <- yup2$observed
yup3 <- yup2[names(yup2) != "observed"]

yup4 <- matrix(rep(obs_ts, 5),dim(yup3))

attempt <- yup3 - yup4

yup4 <- melt(attempt)
yup5 <- subset(weber_climchange, data!="observed")



yup4 <- subset(weber_climchange, data=="observed")
yup5 <- subset(weber_climchange, data!="observed")

yup5$obs_ts <- rep(yup4$monthly_mean, 5)
yup5$diff <- yup5$monthly_mean - yup5$obs_ts

yup5 <- subset(yup5, data!="CTN5")



ggplot(yup5, aes(x=date, y=diff, color=data)) + geom_line() + theme_classic_correct()

p <- ggplot(yup5, aes(x=diff, color=data))
p <- p + geom_density()
p <- p + theme_classic_correct()
p <- p + facet_wrap( ~ month, scales="free")
p

p + scale_colour_manual(values=c("#e31a1c", "#1f78b4", "#fb9a99", "#a6cee3"))

p + scale_colour_brewer(type="seq", palette="paired")





weber_climchange$month_z <- qnorm(weber_climchange$month_perc)

p <- ggplot(weber_climchange, aes(x=month, y=month_z, group=year))
p <- p + geom_line()
p <- p + facet_wrap( ~ data, scales="free")
p <- p + theme_classic_correct()
p


p <- ggplot(yup5, aes(x=month, y=diff, group=year))
p <- p + geom_line()
p <- p + facet_wrap( ~ data)
p <- p + theme_classic_correct()
p


###########################################################################
###  Plot data
###########################################################################

p <- ggplot(weber_climchange, aes(x = monthly_mean, y = data))
p <- p + geom_joy2()
p <- p + theme_classic_correct()
p

ggplot(weber_cc_melt_df, aes(x=date, y=qnorm(month_perc), color=data))+geom_line()



p <- ggplot(weber_climchange, aes(x=date, y = monthly_mean, color=data))
p <- p + geom_line()
p <- p + theme_classic_correct()
p


p <- ggplot(weber_climchange, aes(x=date, y = month_perc, color=data))
p <- p + geom_line()
p <- p + theme_classic_correct()
p



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