
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


###########################################################################
## Set Initial Values
###########################################################################
### Set site data
site_id_list <- c("10128500")
site_name_list <- c("Weber River")
n <- 1

site_id <- site_id_list[n]
site_name <- site_name_list[n]


### Colors to use for climate change scenarios
cc_colors <- c("#D55E00" ,"#56B4E9", "#E69F00" , "#0072B2", "#CC79A7")

### Conversion from m3/s to acft/month
m3s_to_acftmonth <- 2131.97

###########################################################################
###  Read in Data
###########################################################################
### Read in observed flow
#read_location <- file.path(write_output_base_path, paste0(site_id,"_obs_perc_ts.csv"))
#flow_obs <- read.csv(file = read_location)

### Read in climate change
read_location <- file.path(write_output_base_path, paste0(site_id,"_drought_details.csv"))

drought_event_summary <- read.csv(file = read_location)



###########################################################################
###   Distribution of duration
###########################################################################
#data!="base" &
plot_df <- subset(drought_event_summary,  data!="CTN5")
plot_df$data <- factor(plot_df$data, levels=c("paleo", "observed", "base", "WWN5","HWN5","WDN5", "HDN5"), labels= c("Reconst", "Observed",  "Base", "WW", "HW", "WD", "HD"))

p <- ggplot(plot_df, aes(x=dura_months/12, y = ..density.., fill=data))
p <- p + geom_histogram(binwidth=0.5)
p <- p + geom_hline(yintercept=0, colour="grey60")
p <- p + theme_classic_new(12)
p <- p + facet_grid(data ~ .)
p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,1))
p <- p + scale_y_continuous(name="Density", breaks=seq(0,2, 0.5))
p <- p + scale_fill_manual(name=NULL, values= c( cc_colors[c(5)], "black","grey60", cc_colors[c(4,2,3,1)]), labels=c("Reconstructed", "Observed", "Base", "Warm-Wet", "Hot-Wet", "Warm-Dry", "Hot-Dry)"), limits=c("Reconst", "Observed", "Base", "WW", "HW", "WD", "HD"))
p <- p + coord_cartesian(xlim = NULL, ylim = c(0,1.6), expand = FALSE)
p <- p + theme(
  strip.background = element_blank(),
  strip.text = element_blank()
)
#p <- p + theme(legend.position = c(0.85, 0.85))
p <- p+ theme(legend.position="none")

### Save figures
ggsave(file.path(write_output_base_path,"duration_hist.png"),  p, width=5, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"duration_hist.pdf"),  p, width=5, height=6)
ggsave(file.path(write_output_base_path,"duration_hist.svg"),  p, width=5, height=6)

p <- p +theme(legend.position="bottom")
p <- p + guides(fill=guide_legend(nrow=1, byrow=TRUE))
p

### Save figures
ggsave(file.path(write_output_base_path,"duration_hist_legend.svg"),  p, width=10, height=6)


p <- ggplot(plot_df, aes(x=data, y=dura_months/12, fill=data))
p <- p + geom_boxplot(alpha=0.8)
p <- p + theme_classic_new(11)
#p <- p + facet_grid(data ~ .)
#p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
p <- p + scale_x_discrete(name="Scenario", limits=c("Reconst", "Observed",  "Base", "WW", "HW", "WD", "HD"), labels=c("Reconst", "Observed", "Base\nClim Change", "Warm\nWet", "Hot\nWet", "Warm\nDry", "Hot\nDry"))
p <- p + scale_y_continuous(name="Drought Duration (Years)", breaks=seq(0,30,1))
#p <- p + coord_flip()
p <- p + scale_fill_manual(name=NULL, values= c( cc_colors[c(5)], "grey20","grey60", cc_colors[c(4,2,3,1)]), labels=c("Reconstructed", "Observed", "Base", "Warm\nWet", "Hot\nWet", "Warm\nDry", "Hot\nDry"), limits=c("Reconst", "Observed",  "Base", "WW", "HW", "WD", "HD"))
p <- p + theme(legend.position="none")
p

### Save figures
ggsave(file.path(write_output_base_path,"duration_box.png"),  p, width=5, height=4.5, dpi=300)
ggsave(file.path(write_output_base_path,"duration_box.pdf"),  p, width=5, height=4.5)
ggsave(file.path(write_output_base_path,"duration_box.svg"),  p, width=5, height=4.5)



p <- ggplot(plot_df, aes(x=dura_months/12,  colour=data))
p <- p + geom_density(alpha=0.7)
p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,1))
p <- p + scale_y_continuous(name="Density")
p <- p + scale_color_manual(name="Scenario", values= c( cc_colors[c(5)], "grey40",cc_colors[c(4,2,3,1)]), labels=c("Reconstructed", "Observed",  "WW (Clim Change)", "HW (Clim Change)", "WD (Clim Change)", "HD (Clim Change)"), limits=c("Reconst", "Observed",  "WW", "HW", "WD", "HD"))
p <- p + theme_classic_new()
p <- p + theme(legend.position = c(0.85, 0.75))
p <- p + coord_cartesian(xlim = NULL, ylim = c(0,0.9), expand = FALSE)

p

### Save figures
ggsave(file.path(write_output_base_path,"duration_density.png"),  p, width=6, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"duration_density.pdf"),  p, width=6, height=5)
ggsave(file.path(write_output_base_path,"duration_density.svg"),  p, width=6, height=5)

### Could do a Kolmogorov-Smirnov test on distributions, e.g. is future significantly different from current?
### Is current significantly different from paleo?


###########################################################################
###   Distribution of duration
###########################################################################

plot_df <- subset(drought_event_summary,  data!="CTN5")
plot_df$data <- factor(plot_df$data, levels=c("paleo", "observed", "base", "WWN5","HWN5","WDN5", "HDN5"), labels= c("Reconst", "Observed",  "Base", "WW", "HW", "WD", "HD"))

p <- ggplot(plot_df, aes(x=min_flow*m3s_to_acftmonth/1000, y = ..density.., fill=data))
p <- p + geom_histogram()
p <- p + geom_hline(yintercept=0, colour="grey80")
p <- p + theme_classic_new(12)
p <- p + facet_grid(data ~ .)
p <- p + scale_x_continuous(name="Minimum Flow (1,000 ac-ft/month)", breaks=seq(0,30,1))
p <- p + scale_y_continuous(name="Density")#, breaks=seq(0,2, 0.5))
p <- p + scale_fill_manual(name=NULL, values= c( cc_colors[c(5)], "black","grey60", cc_colors[c(4,2,3,1)]), labels=c("Reconstructed", "Observed", "Base", "Warm-Wet", "Hot-Wet", "Warm-Dry", "Hot-Dry)"), limits=c("Reconst", "Observed",  "Base", "WW", "HW", "WD", "HD"))
p <- p + coord_cartesian(xlim = c(1,10), ylim = NULL, expand = FALSE)
p <- p + theme(
  strip.background = element_blank(),
  strip.text = element_blank()
)
#p <- p + theme(legend.position = c(0.85, 0.85))
p <- p+ theme(legend.position="none")
p

### Save figures
ggsave(file.path(write_output_base_path,"min_flow_hist.png"),  p, width=5, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"min_flow_hist.pdf"),  p, width=5, height=6)
ggsave(file.path(write_output_base_path,"min_flow_hist.svg"),  p, width=5, height=6)



p <- ggplot(plot_df, aes(x=data, y=min_flow*m3s_to_acftmonth/1000, fill=data))
p <- p + geom_boxplot(alpha=0.8)
p <- p + theme_classic_new(11)
#p <- p + facet_grid(data ~ .)
#p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
p <- p + scale_x_discrete(name="Scenario", limits=c("Reconst", "Observed",  "Base", "WW", "HW", "WD", "HD"), labels=c("Reconst", "Observed", "Base\nClim Change", "Warm\nWet", "Hot\nWet", "Warm\nDry", "Hot\nDry"))
p <- p + scale_y_continuous(name="Minimum Flow (1,000 ac-ft/month)", breaks=seq(0,30,1))
#p <- p + coord_flip()
p <- p + scale_fill_manual(name=NULL, values= c( cc_colors[c(5)], "grey20","grey60", cc_colors[c(4,2,3,1)]), labels=c("Reconstructed", "Observed", "Base", "Warm\nWet", "Hot\nWet", "Warm\nDry", "Hot\nDry"), limits=c("Reconst", "Observed",  "Base", "WW", "HW", "WD", "HD"))
p <- p+ theme(legend.position="none")
p

### Save figures
ggsave(file.path(write_output_base_path,"min_flow_box.png"),  p, width=5, height=4.5, dpi=300)
ggsave(file.path(write_output_base_path,"min_flow_box.pdf"),  p, width=5, height=4.5)
ggsave(file.path(write_output_base_path,"min_flow_box.svg"),  p, width=5, height=4.5)



p <- ggplot(plot_df, aes(x=min_flow*m3s_to_acftmonth/1000,  colour=data))
p <- p + geom_density(alpha=0.7)
p <- p + scale_x_continuous(name="Minimum Flow (1,000 ac-ft/month))", breaks=seq(0,30,1))
p <- p + scale_y_continuous(name="Density")
p <- p + scale_color_manual(name="Scenario", values= c( cc_colors[c(5)], "grey40",cc_colors[c(4,2,3,1)]), labels=c("Reconstructed", "Observed",  "WW (Clim Change)", "HW (Clim Change)", "WD (Clim Change)", "HD (Clim Change)"), limits=c("Reconst", "Observed",  "WW", "HW", "WD", "HD"))
p <- p + theme_classic_new()
p <- p + theme(legend.position = c(0.85, 0.75))
p <- p + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)

p

### Save figures
ggsave(file.path(write_output_base_path,"min_flow_density.png"),  p, width=6, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"min_flow_density.pdf"),  p, width=6, height=5)
ggsave(file.path(write_output_base_path,"min_flow_density.svg"),  p, width=6, height=5)

### Could do a Kolmogorov-Smirnov test on distributions, e.g. is future significantly different from current?
### Is current significantly different from paleo?




###########################################################################
###  Plot dur vs percentile
###########################################################################
plot_df <- subset(drought_event_summary, data=="observed" | data=="paleo" | data=="HDN5")

	### Some example plots
	p <- ggplot(plot_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
	p <- p + geom_text(aes(colour=data), size=2.7)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Monthly Flow Percentile")
	#p <- p + scale_color_manual(name="Data Source", values=c("black","cadetblue",  "red"))
	p <- p + scale_colour_manual(name="Scenario", values= c(  "grey20","#377eb8", "#e41a1c"), labels=c("Observed",  "Reconstructed", "Climate Change (HD)"), limits=c("observed", "paleo", "HDN5" ))
	p <- p + theme_classic_new(12)
	p <- p + theme(legend.position = c(0.85, 0.85))
	p

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_text.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_text.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_text.svg"),  p, width=6.5, height=5)


obs_plot <- subset(drought_event_summary, data=="observed")
paleo_plot <- subset(drought_event_summary, data=="paleo")
cc_plot <- subset(drought_event_summary, data=="HDN5")
	
	### Some example plots
	p <- ggplot(obs_plot, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
	p <- p + geom_point(aes(colour=data), size=2.5)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,1))
	p <- p + scale_y_continuous(name="Min Monthly Flow Percentile")
	#p <- p + scale_color_manual(name="Data Source", values=c("black","cadetblue",  "red"))
	p <- p + scale_colour_manual(name="Scenario", values= c(  "grey20","#377eb8", "#e41a1c"), labels=c("Observed",  "Reconstructed", "Climate Change (HD)"), limits=c("observed", "paleo", "HDN5" ))
#	p <- p + scale_colour_manual(name="Scenario", values= c(  "grey20","#80b1d3", "#fb8072"), labels=c("Observed",  "Reconstr", "Climate Change (HD)"), limits=c("observed", "paleo", "HDN5" ))
	p <- p + theme_classic_new(12)
	p <- p + theme(legend.position = c(0.85, 0.85))
	p <- p + coord_cartesian(xlim = c(0.35,6.5), ylim = c(0,47), expand = TRUE)
	p

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_points_obsonly.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_points_obsonly.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_points_obsonly.svg"),  p, width=6.5, height=5)


	### Some example plots
	p <- ggplot(obs_plot, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
	p <- p + geom_point(data=paleo_plot, aes(colour=data), size=2.5, alpha=0.4)	
	p <- p + geom_point(aes(colour=data), size=2.5)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,1))
	p <- p + scale_y_continuous(name="Min Monthly Flow Percentile")
	#p <- p + scale_color_manual(name="Data Source", values=c("black","cadetblue",  "red"))
	p <- p + scale_colour_manual(name="Scenario", values= c(  "grey20","#377eb8", "#e41a1c"), labels=c("Observed",  "Reconstructed", "Climate Change (HD)"), limits=c("observed", "paleo", "HDN5" ))
#	p <- p + scale_colour_manual(name="Scenario", values= c(  "grey20","#80b1d3", "#fb8072"), labels=c("Observed",  "Reconstr", "Climate Change (HD)"), limits=c("observed", "paleo", "HDN5" ))
	p <- p + theme_classic_new(12)
	p <- p + theme(legend.position = c(0.85, 0.85))
	p <- p + coord_cartesian(xlim = c(0.35,6.5), ylim = c(0,47), expand = TRUE)
	p

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_points_paleoandobs.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_points_paleoandobs.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_points_paleoandobs.svg"),  p, width=6.5, height=5)



### Some example plots
	p <- ggplot(obs_plot, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
	p <- p + geom_point(data=paleo_plot, aes(colour=data), size=2.5, alpha=0.4)	
	p <- p + geom_point(aes(colour=data), size=2.5)
	p <- p + geom_point(data=cc_plot, aes(colour=data), size=2.5, alpha=0.9)	
	
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,1))
	p <- p + scale_y_continuous(name="Min Monthly Flow Percentile")
	#p <- p + scale_color_manual(name="Data Source", values=c("black","cadetblue",  "red"))
	p <- p + scale_colour_manual(name="Scenario", values= c(  "grey20","#377eb8", "#e41a1c"), labels=c("Observed",  "Reconstructed", "Climate Change (HD)"), limits=c("observed", "paleo", "HDN5" ))
#	p <- p + scale_colour_manual(name="Scenario", values= c(  "grey20","#80b1d3", "#fb8072"), labels=c("Observed",  "Reconstr", "Climate Change (HD)"), limits=c("observed", "paleo", "HDN5" ))
	p <- p + theme_classic_new(12)
	p <- p + theme(legend.position = c(0.85, 0.85))
	p <- p + coord_cartesian(xlim = c(0.35,6.5), ylim = c(0,47), expand = TRUE)
	p

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_points.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_points.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_points.svg"),  p, width=6.5, height=5)




###########################################################################
###  Plot dur vs percentile for 2000 event
###########################################################################

drought_event_summary$year <- substr(drought_event_summary$begin, 1, 4)

background_df <- subset(drought_event_summary, data=="observed" | data=="paleo")
cc_df <- subset(drought_event_summary, data!="observed" & data!="paleo" & data!="CTN5")

only_event <- subset(cc_df, year >= 2000 & year <=2001)

	### Some example plots
	p <- ggplot(background_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	p <- p + geom_point(size=2.7, alpha=0.3, colour="grey70")
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
	#p <- p + geom_text(data=only_event, aes(colour=data), size=2.7)
	p <- p + geom_point(data=only_event, aes(colour=data), size=5)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Monthly Flow Percentile")
	#p <- p + scale_color_manual(name="Data Source", values=c("black","cadetblue",  "red"))
	p <- p + scale_colour_manual(name="Scenario", values= c( "black",cc_colors[c(4,2,3,1)]), labels=c("Observed",  "WW (Clim Change)", "HW (Clim Change)", "WD (Clim Change)", "HD (Clim Change)"), limits=c("base", "WWN5","HWN5","WDN5", "HDN5"))
	p <- p + theme_classic_new(12)
	p <- p + theme(legend.position = c(0.85, 0.85))
	p

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_climchange_2000.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_climchange_2000.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_climchange_2000.svg"),  p, width=6.5, height=5)


###########################################################################
###  Plot dur vs percentile for 1987 event
###########################################################################

background_df <- subset(drought_event_summary, data=="observed" | data=="paleo")
cc_df <- subset(drought_event_summary, data!="observed" & data!="paleo" & data!="CTN5")

only_event <- subset(cc_df, year >= 1987 & year <=1989)

	### Some example plots
	p <- ggplot(background_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	p <- p + geom_point(size=2.7, alpha=0.3, colour="grey70")
	#p <- p + geom_point(data=paleo_drought_event_df, aes(colour=timeperiod))
#	p <- p + geom_text(data=only_event, aes(colour=data), size=2.7)
	p <- p + geom_point(data=only_event, aes(colour=data), size=5)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,2))
	p <- p + scale_y_continuous(name="Min Monthly Flow Percentile")
	p <- p + scale_colour_manual(name="Scenario", values= c( "black",cc_colors[c(4,2,3,1)]), labels=c("Observed",  "WW (Clim Change)", "HW (Clim Change)", "WD (Clim Change)", "HD (Clim Change)"), limits=c("base", "WWN5","HWN5","WDN5", "HDN5"))

	p <- p + theme_classic_new(12)
	p <- p + theme(legend.position = c(0.85, 0.85))
	p

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_climchange_1987.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_climchange_1987.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_climchange_1987.svg"),  p, width=6.5, height=5)



































###########################################################################
###  Prepare data
###########################################################################
### Remove CTN5
flow_cc <- subset(flow_cc, data!="CTN5")

### Add water year column
flow_cc$wy <- usgs_wateryear(year=flow_cc$year, month=flow_cc$month)


### Make months a factor
flow_cc$month <- factor(flow_cc$month, levels=c(seq(10, 12), seq(1,9)))


###########################################################################
###  Make data descriptors
###########################################################################
### Add column for temperature
flow_cc$temp <- NA
flow_cc$temp[flow_cc$data == "HDN5" | flow_cc$data == "HWN5"] <- "Hot"
flow_cc$temp[flow_cc$data == "WDN5" | flow_cc$data == "WWN5"] <- "Warm"
flow_cc$temp[flow_cc$data == "base"] <- "Base"
flow_cc$temp <- factor(flow_cc$temp, levels= c("Base", "Warm", "Hot"))

### Add column for precipitation
flow_cc$precip <- NA
flow_cc$precip[flow_cc$data == "HDN5" | flow_cc$data == "WDN5"] <- "Dry"
flow_cc$precip[flow_cc$data == "HWN5" | flow_cc$data == "WWN5"] <- "Wet"
flow_cc$precip[flow_cc$data == "base"] <- "Base"
flow_cc$precip <- factor(flow_cc$precip, levels= c("Base", "Wet", "Dry"))


###########################################################################
###  Plot Climate Change Difference
###########################################################################
### Separate climate change and base
flow_cc_only <- subset(flow_cc, data!="base")
flow_base <- subset(flow_cc, data=="base")

### Merge base back with climate change
flow_cc_only$base_m3s <- rep(flow_base$flow_m3s, 4)

### Subtract flows
flow_cc_only$flow_diff <- flow_cc_only$flow_m3s - flow_cc_only$base_m3s


p <- ggplot(flow_cc_only, aes(x=month))
p <- p + geom_hline(yintercept=0, colour="black", linetype="dotted")
p <- p + geom_line(aes(y=flow_diff, group=wy, colour=data))
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(values= cc_colors)
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow Change (m3/s)") + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_change_line_m3s.png"),  p, width=6, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_change_line_m3s.pdf"),  p, width=6, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_change_line_m3s.svg"),  p, width=6, height=6)



p <- ggplot(flow_cc_only, aes(x=month))
p <- p + geom_hline(yintercept=0, colour="black", linetype="dotted")
p <- p + geom_line(aes(y=flow_diff*m3s_to_acftmonth/1000, group=wy, colour=data))
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(values= cc_colors)
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow Change (1,000 ac-ft/month)") + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_change_line_acft.png"),  p, width=6, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_change_line_acft.pdf"),  p, width=6, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_change_line_acft.svg"),  p, width=6, height=6)



p <- ggplot(flow_cc_only, aes(x=month))
p <- p + geom_hline(yintercept=0, colour="black", linetype="dotted")
p <- p + geom_boxplot(aes( y=flow_diff, colour=data), fill=NA)
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(values= cc_colors)
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow Change (m3/s)") + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_change_box_m3s.png"),  p, width=6, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_change_box_m3s.pdf"),  p, width=6, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_change_box_m3s.svg"),  p, width=6, height=6)

p <- ggplot(flow_cc_only, aes(x=month))
p <- p + geom_hline(yintercept=0, colour="black", linetype="dotted")
p <- p + geom_boxplot(aes( y=flow_diff*m3s_to_acftmonth/1000, colour=data), fill=NA)
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(values= cc_colors)
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow Change (1,000 ac-ft/month)") + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_change_box_acft.png"),  p, width=6, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_change_box_acft.pdf"),  p, width=6, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_change_box_acft.svg"),  p, width=6, height=6)





###########################################################################
###  Plot Climate Change Flows
###########################################################################
cutoff <- data.frame(subset(flow_cc,wy==2000), cutoff="Median\n1925-2005")

p <- ggplot(flow_cc, aes(x=month))
p <- p + geom_line(aes(y=flow_m3s, group=wy, colour=data), size=0.3)
#p <- p + geom_line(data= aes(y=thresh_m3s, group=1), colour="black", linetype="longdash")
p <- p + geom_line(aes( x=month, y=thresh_m3s, linetype = cutoff, group=wy), cutoff, colour="black")
#p <- p + facet_wrap( ~ data, ncol=3)
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(name="Scenario", values= c("grey60", cc_colors), labels=c("Base", "HD", "HW", "WD", "WW"))
p <- p + scale_linetype_manual(name="Reference", values=c("longdash"))
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow (m3/s)", breaks=seq(0,100,10)) + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_line_m3s.png"),  p, width=9, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_line_m3s.pdf"),  p, width=9, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_line_m3s.svg"),  p, width=9, height=6)


p <- ggplot(flow_cc, aes(x=month))
p <- p + geom_line(aes(y=flow_m3s*m3s_to_acftmonth/1000, group=wy, colour=data), size=0.3)
#p <- p + geom_line(data= aes(y=thresh_m3s, group=1), colour="black", linetype="longdash")
p <- p + geom_line(aes( x=month, y=thresh_m3s*m3s_to_acftmonth/1000, linetype = cutoff, group=wy), cutoff, colour="black")
#p <- p + facet_wrap( ~ data, ncol=3)
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(name="Scenario", values= c("grey60", cc_colors), labels=c("Base", "HD", "HW", "WD", "WW"))
p <- p + scale_linetype_manual(name="Reference", values=c("longdash"))
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow (1,000 ac-ft/month)", breaks=seq(0,160,20)) + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_line_acft.png"),  p, width=9, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_line_acft.pdf"),  p, width=9, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_line_acft.svg"),  p, width=9, height=6)


cutoff_nobase <- data.frame(subset(flow_cc_only,wy==2000), cutoff="Median\n1925-2005")


p <- ggplot(flow_cc_only, aes(x=month))
p <- p + geom_line(aes(y=flow_m3s, group=wy, colour=data), size=0.3)
#p <- p + geom_line(data= aes(y=thresh_m3s, group=1), colour="black", linetype="longdash")
p <- p + geom_line(aes( x=month, y=thresh_m3s, linetype = cutoff, group=wy), cutoff_nobase, colour="black")
#p <- p + facet_wrap( ~ data, ncol=3)
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(name="Scenario", values= cc_colors, labels=c("Base", "HD", "HW", "WD", "WW"))
p <- p + scale_linetype_manual(name="Reference", values=c("longdash"))
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow (m3/s)", breaks=seq(0,100,10)) + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_line_m3s_nobase.png"),  p, width=8, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_line_m3s_nobase.pdf"),  p, width=8, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_line_m3s_nobase.svg"),  p, width=8, height=6)


p <- ggplot(flow_cc_only, aes(x=month))
p <- p + geom_line(aes(y=flow_m3s*m3s_to_acftmonth/1000, group=wy, colour=data), size=0.3)
#p <- p + geom_line(data= aes(y=thresh_m3s, group=1), colour="black", linetype="longdash")
p <- p + geom_line(aes( x=month, y=thresh_m3s*m3s_to_acftmonth/1000, linetype = cutoff, group=wy), cutoff_nobase, colour="black")
#p <- p + facet_wrap( ~ data, ncol=3)
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(name="Scenario", values= cc_colors, labels=c("Base", "HD", "HW", "WD", "WW"))
p <- p + scale_linetype_manual(name="Reference", values=c("longdash"))
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow (1,000 ac-ft/month)", breaks=seq(0,160,20)) + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_line_acft_nobase.png"),  p, width=8, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_line_acft_nobase.pdf"),  p, width=8, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_line_acft_nobase.svg"),  p, width=8, height=6)





###########################################################################
###  Plot Climate Change Flow boxplots
###########################################################################
cutoff <- data.frame(subset(flow_cc,wy==2000), cutoff="Median\n1925-2005")

p <- ggplot(flow_cc, aes(x=month))
p <- p +geom_boxplot(aes( y=flow_m3s, colour=data), fill=NA, size=0.4)
p <- p + geom_line(aes( x=month, y=thresh_m3s, linetype = cutoff, group=wy), cutoff, colour="black")
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(name="Scenario", values= c("grey60", cc_colors), labels=c("Base", "HD", "HW", "WD", "WW"))
p <- p + scale_linetype_manual(name="Reference", values=c("longdash"))
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow (m3/s)", breaks=seq(0,100,10)) + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_box_m3s.png"),  p, width=9, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_box_m3s.pdf"),  p, width=9, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_box_m3s.svg"),  p, width=9, height=6)


p <- ggplot(flow_cc, aes(x=month))
p <- p +geom_boxplot(aes( y=flow_m3s*m3s_to_acftmonth/1000, colour=data), fill=NA, size=0.4)
p <- p + geom_line(aes( x=month, y=thresh_m3s*m3s_to_acftmonth/1000, linetype = cutoff, group=wy), cutoff, colour="black")
#p <- p + facet_wrap( ~ data, ncol=3)
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(name="Scenario", values= c("grey60", cc_colors), labels=c("Base", "HD", "HW", "WD", "WW"))
p <- p + scale_linetype_manual(name="Reference", values=c("longdash"))
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow (1,000 ac-ft/month)", breaks=seq(0,160,20)) + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_box_acft.png"),  p, width=9, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_box_acft.pdf"),  p, width=9, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_box_acft.svg"),  p, width=9, height=6)


cutoff_nobase <- data.frame(subset(flow_cc_only,wy==2000), cutoff="Median\n1925-2005")


p <- ggplot(flow_cc_only, aes(x=month))
p <- p + geom_boxplot(aes( y=flow_m3s, colour=data), fill=NA, size=0.4)
p <- p + geom_line(aes( x=month, y=thresh_m3s, linetype = cutoff, group=wy), cutoff_nobase, colour="black")
#p <- p + facet_wrap( ~ data, ncol=3)
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(name="Scenario", values= cc_colors, labels=c("Base", "HD", "HW", "WD", "WW"))
p <- p + scale_linetype_manual(name="Reference", values=c("longdash"))
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow (m3/s)", breaks=seq(0,100,10)) + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_box_m3s_nobase.png"),  p, width=8, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_box_m3s_nobase.pdf"),  p, width=8, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_box_m3s_nobase.svg"),  p, width=8, height=6)


p <- ggplot(flow_cc_only, aes(x=month))
p <- p + geom_boxplot(aes( y=flow_m3s*m3s_to_acftmonth/1000, colour=data), fill=NA, size=0.4)
p <- p + geom_line(aes( x=month, y=thresh_m3s*m3s_to_acftmonth/1000, linetype = cutoff, group=wy), cutoff_nobase, colour="black")
#p <- p + facet_wrap( ~ data, ncol=3)
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(name="Scenario", values= cc_colors, labels=c("Base", "HD", "HW", "WD", "WW"))
p <- p + scale_linetype_manual(name="Reference", values=c("longdash"))
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow (1,000 ac-ft/month)", breaks=seq(0,160,20)) + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_flow_box_acft_nobase.png"),  p, width=8, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_flow_box_acft_nobase.pdf"),  p, width=8, height=6)
ggsave(file.path(write_output_base_path,"clim_change_flow_box_acft_nobase.svg"),  p, width=8, height=6)





###########################################################################
###  Plot Climate Change Normalized boxplots
###########################################################################
cutoff <- data.frame(subset(flow_cc,wy==2000), cutoff="Median\n1925-2005")
cutoff$norm <- 0

p <- ggplot(flow_cc, aes(x=month))
p <- p +geom_boxplot(aes( y=norm, colour=data), fill=NA, size=0.4)
p <- p + geom_line(aes( x=month, y=norm, linetype = cutoff, group=wy), cutoff, colour="black")
#p <- p + geom_hline(yintercept=0, colour="black", linetype="longdash")
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(name="Scenario", values= c("grey60", cc_colors), labels=c("Base", "HD", "HW", "WD", "WW"))
p <- p + scale_linetype_manual(name="Reference", values=c("longdash"))
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow Change (Standard Dev)", breaks=seq(-10,10,1)) + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_norm_box.png"),  p, width=9, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_norm_box.pdf"),  p, width=9, height=6)
ggsave(file.path(write_output_base_path,"clim_change_norm_box.svg"),  p, width=9, height=6)



cutoff_nobase <- data.frame(subset(flow_cc_only,wy==2000), cutoff="Median\n1925-2005")


p <- ggplot(flow_cc_only, aes(x=month))
p <- p + geom_boxplot(aes( y=flow_m3s, colour=data), fill=NA, size=0.4)
p <- p + geom_line(aes( x=month, y=thresh_m3s, linetype = cutoff, group=wy), cutoff_nobase, colour="black")
#p <- p + facet_wrap( ~ data, ncol=3)
p <- p + facet_grid(precip ~ temp)
p <- p + scale_colour_manual(name="Scenario", values= cc_colors, labels=c("Base", "HD", "HW", "WD", "WW"))
p <- p + scale_linetype_manual(name="Reference", values=c("longdash"))
p <- p + theme_bw(10)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(legend.position="none")
p <- p + scale_y_continuous(name="Flow (m3/s)", breaks=seq(0,100,10)) + scale_x_discrete(name="Month")
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_norm_box_m3s_nobase.png"),  p, width=8, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_norm_box_m3s_nobase.pdf"),  p, width=8, height=6)
ggsave(file.path(write_output_base_path,"clim_change_norm_box_m3s_nobase.svg"),  p, width=8, height=6)




