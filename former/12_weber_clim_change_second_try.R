
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

weber_data_path <- file.path(weber_output_path,"percent_ts")
write_output_base_path <- file.path(weber_output_path, output_name)

dir.create(write_output_base_path)


###########################################################################
###  Load functions
###########################################################################
### Load these functions for all code
require(colorout)
require(assertthat)

### Load these functions for this unique project
#require(ggplot2)
require(tidyverse)
require(lubridate)
require(scales)
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
monthly_distr <- "gamma"

###########################################################################
###  Read in data
###########################################################################

### Read to csv
read_location <- file.path(weber_data_path, paste0(site_id,"_obs_perc_ts.csv"))
weber_obs <- read.csv(read_location)
weber_obs$date <- as.Date(weber_obs$date)

### Read to csv
read_location <- file.path(weber_data_path, paste0(site_id,"_climchange_perc_ts.csv"))
weber_climchange <- read.csv(read_location)
weber_climchange$date <- as.Date(weber_climchange$date)


################################################
### Read in monthly and annual parameters to Transform Percentiles
#################################################
percentile_path <- file.path(file.path(weber_output_path, "paleo_reconst"), "ap_model")


paleo_monthly_weber/output/paleo_reconst/norm_fit/10128500/


percentile_path <- file.path(write_output_base_path, site_id)

monthly_param <- list(param=read.csv(file.path(write_output_base_path, paste0(site_id,"/",site_id,"_param_month_",monthly_distr,".csv"))), distr=monthly_distr)
monthly_param$param <- monthly_param$param[,c(2,3)]


###########################################################################
###  Prepare data for plotting
###########################################################################
weber_climchange <- tbl_df(weber_climchange)
weber_climchange

### Create wide table for monthly mean 
weber_cc_flow <- weber_climchange[,c(5, 6, 8)]
weber_cc_flow <- spread(weber_cc_flow, data, monthly_mean)

### Separate scenarios
weber_cc_scenarios <- select(weber_cc_flow, c(-observed, -CTN5))
weber_cc_obs <- select(weber_cc_flow, c(date, observed))

### Calculate flow difference
weber_cc_diff <- select(weber_cc_scenarios, -date) - data.frame(rep(select(weber_cc_obs, -date), 4))
weber_cc_diff <- data.frame(date=weber_cc_scenarios$date, weber_cc_diff) %>% tbl_df
weber_cc_diff$month <- month(weber_cc_diff$date)
weber_cc_diff$month_name <- month(weber_cc_diff$date, label=TRUE, abbr=TRUE)
weber_cc_diff$year <- year(weber_cc_diff$date)

weber_cc_diff <- gather(weber_cc_diff, model, flow, 2:5)
weber_cc_diff

### Re-label models
weber_cc_diff$model <- factor(weber_cc_diff$model, levels=c("WWN5", "HWN5", "WDN5", "HDN5"), labels=c("WW", "HW", "WD", "HD"))

###########################################################################
###  Custom theme
###########################################################################
theme_classic_bottom <- theme_classic_correct()+ theme(legend.direction = "horizontal", legend.position="bottom", axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))

###########################################################################
###  Plot Distribution of Flows ac-ft
###########################################################################
p <- ggplot(weber_cc_diff, aes(x=flow * 2.13056, color=model))
#p <- p + geom_density(adjust=2)
p <- p + geom_vline(xintercept = 0, colour="grey50", linetype="dashed")
p <- p + stat_density(geom="path", position="identity", adjust=2) 
p <- p + scale_colour_manual(name="Model", values=c("#1f78b4","#a6cee3", "#fb9a99", "#e31a1c" ))
p <- p + theme_classic_bottom
p <- p + guides(colour = guide_legend( legend.position="bottom", legend.direction = "horizontal", title.position = "top", title.hjust=0.5))
p <- p + facet_wrap( ~ month_name, scales="free")
p <- p + scale_x_continuous(name="Flow Change (1,000 ac-ft)", expand = c(0.02, 0.02))
p <- p + scale_y_continuous(name="Density")
p 

### Save figures
write_folder <- file.path(write_output_base_path, "distribution")
dir.create(write_folder)
file_name <- "clim_change_distribution_freescale_acft"
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=9, height=5.5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=9, height=5.5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=9, height=5.5)
	
### Make scales free
p <- p + facet_wrap( ~ month_name, scales="free_y")
	
### Save figures
file_name <- "clim_change_distribution_equalscale_acft"
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=9, height=5.5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=9, height=5.5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=9, height=5.5)
	
###########################################################################
###  Plot Distribution of Flows m3/s
###########################################################################
p <- ggplot(weber_cc_diff, aes(x=flow, color=model))
#p <- p + geom_density(adjust=2)
p <- p + geom_vline(xintercept = 0, colour="grey50", linetype="dashed")
p <- p + stat_density(geom="path", position="identity", adjust=2) 
p <- p + scale_colour_manual(name="Model", values=c("#1f78b4","#a6cee3", "#fb9a99", "#e31a1c" ))
p <- p + theme_classic_bottom
p <- p + guides(colour = guide_legend( legend.position="bottom", legend.direction = "horizontal", title.position = "top", title.hjust=0.5))
p <- p + facet_wrap( ~ month_name, scales="free")
p <- p + scale_x_continuous(name=expression( bold(Flow~Change~~(m^3/s) ) ) , expand = c(0.02, 0.02))
p <- p + scale_y_continuous(name="Density")
p 

### Save figures
file_name <- "clim_change_distribution_freescale_m3s"
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=9, height=5.5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=9, height=5.5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=9, height=5.5)
	
### Make scales free
p <- p + facet_wrap( ~ month_name, scales="free_y")
	
### Save figures
file_name <- "clim_change_distribution_equalscale_m3s"
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=9, height=5.5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=9, height=5.5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=9, height=5.5)
	
	
###########################################################################
###  Plot Box Plot acft
###########################################################################

# Basic box plot
p <- ggplot(weber_cc_diff, aes(x=model, y=flow  * 2.13056))
#p <- p + geom_boxplot()
p <- p + geom_hline(yintercept = 0, colour="grey50", linetype="dashed")
p <- p + geom_boxplot(outlier.colour="grey40", outlier.shape=18)
p <- p + stat_summary(fun.y=mean, geom="point", colour="red", shape=23)
p <- p + theme_classic_correct()
p <- p + facet_wrap( ~ month_name, scales="free_y")
p <- p + scale_y_continuous(name="Flow Change (1,000 ac-ft)")
p <- p + scale_x_discrete(name="Month")
p

### Save figures
write_folder <- file.path(write_output_base_path, "boxplot")
dir.create(write_folder)
file_name <- "clim_change_boxplot_freescale_acft"
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=8, height=5.5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=8, height=5.5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=8, height=5.5)
	
### Make scales free
p <- p + facet_wrap( ~ month_name)
	
### Save figures
file_name <- "clim_change_boxplot_equalscale_acft"
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=8, height=5.5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=8, height=5.5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=8, height=5.5)

	
###########################################################################
###  Plot Box Plot m3s
###########################################################################

# Basic box plot
p <- ggplot(weber_cc_diff, aes(x=model, y=flow))
#p <- p + geom_boxplot()
p <- p + geom_hline(yintercept = 0, colour="grey50", linetype="dashed")
p <- p + geom_boxplot(outlier.colour="grey40", outlier.shape=18)
p <- p + stat_summary(fun.y=mean, geom="point", colour="red", shape=23)
p <- p + theme_classic_correct()
p <- p + facet_wrap( ~ month_name, scales="free_y")
p <- p + scale_y_continuous(name=expression( bold(Flow~Change~~(m^3/s) ) ))
p <- p + scale_x_discrete(name="Month")
p

### Save figures
file_name <- "clim_change_boxplot_freescale_m3s"
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=8, height=5.5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=8, height=5.5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=8, height=5.5)
	
### Make scales free
p <- p + facet_wrap( ~ month_name)
	
### Save figures
file_name <- "clim_change_boxplot_equalscale_m3s"
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=8, height=5.5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=8, height=5.5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=8, height=5.5)


	
###########################################################################
###  Plot Month Diff acft
###########################################################################
### Calculate Monthly mean
weber_cc_mean <- weber_cc_diff %>%
       group_by(month_name, model) %>%
       summarise(mean = mean(flow), n = n())
 
# Plot Time Series

p <- ggplot(weber_cc_diff, aes(x=month_name))
p <- p + geom_line(colour="grey70", aes(y=flow * 2.13056, group=year))
p <- p + geom_hline(yintercept = 0, linetype="dashed")
p <- p + geom_line(data=weber_cc_mean, aes(y=mean * 2.13056, group=model), colour="red", size=1)
p <- p + facet_wrap( ~ model, scales="free_x")
p <- p + scale_x_discrete(name="Month", expand = c(0.02, 0.02))  
p <- p + scale_y_continuous(name="Flow Change (1,000 ac-ft)", labels = comma)       
p <- p + theme_classic_correct()
p

### Save figures
write_folder <- file.path(write_output_base_path, "month_ts_line")
dir.create(write_folder)
file_name <- "clim_change_month_ts_line_acft"
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=7, height=5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=7, height=5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=7, height=5)
	
p <- ggplot(weber_cc_diff, aes(x=month_name, y=flow * 2.13056))
p <- p + geom_hline(yintercept = 0, linetype="dashed", colour="grey50")
p <- p + geom_boxplot(outlier.colour="grey40", outlier.shape=18)
p <- p + stat_summary(fun.y=mean, geom="point", colour="red", shape=23)
p <- p + facet_wrap( ~ model, scales="free_x")
p <- p + theme_classic_correct()
p <- p + scale_x_discrete(name="Month", expand = c(0.02, 0.02))  
p <- p + scale_y_continuous(name="Flow Change (1,000 ac-ft)", labels = comma)    
p

### Save figures
write_folder <- file.path(write_output_base_path, "month_ts_box")
dir.create(write_folder)
file_name <- "clim_change_month_ts_box_acft"
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=7, height=5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=7, height=5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=7, height=5)

###########################################################################
###  Plot Month Diff m3s
###########################################################################
# Plot Time Series
p <- ggplot(weber_cc_diff, aes(x=month_name))
p <- p + geom_line(colour="grey70", aes(y=flow, group=year))
p <- p + geom_hline(yintercept = 0, linetype="dashed")
p <- p + geom_line(data=weber_cc_mean, aes(y=mean, group=model), colour="red", size=1)
p <- p + facet_wrap( ~ model, scales="free_x")
p <- p + scale_x_discrete(name="Month", expand = c(0.02, 0.02))  
p <- p + scale_y_continuous(name=expression( bold(Flow~Change~~(m^3/s) ) ))
p <- p + theme_classic_correct()
p

### Save figures
file_name <- "clim_change_month_ts_line_m3s"
write_folder <- file.path(write_output_base_path, "month_ts_line")
dir.create(write_folder)
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=7, height=5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=7, height=5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=7, height=5)
	
p <- ggplot(weber_cc_diff, aes(x=month_name, y=flow))
p <- p + geom_hline(yintercept = 0, linetype="dashed", colour="grey50")
p <- p + geom_boxplot(outlier.colour="grey40", outlier.shape=18)
p <- p + stat_summary(fun.y=mean, geom="point", colour="red", shape=23)
p <- p + facet_wrap( ~ model, scales="free_x")
p <- p + theme_classic_correct()
p <- p + scale_x_discrete(name="Month", expand = c(0.02, 0.02))  
p <- p + scale_y_continuous(name=expression( bold(Flow~Change~~(m^3/s) ) ))
p

### Save figures
write_folder <- file.path(write_output_base_path, "month_ts_box")
dir.create(write_folder)
file_name <- "clim_change_month_ts_box_m3s"
ggsave(file.path(write_folder, paste0(file_name,".png")),  p, width=7, height=5, dpi=600)
ggsave(file.path(write_folder, paste0(file_name,".pdf")),  p, width=7, height=5)
ggsave(file.path(write_folder, paste0(file_name,".svg")),  p, width=7, height=5)





################################################
### Calculate Thresholds
#################################################

threshold_month <- data.frame(month=seq(1,12), q_1=NA, q_5=NA, q_10=NA, q_25=NA, q_50=NA, q_75=NA)

for (j in 1:12) {

	threshold_month$q_1[j] <- qgamma(0.01, shape= monthly_param$param[j,1], rate=monthly_param$param[j,2])
	threshold_month$q_5[j] <- qgamma(0.05, shape= monthly_param$param[j,1], rate=monthly_param$param[j,2])
	threshold_month$q_10[j] <- qgamma(0.1, shape= monthly_param$param[j,1], rate=monthly_param$param[j,2])
	threshold_month$q_25[j] <- qgamma(0.25, shape= monthly_param$param[j,1], rate=monthly_param$param[j,2])
	threshold_month$q_50[j] <- qgamma(0.5, shape= monthly_param$param[j,1], rate=monthly_param$param[j,2])
	threshold_month$q_75[j] <- qgamma(0.75, shape= monthly_param$param[j,1], rate=monthly_param$param[j,2])
}


###########################################################################
###  Process percentiles
###########################################################################
### Create wide table for monthly percentile 
weber_cc_perc <- weber_climchange[,c(5, 7, 8)]
weber_cc_perc <- spread(weber_cc_perc, data, month_perc)

### Separate scenarios
weber_perc_scenarios <- select(weber_cc_perc, c(-observed, -CTN5))
weber_perc_obs <- select(weber_cc_perc, c(date, observed))


apply(select(weber_perc_scenarios, -date), 2, pnorm)

plot(qnorm(weber_perc_obs$observed), type="l")
lines(qnorm(weber_perc_scenarios$HDN5), col="red")


### Calculate flow difference
weber_cc_diff <- select(weber_cc_scenarios, -date) - data.frame(rep(select(weber_cc_obs, -date), 4))
weber_cc_diff <- data.frame(date=weber_cc_scenarios$date, weber_cc_diff) %>% tbl_df
weber_cc_diff$month <- month(weber_cc_diff$date)
weber_cc_diff$month_name <- month(weber_cc_diff$date, label=TRUE, abbr=TRUE)
weber_cc_diff$year <- year(weber_cc_diff$date)

weber_cc_diff <- gather(weber_cc_diff, model, flow, 2:5)
weber_cc_diff

### Re-label models
weber_cc_diff$model <- factor(weber_cc_diff$model, levels=c("WWN5", "HWN5", "WDN5", "HDN5"), labels=c("WW", "HW", "WD", "HD"))






weber_cc_flow$month <- month(weber_cc_flow$date)
weber_cc_flow$year <- year(weber_cc_flow$date)





