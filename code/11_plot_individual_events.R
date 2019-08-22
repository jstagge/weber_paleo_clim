# *------------------------------------------------------------------
# | PROGRAM NAME: 11_plot_individual_events
# | FILE NAME: 11_plot_individual_events.R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  This code 
# | 
# | 
# *------------------------------------------------------------------

### Load the here package to set up the working path
require(here) 

###########################################################################
## Set the Paths
###########################################################################
### Path for Data and Output	
data_path <- file.path(here(), "data")
output_path <- file.path(here(), "output")
global_path <- "../global_func"
function_path <- file.path(here(), "functions")

### Set output location
output_name <- "clim_change"

weber_output_path <- output_path

write_output_base_path <- file.path(weber_output_path, output_name)

dir.create(write_output_base_path, recursive = TRUE)

###########################################################################
###  Load functions
###########################################################################

### Load these functions for this unique project
require(ggplot2)
require(staggefuncs)
require(lubridate)
require(reshape2)
require(scales)
require(tidyverse)

###########################################################################
## Open Output File
###########################################################################
load(file.path(weber_output_path, "weber_clustering.RData"))

load(file.path(weber_output_path, "weber_storage_output.RData"))
load(file.path(weber_output_path, "weber_delivery_output.RData"))


###########################################################################
## Read in Flow data
###########################################################################
### Read in Observed Flow data
read_location <- file.path(write_output_base_path, paste0(site_id,"_obs_perc_ts.csv"))
weber_obs_ts <- read.csv(read_location)
weber_obs_ts$date <- as.Date(weber_obs_ts$date)

### Read in climate change Flow data
read_location <- file.path(write_output_base_path, paste0(site_id,"_climchange_perc_ts.csv"))
weber_cc_ts <- read.csv(read_location)
weber_cc_ts$date <- as.Date(weber_cc_ts$date)

### Read in paleo Flow data
read_location <- file.path(write_output_base_path, paste0(site_id,"_paleo_perc_ts.csv"))
weber_paleo_ts <- read.csv(read_location)
weber_paleo_ts$date <- as.Date(weber_paleo_ts$date)

###########################################################################
## Create Output Paths
###########################################################################
storage_output_path <- file.path(write_output_base_path, "storage")
dir.create(storage_output_path)


###########################################################################
## Set Colors
###########################################################################
### Colors to use for climate change scenarios
cc_colors <- c("#0072B2", "#56B4E9", "#9BAEBC", "#E69F00" , "#D55E00")
data_colors <- c("#7fc97f", "grey40", "grey80", cc_colors)

data_labels <- c( "Reconstructed", "Observed", "Base", "Climate Change (WW)", "Climate Change (HW)", "Climate Change (Median)", "Climate Change (WD)", "Climate Change (HD)")
data_labels_twolines <- c( "Reconstructed", "Observed", "Base", "Clim Change\n(Warm-Wet)", "Clim Change\n(Hot-Wet)", "Clim Change\n(Median)", "Clim Change\n(Warm-Dry)", "Clim Change\n(Hot-Dry)")

### Colors for triggers
trigger_colors <- c("#f03b20",  "#feb24c", "#ffeda0")

### Set regions
region_levels <- as.factor(c("total_res", "current_res", "upper_weber", "upper_ogden", "lower"))
region_names <- as.factor(c("Total System", "Current Reservoir", "Upper Weber", "Upper Ogden", "Lower Weber"))

region_colors <- cb_pal(pal="wong", 3, sort=FALSE)
region_colors <- region_colors[c(1,3,2)]

theme_ts <- theme_classic_new(10)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))


###########################################################################
## Extract only base scenario
###########################################################################
stor_base <- stor_all %>% 
	filter(response == "Base")

stor_base_percent <- stor_percent %>% 
	filter(response == "Base")

###########################################################################
## Calculate total storage
###########################################################################
### Calculate Current total active storage
current_active_stor <- sum(total_storage$Active[seq(1,8)], na.rm=TRUE)
total_active_stor <- current_active_stor

active_stor_vec <- c(current_active_stor, total_active_stor, sum(total_storage$Active[seq(1,5)], na.rm=TRUE), sum(total_storage$Active[seq(6,7)], na.rm=TRUE), sum(total_storage$Active[8], na.rm=TRUE))
names(active_stor_vec) <- region_levels



###########################################################################
###  Set date range
###########################################################################
### Set date range
comparison_dates <- as.Date(c("1999-06-01", "2007-01-01"))
event_dates <- as.Date(c("1633-06-01", "1641-01-01"))

### Calculate date shift
date_shift <- min(comparison_dates) - min(event_dates)
### Plot range is 6 months to either side
plot_range <- c(comparison_dates[1] %m+% months(-3), comparison_dates[2] %m+% months(3))  

###########################################################################
###  Plot System Storage
###########################################################################
### Separate paleo and observed
comparison_df <- stor_base %>%
	filter(response == "Base" & data=="base") %>%
	filter(date >= (comparison_dates[1] %m+% months(-12)) & date <= (comparison_dates[2] %m+% months(12))) %>%
	arrange(date)

event_df <- stor_base %>%
	filter(response == "Base" & data=="paleo") %>%
	filter(date >= (event_dates[1] %m+% months(-12)) & date <= (event_dates[2] %m+% months(12))) %>%
	arrange(date)

comparison_df$plot_date <- comparison_df$date	
event_df$plot_date <- event_df$date  + date_shift


plot_df <- rbind(comparison_df, event_df)

p <- ggplot(plot_df, aes(x=plot_date, y=current_res/1000, color=data)) 
p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, linetype="dotted")
p <- p+ geom_line()
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1400-01-01"), as.Date("2100-01-01"), by="1 years"), labels = date_format("%Y"))#  -%m sec_axis(name="Date", trans=~  + date_shift))
p <- p + scale_y_continuous(name="Total Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="Total Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 
p <- p + scale_colour_manual(name="Scenario", limits=c("base", "paleo"), labels=c("1999-2007", "1633-1641"), values=c("red", "blue"))
p <- p + theme_classic_new() + theme(legend.position="bottom")
p <- p + expand_limits(y=c(-10, current_active_stor/1000+10))
p <- p + coord_cartesian(xlim=plot_range)
p

### Save figures
ggsave("storage_timeseries_1639.png",  p, width=7.5, height=3.5, dpi=600)


breaks_qtr <- seq(as.Date("1400-1-1"), by = "6 months", length.out = 6000)
labels_year = format(seq(from = min(breaks_qtr), to = max(breaks_qtr), by = "1 years"), "%Y")
labs = c(sapply(labels_year, function(x) {
    c(x, rep("", 1))
    }))    

p <- ggplot(plot_df, aes(x=plot_date, y=current_res/1000, color=data, size=data, alpha=data)) #, 
#p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_rect(aes(xmin=as.Date("2015-04-01"), xmax=as.Date("2015-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
#p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
#p <- p + geom_vline(xintercept=as.Date("2015-10-01"), size=0.6, colour="grey50")
p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, linetype="longdash", size=0.1, alpha=0.9)
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("#FF7F0EFF", "grey20"), name="") ##FF7F0EFF #377eb8
p <- p + scale_size_manual(values=c(0.7, 0.35), name="")
p <- p + scale_alpha_manual(values=c(1, 0.5), name="")
#p <- p  + theme_classic_new(14)
p <- p + theme_ts
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1400-01-01"), as.Date("2100-01-01"), by="1 years"), labels = date_format("%Y"))#  -%m 
p <- p + scale_y_continuous(name="Total Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="Total Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 
p <- p + coord_cartesian(ylim=c(0,550), xlim=plot_range)
p <- p  + theme(legend.position="bottom") #+ theme(legend.position = c(0.90, 0.9))
p
#

### Save figures
ggsave("storage_timeseries_1639_v2.png",  p,  width=7.5, height=3.5, dpi=600)


storage_plot_df <- plot_df %>% 
	select(date, month, year, wy, plot_date, current_res, data, response, trigger_category) 

write.csv(storage_plot_df, file="storage_1639.csv")


###########################################################################
###  Plot System delivery / demand /  shortage
###########################################################################
node_to_plot <- "system"

comparison_df <- demand_deliv_df %>%
	filter(response == "Base" & data=="base" & node==node_to_plot) %>%
	filter(date >= (comparison_dates[1] %m+% months(-12)) & date <= (comparison_dates[2] %m+% months(12))) %>%
	arrange(date)

event_df <- demand_deliv_df %>%
	filter(response == "Base" & data=="paleo" & node==node_to_plot) %>%
	filter(date >= (event_dates[1] %m+% months(-12)) & date <= (event_dates[2] %m+% months(12))) %>%
	arrange(date)

comparison_df$plot_date <- comparison_df$date	
event_df$plot_date <- event_df$date  + date_shift

plot_df <- rbind(comparison_df, event_df)

p <- ggplot(plot_df, aes(x=plot_date, color=data))
p <- p + geom_line(data = plot_df, aes(y=demand/1000), linetype="dashed", color="black")
p <- p+ geom_line(aes(y=delivery/1000 ))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1400-01-01"), as.Date("2100-01-01"), by="1 years"), labels = date_format("%Y"))#  -%m sec_axis(name="Date", trans=~  + date_shift))
p <- p + scale_y_continuous(name="Total System Delivery (1,000 ac-ft)")
p <- p + scale_colour_manual(name="Scenario", limits=c("base", "paleo"), labels=c("1999-2007", "1633-1641"), values=c("red", "blue"))
p <- p + theme_classic_new() + theme(legend.position="bottom")
#p <- p + expand_limits(y=c(-10, current_active_stor/1000+10))
p


### Save figures
ggsave("delivery_timeseries_1639.png",  p, width=7.5, height=3.5, dpi=600)


p <- ggplot(plot_df, aes(x=plot_date, y=demand_shortage/1000, color=data))
p <- p+ geom_line()
#p <- p + scale_y_continuous(name="Total Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="Total Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 
p <- p + scale_y_continuous(name="Total System Shortage (1,000 ac-ft)")
p <- p + scale_colour_manual(name="Scenario", limits=c("base", "paleo"), labels=c("1999-2007", "1633-1641"), values=c("red", "blue"))
p <- p + theme_classic_new() + theme(legend.position="bottom")
#p <- p + expand_limits(y=c(-10, current_active_stor/1000+10))
p

### Save figures
ggsave("shortage_timeseries_1639.png",  p, width=7.5, height=3.5, dpi=600)



p <- ggplot(plot_df, aes(x=plot_date, y=demand_shortage/1000, color=data, size=data, alpha=data)) #, 
#p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_rect(aes(xmin=as.Date("2015-04-01"), xmax=as.Date("2015-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
#p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
#p <- p + geom_vline(xintercept=as.Date("2015-10-01"), size=0.6, colour="grey50")
#p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, linetype="longdash", size=0.1, alpha=0.9)
p <- p + geom_line()
p <- p + scale_colour_manual(values=c("#FF7F0EFF", "grey20"), name="Scenario", limits=c("paleo", "base"), labels=c("1633-1641", "1999-2007")) ##FF7F0EFF #377eb8
p <- p + scale_size_manual(values=c(0.7, 0.35), name="Scenario", limits=c("paleo", "base"), labels=c("1633-1641", "1999-2007")) 
p <- p + scale_alpha_manual(values=c(1, 0.5), name="Scenario", limits=c("paleo", "base"), labels=c("1633-1641", "1999-2007")) 
#p <- p  + theme_classic_new(14)
p <- p + theme_ts
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1400-01-01"), as.Date("2100-01-01"), by="1 years"), labels = date_format("%Y"))#  -%m 
p <- p + scale_y_continuous(name="System Demand Shortage (1,000 ac-ft)", breaks=seq(0,800,20), expand = c(0,0)) 
p <- p + coord_cartesian(ylim=c(0,80), xlim=plot_range)
p <- p +  theme(legend.position="bottom")# theme(legend.position = c(0.90, 0.9))
p
#

### Save figures
ggsave("shortage_timeseries_1639_v2.png",  p, width=7.5, height=3.5, dpi=600)


### Wanted a rolling cumulative per year, but makes it confusing based on non-water months
plot_df$plot_year <- year(plot_df$plot_date)

plot_annual <- plot_df %>%
	group_by(data, response, plot_year) %>%
	summarize(deliv_annual=sum(delivery), demand_annual=sum(demand, na.rm=TRUE), request_annual=sum(request, na.rm=TRUE), shortage_demand=sum(demand_shortage), shortage_request=sum(request_shortage, na.rm=TRUE)) %>%
	arrange(data, response, plot_year)

plot_annual$demand_shortage <- plot_annual$demand_annual - plot_annual$deliv_annual

plot_annual <- plot_annual %>% filter(plot_year > 1998 & plot_year < 2008)

p <- ggplot(plot_annual, aes(x=plot_year, y=demand_shortage/1000, fill=data)) #, 
#p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_rect(aes(xmin=as.Date("2015-04-01"), xmax=as.Date("2015-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
#p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
#p <- p + geom_vline(xintercept=as.Date("2015-10-01"), size=0.6, colour="grey50")
#p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, linetype="longdash", size=0.1, alpha=0.9)
p <- p + geom_bar(stat="identity", position="dodge")
p <- p + scale_fill_manual(values=c("#FF7F0EFF", "grey70"), name="Scenario", limits=c("paleo", "base"), labels=c("1633-1641", "1999-2007")) ##FF7F0EFF #377eb8
p <- p  + theme_classic_new(10)
#p <- p + theme_ts
p <- p + scale_x_continuous(breaks=seq(1900,2018), name = "Year") 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1400-01-01"), as.Date("2100-01-01"), by="1 years"), labels = date_format("%Y"))#  -%m 
p <- p + scale_y_continuous(name="System Demand Shortage (1,000 ac-ft)", breaks=seq(0,800,20), expand = c(0,0)) 
#p <- p + coord_cartesian(ylim=c(0,80), xlim=plot_range)
p <- p +  theme(legend.position="bottom")# theme(legend.position = c(0.90, 0.9))
p

### Save figures
ggsave("shortage_timeseries_1639_bar.png",  p, width=7.5, height=3.5, dpi=600)





delivery_plot_df <- plot_df %>% 
	select(date, month, year, wy, plot_date, delivery, demand, demand_shortage, data, response) 

write.csv(storage_plot_df, file="delivery_1639_monthly.csv")
write.csv(plot_annual, file="delivery_1639_annual.csv")



###########################################################################
###  Plot Flow
###########################################################################
#weber_obs_ts
#weber_cc_ts
#weber_paleo_ts

### Separate paleo and observed
comparison_df <- weber_obs_ts %>%
	#filter(response == "Base" & data=="base") %>%
	filter(date >= (comparison_dates[1] %m+% months(-12)) & date <= (comparison_dates[2] %m+% months(12))) %>%
	arrange(date)

event_df <- weber_paleo_ts %>%
	#filter(response == "Base" & data=="paleo") %>%
	filter(date >= (event_dates[1] %m+% months(-12)) & date <= (event_dates[2] %m+% months(12))) %>%
	arrange(date)

comparison_df$plot_date <- comparison_df$date	
event_df$plot_date <- event_df$date  + date_shift


plot_df <- rbind(comparison_df, event_df)

thresh_df <- data.frame(plot_date=event_df$plot_date, flow_m3s=event_df$thresh_m3s, norm=0, label="50th Percentile")

p <- ggplot(plot_df, aes(x=plot_date, y=flow_m3s*2131.96802/1000)) #, 
#p <- p + geom_line(data=thresh_df, linetype="longdash", size=0.1, alpha=0.9)
p <- p + geom_line(aes(color=data, size=data, alpha=data))
p <- p + scale_colour_manual(values=c("#FF7F0EFF", "grey20"), name="Scenario", limits=c("paleo", "observed"), labels=c("1633-1641", "1999-2007")) ##FF7F0EFF #377eb8
p <- p + scale_size_manual(values=c(0.7, 0.35), name="Scenario", limits=c("paleo", "observed"), labels=c("1633-1641", "1999-2007")) 
p <- p + scale_alpha_manual(values=c(1, 0.5), name="Scenario", limits=c("paleo", "observed"), labels=c("1633-1641", "1999-2007")) 
#p <- p  + theme_classic_new(14)
p <- p + theme_ts
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1400-01-01"), as.Date("2100-01-01"), by="1 years"), labels = date_format("%Y"))#  -%m 
p <- p + scale_y_continuous(name="Flow (1,000 ac-ft)", breaks=seq(0,800,20), expand = c(0,0)) 
p <- p + coord_cartesian(ylim=c(0,100), xlim=plot_range)
p <- p +  theme(legend.position="bottom")# theme(legend.position = c(0.90, 0.9))
p
#
### Save figures
ggsave("flow_timeseries_1639.png",  p, width=7.5, height=3.5, dpi=600)


p <- ggplot(weber_obs_ts, aes(x=date, y=pnorm(norm)))
#p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=-Inf, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_rect(aes(xmin=as.Date("2015-04-01"), xmax=as.Date("2015-10-01"), ymin=-Inf, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
#p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
#p <- p + geom_vline(xintercept=as.Date("2015-10-01"), size=0.6, colour="grey50")
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="#d95f02")
p <- p + geom_line(colour="black")#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
#p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2011-01-01"), as.Date("2017-01-01")))
p
#
### Save figures
ggsave("flow_percentile_timeseries_1639.png",  p, width=7.5, height=3.5, dpi=600)

p <- ggplot(plot_df, aes(x=plot_date, y=pnorm(norm))) #, 
#p <- p + geom_rect(aes(xmin=as.Date("2013-04-01"), xmax=as.Date("2013-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_rect(aes(xmin=as.Date("2014-04-01"), xmax=as.Date("2014-10-01"), ymin=-Inf, ymax=Inf), fill="grey93", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_rect(aes(xmin=as.Date("2015-04-01"), xmax=as.Date("2015-10-01"), ymin=0, ymax=Inf), fill="grey90", colour=NA, alpha=0.1, show.legend = NA)
#p <- p + geom_vline(xintercept=as.Date("2013-10-01"), size=0.6, colour="grey50")
#p <- p + geom_vline(xintercept=as.Date("2014-10-01"), size=0.6, colour="grey50")
#p <- p + geom_vline(xintercept=as.Date("2015-10-01"), size=0.6, colour="grey50")
#p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, linetype="longdash", size=0.1, alpha=0.9)
p <- p + geom_line(data=thresh_df, linetype="longdash", size=0.1, alpha=0.9)
p <- p + geom_line(aes(color=data, size=data, alpha=data))
p <- p + scale_colour_manual(values=c("#FF7F0EFF", "grey20"), name="Scenario", limits=c("paleo", "observed"), labels=c("1633-1641", "1999-2007")) ##FF7F0EFF #377eb8
p <- p + scale_size_manual(values=c(0.7, 0.35), name="Scenario", limits=c("paleo", "observed"), labels=c("1633-1641", "1999-2007")) 
p <- p + scale_alpha_manual(values=c(1, 0.5), name="Scenario", limits=c("paleo", "observed"), labels=c("1633-1641", "1999-2007")) 
#p <- p  + theme_classic_new(14)
p <- p + theme_ts
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1400-01-01"), as.Date("2100-01-01"), by="1 years"), labels = date_format("%Y"))#  -%m 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=plot_range)
p <- p +  theme(legend.position="bottom")# theme(legend.position = c(0.90, 0.9))
p
#
### Save figures
ggsave("flow_percentile_timeseries_1639_v2.png",  p, width=7.5, height=3.5, dpi=600)




flow_plot_df <- plot_df %>% 
	mutate(flow_acft = flow_m3s*2131.96802, percentile = pnorm(norm)) %>%
	select(date, month, year, plot_date, flow_acft, percentile, data)
	

write.csv(storage_plot_df, file="flow_1639.csv")


###########################################################################
###  Extract annual values for Chris
###########################################################################


response == "Base" & data=="paleo" & 
	filter(data=="paleo" | )
	
yup <- demand_deliv_df %>%
	filter(node=="system" & response == "Base") %>%
	group_by(data, response, year) %>%
	summarize(deliv_annual=sum(delivery), demand_annual=sum(demand, na.rm=TRUE), request_annual=sum(request, na.rm=TRUE), shortage_demand=sum(demand_shortage), shortage_request=sum(request_shortage, na.rm=TRUE)) %>%
	arrange(response, data, year)
	
#	ggplot(yup, aes(x=year, y=deliv_annual, group=data, fill=response)) + geom_bar(stat="identity")
& response == "Base"

ggplot(yup, aes(x=year, group=data, colour=data)) + geom_line(aes(y=deliv_annual)) + geom_line(aes(y=demand_annual), col="black")

ggplot(yup, aes(x=data, y= shortage_demand, fill=data)) + geom_boxplot() 

write.csv(yup, "deliver_shortages_by_year.csv")



yup_delete <- stor_base[1,]
yup_delete$data <- "delete"
yup_delete$total_res <- NA
yup <- rbind(stor_base, yup_delete)
yup$data <- factor(yup$data, levels=c(data_levels[seq(1,3)], "delete", data_levels[seq(4,7)]), labels=c(data_labels[seq(1,3)], "", data_labels[seq(4,7)]))

yup <- stor_base
yup$data <- factor(yup$data, levels=data_levels[c(1,2,3,6,4,5,7,8)], labels=data_labels[c(1,2,3,6,4,5,7,8)])

p <- ggplot(yup, aes(x=month))
p <- p + geom_line(aes(y=total_res/1000, group=wy, colour=data), size=0.3)
p <- p + scale_x_discrete(name="Month")
p <- p + scale_y_continuous(name="Total Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="Total Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 
p <- p + facet_wrap(~data, nrow=2, scales="free_x")
p <- p + scale_colour_manual(name="Scenario", values=c(data_colors[c(1,2,3,6,4,5,7,8)]), guide=guide_legend(byrow=TRUE))
p <- p + theme_classic_new() + theme(legend.position="bottom")
p <- p + expand_limits(y=c(-10, current_active_stor/1000+10))
p
















# *------------------------------------------------------------------
# | PROGRAM NAME: plot_ap_apr_model
# | FILE NAME: 05_plot_ap_apr_model.R
# | DATE: 09/15/2017
# | CREATED BY:  James Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  This script plots the model fits, including all coefficients for 
# | 			the Annual Percentile (AP) and Annual Percentile with Regression (APR) 
# |				Models for the Bear and Logan Rivers.
# |				The method is detailed in Stagge et al. (2017), which is included as a reference.     
# |				 
# *------------------------------------------------------------------
# | DATA USED:               
# | Reconstruction models calculated in 04_ap_apr_model.R
# |
# *------------------------------------------------------------------
# | COMMENTS:               
# |
# |  1:  
# |  2:  
# |  3: 
# *------------------------------------------------------------------

### Clear any existing data or functions.
rm(list=ls())

###########################################################################
## Set the Paths
###########################################################################


###########################################################################
###  Load functions
###########################################################################
### Load these functions for all code
#require(colorout)
require(assertthat)
require(testthat)

### Load these functions for this unique project
#require(monthlypaleo)
require(staggefuncs)
require(paleoAPR)
require(ggplot2)
require(ggplot2)
require(svglite)
require(reshape2)

#################################################
### Read 
#################################################
plot_ts <- read.csv("10109001_ap_rec_region_ts.csv")

plot_ts$date <- as.Date(paste0(plot_ts$year,"-",plot_ts$month,"-15"))


#################################################
### Plot annual norm
#################################################


### Create plot with observed in black and predicted in red
p <- ggplot(plot_ts, aes(x=date, y=norm_est))
p <- p + geom_hline(yintercept = 0, size = 0.3, linetype="longdash")
p <- p + geom_line(size=0.6, colour="#386cb0")
p <- p + theme_classic_new(14)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("200-01-01"), as.Date("2030-01-01"), by="1 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Standard Normal", breaks=seq(-5,5,by=0.5))
p_break <- p + coord_cartesian(xlim=c(as.Date("1960-01-01"), as.Date("1974-06-01")), ylim=c(-2.5,2.5))

p_break
ggsave("ap_norm_plot.png", p_break, width=10, height=4, dpi=1200)




#### Run the code in github folder
n <- 8

p <- plot_result + theme_classic_new(14)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("200-01-01"), as.Date("2030-01-01"), by="1 years"), date_labels = "%Y")
#p <- p + scale_y_continuous(name="Standard Normal", breaks=seq(-5,5,by=0.5))
p_break <- p + coord_cartesian(xlim=c(as.Date("1960-01-01"), as.Date("1974-06-01")))

ggsave("ap_flow_plot_legend.png", p_break, width=10, height=4, dpi=1200)

p_break <- p_break + theme(legend.position="none")

p_break
ggsave("ap_flow_plot.png", p_break, width=10, height=4, dpi=1200)
















#################################################
### Read 
#################################################
plot_rec_ts <- read.csv("10109001_apr_rec_region_clim_pca_impute_postproc_ts.csv")

plot_rec_ts$date <- as.Date(paste0(plot_ts$year,"-",plot_ts$month,"-15"))


#################################################
### Plot annual norm
#################################################


### Create plot with observed in black and predicted in red
p <- ggplot(plot_ts, aes(x=date, y=norm_est))
p <- p + geom_hline(yintercept = 0, size = 0.3, linetype="longdash")
p <- p + geom_line(data=plot_ts, size=0.6, colour="#386cb0")
p <- p + geom_line(data=plot_rec_ts, size=0.6, colour="#d95f02")
p <- p + theme_classic_new(14)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("200-01-01"), as.Date("2030-01-01"), by="1 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Standard Normal", breaks=seq(-5,5,by=0.5))
p_break <- p + coord_cartesian(xlim=c(as.Date("1992-01-01"), as.Date("2005-06-01")), ylim=c(-2.5,2.5))

p_break
ggsave("apr_norm_plot.png", p_break, width=10, height=4, dpi=1200)




#### Run the code in github folder
n <- 20

p <- plot_result + theme_classic_new(14)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("200-01-01"), as.Date("2030-01-01"), by="1 years"), date_labels = "%Y")
#p <- p + scale_y_continuous(name="Standard Normal", breaks=seq(-5,5,by=0.5))
p_break <- p + coord_cartesian(xlim=c(as.Date("1992-01-01"), as.Date("2005-06-01")))

ggsave("apr_flow_plot_legend.png", p_break, width=10, height=4, dpi=1200)

p_break <- p_break + theme(legend.position="none")

p_break
ggsave("apr_flow_plot.png", p_break, width=10, height=4, dpi=1200)


















### Create plot with observed in black and predicted in red
p <- ggplot(plot_ts, aes(x=date, y=flow_est))
p <- p + geom_line(size=0.25)
p <- p + theme_classic_new()
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("200-01-01"), as.Date("2030-01-01"), by="1 years"), date_labels = "%Y")
p_break <- p + coord_cartesian(xlim=c(as.Date("1992-01-01"), as.Date("2005-06-01")))

p_break

ggsave("test_plot.png", p_break, width=10, height=4, dpi=600)




p <- p + scale_colour_manual(name= NULL, values = c("black", "blue", "red"))
p <- p + scale_linetype_manual(values=c("solid", "longdash", "solid"), guide=FALSE)
#p <- p + scale_y_continuous(name="Discharge (m3/s)")
p <- p + scale_y_continuous(name=expression(bold(paste("Monthly Mean Discharge  ( ",m^3,"/s )"))))
p <- p +  theme(legend.position="bottom")


### Save full plot
ggsave(paste0(file.path(write_folder,"png/"), write_file, "_full.png"), p, width=8, height=4, dpi=600)
ggsave(paste0(file.path(write_folder,"pdf/"), write_file, "_full.pdf"), p, width=8, height=4)
ggsave(paste0(file.path(write_folder,"svg/"), write_file, "_full.svg"), p, width=8, height=4)


### Loop through 15 year periods and save
plot_breaks <- seq(as.Date("1905-01-01"), as.Date("2030-01-01"), by="15 years")

for (k in seq(1,length(plot_breaks)-1)) {
### Identify the start of each break
start_break <- plot_breaks[k]
end_break <- plot_breaks[k+1]

### Cut to break point and reorganize x axis to 2 years
p_break <- p + coord_cartesian(xlim=c(as.Date(start_break), as.Date(end_break)))
suppressMessages(p_break <- p_break + scale_x_date(name="Date", breaks=seq(as.Date("1800-01-01"), as.Date("2030-01-01"), by="2 years"), date_labels = "%Y"))

### Save zoomed plots
write_file_subset <- paste0(write_file, "_",start_break,"_",end_break)
ggsave(paste0(file.path(write_folder,"png/"), write_file_subset, "_",start_break,"_",end_break,".png"), p_break, width=6, height=4, dpi=600)
ggsave(paste0(file.path(write_folder,"pdf/"), write_file_subset, "_",start_break,"_",end_break,".pdf"), p_break, width=6, height=4)
ggsave(paste0(file.path(write_folder,"svg/"), write_file_subset, "_",start_break,"_",end_break,".svg"), p_break, width=6, height=4)

}





### Melt to get into proper format
plot_data <- melt(data,  id.vars="date", measure.vars=c("flow_est", "flow_obs", "flow_annual"))

### Fix column names and convert method to characters
colnames(plot_data) <- c("Date", "Measure", "Flow")
plot_data$Measure <- as.character(plot_data$Measure)

### Factor to make sure in correct order
#plot_data$Measure <- factor(plot_data$Measure, levels = c( recon_name, obs_name, annual_name))
plot_data$Measure <- factor(plot_data$Measure, levels = c("flow_obs", "flow_annual", "flow_est"), labels = c("Observed      ", "Annual Reconstruction", "Monthly Reconstruction"))



