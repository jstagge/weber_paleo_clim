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

weber_stor_path <- file.path(data_path,"weber_storage")
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
require(lubridate)
require(reshape2)
require(scales)
require(tidyverse)

###########################################################################
## Open Output File
###########################################################################
load(file.path(weber_output_path, "weber_clustering.RData"))

load(file.path(weber_output_path, "weber_storage_output.RData"))

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

###########################################################################
## Plot Triggers
###########################################################################

plot_triggers <- weber_triggers
plot_triggers$total <-  sum(total_storage$Total)
names(plot_triggers)[2] <- "Mean 2013-2017"
plot_triggers <- melt(plot_triggers[,c(seq(1,5), 9)], "Month")
names(plot_triggers)[3] <- "storage"
plot_triggers$perc <- plot_triggers$storage / sum(total_storage$Total)

p <- ggplot(plot_triggers, aes(x=Month, y=storage/1000, colour=variable, group=variable))
#p <- ggplot(yup, aes(x=date, y=value/1000, fill=variable))
p <- p + geom_line(size=0.8)
#p <- p + geom_line(data=trigger_plot, aes(y=res_stor/1000, group=trigger_level, fill=NA), colour="grey30", linetype="longdash", size=0.8)
p <- p + theme_classic_new(11)
p <- p + scale_colour_manual(name="", values= c( "#ffeda0", "#feb24c", "#f03b20", "#8da0cb", "grey50"), limits= c( "moderate", "severe", "extreme", "Mean 2013-2017", "total"), labels=c("Moderate Trigger", "Severe Trigger", "Extreme Trigger", "Mean 2013-2017", "Full Storage"), guide = guide_legend(nrow=2,byrow=TRUE))
p <- p + scale_x_discrete(name="Month")
p <- p + scale_y_continuous(name="Total System Storage (1,000 ac-ft)", breaks=seq(0, 600, 100))
p <- p + coord_cartesian(xlim=c(1,12), ylim=c(0,sum(total_storage$Total/1000)*1.1), expand=FALSE)
p <- p + theme(legend.position="bottom")
p

### Save figures
ggsave(file.path(storage_output_path,"trigger_levels_perc.png"),  p, width=4.5, height=4, dpi=600)
ggsave(file.path(storage_output_path,"trigger_levels_perc.pdf"),  p, width=4.5, height=4)
ggsave(file.path(storage_output_path,"trigger_levels_perc.svg"),  p, width=4.5, height=4)


p <- ggplot(plot_triggers, aes(x=Month, y=perc, colour=variable, group=variable))
#p <- ggplot(yup, aes(x=date, y=value/1000, fill=variable))
p <- p + geom_line(size=0.8)
#p <- p + geom_line(data=trigger_plot, aes(y=res_stor/1000, group=trigger_level, fill=NA), colour="grey30", linetype="longdash", size=0.8)
p <- p + theme_classic_new(11)
p <- p + scale_colour_manual(name="", values= c( "#ffeda0", "#feb24c", "#f03b20", "#8da0cb", "grey50"), limits= c( "moderate", "severe", "extreme", "Mean 2013-2017", "total"), labels=c("Moderate Trigger", "Severe Trigger", "Extreme Trigger", "Mean 2013-2017", "Full Storage"), guide = guide_legend(nrow=2,byrow=TRUE))
p <- p + scale_x_discrete(name="Month")
#p <- p + scale_y_continuous(name="Total System Storage (%)", breaks=seq(0, 1, 0.1), labels=percent)
p <- p + scale_y_continuous(name="Total System Storage (%)", breaks= seq(0, 1, 0.1), labels = c("0", "", "20%", "", "40%", "", "60%", "", "80%", "", "100%"))
p <- p + coord_cartesian(xlim=c(1,12), ylim=c(0,1.02), expand=FALSE)
p <- p + theme(legend.position="bottom")
p

### Save figures
ggsave(file.path(storage_output_path,"trigger_levels_perc.png"),  p, width=4.5, height=4, dpi=600)
ggsave(file.path(storage_output_path,"trigger_levels_perc.pdf"),  p, width=4.5, height=4)
ggsave(file.path(storage_output_path,"trigger_levels_perc.svg"),  p, width=4.5, height=4)


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
###  Test plots
###########################################################################

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

### Save figures
ggsave(file.path(storage_output_path,"storage_bymonth_line.png"),  p, width=7.5, height=5, dpi=600)

p <- ggplot(yup, aes(x=month, y=total_res/1000, group=wy, colour=data))
p <- p + geom_line(size=0.3)
p <- p + scale_colour_manual(name="Scenario", values=c(data_colors[seq(1,3)], "pink", data_colors[seq(4,7)]))
p <- p + theme_classic_new()
p

p <- ggplot(yup, aes(x=month))
p <- p + geom_boxplot(aes(y=total_res/1000, colour=data), size=0.3)
p <- p + scale_x_discrete(name="Month")
p <- p + scale_y_continuous(name="Total Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="Total Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 
p <- p + facet_wrap(~data, nrow=2, scales="free_x")
p <- p + scale_colour_manual(name="Scenario", values=c(data_colors[c(1,2,3,6,4,5,7,8)]), guide=guide_legend(byrow=TRUE))
p <- p + theme_classic_new() + theme(legend.position="bottom")
p <- p + expand_limits(y=c(-10, current_active_stor/1000+10))
p

### Save figures
ggsave(file.path(storage_output_path,"storage_bymonth_boxplot.png"),  p, width=7.5, height=5, dpi=600)


yup <- stor_base
yup$data <- factor(yup$data, levels=data_levels, labels=data_labels)

p <- ggplot(yup, aes(x=month))
p <- p + geom_boxplot(aes(y=total_res/1000, fill=data), size=0.3, colour=NA, alpha=0.6)
p <- p + geom_boxplot(aes(y=total_res/1000, colour=data), size=0.3, fill=NA)
p <- p + scale_colour_manual(name="Scenario", values=c(data_colors))
p <- p + scale_fill_manual(name="Scenario", values=c(data_colors))
p <- p + scale_x_discrete(name="Month")
p <- p + scale_y_continuous(name="Total Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="Total Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 
p <- p + theme_classic_new()
p <- p + theme(legend.position="bottom")
p <- p + guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
p <- p + expand_limits(y=c(-10, current_active_stor/1000+10))
p

### Save figures
ggsave(file.path(storage_output_path,"storage_bymonth_boxplot_merge.png"),  p, width=7.5, height=4.5, dpi=600)


p <- ggplot(subset(yup, month==6), aes(x=total_res/1000, y=..density..,  color=data))
p <- p +  geom_freqpoly(binwidth = 25)
p <- p + scale_colour_manual(name="Scenario", values=c(data_colors))
p <- p + theme_classic_new()
p


p <- ggplot(subset(yup, month==6), aes(x=total_res/1000, y=..density..,  fill=data))
p <- p + geom_histogram(position="dodge",  binwidth = 50)
p <- p + scale_fill_manual(name="Scenario", values=c(data_colors))
p <- p + theme_classic_new()
p

#p <- p + geom_density()
p <- p + facet_wrap(~data, ncol=1)
p

p <- ggplot(subset(yup, month==6), aes(x=total_res/1000, color=data))
p <- p + geom_density()
p <- p + scale_colour_manual(name="Scenario", values=c(data_colors))
p <- p + theme_classic_new()
p

require(ggridges)

yup <- stor_base
yup$data <- factor(yup$data, levels=data_levels[c(4, 5, 1:3, 6:8)], labels=data_labels[c(4, 5, 1:3, 6:8)])

p <- ggplot(subset(yup, month==6), aes(x = total_res/1000, y = data, fill = data, height = ..density..))
p <- p + geom_density_ridges(stat = "density", scale=20)
p <- p + geom_vline(xintercept=unlist(weber_triggers[c(2,3,4)])/1000, linetype="dotted")
p <- p + scale_y_discrete(name="", limits = rev(levels(yup$data)), labels=rev(data_labels_twolines[c(4, 5, 1:3, 6:8)]))
p <- p + scale_x_continuous(name="June Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="June Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0))
p <- p + scale_fill_manual(name="Scenario", values=data_colors[c(4, 5, 1:3, 6:8)])
p <- p + theme_classic_new(9)
#p <- p + theme(legend.position="left")
p <- p + expand_limits(x=c(-10, current_active_stor/1000+10))
p <- p + coord_cartesian(ylim=c(1.2,8.5))
p

### Save figures
ggsave(file.path(storage_output_path,"storage_june_ridge.png"),  p, width=7.5, height=4.5, dpi=600)

yup <- stor_base
yup$data <- factor(yup$data, levels=data_levels[c(4, 5, 1:3, 6:8)], labels=data_labels[c(4, 5, 1:3, 6:8)])

p <- ggplot(subset(yup, month==6), aes(x = total_res/1000, y = data, fill = data, height = ..density..))
p <- p + geom_density_ridges(stat = "binline", bins = 20, scale = 1.8, draw_baseline = FALSE, alpha=0.8)
p <- p + geom_vline(xintercept=unlist(weber_triggers[c(2,3,4)])/1000, linetype="dotted")
p <- p + scale_y_discrete(name="", limits = rev(levels(yup$data)), labels=rev(data_labels_twolines[c(4, 5, 1:3, 6:8)]))
p <- p + scale_x_continuous(name="June Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="June Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0))
p <- p + scale_fill_manual(name="Scenario", values=data_colors[c(4, 5, 1:3, 6:8)])
p <- p + theme_classic_new(9)
#p <- p + theme(legend.position="left")
p <- p + expand_limits(x=c(400, current_active_stor/1000+10))
p <- p + coord_cartesian(ylim=c(1.2,8.5))
p
### Save figures
ggsave(file.path(storage_output_path,"storage_june_ridge_box.png"),  p, width=7.5, height=4.5, dpi=600)

  
  
yup <- stor_base
yup$data <- factor(yup$data, levels=data_levels, labels=data_labels_twolines)

p <- ggplot(subset(yup, month==6), aes(x = data, y = total_res/1000, fill = data))
p <- p + geom_boxplot(size=0.3, colour=NA, alpha=0.6)
p <- p + geom_boxplot(size=0.3, fill=NA)
p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, linetype="dotted")
p <- p + scale_colour_manual(name="Scenario", values=c(data_colors))
p <- p + scale_fill_manual(name="Scenario", values=c(data_colors))
p <- p + scale_x_discrete(name="Scenario")
p <- p + scale_y_continuous(name="June Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="June Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 
p <- p + theme_classic_new()
p <- p + theme(legend.position="bottom")
p <- p + guides(colour = guide_legend(nrow = 2, byrow=TRUE), fill = guide_legend(nrow = 2, byrow=TRUE))
p <- p + expand_limits(y=c(-10, current_active_stor/1000+10))
p

### Save figures
ggsave(file.path(storage_output_path,"storage_june_box.png"),  p, width=6.75, height=4.5, dpi=600)



### Check plot
p <- ggplot(subset(stor_base_percent, data!="base") , aes(colour=data))
p <- p + geom_line(aes(x=date, y=current_res), size=0.3)
p <- p + geom_line(data=subset(stor_base_percent, data=="base"), aes(x=base_date, y=current_res), size=0.3)
p <- p + theme_classic_new()
p

###########################################################################
## Plot Storage Time Series as Area
###########################################################################

area_df <- stor_base %>% 
	filter(data %in% c("paleo", "observed", "hd"))
	
time_test <- area_df$data == "paleo" & area_df$year >= 1905
area_df <- area_df[!time_test,]

p <- ggplot(area_df, aes(x=date, y=total_res/1000, fill=data))
p <- p + geom_area( position = "identity", alpha=0.8)
p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, linetype="dotted")
#p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Scenario", values= data_colors[c(1, 2, 7)], labels=c( "Reconstructed", "Observed", "Climate Change (Hot-Dry)"), limits=c("paleo","observed",  "hd"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1425-01-01"), as.Date("2070-01-01")), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Total Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="Total Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 
p <- p + theme(legend.position="bottom")
p <- p + expand_limits(y=c(0, current_active_stor/1000+10))
p


### Save figures
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_full.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_full.pdf"),  p, width=8, height=3.5)
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_full.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_1900.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_1900.pdf"),  p, width=8, height=3.5)
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_1900.svg"),  p, width=8, height=3.5)


area_df <- stor_base %>% 
	filter(data %in% c("paleo", "observed", "hd")) %>%
	filter(month == 6)
	
time_test <- area_df$data == "paleo" & area_df$year >= 1905
area_df <- area_df[!time_test,]

p <- ggplot(area_df, aes(x=date, y=total_res/1000, fill=data))
p <- p + geom_area( position = "identity", alpha=0.8)
p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, linetype="dotted")
#p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Scenario", values= data_colors[c(1, 2, 7)], labels=c( "Reconstructed", "Observed", "Climate Change (Hot-Dry)"), limits=c("paleo","observed",  "hd"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1425-01-01"), as.Date("2070-01-01")), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="June Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="June Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 
p <- p + theme(legend.position="bottom")
p <- p + expand_limits(y=c(0, current_active_stor/1000+10))
p

### Save figures
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_full_june.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_full_june.pdf"),  p, width=8, height=3.5)
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_full_june.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_1900_june.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_1900_june.pdf"),  p, width=8, height=3.5)
ggsave(file.path(storage_output_path,"paleo_future_stor_area_acft_hd_1900_june.svg"),  p, width=8, height=3.5)



###########################################################################
## Create dataframe for triggers
###########################################################################
### Add a 6 year buffer for plotting
trigger_plot <- data.frame(date = rep(c(min(stor_percent$date) %m+% years(-6), max(stor_percent$date) %m+% years(6)), 3))
trigger_plot$month <- month(trigger_plot$date)
trigger_plot$year <- year(trigger_plot$date)

### Add water year column
trigger_plot$wy <- usgs_wateryear(year=trigger_plot$year, month=trigger_plot$month)

### Add triggers
trigger_plot$trigger_level <- rep(c("Moderate", "Severe", "Extreme"), each=2)
trigger_plot$trigger_level <- factor(trigger_plot$trigger_level, levels=c("Moderate", "Severe", "Extreme"))
trigger_plot$res_perc <- rep(c(mod_perc_annual-sev_perc_annual, sev_perc_annual-ext_perc_annual, ext_perc_annual), each=2)
### This line wasn't previously in here
trigger_plot$res_stor <- rep(c(mod_annual, sev_annual, ext_annual), each=2)


###########################################################################
## Plot Storage Time Series as Line
###########################################################################
line_df <- stor_all  %>% 
	filter(data %in% c("paleo", "observed", "hd"))

#line_df2 <- stor_base  %>% #
#	filter(data %in% c("base")) %>%
#	mutate(date = base_date)

#line_df <- rbind(line_df, line_df2)

for (response_i in response_levels){
	line_plot <- line_df %>% 
		filter(response == response_i)
		
for (i in seq(1, length(region_levels))) {
region_i <- as.character(region_levels[i])
region_plot_name <- region_names[i]

y_name_left <- paste0(region_plot_name, " Active Storage (1,000 ac-ft)")
y_name_right <- paste0(region_plot_name, " Active Storage (%)")

if (i <= 2){
	y_name_left <- "Active Storage (1,000 ac-ft)"
	y_name_right <- "Active Storage (%)"
}

region_stor_i <- active_stor_vec[i]
trigger_plot$res_stor_region <- trigger_plot$res_perc * region_stor_i

p <- ggplot(line_plot, aes(x=date))
p <- p + geom_area(data=trigger_plot, aes(y=res_stor_region/1000, fill=trigger_level), alpha=0.85)
p <- p + geom_line(size=0.18, aes(y=get(region_i)/1000, group=data))#, alpha=0.8)
#p <- p + geom_hline(yintercept=0, size=0.15)
p <- p + theme_classic_new()
p <- p + scale_colour_manual(name="Scenario", values= data_colors[c(1, 2,3, 7)], labels=c( "Reconstructed", "Observed", "Base", "Climate Change (HD)"), limits=c("paleo","observed", "base", "hd"), guide = guide_legend())
p <- p + scale_fill_manual(name="Storage Trigger", values= trigger_colors, limits= c("Extreme", "Severe", "Moderate"), labels=c("Extreme", "Severe", "Moderate"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1425-01-01"), as.Date("2070-01-01")), expand=FALSE) + expand_limits(y=0)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name=y_name_left, sec.axis = sec_axis(name=y_name_right, trans=~ ((.*1000)/region_stor_i)*100, breaks=seq(0,100,25)),  expand = c(0,0)) 
p <- p + theme(legend.position="bottom")
p <- p + expand_limits(y=c(0, region_stor_i/1000+10))
p

plot_location <- file.path(storage_output_path,paste0("time_series/", response_i))
dir.create(plot_location, recursive=TRUE, showWarnings=FALSE)

### Save figures
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,".png")),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,".pdf")),  p, width=8, height=3.5)
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,".svg")),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")), expand=FALSE) + expand_limits(y=0)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,"_1900.png")),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,"_1900.pdf")),  p, width=8, height=3.5)
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,"_1900.svg")),  p, width=8, height=3.5)

### Cut to only June
june_plot <- line_plot %>%
	filter(month == 6)

y_name_left <- paste0(region_plot_name, " June Active Storage (1,000 ac-ft)")
y_name_right <- paste0(region_plot_name, " June Active Storage (%)")

if (i <= 2){
	y_name_left <- "June Active Storage (1,000 ac-ft)"
	y_name_right <- "June Active Storage (%)"
}
	
	
p <- ggplot(june_plot, aes(x=date))
p <- p + geom_area(data=trigger_plot, aes(y=res_stor_region/1000, fill=trigger_level), alpha=0.85)
p <- p + geom_line(size=0.35, aes(y=get(region_i)/1000, group=data))#, alpha=0.8)
p <- p + theme_classic_new()
p <- p + scale_colour_manual(name="Scenario", values= data_colors[c(1, 2,3, 7)], labels=c( "Reconstructed", "Observed", "Base", "Climate Change (HD)"), limits=c("paleo","observed", "base", "hd"), guide = guide_legend())
p <- p + scale_fill_manual(name="Storage Trigger", values= trigger_colors, limits= c("Extreme", "Severe", "Moderate"), labels=c("Extreme", "Severe", "Moderate"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1425-01-01"), as.Date("2070-01-01")), expand=FALSE) + expand_limits(y=0)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name=y_name_left, sec.axis = sec_axis(name=y_name_right, trans=~ ((.*1000)/region_stor_i)*100, breaks=seq(0,100,25)),  expand = c(0,0)) 
p <- p + theme(legend.position="bottom")
p <- p + expand_limits(y=c(0, region_stor_i/1000+10))
p

plot_location <- file.path(storage_output_path,paste0("time_series/", response_i))
dir.create(plot_location, recursive=TRUE, showWarnings=FALSE)

### Save figures
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,"_june.png")),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,"_june.pdf")),  p, width=8, height=3.5)
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,"_june.svg")),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")), expand=FALSE) + expand_limits(y=0)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,"_1900_june.png")),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,"_1900_june.pdf")),  p, width=8, height=3.5)
ggsave(file.path(plot_location,paste0("stor_timeseries_acft_hd_", region_i, "_", response_i,"_1900_june.svg")),  p, width=8, height=3.5)


}
}



###########################################################################
## Plot Shortage Time Series by region
###########################################################################
for (response_i in response_levels){

area_region <- stor_all %>% 
	filter(data %in% c("paleo", "observed", "base", "response", "hd")) %>%
	filter(response == response_i) %>%
	select(date, base_date, data, upper_ogden, upper_weber, lower) %>%
	gather(res, value, c("upper_ogden", "upper_weber", "lower"))

area_region$data <- factor(area_region$data)
area_region$res <- factor(area_region$res, levels=c("upper_ogden", "upper_weber", "lower"))
	
p <- ggplot(subset(area_region, data=="observed"), aes(x=date, y=value/1000, fill=res))
p <- p + geom_area()
p <- p + geom_area(data=subset(area_region, data=="paleo" & date < as.Date("1905-01-01")))
p <- p + geom_area(data=subset(area_region, data=="observed" & date >= as.Date("1905-01-01")))
p <- p + geom_area(data=subset(area_region, data=="hd"))
#p <- p + geom_area(data=subset(area_region, data=="base"), aes(x=base_date))
#p <- p + geom_line(data=trigger_plot, aes(y=res_stor/1000, group=trigger_level, fill=NA), colour="grey30", linetype="longdash", size=0.35)
p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, colour="grey30", linetype="longdash", size=0.35)
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Region", values= res_colors[c(1,3,2)], labels=c("Upper Ogden", "Upper Weber", "Lower Weber"), breaks=c("upper_ogden", "upper_weber", "lower"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1425-01-01"), as.Date("2070-01-01")), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Total Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="Total Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 
p <- p + theme(legend.position="bottom")
p <- p + expand_limits(y=c(0, current_active_stor/1000+10))
p


plot_location <- file.path(storage_output_path,paste0("time_series_regions/", response_i))
dir.create(plot_location, recursive=TRUE, showWarnings=FALSE)

### Save figures
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,".png")),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,".pdf")),  p, width=8, height=3.5)
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,".svg")),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")), expand=FALSE) + expand_limits(y=0)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,"_1900.png")),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,"_1900.pdf")),  p, width=8, height=3.5)
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,"_1900.svg")),  p, width=8, height=3.5)

area_region <- area_region %>% 
	mutate(month = month(date)) %>%
	filter(month == 6)
	
p <- ggplot(subset(area_region, data=="observed"), aes(x=date, y=value/1000, fill=res))
p <- p + geom_area()
p <- p + geom_area(data=subset(area_region, data=="paleo" & date < as.Date("1905-01-01")))
p <- p + geom_area(data=subset(area_region, data=="observed" & date >= as.Date("1905-01-01")))
p <- p + geom_area(data=subset(area_region, data=="hd"))
#p <- p + geom_area(data=subset(area_region, data=="base"), aes(x=base_date))
#p <- p + geom_line(data=trigger_plot, aes(y=res_stor/1000, group=trigger_level, fill=NA), colour="grey30", linetype="longdash", size=0.35)
p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, colour="grey30", linetype="longdash", size=0.35)
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Region", values= res_colors[c(1,3,2)], labels=c("Upper Ogden", "Upper Weber", "Lower Weber"), breaks=c("upper_ogden", "upper_weber", "lower"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1425-01-01"), as.Date("2070-01-01")), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="June Active Storage (1,000 ac-ft)", breaks=seq(0,800,100), sec.axis = sec_axis(name="June Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 
p <- p + theme(legend.position="bottom")
p <- p + expand_limits(y=c(0, current_active_stor/1000+10))
p


plot_location <- file.path(storage_output_path,paste0("time_series_regions/", response_i))
dir.create(plot_location, recursive=TRUE, showWarnings=FALSE)

### Save figures
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,"_june.png")),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,"_june.pdf")),  p, width=8, height=3.5)
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,"_june.svg")),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")), expand=FALSE) + expand_limits(y=0)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,"_1900_june.png")),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,"_1900_june.pdf")),  p, width=8, height=3.5)
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_", response_i,"_1900_june.svg")),  p, width=8, height=3.5)

}


###########################################################################
## Plot Shortage Time Series by region using a facet
###########################################################################
area_region <- stor_all %>% 
	filter(data %in% c("paleo", "observed", "base", "response", "hd")) %>%
	select(date, base_date, data, response, upper_ogden, upper_weber, lower) %>%
	gather(res, value, c("upper_ogden", "upper_weber", "lower"))

area_region$data <- factor(area_region$data)
area_region$res <- factor(area_region$res, levels=c("upper_ogden", "upper_weber", "lower"))
	
p <- ggplot(subset(area_region, data=="observed"), aes(x=date, y=value/1000, fill=res))
#p <- p + geom_area()
p <- p + stat_smooth(geom = 'area', method = 'loess', span = 1/30, position='stack')
      
p <- p + geom_area(data=subset(area_region, data=="paleo" & date < as.Date("1980-10-01")))
p <- p + geom_area(data=subset(area_region, data=="hd"))
p <- p + geom_area(data=subset(area_region, data=="base"), aes(x=base_date))
p <- p + geom_line(data=trigger_plot, aes(y=res_stor/1000, group=trigger_level, fill=NA), colour="grey30", linetype="longdash", size=0.35)
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Region", values= res_colors[c(1,3,2)], labels=c("Upper Ogden", "Upper Weber", "Lower Weber"), breaks=c("upper_ogden", "upper_weber", "lower"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1425-01-01"), as.Date("2070-01-01")), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="System Storage (1,000 ac-ft)", breaks=seq(0,800,100))
p <- p + theme(legend.position="bottom")
p <- p + facet_grid(response ~ .)
p

plot_location <- file.path(storage_output_path,"time_series_regions")
dir.create(plot_location, recursive=TRUE, showWarnings=FALSE)

### Save figures
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_facet.png")),  p, width=8, height=8, dpi=300)
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_facet.pdf")),  p, width=8, height=8)
ggsave(file.path(plot_location,paste0("stor_byregion_timeseries_percent_hd_facet.svg")),  p, width=8, height=8)




### THere were differences in this figure, here's how it used to look

#area_region <- area_df[, names(area_df) %in% c("date", "data", "upper_ogden", "upper_weber", "lower")]
#area_region <- melt(area_region, id.vars=c("date", "data"))#, measure.vars=c("upper_ogden", "upper_weber", "lower"))
#area_region$data <- factor(area_region$data)
#area_region$variable <- factor(area_region$variable)

#p <- ggplot(subset(area_region, data=="observed"), aes(x=date, y=value/1000, fill=variable))
#p <- ggplot(yup, aes(x=date, y=value/1000, fill=variable))
#p <- p + geom_area()
#p <- p + geom_area(data=subset(area_region, data=="paleo"))
#p <- p + geom_area(data=subset(area_region, data=="hd"))
#p <- p + theme_classic_new()
#p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "grey30", cc_colors), labels=c("Observed", "Base", "HD", "HW", "WD", "WW", "Reconst" ))
#p <- p + scale_fill_manual(name="Region", values= res_colors, labels=c("Upper Ogden", "Upper Weber", "Lower Weber"), guide = guide_legend())
#p <- p + coord_cartesian(xlim=c(as.Date("1425-01-01"), as.Date("2070-01-01")), expand=FALSE)
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
#p <- p + scale_y_continuous(name="System Storage (1,000 ac-ft)", breaks=seq(0,800,100))
#p <- p + theme(legend.position="bottom")
#p


###########################################################################
## Plot Shortage Time Series
###########################################################################
### Need to fix the climate change run

area_df <- stor_all %>% 
	filter(data %in% c("paleo", "observed", "base", "response", "hd")) %>%
	select(date, base_date, data, response, upper_ogden, upper_weber, lower) %>%
	gather(res, value, c("upper_ogden", "upper_weber", "lower"))

area_df <- stor_all  %>% 
	filter(data %in% c("paleo", "observed", "hd"))

area_df2 <- stor_all  %>% 
	filter(data %in% c("base")) %>%
	mutate(date = base_date)

area_df <- rbind(area_df, area_df)
area_plot <- area_df %>% 
		filter(response == "Base") %>%
		filter(month == 6)
		
#plot_test <- stor_all$data %in% c("observed", "paleo")
#cc_test <- stor_all$data %in% c("base", "hd", "hw", "wd", "ww")
#cc_df <- stor_all[cc_test,]

#plot_df <- rbind(stor_all[plot_test,], cc_df[cc_df$data %in% c("base", "HDN5"),])

#base_df <- cc_df[cc_df$data %in% c("base"),]
#base_df$data <- factor("observed", levels="observed")

#area_df <- rbind(cc_df[cc_df$data %in% c("hd", "hw"),], stor_all[stor_all$data %in% c("paleo", "observed"),])

p <- ggplot(area_plot, aes(x=date, y=-trigger_deficit/1000, fill=data))
p <- p + geom_area( position = "identity", alpha=0.8)
p <- p + geom_hline(yintercept=0, size=0.2)
p <- p + geom_hline(yintercept=-current_active_stor/1000, size=0.4, colour="black", linetype="longdash")

#p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
#p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "grey30", cc_colors), labels=c("Observed", "Base", "HD", "HW", "WD", "WW", "Reconst" ))
p <- p + scale_fill_manual(name="Scenario", values= c("#D55E00", "#56B4E9", "grey30", "#CC79A7"), labels=c("HD", "HW", "Observed", "Reconstr"), guide = guide_legend())
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="System Storage Volume Below Moderate Trigger (1,000 ac-ft)", breaks=seq(-10000,500,50))


p <- p + scale_y_continuous(name="June Active Storage (1,000 ac-ft)", breaks=seq(-10000,500,50), sec.axis = sec_axis(name="June Active Storage (%)", trans=~ (((.*-1000)+380000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)) 

-weber_triggers[2]


p <- p + theme(legend.position="bottom")
p

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_hot_full.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_hot_full.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_hot_full.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_hot_1900.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_hot_1900.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_hot_1900.svg"),  p, width=8, height=3.5)



base_df <- cc_df[cc_df$data %in% c("base"),]
base_df$data <- factor("observed", levels="observed")

area_df <- rbind(cc_df[cc_df$data %in% c("wd", "ww"),], stor_all[stor_all$data %in% c("paleo", "observed"),])

p <- ggplot(area_df, aes(x=date, y=-system_def/1000, fill=data))
p <- p + geom_area( position = "identity", alpha=0.8)
p <- p + geom_hline(yintercept=0, size=0.2)
p <- p + geom_hline(yintercept=-current_active_stor/1000, size=0.4, colour="black", linetype="longdash")
p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "#CC79A7", "#E69F00", "#0072B2"), labels=c("Observed", "Reconstr", "WD", "WW"), guide = guide_legend())
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="System Storage Deficit (1,000 ac-ft)", breaks=seq(-10000,500,50))
p <- p + theme(legend.position="bottom")
p

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_warm_full.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_warm_full.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_warm_full.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_warm_1900.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_warm_1900.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_deficit_area_acft_warm_1900.svg"),  p, width=8, height=3.5)



###########################################################################
## Plot Area Time Series for climate change
###########################################################################

plot_test <- stor_base$data %in% c( "hd", "hw", "ct", "wd", "ww")

area_df <- stor_base[plot_test,]
area_df$data <- factor(area_df$data, levels=c( "ww", "hw","ct", "wd", "hd"), labels=c("Warm-Wet", "Hot-Wet", "Median", "Warm-Dry", "Hot-Dry"))

line_df <- stor_base[stor_base$data %in% "base", ]
line_df$date <- line_df$date %m+% years(55)
line_df <- do.call("rbind", replicate(5, line_df, simplify = FALSE))
line_df$data <- area_df$data
line_df$precip <- area_df$precip
line_df$temp <- area_df$temp

p <- ggplot(area_df, aes(x=date, y=total_res/1000))
p <- p + geom_area(aes(fill=data))
p <- p + geom_hline(yintercept=0)
p <- p + geom_line(data=line_df, colour="black")#, linetype="longdash")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Scenario", values= cc_colors, labels=c("WW", "HW", "CT", "WD", "HD"), breaks=c("ww", "hw", "ct", "wd", "hd"))
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="5 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Active Storage (1,000 ac-ft)", breaks=seq(-5000,5000,200))
p <- p + facet_grid(data ~ .)
p

#, sec.axis = sec_axis(name="Total Active Storage (%)", trans=~ ((.*1000)/current_active_stor)*100, breaks=seq(0,100,25)), expand = c(0,0)


### Save figures
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_vert.png"),  p, width=7.5, height=7, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_vert.pdf"),  p, width=7.5, height=7)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_vert.svg"),  p, width=7.5, height=7)

p <- ggplot(subset(area_df, temp != "Median"), aes(x=date, y=total_res/1000))
p <- p + geom_area(aes(fill=data))
p <- p + geom_hline(yintercept=0)
p <- p + geom_line(data=subset(line_df, temp != "Median"), colour="black")#, linetype="longdash")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Scenario", values= cc_colors[c(1, 2, 4, 5)], labels=c("WW", "HW", "WD", "HD"), breaks=c("ww", "hw", "wd", "hd"))
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="5 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Active Storage (1,000 ac-ft)", breaks=seq(-5000,5000,200))
p <- p + facet_grid(precip ~ temp)

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_square.png"),  p, width=8, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_square.pdf"),  p, width=8, height=6)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_square.svg"),  p, width=8, height=6)


###########################################################################
## Plot Area Time Series for climate change with June Only
###########################################################################

p <- ggplot(subset(area_df, month==6), aes(x=date, y=total_res/1000))
p <- p + geom_area(aes(fill=data))
#p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, colour="grey30", size=0.35) #linetype="longdash", 
#p <- p + geom_hline(yintercept=unlist(weber_triggers[2])/1000, colour="grey30", linetype="longdash", size=0.35)
p <- p + geom_hline(yintercept=0)
p <- p + geom_line(data=subset(line_df, month==6), colour="black")#, linetype="longdash")
#p <- p + geom_hline(yintercept=unlist(weber_triggers[c(2,3,4)])/1000, colour="grey30", linetype="longdash", size=0.35)
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Scenario", values= cc_colors, labels=c("WW", "HW", "CT", "WD", "HD"), breaks=c("ww", "hw", "ct", "wd", "hd"))
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="5 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="June Active Storage (1,000 ac-ft)", breaks=seq(-5000,5000,200))
p <- p + facet_grid(data ~ .)
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_vert_june.png"),  p, width=7, height=7, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_vert_june.pdf"),  p, width=7, height=7)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_vert_june.svg"),  p, width=7, height=7)

p <- ggplot(subset(area_df, temp != "Median" & month==6), aes(x=date, y=total_res/1000))
p <- p + geom_area(aes(fill=data))
p <- p + geom_hline(yintercept=0)
p <- p + geom_line(data=subset(line_df, temp != "Median" & month==6), colour="black")#, linetype="longdash")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Scenario", values= cc_colors[c(1, 2, 4, 5)], labels=c("WW", "HW", "WD", "HD"), breaks=c("ww", "hw", "wd", "hd"))
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="5 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="June Active Storage (1,000 ac-ft)", breaks=seq(-5000,5000,200))
p <- p + facet_grid(precip ~ temp)

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_square_june.png"),  p, width=8, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_square_june.pdf"),  p, width=8, height=6)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_square_june.svg"),  p, width=8, height=6)



###########################################################################
## Plot Area Time Series for climate change By Region
###########################################################################

area_region <- area_df[, names(area_df) %in% c("date", "data", "upper_ogden", "upper_weber", "lower")]
area_region <- melt(area_region, id.vars=c("date", "data"))#, measure.vars=c("upper_ogden", "upper_weber", "lower"))
area_region$data <- factor(area_region$data)
area_region$variable <- factor(area_region$variable)
line_df$value <- line_df$total_res

p <- ggplot(area_region, aes(x=date, y=value/1000))
p <- p + geom_area(aes(fill=variable))
#p <- p + geom_hline(yintercept=0)
#p <- p + geom_line(data=line_df, colour="black")#, linetype="longdash")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Region", values= res_colors, labels=c("Upper Ogden", "Upper Weber", "Lower Weber"), guide = guide_legend())
p <- p + theme(legend.position="bottom")
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="5 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="System Storage (1,000 ac-ft)", breaks=seq(-5000,5000,100))
p <- p + facet_grid(data ~ .)
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_storage_byregion_acft_vert.png"),  p, width=7, height=7, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_storage_byregion_acft_vert.pdf"),  p, width=7, height=7)
ggsave(file.path(write_output_base_path,"clim_change_storage_byregion_acft_vert.svg"),  p, width=7, height=7)


p <- p + facet_grid(precip ~ temp)

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_storage_byregion_acft_square.png"),  p, width=8, height=7, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_storage_byregion_acft_square.pdf"),  p, width=8, height=7)
ggsave(file.path(write_output_base_path,"clim_change_storage_byregion_acft_square.svg"),  p, width=8, height=7)













###################################
###  Silly attempts
###################################
require(ggridges)

p <- ggplot(yup %>% filter(data == "paleo"), aes(x=total_res, y=month))
#p <- p + geom_density_ridges()
p <- p + geom_density_ridges(scale = 3)#, rel_min_height = 0.03)scale = 10, 
#p <- p + scale_colour_manual(name="Scenario", values=data_colors)
p <- p +  theme_ridges()
p

p <- ggplot(yup %>% filter(data %in% c("paleo", "hd")), aes(x=total_res, y=month, fill=data))
#p <- p + geom_density_ridges()
p <- p + geom_density_ridges(scale = 2, alpha=0.3)#, rel_min_height = 0.03)scale = 10, 
#p <- p + scale_colour_manual(name="Scenario", values=data_colors)
p <- p +  theme_ridges()
p

p <- ggplot(yup %>% filter(data == "paleo"), aes(x=total_res, y=month, fill= ..x..))
#p <- p + geom_density_ridges()
p <- p + geom_density_ridges_gradient(scale = 3)#, rel_min_height = 0.03)
#p <- p + scale_colour_manual(name="Scenario", values=data_colors)
p <- p +  theme_ridges()
p

















###########################################################################
###  Plot against triggers
###########################################################################


plot_triggers <- merge(plot_df, weber_triggers, by.x="month", by.y="Month")

plot_df$trigger <- "None"
plot_df$trigger[plot_df$system_min < weber_triggers[6,2] * 0.25] <- "Extreme"
plot_df$trigger[plot_df$system_min < weber_triggers[6,2] * 0.5 & plot_df$trigger == "None"] <- "Severe"
plot_df$trigger[plot_df$system_min < weber_triggers[6,2] * 0.7 & plot_df$trigger == "None"] <- "Moderate"




p <- ggplot(plot_df, aes(x=dura_months/12, y=min_perc*100, colour=trigger))
p <- p + geom_point(data=subset(plot_df, trigger=="None"), size=4, alpha=0.7)#, aes(shape=data))
p <- p + geom_point(data=subset(plot_df, trigger!="None"), size=4)#, 
#p <- p + scale_colour_viridis(direction=-1, option="plasma", name="Minimum\nStorage\n(1,000 ac-ft)")
p <- p + theme_classic_new(12)
p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,1))
p <- p + scale_y_continuous(name="Min Monthly Flow Percentile")
p <- p + scale_colour_manual(name="Storage Trigger", values= c("grey80", "#4daf4a", "#e5c494", "#e41a1c"), limits= c("None", "Moderate", "Severe", "Extreme"), guide = guide_legend())
p <- p + theme(legend.position = c(0.85, 0.75))
p <- p + coord_cartesian(xlim = c(0.35,6.5), ylim = c(0,47), expand = TRUE)
p



### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_trigger.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_trigger.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_trigger.svg"),  p, width=6.5, height=5)




###########################################################################
## Plot Shortage Time Series
###########################################################################
trigger_subset <- subset(stor_all, select=c("date", "data", "moderate", "severe", "extreme"))

trigger_subset <- melt(trigger_subset, id.vars=c("date", "data"), measure.vars=c("moderate", "severe", "extreme"))
trigger_subset$variable <- factor(trigger_subset$variable, levels=c("moderate", "severe", "extreme"))
trigger_subset$data_variable <- paste0(trigger_subset$data,"_",trigger_subset$variable)

plot_test <- trigger_subset$data %in% c("observed", "paleo", "hd")
area_df <- trigger_subset[plot_test,]
area_df$variable <- factor(area_df$variable, levels=c("extreme", "severe", "moderate"))

#cc_test <- trigger_subset$data %in% c("base", "hd", "hw", "wd", "ww")
#cc_df <- trigger_subset[cc_test,]

#base_df <- cc_df[cc_df$data %in% c("base"),]
#base_df$data <- factor("observed", levels="observed")

#area_df <- rbind(cc_df[cc_df$data %in% c("hd", "hw"),], stor_all[stor_all$data %in% c("paleo", "observed"),])

p <- ggplot(area_df, aes(x=date, y=-value/1000, fill=variable))#, group=data_variable))
p <- p + geom_rect(aes(xmin=as.Date("00-12-01"),xmax=as.Date("1903-12-01"),ymin=-Inf,ymax=Inf),alpha=0.1,fill="grey90")
p <- p + geom_rect(aes(xmin=as.Date("2035-01-01"),xmax=as.Date("2065-12-01"),ymin=-Inf,ymax=Inf),alpha=0.1,fill="grey90")
p <- p + geom_area( data=subset(area_df, data=="observed"), position = "stack")
p <- p + geom_area( data=subset(area_df, data=="paleo"), position = "stack")
p <- p + geom_area( data=subset(area_df, data=="hd"), position = "stack")
p <- p + geom_hline(yintercept=0, size=0.2)
#p <- p + geom_hline(yintercept=-max_stor/1000, size=0.4, colour="black", linetype="longdash")
#p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Storage Trigger", values= c("#4daf4a", "#e5c494", "#e41a1c"), limits= c("moderate", "severe", "extreme"), labels=c("Moderate", "Severe", "Extreme"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1400-01-01"), as.Date("2070-01-01")), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Storage Trigger Deficit (1,000 ac-ft)", breaks=seq(-10000,500,50))
p <- p + theme(legend.position="bottom")
#p <- p + coord_cartesian(xlim=c(as.Date("1400-01-01", as.Date("2070-12-01"))
p


### Save figures
ggsave(file.path(write_output_base_path,"trigger_storage.png"),  p, width=8, height=4, dpi=300)
ggsave(file.path(write_output_base_path,"trigger_storage.pdf"),  p, width=8, height=4)
ggsave(file.path(write_output_base_path,"trigger_storage.svg"),  p, width=8, height=4)


p <- ggplot(area_df, aes(x=date, y=-value/1000, fill=variable))#, group=data_variable))
#p <- p + geom_rect(aes(xmin=as.Date("00-12-01"),xmax=as.Date("1903-12-01"),ymin=-Inf,ymax=Inf),alpha=0.1,fill="grey90")
#p <- p + geom_rect(aes(xmin=as.Date("2035-01-01"),xmax=as.Date("2065-12-01"),ymin=-Inf,ymax=Inf),alpha=0.1,fill="grey90")
p <- p + geom_vline(xintercept=as.Date("1904-01-01"), linetype="dotted")
#p <- p + geom_vline(xintercept=as.Date("2002-09-01"), linetype="dotted")
p <- p + geom_vline(xintercept=as.Date("2035-01-01"), linetype="dotted")

p <- p + geom_area( data=subset(area_df, data=="observed"), position = "stack")
p <- p + geom_area( data=subset(area_df, data=="paleo"), position = "stack")
p <- p + geom_area( data=subset(area_df, data=="hd"), position = "stack")
p <- p + geom_hline(yintercept=0, size=0.2)
#p <- p + geom_hline(yintercept=-max_stor/1000, size=0.4, colour="black", linetype="longdash")
#p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Storage Trigger", values= c("#f03b20",  "#feb24c", "#ffeda0"), limits= c("extreme", "severe", "moderate"), labels=c("Extreme", "Severe", "Moderate"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1400-01-01"), as.Date("2070-01-01")), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Storage Trigger Deficit (1,000 ac-ft)", breaks=seq(-10000,500,50))
p <- p + theme(legend.position="bottom")
#p <- p + coord_cartesian(xlim=c(as.Date("1400-01-01", as.Date("2070-12-01"))
p


### Save figures
ggsave(file.path(write_output_base_path,"trigger_storage_red.png"),  p, width=8, height=4, dpi=300)
ggsave(file.path(write_output_base_path,"trigger_storage_red.pdf"),  p, width=8, height=4)
ggsave(file.path(write_output_base_path,"trigger_storage_red.svg"),  p, width=8, height=4)





###########################################################################
## Plot 
###########################################################################
stor_bar <- stor_all
stor_bar$trigger <- "None"
stor_bar$trigger[stor_bar$extreme > 0] <- "Extreme"
stor_bar$trigger[stor_bar$severe > 0 & stor_bar$trigger == "None"] <- "Severe"
stor_bar$trigger[stor_bar$moderate > 0 & stor_bar$trigger == "None"] <- "Moderate"

stor_bar$data <- factor(stor_bar$data, levels=c( "paleo", "observed", "base" ,"ww", "hw", "wd", "hd"), labels=c("Paleo", "Observed", "Base", "Warm-Wet", "Hot-Wet", "Warm-Dry", "Hot-Dry"))
 
stor_bar$trigger <- factor(stor_bar$trigger, levels=c( "None","Extreme", "Severe", "Moderate"))

 
p <- ggplot(stor_bar, aes(x=data, fill=trigger))
p <- p + geom_bar(position = "fill")
#p <- p + geom_text(size = 3, position = position_stack(vjust = 0.5))
p <- p + coord_cartesian(ylim=c(0,0.4), expand=FALSE)
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Storage Trigger", values= c( "#e41a1c", "#e5c494", "#4daf4a"), limits= c("Extreme", "Severe", "Moderate"), guide = guide_legend())
p




p <- ggplot(stor_bar, aes(x=data, fill=trigger))
p <- p + geom_bar(position = "fill")
#p <- p + geom_text(size = 3, position = position_stack(vjust = 0.5))
p <- p + coord_cartesian(ylim=c(0,0.25), expand=FALSE)
p <- p + theme_classic_new()
p <- p + scale_y_continuous(name="Proportion of Months", breaks=seq(0,1,0.05), labels = scales::percent)
p <- p + scale_x_discrete(name="Data")
p <- p + scale_fill_manual(name="Storage Trigger", values= c("#f03b20",  "#feb24c", "#ffeda0"), limits= c("Extreme", "Severe", "Moderate"))
p <- p + theme(legend.position = c(0.2, 0.8))
p <- p + theme(panel.grid.major.y = element_line(size = 0.2, linetype = 'solid', colour = "grey85"))
p
                             

### Save figures
ggsave(file.path(write_output_base_path,"trigger_proportion.png"),  p, width=5, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"trigger_proportion.pdf"),  p, width=5, height=3.5)
ggsave(file.path(write_output_base_path,"trigger_proportion.svg"),  p, width=5, height=3.5)


### Create Table of occurences
trig_occur <- subset(stor_bar, select=c("data", "trigger"))
trig_occur <- table(trig_occur)
trig_occur

### Calculate total number of time steps and then divide to create proportions
total_time_steps <- apply(trig_occur, 1, sum)
#
trig_prop <- sweep(trig_occur, 1, total_time_steps, FUN="/")


### Calculate chi-square pairs
chi_p <- matrix(NA, dim(trig_occur)[1], dim(trig_occur)[1])
colnames(chi_p) <- rownames(trig_occur)
rownames(chi_p) <- rownames(trig_occur)

for(i in seq(1, dim(trig_occur)[1])){
	for(j in seq(1, dim(trig_occur)[1])){
		if (i > j) {
			chi_result <- chisq.test(trig_occur[c(i,j),], simulate.p.value=TRUE)
			chi_p[i,j] <- chi_result$p.value
		}	
	}
}


### Write data
trig_occur <- cbind(trig_occur, total_time_steps)

write.csv(trig_prop , file.path(write_output_base_path,"trigger_proportion.csv") )
write.csv(trig_occur , file.path(write_output_base_path,"trigger_occur.csv") )
write.csv(chi_p , file.path(write_output_base_path,"trigger_chi_p.csv") )






##################################################################
### Hashimoto figures
##################################################################

p <- ggplot(subset(hash_stor, response=="Base"), aes(x=reliability, y=resilience, color=vulnerability, shape=data))
p <- p + geom_point(size=6)
p <- p + scale_shape_manual(values=c(16,17, 15, 0, 3, 2, 7))
#p <- p + scale_colour_viridis(direction=-1, option="plasma")
#p <- p + scale_colour_distiller(type = "seq", palette = "YlOrRd", direction=1)
p <- p + scale_colour_gradient(low="#feb24c", high="#800026")
p <- p + theme_classic_correct()
p

















###########################################################################
## Plot drought events duration vs severity
###########################################################################


plot_df <- drought_event_summary[!is.na(drought_event_summary$system_min), ]
plot_df_subset <- plot_df[plot_df$data %in% c("paleo", "observed", "hd"),]

require(viridis)

p <- ggplot(plot_df, aes(x=dura_months/12, y=min_perc*100, colour=system_min/1000))
p <- p + geom_point(size=4)#, aes(shape=data))
p <- p + scale_colour_viridis(direction=-1, option="plasma", name="Minimum\nStorage\n(1,000 ac-ft)")
p <- p + theme_classic_new(12)
p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,1))
p <- p + scale_y_continuous(name="Min Monthly Flow Percentile")
p <- p + theme(legend.position = c(0.85, 0.75))
p <- p + coord_cartesian(xlim = c(0.35,6.5), ylim = c(0,47), expand = TRUE)
p

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_system_storage.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_system_storage.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_system_storage.svg"),  p, width=6.5, height=5)

p <- p + scale_colour_distiller(name="Minimum\nStorage\n(1,000 ac-ft)", type = "seq", palette = "Oranges", direction = -1)


#p + scale_colour_distiller(name="Minimum\nStorage\n(1,000 ac-ft)", type = "seq", palette = "OrRd", direction = -1)
#p + scale_colour_distiller(name="Minimum\nStorage\n(1,000 ac-ft)", type = "seq", palette = "OrRd", direction = -1)
#p + scale_colour_distiller(name="Minimum\nStorage\n(1,000 ac-ft)", type = "seq", palette = "YlOrBr", direction = -1)
#p + scale_colour_distiller(name="Minimum\nStorage\n(1,000 ac-ft)", type = "seq", palette = "YlOrRd", direction = -1)
#p + scale_colour_distiller(name="Minimum\nStorage\n(1,000 ac-ft)", type = "seq", palette = "Oranges", direction = -1)
#p + scale_colour_distiller(name="Minimum\nStorage\n(1,000 ac-ft)", type = "seq", palette = "Reds", direction = -1)

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_system_storage_orange.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_system_storage_orange.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_system_storage_orange.svg"),  p, width=6.5, height=5)




p <- ggplot(plot_df, aes(x=dura_months/12, y=min_perc*100, colour=upper_ogden_min/1000))
p <- p + geom_point(size=4)#, aes(shape=data))
p <- p + scale_colour_viridis(direction=-1, option="plasma", name="Minimum\nStorage\n(1,000 ac-ft)")
p <- p + theme_classic_new(12)
p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,1))
p <- p + scale_y_continuous(name="Min Monthly Flow Percentile")
p <- p + theme(legend.position = c(0.85, 0.75))
p <- p + coord_cartesian(xlim = c(0.35,6.5), ylim = c(0,47), expand = TRUE)
p

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_ogden_storage.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_ogden_storage.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_ogden_storage.svg"),  p, width=6.5, height=5)

p <- p + scale_colour_distiller(name="Minimum\nStorage\n(1,000 ac-ft)", type = "seq", palette = "Oranges", direction = -1)

ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_ogden_storage_orange.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_ogden_storage_orange.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_ogden_storage_orange.svg"),  p, width=6.5, height=5)



p <- ggplot(plot_df, aes(x=dura_months/12, y=min_perc*100, colour=upper_weber_min/1000))
p <- p + geom_point(size=4)#, aes(shape=data))
p <- p + scale_colour_viridis(direction=-1, option="plasma", name="Minimum\nStorage\n(1,000 ac-ft)")
p <- p + theme_classic_new(12)
p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,1))
p <- p + scale_y_continuous(name="Min Monthly Flow Percentile")
p <- p + theme(legend.position = c(0.85, 0.75))
p <- p + coord_cartesian(xlim = c(0.35,6.5), ylim = c(0,47), expand = TRUE)
p

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_weber_storage.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_weber_storage.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_weber_storage.svg"),  p, width=6.5, height=5)

p <- p + scale_colour_distiller(name="Minimum\nStorage\n(1,000 ac-ft)", type = "seq", palette = "Oranges", direction = -1)

ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_weber_storage_orange.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_weber_storage_orange.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_upper_weber_storage_orange.svg"),  p, width=6.5, height=5)




p <- ggplot(plot_df, aes(x=dura_months/12, y=min_perc*100, colour=lower_min/1000))
p <- p + geom_point(size=4)#, aes(shape=data))
p <- p + scale_colour_viridis(direction=-1, option="plasma", name="Minimum\nStorage\n(1,000 ac-ft)")
p <- p + theme_classic_new(12)
p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,30,1))
p <- p + scale_y_continuous(name="Min Monthly Flow Percentile")
p <- p + theme(legend.position = c(0.85, 0.75))
p <- p + coord_cartesian(xlim = c(0.35,6.5), ylim = c(0,47), expand = TRUE)
p

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_lower_storage.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_lower_storage.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_lower_storage.svg"),  p, width=6.5, height=5)

p <- p + scale_colour_distiller(name="Minimum\nStorage\n(1,000 ac-ft)", type = "seq", palette = "Oranges", direction = -1)

### Save figures
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_lower_storage_orange.png"),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_lower_storage_orange.pdf"),  p, width=6.5, height=5)
ggsave(file.path(write_output_base_path,"drought_perc_vs_duration_lower_storage_orange.svg"),  p, width=6.5, height=5)



