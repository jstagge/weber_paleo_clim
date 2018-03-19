
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
read_location <- file.path(write_output_base_path, paste0(site_id,"_climchange_perc_ts.csv"))
flow_cc <- read.csv(file = read_location)

### Read in paleo flow
#read_location <- file.path(write_output_base_path, paste0(site_id,"_paleo_perc_ts.csv"))
#flow_paleo <- read.csv(file = read_location)

### Combine all data
#flow_all <- rbind(flow_obs, flow_cc, flow_paleo)
#flow_all$data <- as.factor(flow_all$data)
#flow_all$date <- as.Date(flow_all$date)

#data_levels <- levels(flow_all$data)

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
#p <- p + scale_colour_manual(name="Scenario", values= cc_colors, labels=c("Hot-Dry", "Hot-Wet", "Warm-Dry", "Warm-Wet"))
p <- p + scale_colour_manual(name="Scenario", values= cc_colors[c(4,2,3,1)], labels=c("Warm-Wet", "Hot-Wet", "Warm-Dry", "Hot-Dry"), limits=c("WWN5", "HWN5", "WDN5", "HDN5"))
p <- p + scale_linetype_manual(name=NULL, values=c("longdash"))
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








###########################################################################
###  Could calculate paired stats
###########################################################################

x_data <- flow_cc$flow_m3s[flow_cc$data=="base" & flow_cc$month==12] 
y_data <- flow_cc$flow_m3s[flow_cc$data=="HWN5" & flow_cc$month==12] 

ks.test(x_data,y_data)

wilcox.test(x_data, y_data, paired=TRUE)



x_data <- flow_obs$flow_m3s[flow_obs$data=="observed" & flow_obs$month==3] 
y_data <- flow_paleo$flow_m3s[flow_paleo$data=="paleo" & flow_paleo$month==3] 

ks.test(x_data,y_data)

wilcox.test(x_data, y_data, paired=FALSE)




