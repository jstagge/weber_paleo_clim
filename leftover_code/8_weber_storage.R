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

res_colors <- cb_pal(pal="wong", 3, sort=FALSE)

#9E67AB  - pinky, purply, about the same, a little darker
#CC79A7  - pink, from the correct palette, looks pretty good
#F0E442   - too light
#009E73   - green, another option

### Conversion from m3/s to acft/month
m3s_to_acftmonth <- 2131.97

m3s_to_m3month <- 2629746

m3_to_acft <- 1233.4818


###########################################################################
###  Read in Data
###########################################################################
### Read in reservoir storage
read_location <- file.path(weber_stor_path, "weber_storage.csv")
total_storage <- read.csv(file = read_location)

### Read in reservoir storage trigger points
weber_triggers <- file.path(weber_stor_path, "mean_weber_stor_levels.csv")
weber_triggers <- read.csv(file = weber_triggers)


### Read in climate change
read_location <- file.path(weber_stor_path, "CMIP5/Base/CMIP5WeberOutput-Base_hist.csv")
stor_base <- read.csv(file = read_location)
### Cut to only reservoirs and save reservoir names
stor_base <- stor_base[,seq(1,11)]
res_names <- names(stor_base)[2:11]
### Process date
stor_base[,1] <- as.Date(stor_base[,1], format="%m/%d/%Y")
### Rename columns
names(stor_base) <- c("date", paste0("res_",seq(1,10)))
stor_base$data <- "base"
stor_cc <- stor_base

for (name in c("hd", "hw", "wd", "ww")) {
### Read in climate change
read_location <- file.path(weber_stor_path, paste0("CMIP5/Base/CMIP5WeberOutput-Base_",name,".csv"))
stor_temp <- read.csv(file = read_location)
### Cut to only reservoirs and save reservoir names
stor_temp <- stor_temp[,seq(1,11)]
### Process date
stor_temp[,1] <- as.Date(stor_temp[,1], format="%m/%d/%Y")
### Rename columns
names(stor_temp) <- c("date", paste0("res_",seq(1,10)))
stor_temp$data <- name
stor_cc <- rbind(stor_cc, stor_temp)
}


#### Shift climate change forward in time (55 years)
stor_cc$base_date <- stor_cc$date
stor_cc$date <- stor_cc$date %m+% years(55)


### Read in paleo flow
read_location <- file.path(weber_stor_path, "Paleo/Base/PaleoWeberOutput-Base.csv")
stor_paleo <- read.csv(file = read_location)
### Cut to only reservoirs
stor_paleo <- stor_paleo[,seq(2,12)]
### Process date
stor_paleo[,1] <- as.Date(paste0(substr(stor_paleo[,1],4,9), "-", substr(stor_paleo[,1],1,2), "-01"))
### Rename columns
names(stor_paleo) <- c("date", paste0("res_",seq(1,10)))
stor_paleo$data <- "paleo"
### Observed starts in 1904
stor_paleo$data[year(stor_paleo$date) >= 1904] <- "observed"




###########################################################################
###  Read in Drought Events
###########################################################################
read_location <- file.path(write_output_base_path, paste0(site_id,"_drought_details.csv"))

drought_event_summary <- read.csv(file = read_location)


###########################################################################
###  Prepare data
###########################################################################
### Combine all data
stor_paleo$base_date <- stor_paleo$date
stor_all <- rbind(stor_paleo, stor_cc)
stor_all$data <- as.factor(stor_all$data)
data_levels <- levels(stor_all$data)

### Cut to records with dates and save month/year
stor_all <- stor_all[!is.na(stor_all$date),]
stor_all$month <- month(stor_all$date)
stor_all$year <- year(stor_all$date)
### Add water year column
stor_all$wy <- usgs_wateryear(year=stor_all$year, month=stor_all$month)
### Make months a factor
stor_all$month <- factor(stor_all$month, levels=c(seq(10, 12), seq(1,9)))

### Sum all reservoirs
stor_all$total_res <- apply(stor_all[,2:11],1,sum, na.rm=TRUE)

###########################################################################
###  Calculate Percent Storage
###########################################################################
stor_percent <- stor_all
res_test <- names(stor_percent) %in% paste0("res_", seq(1,10))

### Divide each reservoir by its total storage
stor_percent[,res_test] <- sweep(stor_percent[,res_test], 2, c(total_storage$Total, NA, NA), "/")

### Divide for total system storage
stor_percent$total_res <- stor_percent$total_res / sum(total_storage$Total, na.rm=TRUE)

head(stor_percent)


###########################################################################
###  Calculate Triggers
###########################################################################
weber_triggers$moderate <- weber_triggers[,2] * 0.7
weber_triggers$severe <- weber_triggers[,2] * 0.5
weber_triggers$extreme <- weber_triggers[,2] * 0.25

### Estimate annual trigger for plots based on May storage
#mod_annual <- mean(weber_triggers$moderate)
#sev_annual <- mean(weber_triggers$severe)
#ext_annual <- mean(weber_triggers$extreme)
mod_annual <- weber_triggers$moderate[5]
sev_annual <- weber_triggers$severe[5]
ext_annual <- weber_triggers$extreme[5]

### Calculate trigger percents based on storage
total_system_stor <- sum(total_storage$Total, na.rm=TRUE)
weber_triggers$mod_perc <- weber_triggers$moderate/total_system_stor
weber_triggers$sev_perc <- weber_triggers$severe/total_system_stor
weber_triggers$ext_perc <- weber_triggers$extreme/total_system_stor

### Estimate annual trigger for plots based on May storage
mod_perc_annual <- weber_triggers$mod_perc[5]
sev_perc_annual <- weber_triggers$sev_perc[5]
ext_perc_annual <- weber_triggers$ext_perc[5]



###########################################################################
###  Calculate Regions
###########################################################################

#Upper - Ogden
# "RES1.Smith.And.Morehouse.Storage"  
# "RES2.Rockport.Storage" 
#  "RES3.Echo.Storage"
# "RES4.Lost.Creek.Storage"   
# "RES5.East.Canyon.Storage"   second row

#Upper - Weber
# "RES6.Causey.Storage"  first row
#  "RES7.Pineview.Storage"  

#Lower
#"RES8.Willard.Bay.Storage"
#"RES9.Gravel.Pit.Storage"  
#"RES10.Chalk.Creek.Storage" 

#Page 14 of drought contingency plan for active storage
upper_ogden <- paste0("res_",seq(1,5)) 
upper_weber <- paste0("res_",seq(6,7)) 
lower <- paste0("res_",seq(8,10)) 

### Calculate storage by region
stor_all$upper_ogden <- apply(stor_all[,seq(2,6)],1,sum, na.rm=TRUE)
stor_all$upper_weber <- apply(stor_all[,seq(7,8)],1,sum, na.rm=TRUE)
stor_all$lower <- apply(stor_all[,seq(9,11)],1,sum, na.rm=TRUE)
#stor_all$system <- stor_all$upper_ogden + stor_all$upper_weber + stor_all$lower

max_stor <- max(stor_all$total_res, na.rm=TRUE)
stor_all$system_def <- max_stor - stor_all$total_res


### Calculate storage percent by region
stor_percent$upper_ogden <- stor_all$upper_ogden / sum(total_storage$Total[seq(1,5)], na.rm=TRUE)
stor_percent$upper_weber <- stor_all$upper_weber / sum(total_storage$Total[seq(6,7)], na.rm=TRUE)
stor_percent$lower <- stor_all$lower  / sum(total_storage$Total[seq(8,10)], na.rm=TRUE)


###########################################################################
###  Create trigger column and calculate volume below trigger
###########################################################################
trigger_df <- data.frame(row_index=as.numeric(row.names(stor_all)),  month=stor_all$month)

### Merge with original data
trigger_df <- merge(trigger_df, weber_triggers, by.x="month", by.y="Month")
### Sort and rename
trigger_df <- trigger_df[ with(trigger_df, order(row_index)),]
names(trigger_df)[3] <- "trigger"
head(trigger_df)

stor_all$trigger <- trigger_df$trigger

### Moderate 50-70% of mean
### Severe 25-50% of mean
### Extreme 0-25% of mean

moderate_trigger <- trigger_df$trigger * 0.7
severe_trigger <- trigger_df$trigger * 0.5
extreme_trigger <- trigger_df$trigger * 0.25


stor_all$moderate <- moderate_trigger - stor_all$system
stor_all$severe <- severe_trigger - stor_all$system
stor_all$extreme <- extreme_trigger - stor_all$system

### Clear all negative deficits (above threshold) 
stor_all$moderate[stor_all$moderate < 0] <- 0
stor_all$severe[stor_all$severe < 0] <- 0
stor_all$extreme[stor_all$extreme < 0] <- 0

### Clear all negative deficits (above threshold) 
stor_all$moderate[stor_all$moderate > (moderate_trigger -severe_trigger)] <- (moderate_trigger - severe_trigger)[stor_all$moderate > (moderate_trigger -severe_trigger)]
stor_all$severe[stor_all$severe > (severe_trigger - extreme_trigger)] <- (severe_trigger - extreme_trigger)[stor_all$severe > (severe_trigger - extreme_trigger)] 


#### Plot triggers
plot_triggers <- weber_triggers
plot_triggers$total <-  sum(total_storage$Total)
names(plot_triggers)[2] <- "Mean 2013-2017"
plot_triggers <- melt(plot_triggers[,c(seq(1,5), 9)], "Month")
names(plot_triggers)[3] <- "storage"
plot_triggers$perc <- plot_triggers$storage / sum(total_storage$Total)

p <- ggplot(plot_triggers, aes(x=Month, y=storage/1000, colour=variable))
#p <- ggplot(yup, aes(x=date, y=value/1000, fill=variable))
p <- p + geom_line(size=0.8)
#p <- p + geom_line(data=trigger_plot, aes(y=res_stor/1000, group=trigger_level, fill=NA), colour="grey30", linetype="longdash", size=0.8)
p <- p + theme_classic_new(11)
p <- p + scale_colour_manual(name="", values= c( "#ffeda0", "#feb24c", "#f03b20", "#8da0cb", "grey50"), limits= c( "moderate", "severe", "extreme", "Mean 2013-2017", "total"), labels=c("Moderate Trigger", "Severe Trigger", "Extreme Trigger", "Mean 2013-2017", "Full Storage"), guide = guide_legend(nrow=2,byrow=TRUE))
p <- p + scale_x_continuous(name="Month", breaks=seq(1,12,1))
p <- p + scale_y_continuous(name="Total System Storage (1,000 ac-ft)", breaks=seq(0, 600, 100))
p <- p + coord_cartesian(xlim=c(1,12), ylim=c(0,sum(total_storage$Total/1000)*1.1), expand=FALSE)
p <- p + theme(legend.position="bottom")
p

### Save figures
ggsave(file.path(write_output_base_path,"trigger_levels_perc.png"),  p, width=4.5, height=4, dpi=600)
ggsave(file.path(write_output_base_path,"trigger_levels_perc.pdf"),  p, width=4.5, height=4)
ggsave(file.path(write_output_base_path,"trigger_levels_perc.svg"),  p, width=4.5, height=4)


p <- ggplot(plot_triggers, aes(x=Month, y=perc, colour=variable))
#p <- ggplot(yup, aes(x=date, y=value/1000, fill=variable))
p <- p + geom_line(size=0.8)
#p <- p + geom_line(data=trigger_plot, aes(y=res_stor/1000, group=trigger_level, fill=NA), colour="grey30", linetype="longdash", size=0.8)
p <- p + theme_classic_new(11)
p <- p + scale_colour_manual(name="", values= c( "#ffeda0", "#feb24c", "#f03b20", "#8da0cb", "grey50"), limits= c( "moderate", "severe", "extreme", "Mean 2013-2017", "total"), labels=c("Moderate Trigger", "Severe Trigger", "Extreme Trigger", "Mean 2013-2017", "Full Storage"), guide = guide_legend(nrow=2,byrow=TRUE))
p <- p + scale_x_continuous(name="Month", breaks=seq(1,12,1))
#p <- p + scale_y_continuous(name="Total System Storage (%)", breaks=seq(0, 1, 0.1), labels=percent)
p <- p + scale_y_continuous(name="Total System Storage (%)", breaks= seq(0, 1, 0.1), labels = c("0", "", "20%", "", "40%", "", "60%", "", "80%", "", "100%"))
p <- p + coord_cartesian(xlim=c(1,12), ylim=c(0,1.02), expand=FALSE)
p <- p + theme(legend.position="bottom")
p

### Save figures
ggsave(file.path(write_output_base_path,"trigger_levels_perc.png"),  p, width=4.5, height=4, dpi=600)
ggsave(file.path(write_output_base_path,"trigger_levels_perc.pdf"),  p, width=4.5, height=4)
ggsave(file.path(write_output_base_path,"trigger_levels_perc.svg"),  p, width=4.5, height=4)

###########################################################################
###  Test plots
###########################################################################
p <- ggplot(stor_all, aes(x=month))
p <- p + geom_line(aes(y=total_res, group=wy, colour=data), size=0.3)
p <- p + theme_classic_new()
p


p <- ggplot(stor_all, aes(x=date, fill=data))
p <- p + geom_area(aes(y=total_res), size=0.3, position="identity")
p <- p + theme_classic_new()
p

### Check plot
p <- ggplot(subset(stor_percent, data!="base") , aes(colour=data))
p <- p + geom_line(aes(x=date, y=total_res), size=0.3)
p <- p + geom_line(data=subset(stor_percent, data=="base"), aes(x=base_date, y=total_res), size=0.3)
p <- p + theme_classic_new()
p


###########################################################################
###  Make data descriptors
###########################################################################
### Add column for temperature
stor_all$temp <- NA
stor_all$temp[stor_all$data == "hd" | stor_all$data == "hw"] <- "Hot"
stor_all$temp[stor_all$data == "wd" | stor_all$data == "ww"] <- "Warm"
stor_all$temp[stor_all$data == "base"] <- "Base"
stor_all$temp <- factor(stor_all$temp, levels= c("Base", "Warm", "Hot"))

### Add column for precipitation
stor_all$precip <- NA
stor_all$precip[stor_all$data == "hd" | stor_all$data == "wd"] <- "Dry"
stor_all$precip[stor_all$data == "hw" | stor_all$data == "ww"] <- "Wet"
stor_all$precip[stor_all$data == "base"] <- "Base"
stor_all$precip <- factor(stor_all$precip, levels= c("Base", "Wet", "Dry"))

### Add to Percent
stor_percent$temp <- stor_all$temp 
stor_percent$precip <- stor_all$precip 

###########################################################################
## Plot Storage Time Series as Area
###########################################################################
area_df <- stor_all[stor_all$data %in% c("paleo", "observed", "hd"),]

p <- ggplot(area_df, aes(x=date, y=total_res/1000, fill=data))
p <- p + geom_area( position = "identity", alpha=0.8)
p <- p + geom_hline(yintercept=0, size=0.2)
#p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
#p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "grey30", cc_colors), labels=c("Observed", "Base", "HD", "HW", "WD", "WW", "Reconst" ))
p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "#CC79A7", "#D55E00"), labels=c( "Observed", "Reconstr", "Climate Change (HD)"), limits=c("observed", "paleo", "hd"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1425-01-01"), as.Date("2070-01-01")), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="System Storage (1,000 ac-ft)", breaks=seq(0,800,100))
p <- p + theme(legend.position="bottom")
p

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_full.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_full.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_full.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_1900.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_1900.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_1900.svg"),  p, width=8, height=3.5)

###########################################################################
## Plot Storage Time Series as Line
###########################################################################
line_df <- stor_all[stor_all$data %in% c("paleo", "observed", "hd"),]
line_df2 <- stor_all[stor_all$data %in% c("base"),]
line_df2$date <- line_df2$base_date
line_df <- rbind(line_df, line_df2)

p <- ggplot(line_df, aes(x=date, y=total_res/1000, colour=data))
p <- p + geom_line()#, alpha=0.8)
p <- p + geom_hline(yintercept=0, size=0.2)
#p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
#p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "grey30", cc_colors), labels=c("Observed", "Base", "HD", "HW", "WD", "WW", "Reconst" ))
p <- p + scale_colour_manual(name="Scenario", values= c("grey30", "#1f78b4", "#CC79A7", "#D55E00"), labels=c( "Observed", "Base", "Reconstr", "Climate Change (HD)"), limits=c("observed", "base", "paleo", "hd"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1425-01-01"), as.Date("2070-01-01")), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="System Storage (1,000 ac-ft)", breaks=seq(0,800,100))
p <- p + theme(legend.position="bottom")
p

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_full.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_full.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_full.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_1900.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_1900.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_area_acft_hd_1900.svg"),  p, width=8, height=3.5)


###########################################################################
## Plot Storage percent Time Series
###########################################################################
### Create percentile plot dataframe
perc_plot <- stor_percent[stor_all$data %in% c("base"),]
perc_plot$date <- perc_plot$base_date

perc_temp <- stor_percent[stor_all$data %in% c("paleo", "observed", "hd"),]
perc_test <- perc_temp$date < min(perc_plot$date) | perc_temp$date > max(perc_plot$date)
perc_plot <- rbind(perc_temp[perc_test, ], perc_plot)

### Create dataframe for triggers
trigger_plot <- data.frame(date = rep(c(min(perc_plot$date), max(perc_plot$date)), 3))
trigger_plot$month <- month(trigger_plot$date)
trigger_plot$year <- year(trigger_plot$date)

### Add water year column
trigger_plot$wy <- usgs_wateryear(year=trigger_plot$year, month=trigger_plot$month)

### Add triggers
trigger_plot$trigger_level <- rep(c("Moderate", "Severe", "Extreme"), each=2)
trigger_plot$trigger_level <- factor(trigger_plot$trigger_level, levels=c("Moderate", "Severe", "Extreme"))
trigger_plot$res_perc <- rep(c(mod_perc_annual-sev_perc_annual, sev_perc_annual-ext_perc_annual, ext_perc_annual), each=2)
trigger_plot$res_stor <- rep(c(mod_annual, sev_annual, ext_annual), each=2)



### To plot
p <- ggplot(perc_plot, aes(x=date, y=total_res))
p <- p + geom_area(data=trigger_plot, aes(y=res_perc, fill=trigger_level))
p <- p + geom_line(aes(group=data), size=0.12)
#p <- p + geom_line(data=trigger_plot, aes(y=mod_perc_annual), col="red")
#p <- p + geom_line(data=trigger_plot, aes(y=ext_perc_annual), col="green")
#p <- p + geom_line(data=trigger_plot, aes(y=sev_perc_annual), col="blue")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Storage Trigger", values= c("#f03b20",  "#feb24c", "#ffeda0"), limits= c("Extreme", "Severe", "Moderate"), labels=c("Extreme", "Severe", "Moderate"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), ylim=c(0,1), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Percent Storage", breaks=seq(0,1,0.1), labels=percent)
p <- p + theme(legend.position="bottom")
p


### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900.svg"),  p, width=8, height=3.5)




### To plot
p <- ggplot(perc_plot, aes(x=date, y=upper_ogden))
p <- p + geom_area(data=trigger_plot, aes(y=res_perc, fill=trigger_level))
p <- p + geom_line(aes(group=data), size=0.12)
#p <- p + geom_line(data=trigger_plot, aes(y=mod_perc_annual), col="red")
#p <- p + geom_line(data=trigger_plot, aes(y=ext_perc_annual), col="green")
#p <- p + geom_line(data=trigger_plot, aes(y=sev_perc_annual), col="blue")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Storage Trigger", values= c("#f03b20",  "#feb24c", "#ffeda0"), limits= c("Extreme", "Severe", "Moderate"), labels=c("Extreme", "Severe", "Moderate"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), ylim=c(0,1), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Percent Storage", breaks=seq(0,1,0.1), labels=percent)
p <- p + theme(legend.position="bottom")
p


### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full_upper_ogden.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full_upper_ogden.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full_upper_ogden.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900_upper_ogden.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900_upper_ogden.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900_upper_ogden.svg"),  p, width=8, height=3.5)




### To plot
p <- ggplot(perc_plot, aes(x=date, y=upper_weber))
p <- p + geom_area(data=trigger_plot, aes(y=res_perc, fill=trigger_level))
p <- p + geom_line(aes(group=data), size=0.12)
#p <- p + geom_line(data=trigger_plot, aes(y=mod_perc_annual), col="red")
#p <- p + geom_line(data=trigger_plot, aes(y=ext_perc_annual), col="green")
#p <- p + geom_line(data=trigger_plot, aes(y=sev_perc_annual), col="blue")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Storage Trigger", values= c("#f03b20",  "#feb24c", "#ffeda0"), limits= c("Extreme", "Severe", "Moderate"), labels=c("Extreme", "Severe", "Moderate"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), ylim=c(0,1), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Percent Storage", breaks=seq(0,1,0.1), labels=percent)
p <- p + theme(legend.position="bottom")
p


### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full_upper_weber.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full_upper_weber.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full_upper_weber.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900_upper_weber.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900_upper_weber.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900_upper_weber.svg"),  p, width=8, height=3.5)


### To plot
p <- ggplot(perc_plot, aes(x=date, y=lower))
p <- p + geom_area(data=trigger_plot, aes(y=res_perc, fill=trigger_level))
p <- p + geom_line(aes(group=data), size=0.12)
#p <- p + geom_line(data=trigger_plot, aes(y=mod_perc_annual), col="red")
#p <- p + geom_line(data=trigger_plot, aes(y=ext_perc_annual), col="green")
#p <- p + geom_line(data=trigger_plot, aes(y=sev_perc_annual), col="blue")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Storage Trigger", values= c("#f03b20",  "#feb24c", "#ffeda0"), limits= c("Extreme", "Severe", "Moderate"), labels=c("Extreme", "Severe", "Moderate"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), ylim=c(0,1), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Percent Storage", breaks=seq(0,1,0.1), labels=percent)
p <- p + theme(legend.position="bottom")
p


### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full_lower.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full_lower.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_full_lower.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900_lower.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900_lower.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_line_perc_hd_1900_lower.svg"),  p, width=8, height=3.5)





names(perc_plot)[2:9] <- as.character(total_storage$Name)
yup <-  melt(perc_plot, id.vars = c("date", "month", "year", "wy", "data"))
res_test <- yup$variable %in% total_storage$Name
yup <- yup[res_test,]
 
### To plot
p <- ggplot(yup, aes(x=date, y=value))
p <- p + geom_area(data=trigger_plot, aes(y=res_perc, fill=trigger_level), alpha=0.85)
p <- p + geom_line(aes(group=data), size=0.12)
#p <- p + geom_line(data=trigger_plot, aes(y=mod_perc_annual), col="red")
#p <- p + geom_line(data=trigger_plot, aes(y=ext_perc_annual), col="green")
#p <- p + geom_line(data=trigger_plot, aes(y=sev_perc_annual), col="blue")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Storage Trigger", values= c("#f03b20",  "#feb24c", "#ffeda0"), limits= c("Extreme", "Severe", "Moderate"), labels=c("Extreme", "Severe", "Moderate"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), ylim=c(0,1), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Percent of Total Storage", breaks=seq(0,1,0.1), labels=percent)
p <- p + theme(legend.position="bottom")
p <- p + facet_wrap(~variable, nrow = 4)
p


### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_perc.png"),  p, width=11, height=9, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_perc.pdf"),  p, width=11, height=9)




###########################################################################
## Plot Shortage Time Series by region
###########################################################################
area_df <- stor_all[stor_all$data %in% c("paleo", "observed", "base", "hd"),]
area_region <- area_df[, names(area_df) %in% c("date", "base_date", "data", "upper_ogden", "upper_weber", "lower")]
area_region <- melt(area_region, id.vars=c("date", "base_date", "data"))#, measure.vars=c("upper_ogden", "upper_weber", "lower"))
area_region$data <- factor(area_region$data)
area_region$variable <- factor(area_region$variable, levels=c("upper_weber", "upper_ogden", "lower"))

p <- ggplot(subset(area_region, data=="observed"), aes(x=date, y=value/1000, fill=variable))
#p <- ggplot(yup, aes(x=date, y=value/1000, fill=variable))
p <- p + geom_area()
p <- p + geom_area(data=subset(area_region, data=="paleo" & date < as.Date("1980-10-01")))
p <- p + geom_area(data=subset(area_region, data=="hd"))
p <- p + geom_area(data=subset(area_region, data=="base"), aes(x=base_date))
p <- p + geom_line(data=trigger_plot, aes(y=res_stor/1000, group=trigger_level, fill=NA), colour="grey30", linetype="longdash", size=0.5)
p <- p + theme_classic_new()
#p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "grey30", cc_colors), labels=c("Observed", "Base", "HD", "HW", "WD", "WW", "Reconst" ))
p <- p + scale_fill_manual(name="Region", values= res_colors, labels=c("Upper Weber", "Upper Ogden", "Lower Weber"), breaks=c("upper_weber", "upper_ogden", "lower"), guide = guide_legend())
p <- p + coord_cartesian(xlim=c(as.Date("1425-01-01"), as.Date("2070-01-01")), expand=FALSE)
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="System Storage (1,000 ac-ft)", breaks=seq(0,800,100))
p <- p + theme(legend.position="bottom")
p


### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_byregion_acft_hd_full.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_byregion_acft_hd_full.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_byregion_acft_hd_full.svg"),  p, width=8, height=3.5)


### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_stor_byregion_acft_hd_1900.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_stor_byregion_acft_hd_1900.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_stor_byregion_acft_hd_1900.svg"),  p, width=8, height=3.5)




###########################################################################
## Plot Shortage Time Series
###########################################################################
plot_test <- stor_all$data %in% c("observed", "paleo")

cc_test <- stor_all$data %in% c("base", "hd", "hw", "wd", "ww")
cc_df <- stor_all[cc_test,]

#plot_df <- rbind(stor_all[plot_test,], cc_df[cc_df$data %in% c("base", "HDN5"),])

base_df <- cc_df[cc_df$data %in% c("base"),]
base_df$data <- factor("observed", levels="observed")

area_df <- rbind(cc_df[cc_df$data %in% c("hd", "hw"),], stor_all[stor_all$data %in% c("paleo", "observed"),])

p <- ggplot(area_df, aes(x=date, y=-system_def/1000, fill=data))
p <- p + geom_area( position = "identity", alpha=0.8)
p <- p + geom_hline(yintercept=0, size=0.2)
p <- p + geom_hline(yintercept=-max_stor/1000, size=0.4, colour="black", linetype="longdash")
p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
#p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "grey30", cc_colors), labels=c("Observed", "Base", "HD", "HW", "WD", "WW", "Reconst" ))
p <- p + scale_fill_manual(name="Scenario", values= c("#D55E00", "#56B4E9", "grey30", "#CC79A7"), labels=c("HD", "HW", "Observed", "Reconstr"), guide = guide_legend())
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="System Storage Deficit (1,000 ac-ft)", breaks=seq(-10000,500,50))
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
p <- p + geom_hline(yintercept=-max_stor/1000, size=0.4, colour="black", linetype="longdash")
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

plot_test <- stor_all$data %in% c( "hd", "hw", "wd", "ww")

area_df <- stor_all[plot_test,]
area_df$data <- factor(area_df$data, levels=c( "ww", "hw", "wd", "hd"), labels=c("Warm-Wet", "Hot-Wet", "Warm-Dry", "Hot-Dry"))

line_df <- stor_all[stor_all$data %in% "base", ]
line_df <- do.call("rbind", replicate(4, line_df, simplify = FALSE))
line_df$data <- area_df$data
line_df$precip <- area_df$precip
line_df$temp <- area_df$temp

p <- ggplot(area_df, aes(x=date, y=total_res/1000))
p <- p + geom_area(aes(fill=data))
p <- p + geom_hline(yintercept=0)
p <- p + geom_line(data=line_df, colour="black")#, linetype="longdash")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Scenario", values= cc_colors[c(4,2,3,1)], labels=c("WW", "HW", "WD", "HD"), breaks=c("ww", "hw", "wd", "hd"))
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="5 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="System Storage (1,000 ac-ft)", breaks=seq(-5000,5000,100))
p <- p + facet_grid(data ~ .)
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_vert.png"),  p, width=7, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_vert.pdf"),  p, width=7, height=6)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_vert.svg"),  p, width=7, height=6)


p <- p + facet_grid(precip ~ temp)

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_square.png"),  p, width=8, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_square.pdf"),  p, width=8, height=6)
ggsave(file.path(write_output_base_path,"clim_change_storage_area_acft_square.svg"),  p, width=8, height=6)




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
p <- p + geom_line(data=line_df, colour="black")#, linetype="longdash")
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


###########################################################################
###  Calculate drought events
###########################################################################

drought_event_summary <- subset(drought_event_summary, data!="CTN5")

drought_event_summary$data <- as.character(drought_event_summary$data)
drought_event_summary$data[drought_event_summary$data == "WDN5"] <- "wd"
drought_event_summary$data[drought_event_summary$data == "WWN5"] <- "ww"
drought_event_summary$data[drought_event_summary$data == "HDN5"] <- "hd"
drought_event_summary$data[drought_event_summary$data == "HWN5"] <- "hw"
#drought_event_summary$data[drought_event_summary$data == "WDN5"] <= "wd"

#### Shift climate change forward in time (55 years)
drought_event_summary$begin <- as.Date(drought_event_summary$begin)
drought_event_summary$end <- as.Date(drought_event_summary$end)

cc_test <- drought_event_summary$data %in% c("wd", "ww", "hd", "hw")
drought_event_summary$begin[cc_test] <- drought_event_summary$begin[cc_test] %m+% years(55)
drought_event_summary$end[cc_test] <- drought_event_summary$end[cc_test] %m+% years(55)

drought_event_summary$upper_ogden_min <- NA
drought_event_summary$upper_weber_min <- NA
drought_event_summary$lower_min <- NA
drought_event_summary$system_min <- NA

for (i in seq(1, dim(drought_event_summary)[1])){

	event_i <- drought_event_summary[i,]
	begin_i <- as.Date(event_i$begin)
	end_i <- as.Date(event_i$end)
	
	data_i <- as.character(event_i$data)
	
	stor_i <- subset(stor_all, date >= begin_i-60 & date <=end_i+60 & data==data_i)
	
	if (length(stor_i$system)>0){
	drought_event_summary$upper_ogden_min[i] <- min(stor_i$upper_ogden, na.rm=TRUE)
	drought_event_summary$upper_weber_min[i] <- min(stor_i$upper_weber, na.rm=TRUE)
	drought_event_summary$lower_min[i] <- min(stor_i$lower, na.rm=TRUE)
	drought_event_summary$system_min[i] <- min(stor_i$system, na.rm=TRUE)
	}
}

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




inferno
viridis ### no 
plasma
magma### no

###########################################################################
###  Plot against triggers
###########################################################################


plot_triggers <- merge(plot_df, weber_triggers, by.x="month", by.y="Month")

plot_df$trigger <- "None"
plot_df$trigger[plot_df$system_min < weber_triggers[8,2] * 0.25] <- "Extreme"
plot_df$trigger[plot_df$system_min < weber_triggers[8,2] * 0.5 & plot_df$trigger == "None"] <- "Severe"
plot_df$trigger[plot_df$system_min < weber_triggers[8,2] * 0.7 & plot_df$trigger == "None"] <- "Moderate"




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
















"green", "yellow", "red"
"#4daf4a", "#ffff33", "#e41a1c"

"#000004FF" "#330A5FFF" "#781C6DFF" "#BB3754FF" "#ED6925FF" "#FCB519FF"
[7] "#FCFFA4FF"

"100-01-01", "1903-12-01"
viridis(6, option="inferno")

p + scale_fill_viridis(option="inferno", discrete=TRUE)





