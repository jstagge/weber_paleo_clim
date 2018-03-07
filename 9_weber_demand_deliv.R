
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
### Read in Node Demands
read_location <- file.path(weber_stor_path, "weber_demand.csv")
demand_df <- read.csv(file = read_location)

### Read in reservoir storage
#read_location <- file.path(weber_stor_path, "weber_storage.csv")
#total_storage <- read.csv(file = read_location)

### Read in reservoir storage trigger points
#weber_triggers <- file.path(weber_stor_path, "mean_weber_stor_levels.csv")
#weber_triggers <- read.csv(file = weber_triggers)

### Read in climate change
read_location <- file.path(weber_stor_path, "CMIP5/Base/CMIP5WeberOutput-Base_hist.csv")
deliv_base <- read.csv(file = read_location)
### Cut to only nodes and save node names
deliv_base <- deliv_base[,c(1,seq(14,33))]
node_names <- names(deliv_base)[2:21]
### Process date
deliv_base[,1] <- as.Date(deliv_base[,1], format="%m/%d/%Y")
### Rename columns
names(deliv_base) <- c("date", paste0("SA",seq(1,20)))
deliv_base$data <- "base"
deliv_cc <- deliv_base

for (name in c("hd", "hw", "wd", "ww")) {
### Read in climate change
read_location <- file.path(weber_stor_path, paste0("CMIP5/Base/CMIP5WeberOutput-Base_",name,".csv"))
deliv_temp <- read.csv(file = read_location)
### Cut to only reservoirs and save reservoir names
deliv_temp <- deliv_temp[,c(1,seq(14,33))]
### Process date
deliv_temp[,1] <- as.Date(deliv_temp[,1], format="%m/%d/%Y")
### Rename columns
names(deliv_temp) <- c("date", paste0("SA",seq(1,20)))
deliv_temp$data <- name
deliv_cc <- rbind(deliv_cc, deliv_temp)
}


#### Shift climate change forward in time (55 years)
deliv_cc$date <- deliv_cc$date %m+% years(55)


### Read in paleo flow
read_location <- file.path(weber_stor_path, "Paleo/Base/PaleoWeberOutput-Base.csv")
deliv_paleo <- read.csv(file = read_location)
### Cut to only reservoirs
deliv_paleo <- deliv_paleo[,c(2,seq(15,34))]
### Process date
deliv_paleo[,1] <- as.Date(paste0(substr(deliv_paleo[,1],4,9), "-", substr(deliv_paleo[,1],1,2), "-01"))
### Rename columns
names(deliv_paleo) <- c("date", paste0("SA",seq(1,20)))
deliv_paleo$data <- "paleo"
### Observed starts in 1904
deliv_paleo$data[year(deliv_paleo$date) >= 1904] <- "observed"

###########################################################################
###  Prepare data
###########################################################################
### Combine all data
deliv_all <- rbind(deliv_paleo, deliv_cc)
deliv_all$data <- as.factor(deliv_all$data)
data_levels <- levels(deliv_all$data)

### Cut to records with dates and save month/year
deliv_all <- deliv_all[!is.na(deliv_all$date),]
deliv_all$month <- month(deliv_all$date)
deliv_all$year <- year(deliv_all$date)
### Add water year column
deliv_all$wy <- usgs_wateryear(year=deliv_all$year, month=deliv_all$month)
### Make months a factor
deliv_all$month <- factor(deliv_all$month, levels=c(seq(10, 12), seq(1,9)))

### Create an ID column for sorting
deliv_all <- data.frame(id=seq(1, dim(deliv_all)[1]), deliv_all)

###########################################################################
###  Read in Request Data
###########################################################################
### Read in climate change
read_location <- file.path(weber_stor_path, "CMIP5/Base/CMIP5WeberOutput-Base_hist.csv")
request_base <- read.csv(file = read_location)
### Cut to only nodes and save node names
request_base <- request_base[,c(1,seq(34,53))]
request_names <- names(request_base)[2:21]
### Process date
request_base[,1] <- as.Date(request_base[,1], format="%m/%d/%Y")
### Rename columns
names(request_base) <- c("date", paste0("SA",seq(1,20)))
request_base$data <- "base"
request_cc <- request_base

for (name in c("hd", "hw", "wd", "ww")) {
### Read in climate change
read_location <- file.path(weber_stor_path, paste0("CMIP5/Base/CMIP5WeberOutput-Base_",name,".csv"))
request_temp <- read.csv(file = read_location)
### Cut to only reservoirs and save reservoir names
request_temp <- request_temp[,c(1,seq(34,53))]
### Process date
request_temp[,1] <- as.Date(request_temp[,1], format="%m/%d/%Y")
### Rename columns
names(request_temp) <- c("date", paste0("SA",seq(1,20)))
request_temp$data <- name
request_cc <- rbind(request_cc, request_temp)
}


#### Shift climate change forward in time (55 years)
request_cc$date <- request_cc$date %m+% years(55)


### Read in paleo flow
read_location <- file.path(weber_stor_path, "Paleo/Base/PaleoWeberOutput-Base.csv")
request_paleo <- read.csv(file = read_location)
### Cut to only reservoirs
request_paleo <- request_paleo[,c(2,seq(35,54))]
### Process date
request_paleo[,1] <- as.Date(paste0(substr(request_paleo[,1],4,9), "-", substr(request_paleo[,1],1,2), "-01"))
### Rename columns
names(request_paleo) <- c("date", paste0("SA",seq(1,20)))
request_paleo$data <- "paleo"
### Observed starts in 1904
request_paleo$data[year(request_paleo$date) >= 1904] <- "observed"


###########################################################################
###  Prepare data
###########################################################################
### Combine all data
request_all <- rbind(request_paleo, request_cc)
request_all$data <- as.factor(request_all$data)
data_levels <- levels(request_all$data)

### Cut to records with dates and save month/year
request_all <- request_all[!is.na(request_all$date),]
request_all$month <- month(request_all$date)
request_all$year <- year(request_all$date)
### Add water year column
request_all$wy <- usgs_wateryear(year=request_all$year, month=request_all$month)
### Make months a factor
request_all$month <- factor(request_all$month, levels=c(seq(10, 12), seq(1,9)))

### Create an ID column for sorting
request_all <- data.frame(id=seq(1, dim(request_all)[1]), request_all)

#############################################################
###  Create monthly demand matrix
#############################################################
deliv_melt <- melt(deliv_all, id.vars=c("id", "date", "data", "month", "year", "wy"), value.name="deliv")
request_melt <- melt(request_all, id.vars=c("id", "date", "data", "month", "year", "wy"), value.name="request")
demand_melt <- melt(demand_df, id.vars=c("Month"), value.name="demand")

### Merge demand and delivery
demand_deliv_df <- merge(deliv_melt, request_melt)
demand_deliv_df <- merge(demand_deliv_df, demand_melt, by.x=c("month", "variable"), by.y=c("Month", "variable"))

### Sort by node and then id
demand_deliv_df <- demand_deliv_df[order(demand_deliv_df$variable, demand_deliv_df$id),] 

#############################################################
###  Calculate Shortage
#############################################################
head(demand_deliv_df)
### Create a column for shortage events
demand_deliv_df$shortage_event <- demand_deliv_df$demand > demand_deliv_df$deliv

### Create a column for shortages
demand_deliv_df$shortage <- demand_deliv_df$demand - demand_deliv_df$deliv
demand_deliv_df$shortage[demand_deliv_df$shortage_event == FALSE] <- 0

### Make shortage NA where demand is zero
demand_deliv_df$shortage_event[demand_deliv_df$demand == 0] <- NA
demand_deliv_df$shortage[demand_deliv_df$demand == 0] <- NA

### Calculate percent shortage
demand_deliv_df$shortage_perc <- demand_deliv_df$shortage/demand_deliv_df$demand

#############################################################
###  Quick Check Plots
#############################################################
demand_deliv_plot <- demand_deliv_df[demand_deliv_df$data %in% c("paleo", "observed", "hd"),]

### To plot
p <- ggplot(demand_deliv_plot, aes(x=date, y=shortage))
p <- p + geom_line(aes(group=data), size=0.12)
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (ac-ft)")
p <- p + theme(legend.position="bottom")
p <- p + theme(axis.text.x = element_text(angle = 40, hjust = 1))
p <- p + facet_wrap(~variable, nrow = 4)
p


### Save figures
ggsave(file.path(write_output_base_path,"Shortage_all_facets.png"),  p, width=14, height=6.5, dpi=300)
ggsave(file.path(write_output_base_path,"Shortage_all_facets.pdf"),  p, width=14, height=6.5)
ggsave(file.path(write_output_base_path,"Shortage_all_facets.svg"),  p, width=14, height=6.5)


### To plot
p <- ggplot(demand_deliv_plot, aes(x=date, y=shortage_perc))
p <- p + geom_line(aes(group=data), size=0.12)
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shorage (% of Demand)", breaks=seq(0,1,0.1), labels=percent)
p <- p + theme(legend.position="bottom")
p <- p + facet_wrap(~variable, nrow = 4)
p


### Save figures
ggsave(file.path(write_output_base_path,"Shortage_perc_all_facets.png"),  p, width=14, height=6.5, dpi=300)
ggsave(file.path(write_output_base_path,"Shortage_perc_all_facets.pdf"),  p, width=14, height=6.5)
ggsave(file.path(write_output_base_path,"Shortage_perc_all_facets.svg"),  p, width=14, height=6.5)


### To plot
p <- ggplot(demand_deliv_plot, aes(x=date))
p <- p + geom_line(aes(group=data, y=demand), size=0.12, colour="grey20")
p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
#p <- p + geom_line(aes(group=data, y=request), size=0.12, colour="red")
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
p <- p + facet_wrap(~variable, nrow = 4)
p



#############################################################
###  Check SA1
#############################################################
sa1_plot <- subset(demand_deliv_plot, variable=="SA1")

### To plot
p <- ggplot(sa1_plot, aes(x=date))
p <- p + geom_line(aes(group=data, y=demand), size=0.12, colour="grey20")
p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p




#############################################################
###  Check SA12
#############################################################
sa1_plot <- subset(demand_deliv_plot, variable=="SA12")

### To plot
p <- ggplot(sa1_plot, aes(x=date))
p <- p + geom_line(aes(group=data, y=demand), size=0.12, colour="grey20")
p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p


### To plot
p <- ggplot(sa1_plot, aes(x=month, group=year))
p <- p + geom_line(aes(y=shortage), size=0.12, colour="grey20")
p <- p + theme_classic_new()
#p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p



#############################################################
###  Check SA14
#############################################################
sa1_plot <- subset(demand_deliv_plot, variable=="SA14")

### To plot
p <- ggplot(sa1_plot, aes(x=date))
p <- p + geom_line(aes(group=data, y=demand), size=0.12, colour="grey20")
p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p


### To plot
p <- ggplot(sa1_plot, aes(x=month, group=year))
p <- p + geom_line(aes(y=shortage), size=0.12, colour="grey20")
p <- p + theme_classic_new()
#p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p




#############################################################
###  Check SA10
#############################################################
sa10_plot <- subset(demand_deliv_plot, variable=="SA10")
#sa10_plot <- sa10_plot[,c(4, 5, 8, 10)] #9, 
#sa10_plot <- melt(sa10_plot, c("date", "data"))

### To plot
p <- ggplot(sa10_plot, aes(x=date))
p <- p + geom_line(aes(group=data, y=demand), size=0.12, colour="grey20")
p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Delivery (Ac-ft)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p

### Save figures
ggsave(file.path(write_output_base_path,"SA10_deliv_demand.png"),  p, width=8, height=4.5, dpi=300)
ggsave(file.path(write_output_base_path,"SA10_deliv_demand.pdf"),  p, width=8, height=4.5)
ggsave(file.path(write_output_base_path,"SA10_deliv_demand.svg"),  p, width=8, height=4.5)



sa10_plot <- subset(demand_deliv_plot, variable=="SA10")
sa10_plot <- sa10_plot[,c(4, 5, 8, 10)] #9, 
sa10_plot <- melt(sa10_plot, c("date", "data"))
### Sort by node and then id
sa10_plot <- sa10_plot[order(sa10_plot$variable, sa10_plot$data ,sa10_plot$date),] 

### To plot
p <- ggplot(sa10_plot, aes(x=date))
p <- p + geom_line(aes(group=data, colour=variable, y=value), size=0.12)
#p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
p <- p + theme_classic_new()
p <- p + scale_colour_manual(values=c("black", "red"))
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Delivery (ac-ft)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p

### Save figures
ggsave(file.path(write_output_base_path,"Shortage_all_facets.png"),  p, width=14, height=6.5, dpi=300)
ggsave(file.path(write_output_base_path,"Shortage_all_facets.pdf"),  p, width=14, height=6.5)
ggsave(file.path(write_output_base_path,"Shortage_all_facets.svg"),  p, width=14, height=6.5)


#############################################################
###  Plot requests
#############################################################
### It looks like there is a significant difference between paleo and climate change requests at SA 1 and SA10

p <- ggplot(request_all, aes(x=month, group=wy))
p <- p + geom_line(aes(y=SA1, colour=data), size=0.12)
p <- p + theme_classic_new()
#p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Request")
p <- p + theme(legend.position="bottom")
p <- p + facet_wrap(~data, nrow = 4)
p


### Make months a factor
demand_df$month_adj_forplot <- c(seq(4, 12), seq(1,3))

### Make months a factor
demand_df$month_adj_forplot <- c(seq(4, 12), seq(1,3))

p <- ggplot(request_all, aes(x=month))
p <- p + geom_boxplot(aes(y=SA1, colour=data))
#p <- p + geom_line(data=demand_df, aes(x=month_adj_forplot, y=SA1), colour="black")
p <- p + theme_classic_new()
#p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Request")
p <- p + theme(legend.position="bottom")
p

### Funky - April requests are always exactly the same for entire Paleo run.
### Climate change runs vary across scenarios, but each climate change run has the same number, but different for each scenario.


#############################################################
###  Do calculations regionally
#############################################################



###########################################################################
###  Calculate Percent Storage
###########################################################################
stor_percent <- stor_all
res_test <- names(stor_percent) %in% paste0("res_", seq(1,10))

### Divide each reservoir by its total storage
stor_percent[,res_test] <- sweep(stor_percent[,res_test], 2, c(total_storage$Total, NA, NA), "/")


head(stor_percent)


###########################################################################
###  Calculate Triggers
###########################################################################
weber_triggers$moderate <- weber_triggers[,2] * 0.7
weber_triggers$severe <- weber_triggers[,2] * 0.5
weber_triggers$extreme <- weber_triggers[,2] * 0.25

### Calculate trigger percents based on storage
total_system_stor <- sum(total_storage$Total, na.rm=TRUE)
weber_triggers$mod_perc <- weber_triggers$moderate/total_system_stor
weber_triggers$sev_perc <- weber_triggers$severe/total_system_stor
weber_triggers$ext_perc <- weber_triggers$extreme/total_system_stor

### Calculate annual trigger percents based on storage
mod_perc_annual <- mean(weber_triggers$mod_perc)
sev_perc_annual <- mean(weber_triggers$sev_perc)
ext_perc_annual <- mean(weber_triggers$ext_perc)



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

###
stor_all$upper_ogden <- apply(stor_all[,seq(2,6)],1,sum, na.rm=TRUE)
stor_all$upper_weber <- apply(stor_all[,seq(7,8)],1,sum, na.rm=TRUE)
stor_all$lower <- apply(stor_all[,seq(9,11)],1,sum, na.rm=TRUE)
stor_all$system <- stor_all$upper_ogden + stor_all$upper_weber + stor_all$lower

max_stor <- max(stor_all$system, na.rm=TRUE)
stor_all$system_def <- max_stor - stor_all$system


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


###########################################################################
###  Test plots
###########################################################################
p <- ggplot(stor_all, aes(x=month))
p <- p + geom_line(aes(y=system, group=wy, colour=data), size=0.3)
p <- p + theme_classic_new()
p


p <- ggplot(stor_all, aes(x=date, fill=data))
p <- p + geom_area(aes(y=system), size=0.3, position="identity")
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


###########################################################################
## Plot Storage Time Series
###########################################################################
area_df <- stor_all[stor_all$data %in% c("paleo", "observed", "hd"),]

p <- ggplot(area_df, aes(x=date, y=system/1000, fill=data))
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
## Plot Storage percent Time Series
###########################################################################
### Create percentile plot dataframe
perc_plot <- stor_percent[stor_all$data %in% c("paleo", "observed", "hd"),]

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


### To plot
p <- ggplot(perc_plot, aes(x=date, y=res_1))
p <- p + geom_area(data=trigger_plot, aes(y=res_perc, fill=trigger_level))
p <- p + geom_line(aes(y=res_1, group=data))
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
area_region <- area_df[, names(area_df) %in% c("date", "data", "upper_ogden", "upper_weber", "lower")]
area_region <- melt(area_region, id.vars=c("date", "data"))#, measure.vars=c("upper_ogden", "upper_weber", "lower"))
area_region$data <- factor(area_region$data)
area_region$variable <- factor(area_region$variable)

p <- ggplot(subset(area_region, data=="observed"), aes(x=date, y=value/1000, fill=variable))
#p <- ggplot(yup, aes(x=date, y=value/1000, fill=variable))
p <- p + geom_area()
p <- p + geom_area(data=subset(area_region, data=="paleo"))
p <- p + geom_area(data=subset(area_region, data=="hd"))
p <- p + theme_classic_new()
#p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "grey30", cc_colors), labels=c("Observed", "Base", "HD", "HW", "WD", "WW", "Reconst" ))
p <- p + scale_fill_manual(name="Region", values= res_colors, labels=c("Upper Ogden", "Upper Weber", "Lower Weber"), guide = guide_legend())
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

p <- ggplot(area_df, aes(x=date, y=system/1000))
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
line_df$value <- line_df$system

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
p <- p + coord_cartesian(ylim=c(0,0.4), expand=FALSE)
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





