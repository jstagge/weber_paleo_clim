
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
###  Hashimoto Vulnerability Function
###########################################################################
hash_perform <- function(site, demand,  delivery){
### Follows Hashimoto, T., Stedinger, J.R., Loucks, D.P., 1982. Reliability, resiliency, and vulnerability criteria for water resource system performance evaluation. Water Resour. Res. 18, 14–20. https://doi.org/10.1029/WR018i001p00014
### Also used https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2002WR001778
demand <- demand[,names(demand) %in% site]
delivery <- delivery[,names(delivery) %in% site]

shortage <- demand - delivery

sat <- shortage <= 0
unsat <- shortage > 0

sat_lag <- c(sat[seq(2, length(sat))], NA)
w_trans <- unsat == TRUE & sat_lag == TRUE

time_length <- sum(!is.na(shortage))

### Reliability
### Fraction of time that demand is met
reliability <- sum(sat, na.rm=TRUE)/time_length

### Resilience
### Average probability of recovery in any given failure time step
resilience <- sum(w_trans, na.rm=TRUE)/ (time_length - sum(sat, na.rm=TRUE))

### Vulnerability
### Maximum shortage
vulnerability <- max(shortage, na.rm=TRUE)

return(list(reliability=reliability, resilience=resilience, vulnerability = vulnerability))
}


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
###  Check System Shortages
#############################################################
deliv_wide <-    dcast(demand_deliv_df, date + data ~ variable, value.var="deliv", fun.aggregate = sum, na.rm=TRUE)
demand_wide <-    dcast(demand_deliv_df, date + data ~ variable, value.var="demand", fun.aggregate = sum, na.rm=TRUE)
shortage_wide <-    dcast(demand_deliv_df, date + data ~ variable, value.var="shortage", fun.aggregate = sum, na.rm=TRUE)

shortage_system <-    dcast(demand_deliv_df, date + data ~ ., value.var="shortage", fun.aggregate = sum, na.rm=TRUE)

### Add column for system
deliv_wide$system <- apply(deliv_wide[,!names(deliv_wide) %in% c("date", "data")], 1, sum, na.rm=TRUE)
demand_wide$system <- apply(demand_wide[,!names(demand_wide) %in% c("date", "data")], 1, sum, na.rm=TRUE)
shortage_wide$system <- apply(shortage_wide[,!names(shortage_wide) %in% c("date", "data")], 1, sum, na.rm=TRUE)



#############################################################
###  Run Hashimoto calculations
#############################################################
site_names <- c(paste0("SA", seq(1,20)), "system")
data_names <- levels(demand_wide$data)

### Loop through all data sites and calculate Hashimoto vulnerabilities
for (i in seq(1,length(data_names))){
	data_i <- data_names[[i]]
	
	hash_temp <- sapply(site_names, hash_perform, demand=demand_wide[demand_wide$data==data_i,], delivery=deliv_wide[deliv_wide$data==data_i,])
	hash_temp <- t(hash_temp)
	
	hash_temp <- data.frame(data=data_i, site=site_names, hash_temp)

	if (i == 1){
		hash_df <- hash_temp
	} else {
		hash_df <- rbind(hash_df, hash_temp)
	}
}

### Reorganize results
hash_df$reliability <- as.numeric(hash_df$reliability)
hash_df$resilience <- as.numeric(hash_df$resilience)
hash_df$vulnerability <- as.numeric(hash_df$vulnerability)
hash_df$data <- factor(hash_df$data, levels=c("paleo", "observed", "base", "ww", "wd", "hw", "hd"))
hash_df$site <- factor(hash_df$site, levels=c(paste0("SA", seq(1,20)),"system"))




###########################################################################
## Save Workspace
###########################################################################
save.image(file.path(weber_output_path, "weber_delivery_output.RData"))









