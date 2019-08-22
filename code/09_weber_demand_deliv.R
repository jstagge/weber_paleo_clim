# *------------------------------------------------------------------
# | PROGRAM NAME: 09_weber_demand_deliv
# | FILE NAME: 09_weber_demand_deliv.R
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


### Load project specific functions
file.sources = list.files(function_path, pattern="*.R", recursive=TRUE)
sapply(file.path(function_path, file.sources),source)

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
cc_colors <- c("#0072B2", "#56B4E9", "#E69F00" , "#D55E00")
#, "#CC79A7"

res_colors <- cb_pal(pal="wong", 3, sort=FALSE)

#9E67AB  - pinky, purply, about the same, a little darker
#CC79A7  - pink, from the correct palette, looks pretty good
#F0E442   - too light
#009E73   - green, another option

### Conversion from m3/s to acft/month
m3s_to_acftmonth <- 2131.97
m3s_to_m3month <- 2629746
m3_to_acft <- 1233.4818

### List of drought responses and data sources
drought_response_list <- c("Base", "ChalkCreekRes", "DemandManagement", "GravelPitRes")
data_list <- c("paleo", "observed", "base", "ww", "hw", "ct", "wd", "hd")
node_levels <- c(paste0("SA", seq(1,20)), "upper_weber", "upper_ogden", "lower", "system")

###########################################################################
###  Read in Demand Data
###########################################################################
### Read in Node Demands
read_location <- file.path(weber_stor_path, "weber_demand.csv")
demand_df <- read.csv(file = read_location)
### Convert demand month column into factor
demand_df$Month <- factor(demand_df$Month, levels=c(seq(10,12), seq(1,9)))


###########################################################################
###  Read in Drought Events
###########################################################################
read_location <- file.path(write_output_base_path, paste0(site_id,"_drought_details.csv"))

drought_event_summary <- read.csv(file = read_location)



###########################################################################
###  Prepare all scenarios for read in
###########################################################################
scenario_df <- expand.grid(response=drought_response_list, data=data_list)

### Create a column for base folders
scenario_df <- scenario_df %>%
       mutate(base_folder = case_when(
           data == "paleo" ~ "Paleo",
           TRUE ~ "CMIP5") )

### Create a column for response folder
scenario_df <- scenario_df %>%
		mutate(response_folder = case_when(
  			response == "Base"  ~ file.path(weber_stor_path, paste0(base_folder,"/",response)),
  			TRUE ~ file.path(weber_stor_path, paste0(base_folder,"/DroughtResponse/",response)))
  		)
  		
### Create a column for file location
scenario_df <- scenario_df %>%
		mutate(read_location = case_when(
  			data == "paleo"  ~ file.path(response_folder, paste0(base_folder,"WeberOutput-",response,".csv")),
  			data == "observed" ~ file.path(response_folder, paste0(base_folder,"WeberOutput-",response,"_hist.csv")),
  			data == "base" ~ file.path(response_folder, paste0(base_folder,"WeberOutput-",response,"_hist.csv")),
  			TRUE ~ file.path(response_folder, paste0(base_folder,"WeberOutput-",response,"_",data,".csv")))
  		)
				

###########################################################################
###  For now only consider base scenario
###########################################################################



###########################################################################
###  Read in Loop for delivery
###########################################################################

### Loop through all scenarios to read in data
for (i in seq(1,dim(scenario_df)[1])){
	cat(paste0(i, " \n"))

	### Read in data
	deliv_temp <- try(read.csv(file = scenario_df$read_location[i]), silent=TRUE)

	if(class(deliv_temp) != "try-error"){

	### Cut to only delivery nodes
	deliv_and_request <- deliv_temp %>% select(dplyr::contains("SA"))

	### Separate delivery and request
	deliv_nodes <- deliv_and_request %>% select(-dplyr::contains("Requested"))
	request_nodes <- deliv_and_request %>% select(dplyr::contains("Requested"))

	### Save node names
	node_names <- names(deliv_nodes)
	
	### Create date column
	if (scenario_df$data[i]=="paleo"){
		date_vec <- deliv_temp$Actual.Date
		date_vec <- as.Date(paste0(substr(date_vec,4,9), "-", substr(date_vec,1,2), "-01"))
	} else {
		date_vec <- deliv_temp$X
		date_vec <- as.Date(date_vec, format="%m/%d/%Y")
	}
	
	### Combine date with nodes  and rename columns
	deliv_temp <- data.frame(date=date_vec, deliv_nodes, variable="delivery")
	request_temp <- data.frame(date=date_vec, request_nodes, variable="request")

	names(deliv_temp) <- c("date", paste0("SA",seq(1,20)), "variable")
	names(request_temp) <-	names(deliv_temp) 
	
	### Combine deliver and requests and rename columns
	deliv_request_temp <- rbind(deliv_temp, request_temp)

	### Add data and response columns
	deliv_request_temp$data <- scenario_df$data[i]
	deliv_request_temp$response <- scenario_df$response[i]
		
	### Add column to shift climate change forward in time
	deliv_request_temp$base_date <- deliv_request_temp$date

	### Add column to shift climate change forward in time
	deliv_request_temp$base_date <- deliv_request_temp$date
	if (scenario_df$data[i] != "paleo" & scenario_df$data[i] != "base" & scenario_df$data[i] != "observed"){
		deliv_request_temp$date <- deliv_request_temp$date %m+% years(55)
	}

	### Assume Paleo becomes observed in 1904 
	if (scenario_df$data[i] == "paleo"){
		deliv_request_temp <- deliv_request_temp[deliv_request_temp$date < as.Date("1980-10-01"),]
		deliv_request_temp$data[deliv_request_temp$date >= as.Date("1904-10-01")] <- "observed"
	}	
	
	### Combine data
	if (i == 1){
		deliv_request_all <- deliv_request_temp
	} else {
		deliv_request_all <- rbind(deliv_request_all, deliv_request_temp)
	}
	
	}
	
}

###########################################################################
###  No need to adjust dates forward by 1 month
###########################################################################
### What is listed as June 1 is what was delivered from June 1 until June 30.

###########################################################################
###  Prepare data
###########################################################################
### Define order for data and response factors
deliv_request_all$data <- factor(deliv_request_all$data, levels= data_list)
data_levels <- levels(deliv_request_all$data)

deliv_request_all$response <- factor(deliv_request_all$response, levels= drought_response_list)
response_levels <- levels(deliv_request_all$response)

### Cut to records with dates and save month/year
deliv_request_all <- deliv_request_all[!is.na(deliv_request_all$date),]
deliv_request_all$month <- month(deliv_request_all$date)
deliv_request_all$year <- year(deliv_request_all$date)
### Add water year column
deliv_request_all$wy <- usgs_wateryear(year=deliv_request_all$year, month=deliv_request_all$month)
### Make months a factor
deliv_request_all$month <- factor(deliv_request_all$month, levels=c(seq(10, 12), seq(1,9)))


#############################################################
###  Create monthly demand matrix
#############################################################
### Join the demands to the matrix
demand_temp <- deliv_request_all %>% filter(variable == "delivery" ) %>%
	left_join(., demand_df, by=c("month" = "Month"), suffix = c("_old", "_new"))

### Identify old and new columns
old_columns <- stringr::str_detect(names(demand_temp), "_old")
new_columns <- stringr::str_detect(names(demand_temp), "_new")

### Replace old columns with new
demand_temp[,old_columns] <- demand_temp[,new_columns]
### Remove new columns
demand_temp <- demand_temp[,!new_columns]
### Change names
names(demand_temp) <- names(deliv_request_all)

### Swap variable name
demand_temp$variable <- "demand"

### Merge with delivery and request
deliv_request_all <- rbind(deliv_request_all, demand_temp)

#############################################################
###  Calculate System and Regional demand/deliv
#############################################################
### Don't consider Weber-Provo Diversion Canal (SA1) as part of the Weber System
### Define System delivery, demand, request, shortage based on SA2-20

### Sum all nodes from SA2 - SA20
### Any NA returns NA to prevent mistaken comparisons
deliv_request_all$system <- deliv_request_all %>%
	select(dplyr::num_range("SA", 2:20)) %>% 
    apply(., 1, sum)
    
#Upper - Weber
#### Remove SA1 Weber - Provo Diversion Canal                                       
#"Oakley.Diversion.SA2.Oakley.to.Wanship.Diversion"                                                  
#"Wanship.Diversion.SA3.Wanship.to.Echo.Diversion"                                                   
#"Echo.to.Devils.Slide.Diversion.SA4.Echo.to.Devils.Slide.Diversion"                                 
#"Lost.Creek.Diversion.SA5.Lost.Creek.Diversion"                                                     
#"Devils.Slide.to.Stoddard.Diversion.SA6.Devils.Slide.to.Stodard.Diversion"                          
#"Park.City.Diversion.SA7.Park.City.Diversion"                                                       
#"East.Canyon.Diversion.SA8.East.Canyon.Diversion"                                                   
#"Stoddard.To.Gateway.Diversion.SA9.Stoddard.To.Gateway.Diversion"                                   
#"Gateway.Canal.Diversion.SA10.Gateway.Canal.Diversion"                                              
#"Davis.Weber.Canal.Diversion.SA11.Davis.Weber.Canal.Diversion"  

#Upper - Ogden
# Weber.Basin.Project.Ogden.Valley.SA12.Weber.Basin.Project.Ogden.Valley.Diversion
# Ogden.Brigham.and.S.Ogden.Highline.Canals.SA13.Ogden.Brigham.and.S.Ogden.Highline.Canals.Diversion
#  Ogden.River.Below.Pineview.Diversion.SA14.Ogden.River.Below.Pineview.Diversion

#Lower - Weber
#"Slaterville.Diversion.SA15.Slaterville.Diversion"                                                  
#"Warren.Canal.Diversion.SA16.Warren.Canal.Diversion"                                                
#"Ogden.Bay.Bird.Refuge.Diversion.SA17.Ogden.Bay.Bird.Refuge.Diversion"                              
#"GSL.Minerals.Diversion.SA18.GSL.Minerals.Diversion"                                                
#"Gateway.To.Slaterville.Diversion.SA19.Gateway.To.Slaterville.Diversion"                            
#"Additional.WB.Demand.Diversion.SA20.Additional.WB.Demand.Diversion"    

#Page 14 of drought contingency plan for active storage
system_nodes <- paste0("SA",seq(2,20)) 
upper_weber <- paste0("SA",seq(2,11)) 
upper_ogden <- paste0("SA",seq(12,14)) 
lower <- paste0("SA",seq(15,20)) 

### SA20 is all NAs for Base scenario, have to come back to this ultimately
### In the meantime, use na.rm to sum

### Calculate storage by region
deliv_request_all$system <- deliv_request_all %>%
	select(system_nodes) %>% 
    apply(., 1, sum, na.rm=TRUE) 

deliv_request_all$upper_weber <- deliv_request_all %>%
	select(upper_weber) %>% 
    apply(., 1, sum, na.rm=TRUE) 

deliv_request_all$upper_ogden <- deliv_request_all %>%
	select(upper_ogden) %>% 
    apply(., 1, sum, na.rm=TRUE)  
        
deliv_request_all$lower <- deliv_request_all %>%
	select(lower) %>% 
    apply(., 1, sum, na.rm=TRUE) 

###########################################################################
###  Generate Full Node Name List
###########################################################################
node_names_full <- strsplit(node_names, ".", fixed=TRUE)

for (j in seq(1, length(node_names_full))){
	name_j <- node_names_full[[j]]
	sa_location <- which(stringr::str_detect(name_j, "SA"))
	name_j <- name_j[seq(sa_location, length(name_j) - 1)]
	name_j[1] <- paste0(name_j[1], ":")
	node_names_full[[j]] <- paste(name_j, collapse=" ")
}

node_names_full <- unlist(node_names_full)

node_names_full <- c(node_names_full, "Upper Weber", "Upper Ogden", "Lower Weber", "Total System")

#############################################################
###  Calculate Shortage
#############################################################
### Create long data form
 deliv_request_long <-  deliv_request_all %>% 
 	gather(node, flow, -date, -variable, -data, -response, -base_date, -month, -year, -wy)
### Define order for nodes
 deliv_request_long$node <- factor(deliv_request_long$node, levels=node_levels)#, labels=node_names_full)

### Create format with variables separated  
 demand_deliv_df <-  deliv_request_long %>% spread(variable, flow)
 
### Create a column for shortage events
demand_deliv_df$demand_shortage_event <- demand_deliv_df$demand > demand_deliv_df$delivery
demand_deliv_df$request_shortage_event <- demand_deliv_df$request > demand_deliv_df$delivery

### Create a column for shortages
demand_deliv_df$demand_shortage <- demand_deliv_df$demand - demand_deliv_df$delivery
demand_deliv_df$demand_shortage[demand_deliv_df$demand_shortage_event == FALSE] <- 0

demand_deliv_df$request_shortage <- demand_deliv_df$request - demand_deliv_df$delivery
demand_deliv_df$request_shortage[demand_deliv_df$request_shortage_event == FALSE] <- 0

### Make shortage NA where demand is zero
demand_deliv_df$demand_shortage_event[demand_deliv_df$demand == 0] <- NA
demand_deliv_df$demand_shortage[demand_deliv_df$demand == 0] <- NA

demand_deliv_df$request_shortage_event[demand_deliv_df$request == 0] <- NA
demand_deliv_df$request_shortage[demand_deliv_df$request == 0] <- NA

### Calculate percent shortage
demand_deliv_df$demand_shortage_perc <- demand_deliv_df$demand_shortage/demand_deliv_df$demand
demand_deliv_df$request_shortage_perc <- demand_deliv_df$request_shortage/demand_deliv_df$request



###########################################################################
###  Hashimoto Vulnerability Calculation for all combinations
###########################################################################
### Create a dataframe with all possible combinations
### Had to do it this way because observed is not in original scenario_df
all_runs <- expand.grid(response=levels(demand_deliv_df$response), data=levels(demand_deliv_df$data), node=node_levels)

### Loop through all possible combinations and run Hashimoto vulnerability
for (i in seq(1,dim(all_runs)[1])){
	### Extract only the correct data and sort
	demand_i <- demand_deliv_df %>% 
		filter(data == all_runs$data[i] & response == all_runs$response[i] & node == all_runs$node[i]) %>%
		dplyr::arrange(-dplyr::desc(date))

	### Run Hashimoto function for demand and request shortages
	hash_temp_demand <- hash_perform(demand_i, shortage_col="demand_shortage")
	hash_temp_demand <- data.frame(all_runs[i,], short_type="demand", hash_temp_demand)
	
	hash_temp_request <- hash_perform(demand_i, shortage_col="request_shortage")
	hash_temp_request <- data.frame(all_runs[i,], short_type="request", hash_temp_request)
	
	### Combine the results
	if (i == 1){
		hash_full <- rbind(hash_temp_demand, hash_temp_request)
	} else {
		hash_full <- rbind(hash_full, hash_temp_demand, hash_temp_request)
	}
}



		
###########################################################################
## Save Workspace
###########################################################################
save.image(file.path(weber_output_path, "weber_delivery_output.RData"))


