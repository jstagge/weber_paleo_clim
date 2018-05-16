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

### List of drought responses and data sources
drought_response_list <- c("Base", "ChalkCreekRes", "DemandManagement", "GravelPitRes")
data_list <- c("paleo", "observed", "base", "ww", "hw", "ct", "wd", "hd")

###########################################################################
###  Read in Data
###########################################################################
### Read in reservoir storage
read_location <- file.path(weber_stor_path, "weber_storage.csv")
total_storage <- read.csv(file = read_location)

### Read in reservoir storage trigger points
#weber_triggers <- file.path(weber_stor_path, "mean_weber_stor_levels.csv")
#weber_triggers <- read.csv(file = weber_triggers)

### Weber annual triggers
weber_triggers <- data.frame(Month=6, moderate=380000, severe=340000, extreme=280000)

### Convert weber_triggers month column into factor
weber_triggers$Month <- factor(weber_triggers$Month, levels=c(seq(10,12), seq(1,9)))


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
###  Read in Loop
###########################################################################

### Loop through all scenarios to read in data
for (i in seq(1,dim(scenario_df)[1])){
	cat(paste0(i, " \n"))
	
	### Read in data
	#stor_temp <- read.csv(file = scenario_df$read_location[i])
	stor_temp <- try(read.csv(file = scenario_df$read_location[i]), silent=TRUE)

	if(class(stor_temp) != "try-error"){
	### Cut to only reservoirs and save reservoir names
	stor_res <- stor_temp %>% select(dplyr::contains("RES", ignore.case=FALSE))
	res_names <- names(stor_res)
	
	### Create date column
	if (scenario_df$data[i]=="paleo"){
		date_vec <- stor_temp$Actual.Date
		date_vec <- as.Date(paste0(substr(date_vec,4,9), "-", substr(date_vec,1,2), "-01"))
	} else {
		date_vec <- stor_temp$X
		date_vec <- as.Date(date_vec, format="%m/%d/%Y")
	}
	
	### Combine date with reservoirs and rename columns
	stor_temp <- data.frame(date=date_vec, stor_res)
	names(stor_temp) <- c("date", paste0("res_",seq(1,10)))
	
	### Add data and response columns
	stor_temp$data <- scenario_df$data[i]
	stor_temp$response <- scenario_df$response[i]
	
	### Add column to shift climate change forward in time
	stor_temp$base_date <- stor_temp$date
	if (scenario_df$data[i] != "paleo" & scenario_df$data[i] != "base" & scenario_df$data[i] != "observed"){
		stor_temp$date <- stor_temp$date %m+% years(55)
	}

	#if (scenario_df$data[i] == "observed"){
	#	stor_temp$date <- stor_temp$date %m+% years(-1)
	#}
	
	### Assume Observed based on Base run after 1980-10-01
	#if (scenario_df$data[i] == "observed"){
	#	stor_temp <- stor_temp[year(stor_temp$date) >= 1970,]
	#}			
		
	### Assume Paleo becomes observed in 1904 
	if (scenario_df$data[i] == "paleo"){
		stor_temp <- stor_temp[stor_temp$date < as.Date("1980-10-01"),]
		stor_temp$data[stor_temp$date >= as.Date("1904-10-01")] <- "observed"
	}			

	### Combine data
	if (i == 1){
		stor_all <- stor_temp
	} else {
		stor_all <- rbind(stor_all, stor_temp)
	}
	
	}
	
}

###########################################################################
###  Adjust dates forward by 1 month
###########################################################################
### What is listed as June 1 is actually the end of June, or July 1

stor_all$date <- stor_all$date %m+% months(1)
stor_all$base_date <- stor_all$base_date %m+% months(1)

###########################################################################
###  Prepare data
###########################################################################
### Define order for data and response factors
stor_all$data <- factor(stor_all$data, levels= data_list)
data_levels <- levels(stor_all$data)

stor_all$response <- factor(stor_all$response, levels= drought_response_list)
response_levels <- levels(stor_all$response)

### Cut to records with dates and save month/year
stor_all <- stor_all[!is.na(stor_all$date),]
stor_all$month <- month(stor_all$date)
stor_all$year <- year(stor_all$date)
### Add water year column
stor_all$wy <- usgs_wateryear(year=stor_all$year, month=stor_all$month)
### Make months a factor
stor_all$month <- factor(stor_all$month, levels=c(seq(10, 12), seq(1,9)))

### Remove Observed before 1970
#remove_test <- stor_all$year <= 1970 & stor_all$data=="observed"
#stor_all <- stor_all[!remove_test,]

### Assume Paleo after 1904 is "Observed"
#stor_all$data[stor_all$year >= 1904 & stor_all$data=="paleo"] <- "observed"

### Assume Observed is Paleo after 1904, True observed from 1970 until present
#remove_test <- stor_all$year < 1970 & stor_all$data=="observed"
#stor_all <- stor_all[!remove_test,]

#remove_test <- stor_all$year > 1970 & stor_all$data=="paleo"
#stor_all <- stor_all[!remove_test,]
#stor_all$data[stor_all$year >= 1904 & stor_all$year < 1970 & stor_all$data=="paleo"] <- "observed"

### Re-order by response, data, and then date
stor_all <- stor_all %>% 
	arrange(response, data, date)

###########################################################################
###  Only work with Active Storage
###########################################################################
res_test <- names(stor_all) %in% paste0("res_", seq(1,10))

### Subtract each reservoir by its active bottom
stor_all[,res_test] <- sweep(stor_all[,res_test], 2,total_storage$Active_bottom, "-")


###########################################################################
###  Calculate System Total storage
###########################################################################
### Sum all reservoirs
stor_all$total_res <- stor_all %>%
	select(dplyr::contains("res_")) %>% 
    apply(., 1, sum, na.rm=TRUE )
    
### Sum all current reservoirs (exclude Res 9 and 10)
stor_all$current_res <- stor_all %>%
	select(paste0("res_", seq(1,8))) %>% 
    apply(., 1, sum, na.rm=TRUE )   
 
###########################################################################
###  Calculate Percent Storage
###########################################################################
stor_percent <- stor_all
res_test <- names(stor_percent) %in% paste0("res_", seq(1,10))

### Divide each reservoir by its total storage
stor_percent[,res_test] <- sweep(stor_all[,res_test], 2, total_storage$Active, "/")

### Divide for total system storage
stor_percent$total_res <- stor_percent$total_res / sum(total_storage$Active, na.rm=TRUE)
stor_percent$current_res <- stor_percent$current_res / sum(total_storage$Active[seq(1,8)], na.rm=TRUE)

head(stor_percent)

###########################################################################
###  Calculate Triggers
###########################################################################
#weber_triggers$moderate <- weber_triggers[,2] * 0.7
#weber_triggers$severe <- weber_triggers[,2] * 0.5
#weber_triggers$extreme <- weber_triggers[,2] * 0.25

### Estimate annual trigger for plots based on June storage
mod_annual <- weber_triggers$moderate
sev_annual <- weber_triggers$severe
ext_annual <- weber_triggers$extreme

### Calculate trigger percents based on active storage currently in the system
total_system_stor <- sum(total_storage$Active[seq(1,8)], na.rm=TRUE)
weber_triggers$mod_perc <- weber_triggers$moderate/total_system_stor
weber_triggers$sev_perc <- weber_triggers$severe/total_system_stor
weber_triggers$ext_perc <- weber_triggers$extreme/total_system_stor

### Estimate annual trigger for plots based on June storage
mod_perc_annual <- weber_triggers$mod_perc
sev_perc_annual <- weber_triggers$sev_perc
ext_perc_annual <- weber_triggers$ext_perc

###########################################################################
###  Calculate Regions
###########################################################################
#Upper - Weber
# "RES1.Smith.And.Morehouse.Storage"  
# "RES2.Rockport.Storage" 
#  "RES3.Echo.Storage"
# "RES4.Lost.Creek.Storage"   
# "RES5.East.Canyon.Storage"   second row

#Upper - Ogden
# "RES6.Causey.Storage"  first row
#  "RES7.Pineview.Storage"  

#Lower
#"RES8.Willard.Bay.Storage"
#"RES9.Gravel.Pit.Storage"  
#"RES10.Chalk.Creek.Storage" 

#Page 14 of drought contingency plan for active storage
upper_weber <- paste0("res_",seq(1,5)) 
upper_ogden <- paste0("res_",seq(6,7)) 
lower <- paste0("res_",seq(8,10)) 

### Calculate storage by region
stor_all$upper_weber <- stor_all %>%
	select(upper_weber) %>% 
    apply(., 1, sum, na.rm=TRUE ) 
    
stor_all$upper_ogden <- stor_all %>%
	select(upper_ogden) %>% 
    apply(., 1, sum, na.rm=TRUE )   
    
stor_all$lower <- stor_all %>%
	select(lower) %>% 
    apply(., 1, sum, na.rm=TRUE ) 

### I don't think we need this anymore
#max_stor <- max(stor_all$total_res, na.rm=TRUE)
#stor_all$system_def <- max_stor - stor_all$total_res

### Calculate storage percent by region
stor_percent$upper_ogden <- stor_all$upper_ogden / sum(total_storage$Total[seq(6,7)], na.rm=TRUE)
stor_percent$upper_weber <- stor_all$upper_weber / sum(total_storage$Total[seq(1,5)], na.rm=TRUE)
### Ignore new reservoirs in lower, i.e. calculate percent of current storage available
stor_percent$lower <- stor_all$lower  / sum(total_storage$Total[8], na.rm=TRUE)


###########################################################################
###  Create trigger column and calculate volume below trigger
###########################################################################
### Create an ID column
trigger_df <- data.frame(row_index=seq_len(nrow(stor_all)),  month=stor_all$month)

### Merge with original data
trigger_df <- left_join(trigger_df, weber_triggers, by = c("month" = "Month"))

### Sort and rename
trigger_df <- trigger_df[ with(trigger_df, order(row_index)),]
#names(trigger_df)[3] <- "trigger"
head(trigger_df)

#stor_all$trigger <- trigger_df$trigger

moderate_trigger <- trigger_df$moderate
severe_trigger <- trigger_df$severe
extreme_trigger <- trigger_df$extreme

### Calculate deficits 
stor_all$trigger_deficit <- moderate_trigger - stor_all$total_res
stor_all$trigger_deficit[stor_all$trigger_deficit < 0] <- 0
stor_all$moderate <- moderate_trigger - stor_all$total_res
stor_all$severe <- severe_trigger - stor_all$total_res
stor_all$extreme <- extreme_trigger - stor_all$total_res

### Calculate deficit within each area (so they sum properly) 
### Clear all negative deficits (above threshold) 
### Limit maximum to the difference between thresholds
stor_all <- stor_all %>%
	mutate(moderate = case_when(
		is.na(moderate) ~ NA_real_,
		moderate < 0 ~ 0,
		moderate > (mod_annual - sev_annual) ~ (mod_annual - sev_annual), 
		TRUE ~ moderate)) 

stor_all <- stor_all %>%
	mutate(severe = case_when(
		is.na(severe) ~ NA_real_,
		severe < 0 ~ 0,
		severe > (sev_annual - ext_annual) ~ (sev_annual - ext_annual), 
		TRUE ~ severe)) 
		
stor_all <- stor_all %>%
	mutate(extreme = case_when(
		is.na(extreme) ~ NA_real_,
		extreme < 0 ~ 0, 
		TRUE ~ extreme)) 	

### Create column for categories
stor_all <- stor_all %>%
	mutate(trigger_category = case_when(
		month != 6 ~ NA_character_,
		extreme > 0 ~ "Extreme",
		severe > 0 ~ "Severe",
		moderate > 0 ~ "Moderate",
		TRUE ~ "None")) 
		
### Fill in each year based on the previous June, i.e. triggers remain in effect for a year until next June
stor_all <- stor_all %>%
	fill(trigger_category, .direction="down")		
		
###########################################################################
###  Make data descriptors
###########################################################################
### Add column for temperature
stor_all$temp <- NA
stor_all$temp[stor_all$data == "hd" | stor_all$data == "hw"] <- "Hot"
stor_all$temp[stor_all$data == "wd" | stor_all$data == "ww"] <- "Warm"
stor_all$temp[stor_all$data == "ct"] <- "Median"
stor_all$temp[stor_all$data == "base"] <- "Base"
stor_all$temp <- factor(stor_all$temp, levels= c("Base", "Warm", "Median", "Hot"))

### Add column for precipitation
stor_all$precip <- NA
stor_all$precip[stor_all$data == "hd" | stor_all$data == "wd"] <- "Dry"
stor_all$precip[stor_all$data == "hw" | stor_all$data == "ww"] <- "Wet"
stor_all$precip[stor_all$data == "ct"] <- "Median"
stor_all$precip[stor_all$data == "base"] <- "Base"
stor_all$precip <- factor(stor_all$precip, levels= c("Base", "Wet", "Median", "Dry"))

### Add to Percent
stor_percent$temp <- stor_all$temp 
stor_percent$precip <- stor_all$precip 



###########################################################################
###  Hashimoto Vulnerability Calculation for all combinations
###########################################################################
### Create a dataframe with all possible combinations
### Had to do it this way because observed is not in original scenario_df
all_runs <- expand.grid(response=levels(stor_all$response), data=levels(stor_all$data))

### Remove observed, not Base
#obs_test <- all_runs$data == "observed"
#obs_test <- obs_test == TRUE & all_runs$response != "Base"
#all_runs <- all_runs[!obs_test, ]

### Loop through all possible combinations and run Hashimoto vulnerability
for (i in seq(1,dim(all_runs)[1])){
	### Extract only the correct data and sort
	stor_i <- stor_all %>% 
		filter(data == all_runs$data[i] & response == all_runs$response[i], month == 6) %>%
		dplyr::arrange(-dplyr::desc(date))

	if (dim(stor_i)[1] != 0) {
	### Run Hashimoto function for demand and request storage
	hash_temp_stor <- hash_perform(stor_i, shortage_col="trigger_deficit")
	hash_temp_stor <- data.frame(all_runs[i,], short_type="stor_moderate_trigger", hash_temp_stor)
	
	### Combine the results
	if (i == 1){
		hash_stor <- hash_temp_stor
	} else {
		hash_stor <- rbind(hash_stor, hash_temp_stor)
	}
	}
}

### Catch if there are no failures
hash_stor$resilience[hash_stor$reliability == 1] <- 1

###########################################################################
## Save Workspace
###########################################################################
save.image(file.path(weber_output_path, "weber_storage_output.RData"))





  
###########################################################################
###   Extract just June Storage for Chris
###########################################################################

chris_stor <- stor_all %>% 
	filter(response == "Base" & month == 6) %>%
	select(date, base_date:year, data, res_1:res_10, total_res:trigger_category) %>% 
	arrange(data, year)

write.csv(chris_stor, file.path(weber_output_path, "weber_active_storage_triggers_annual.csv"))





  
###########################################################################
###  Calculate drought events
###########################################################################
### Come back here and calculate based on dates

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
	drought_event_summary$system_min[i] <- min(stor_i$total_res, na.rm=TRUE)
	}
}




