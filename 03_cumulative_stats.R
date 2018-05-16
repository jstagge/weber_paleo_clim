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
require(paleoAPR)
#require(data.table)
require(ggrepel)

### Load project specific functions
file.sources = list.files(function_path, pattern="*.R", recursive=TRUE)
sapply(file.path(function_path, file.sources),source)

### Load global functions
file.sources = list.files(global_path, pattern="*.R", recursive=TRUE)
sapply(file.path(global_path, file.sources),source)


### Colors to use for climate change scenarios
cc_colors <- c("#D55E00" ,"#56B4E9", "#E69F00" , "#0072B2", "#CC79A7")

#9E67AB  - pinky, purply, about the same, a little darker
#CC79A7  - pink, from the correct palette, looks pretty good
#F0E442   - too light
#009E73   - green, another option

### Conversion from m3/s to acft/month
m3s_to_acftmonth <- 2131.97

m3s_to_m3month <- 2629746

m3_to_acft <- 1233.4818

###########################################################################
## Set Initial Values
###########################################################################
### Set site data
site_id_list <- c("10128500")
site_name_list <- c("Weber River")

thresh_level <- 0.5
gap_months <- 2
min_duration <- 3

###########################################################################
###  
###########################################################################
n <- 1


site_id <- site_id_list[n]
site_name <- site_name_list[n]

### Create output folders


###########################################################################
###  Read in Data
###########################################################################
### Read in observed flow
read_location <- file.path(write_output_base_path, paste0(site_id,"_obs_perc_ts.csv"))
flow_obs <- read.csv(file = read_location)

### Read in climate change
read_location <- file.path(write_output_base_path, paste0(site_id,"_climchange_perc_ts.csv"))
flow_cc <- read.csv(file = read_location)

### Read in paleo flow
read_location <- file.path(write_output_base_path, paste0(site_id,"_paleo_perc_ts.csv"))
flow_paleo <- read.csv(file = read_location)

### Combine all data
flow_all <- rbind(flow_obs, flow_cc, flow_paleo)
flow_all$data <- as.factor(flow_all$data)
flow_all$date <- as.Date(flow_all$date)

data_levels <- levels(flow_all$data)



###########################################################################
## Calculate cumulative loss
###########################################################################

two_year_limit <- 2*sum(flow_all$thresh_m3s[1:12])

### Set up cumulative loss column
flow_all$cum_loss <- NA
	
for (i in data_levels){
	### Extract data
	i_test <- flow_all$data == i
	flow_i <- flow_all[i_test, ]

	### Loop through each month
	for (j in 2:dim(flow_i)[1]) {

		loss_j <- flow_i$flow_m3s[j] - flow_i$thresh_m3s[j]

		if (!is.na(loss_j)) {
			### If previously is NA and there is loss, start a new drought
			if (is.na(flow_i$cum_loss[j-1])){
				if (loss_j < 0) {
					flow_i$cum_loss[j] <- loss_j
					flow_i$cum_loss[j-1] <- NA				
				}
			### If previously in drought, does it end or continue
			} else {
			### Cumulative loss is previous plus current loss
				cum_loss_j <- flow_i$cum_loss[j-1] + loss_j

				## If current cumulative loss is less than zero, add to column
				if (cum_loss_j < 0) {
					flow_i$cum_loss[j] <- cum_loss_j
				} else if (flow_i$cum_loss[j-1] == 0){
					flow_i$cum_loss[j] <- NA
				} else {
					flow_i$cum_loss[j] <- 0
				}
		}}
		
		if(flow_i$cum_loss[j] < -two_year_limit && is.na(flow_i$cum_loss[j])!=TRUE) {flow_i$cum_loss[j] <- -two_year_limit}
	}
	### Re-insert cumulative loss into the main dataframe
	flow_all$cum_loss[i_test] <- flow_i$cum_loss
}



###########################################################################
## Save Drought Statistics
###########################################################################
save_file <- file.path(write_output_base_path, paste0(site_id,"_cum_loss_ts.csv"))

write.csv(flow_all, save_file, row.names = FALSE)




###########################################################################
###  Make data descriptors
###########################################################################
### Add column for temperature
flow_all$temp <- NA
flow_all$temp[flow_all$data == "HDN5" | flow_all$data == "HWN5"] <- "Hot"
flow_all$temp[flow_all$data == "WDN5" | flow_all$data == "WWN5"] <- "Warm"
flow_all$temp[flow_all$data == "base"] <- "Base"
flow_all$temp <- factor(flow_all$temp, levels= c("Base", "Warm", "Hot"))

### Add column for precipitation
flow_all$precip <- NA
flow_all$precip[flow_all$data == "HDN5" | flow_all$data == "WDN5"] <- "Dry"
flow_all$precip[flow_all$data == "HWN5" | flow_all$data == "WWN5"] <- "Wet"
flow_all$precip[flow_all$data == "base"] <- "Base"
flow_all$precip <- factor(flow_all$precip, levels= c("Base", "Wet", "Dry"))


flow_all$cum_loss[is.na(flow_all$cum_loss)] <- 0

###########################################################################
## Plot Time Series
###########################################################################
plot_test <- flow_all$data %in% c("observed", "paleo")

cc_test <- flow_all$data %in% c("base", "HDN5", "HWN5", "WDN5", "WWN5")
cc_df <- flow_all[cc_test,]
cc_df$date <- cc_df$date + (as.Date("2050-01-01") - mean(cc_df$date))


#plot_df <- rbind(flow_all[plot_test,], cc_df[cc_df$data %in% c("base", "HDN5"),])

base_df <- rbind(flow_all[flow_all$data %in% c("observed"),], cc_df[cc_df$data %in% c("base"),])
base_df$data <- factor("observed", levels="observed")

area_df <- rbind(cc_df[cc_df$data %in% c("HDN5", "HWN5"),], flow_all[flow_all$data %in% c("paleo"),])

p <- ggplot(area_df, aes(x=date, y=cum_loss * m3s_to_acftmonth / 1000, fill=data))
p <- p + geom_area( position = "identity", alpha=0.8)
p <- p + geom_hline(yintercept=0, size=0.2)
p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
#p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "grey30", cc_colors), labels=c("Observed", "Base", "HD", "HW", "WD", "WW", "Reconst" ))
p <- p + scale_fill_manual(name="", values= c("#D55E00", "#56B4E9", "grey30", "#CC79A7"), labels=c("Hot-Dry ", "Hot-Wet ", "Observed", "Reconstr"), guide = guide_legend())
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Cumulative Deficit (1,000 ac-ft)", breaks=seq(-500,500,50))
p <- p + theme(legend.position="bottom")
p

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_hot_full.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_hot_full.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_hot_full.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_hot_1900.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_hot_1900.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_hot_1900.svg"),  p, width=8, height=3.5)



plot_test <- flow_all$data %in% c("observed", "paleo")

cc_test <- flow_all$data %in% c("base", "HDN5", "HWN5", "WDN5", "WWN5")
cc_df <- flow_all[cc_test,]
cc_df$date <- cc_df$date + (as.Date("2050-01-01") - mean(cc_df$date))


#plot_df <- rbind(flow_all[plot_test,], cc_df[cc_df$data %in% c("base", "HDN5"),])

base_df <- rbind(flow_all[flow_all$data %in% c("observed"),], cc_df[cc_df$data %in% c("base"),])
base_df$data <- factor("observed", levels="observed")

area_df <- rbind(cc_df[cc_df$data %in% c("WDN5", "WWN5"),], flow_all[flow_all$data %in% c("paleo"),])

p <- ggplot(area_df, aes(x=date, y=cum_loss * m3s_to_acftmonth / 1000, fill=data))
p <- p + geom_area( position = "identity", alpha=0.8)
p <- p + geom_hline(yintercept=0, size=0.2)
p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
#p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "grey30", cc_colors), labels=c("Observed", "Base", "HD", "HW", "WD", "WW", "Reconst" ))
p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "#CC79A7", "#E69F00", "#0072B2"), labels=c("Observed", "Reconstr", "Warm-Dry ", "Warm-Wet "), guide = guide_legend())
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Cumulative Deficit (1,000 ac-ft)", breaks=seq(-500,500,50))
p <- p + theme(legend.position="bottom")
p

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_warm_full.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_warm_full.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_warm_full.svg"),  p, width=8, height=3.5)

### Cut to 1900s
p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2066-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="10 years"), date_labels = "%Y")

### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_warm_1900.png"),  p, width=8, height=3.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_warm_1900.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_cum_deficit_area_acft_warm_1900.svg"),  p, width=8, height=3.5)


###########################################################################
## Plot Area Time Series for climate change
###########################################################################

plot_test <- flow_all$data %in% c( "HDN5", "HWN5", "WDN5", "WWN5")

area_df <- flow_all[plot_test,]
area_df$data <- factor(area_df$data, levels=c( "WWN5", "HWN5", "WDN5", "HDN5"))

line_df <- flow_all[flow_all$data %in% "base", ]
line_df <- do.call("rbind", replicate(4, line_df, simplify = FALSE))
line_df$data <- area_df$data
line_df$precip <- area_df$precip
line_df$temp <- area_df$temp

p <- ggplot(area_df, aes(x=date, y=cum_loss * m3s_to_acftmonth / 1000))
p <- p + geom_area(aes(fill=data))
p <- p + geom_hline(yintercept=0)
p <- p + geom_line(data=line_df, colour="black", linetype="longdash")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Scenario", values= cc_colors[c(4,2,3,1)], labels=c("Warm-Wet", "Hot-Wet", "Warm-Dry", "Hot-Dry"))
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="5 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Cumulative Deficit (1,000 ac-ft)", breaks=seq(-500,500,50))
p <- p + facet_grid(data ~ .)
p <- p + theme(
  strip.background = element_blank(),
  strip.text = element_blank()
)
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_acft_vert.png"),  p, width=7, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_acft_vert.pdf"),  p, width=7, height=6)
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_acft_vert.svg"),  p, width=7, height=6)


p <- p + facet_grid(precip ~ temp)

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_acft_square.png"),  p, width=8, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_acft_square.pdf"),  p, width=8, height=6)
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_acft_square.svg"),  p, width=8, height=6)



plot_test <- flow_all$data %in% c( "HDN5", "HWN5", "WDN5", "WWN5")

area_df <- flow_all[plot_test,]
area_df$data <- factor(area_df$data, levels=c( "WWN5", "HWN5", "WDN5", "HDN5"), labels=())

line_df <- flow_all[flow_all$data %in% "base", ]
line_df <- do.call("rbind", replicate(4, line_df, simplify = FALSE))
line_df$data <- area_df$data
line_df$precip <- area_df$precip
line_df$temp <- area_df$temp

p <- ggplot(area_df, aes(x=date, y=cum_loss * m3s_to_m3month / 1E6))
p <- p + geom_area(aes(fill=data))
p <- p + geom_hline(yintercept=0)
p <- p + geom_line(data=line_df, colour="black", linetype="longdash")
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Scenario", values= cc_colors[c(4,2,3,1)], labels=c("Warm-Wet", "Hot-Wet", "Warm-Dry", "Hot-Dry"))
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="5 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Cumulative Deficit (1x10^6 m3)", breaks=seq(-500,500,50))
p <- p + facet_grid(data ~ .)
p <- p + theme(
  strip.background = element_blank(),
  strip.text = element_blank()
)
p

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_m3_vert.png"),  p, width=7, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_m3_vert.pdf"),  p, width=7, height=6)
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_m3_vert.svg"),  p, width=7, height=6)


p <- p + facet_grid(precip ~ temp)

### Save figures
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_m3_square.png"),  p, width=8, height=6, dpi=300)
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_m3_square.pdf"),  p, width=8, height=6)
ggsave(file.path(write_output_base_path,"clim_change_cum_deficit_area_m3_square.svg"),  p, width=8, height=6)











