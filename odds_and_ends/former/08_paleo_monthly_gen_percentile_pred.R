# *------------------------------------------------------------------
# | PROGRAM NAME: Data Music
# | FILE NAME: paleo_monthly_gen_null.R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  This is a code wrapper to generate a midicsv file from data.
# | The resulting file can be processed into a midi file using the program csvmidi.
# | This midi file can then be played using timidity.
# | Check the ToRun.txt file for detailed code             
# |
# *------------------------------------------------------------------
# | COMMENTS:               
# |
# |  1:  
# |  2: 
# |  3: 
# |*------------------------------------------------------------------
# | DATA USED:               
# | This is a test instance using reconstructed climate indices ENSO and PDO
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
output_name <- "apr_model"

weber_output_path <- file.path(output_path,"paleo_weber")
write_output_path <- file.path(weber_output_path,output_name)

write_figures_path <- file.path(weber_output_path, "figures")
write_figures_path <- file.path(write_figures_path, output_name)

dir.create(write_output_path)
dir.create(write_figures_path)


###########################################################################
###  Load functions
###########################################################################
### Load these functions for all code
require(colorout)
require(assertthat)

### Load these functions for this unique project
require(ggplot2)
require(monthlypaleo)
require(staggefuncs)

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
site_id_list <- c("10128500")
site_name_list <- c("Weber River")
recons_file_name_list <- c("weber2014flow.txt")

first_month_wy <- 10 ### Water Year starts on Oct 1
param_cd <- "00060"

monthly_distr <- "gamma"
annual_distr <- "logis"

ref_period <- c(1900,2005)

### Number of PCs to consider as predictors
pc_pred <- 8
### Lags to consider
lags_pred <- seq(-2,2)


###########################################################################
## Read in reconstructed flows, tree rings and climate indices
###########################################################################
### Read in PC scores
pc_score_impute_res <- read.csv(file.path(file.path(weber_output_path,"pca_chronol"), "PC_Score_impute_res.csv"))
pc_score_impute_std <- read.csv(file.path(file.path(weber_output_path,"pca_chronol"), "PC_Score_impute_std.csv"))


### Read in climate indices
clim_ind <- read.csv(file.path(file.path(data_path,"paleo_clim_ind"), "clim_ind.csv"))
### Remove ENSO_var, this is simply a running variance calculation - not useful as a predictor
clim_ind <- subset(clim_ind, select=-c(ENSO_var))

################################################
### Create matrix of potential predictor values
#################################################
## Cut PCs to only number in initial settings
pc_cut <- seq(1,pc_pred+1)
pc_score_impute_std <- pc_score_impute_std[,pc_cut]
pc_score_impute_res <- pc_score_impute_res[,pc_cut]


### Create only climate index predictors
pred_clim_only <- clim_ind
pred_clim_only_lag <- shift_indices_df(pred_clim_only, lags=lags_pred)

pred_enso_ind <- clim_ind[,seq(1,3)]
pred_enso_ind_lag <- shift_indices_df(enso_ind, lags=lags_pred)

### Create concurrent predictors (same calendar year), with climate indices and first 8 PCs
pred_clim_pca_impute_std_concur <- merge(clim_ind, pc_score_impute_std, by.x="Year", by.y="X", all=TRUE)
pred_clim_pca_impute_res_concur <- merge(clim_ind, pc_score_impute_res, by.x="Year", by.y="X", all=TRUE)

pred_enso_pca_impute_std_concur <- merge(enso_ind, pc_score_impute_std, by.x="Year", by.y="X", all=TRUE)
pred_enso_pca_impute_res_concur <- merge(enso_ind, pc_score_impute_res, by.x="Year", by.y="X", all=TRUE)

### Lag the predictors
pred_clim_pca_impute_std_lag <- shift_indices_df(pred_clim_pca_impute_std_concur, lags=lags_pred)
pred_clim_pca_impute_res_lag <- shift_indices_df(pred_clim_pca_impute_res_concur, lags=lags_pred)

pred_enso_pca_impute_std_lag <- shift_indices_df(pred_enso_pca_impute_std_concur, lags=lags_pred)
pred_enso_pca_impute_res_lag <- shift_indices_df(pred_enso_pca_impute_res_concur, lags=lags_pred)

###########################################################################
###  Set up a loop to run through all site_ides and Transform based on the Percentile Model
###########################################################################
for (n in seq(1,length(site_id_list))) {

site_id <- site_id_list[n]
site_name <- site_name_list[n]
recons_file_name <- recons_file_name_list[n]

###########################################################################
###  Read in Flow Data
###########################################################################

### Read in observed flow and fix data type
obs_file_name <- paste0(site_id,"_",param_cd,"_mon_wy.csv")
flow_obs <- read.csv(file.path(weber_output_path,paste0("observed_flow/",obs_file_name)))
flow_obs$date <- as.Date(flow_obs$date)  
#head(flow_obs) # Review data frame

### Read in reconst flows (use fread because of large header)
flow_recon <- read_table_wheaders(file.path(data_path,paste0("paleo_flow_annual/",recons_file_name)), sep="\t", na.string="-9999")


flow_recon <- merge(flow_recon, data.frame(age_AD=flow_obs$water_year, flow.obs.m3s=flow_obs$annual_mean), by="age_AD", all.x=TRUE)

################################################
### Read in monthly and annual parameters to Transform Percentiles
#################################################
percentile_path <- file.path(file.path(weber_output_path, "paleo_reconst"), "ap_model")

monthly_param <- list(param=read.csv(file.path(percentile_path, paste0(site_id,"/",site_id,"_param_month_",monthly_distr,".csv"))), distr=monthly_distr)
monthly_param$param <- monthly_param$param[,c(2,3)]

annual_param <- list(param=read.csv(file.path(percentile_path, paste0(site_id,"/",site_id,"_param_annual_",annual_distr,".csv"))), distr=annual_distr)
annual_param$param <- annual_param$param[,seq(1,3)]

#################################################
### Apply the Percentile Model with Lagged years (1 year previous, 1 year following)
#################################################
### Set up list of data combinations to fit
### Only use the imputed values
#pred_list <- c("enso_ind", "clim_only", "enso_ind_lag", "clim_only_lag",  "enso_pca_impute_std_concur", "enso_pca_impute_res_concur", "clim_pca_impute_std_concur", "clim_pca_impute_lag")

pred_list <- c("enso_ind", "clim_only", "enso_pca_impute_std_concur", "enso_pca_impute_res_concur", "clim_pca_impute_std_concur", "clim_pca_impute_res_concur")


data_list <- c("observ_annual", "rec_local_m3s")

run_combinations <- expand.grid(predictors=pred_list, data=data_list)	

### Loop through all combations, fitting, and saving
for (c in seq(1,dim(run_combinations)[1])){

pred_c <- run_combinations$predictors[c]
data_c <- run_combinations$data[c]

if(data_c != "observ_annual"){
	data_name <- as.character(paste0("annual_flow_",data_c))
} else {
	data_name <- as.character(data_c)
}

### Fit model
lag_fit <- perc_fit_pred(data=flow_obs, predictors=get(paste0("pred_",pred_c
)), monthly_param, annual_param, data_name= data_name)

### Write to csv file
write_location <- file.path(write_output_path, paste0(site_id,"_",output_name, "_coef_",data_c,"_",pred_c,".csv"))
write.csv(lag_fit$coef, file = write_location,row.names=TRUE)
write_location <- file.path(write_output_path, paste0(site_id,"_",output_name, "_alpha_lambda_",data_c,"_",pred_c,".csv"))
write.csv(lag_fit$alpha_lambda, file = write_location,row.names=TRUE)
}
	
	
}




###########################################################################
###  Set up a loop to run through all site_ides and Plot results
###########################################################################
dir.create(file.path(write_figures_path,"png/"), showWarnings=FALSE)
dir.create(file.path(write_figures_path,"svg/"), showWarnings=FALSE)
dir.create(file.path(write_figures_path,"pdf/"), showWarnings=FALSE)


for (n in seq(1,length(site_id_list))) {

site_id <- site_id_list[n]
site_name <- site_name_list[n]
recons_file_name <- recons_file_name_list[n]

### Create table of all combinations
### Chose not to plot all the lags (too much) 
#pred_list <- c("clim_only","clim_only_lag","clim_pca_impute_concur", "clim_pca_impute_lag")
pred_list <- c("enso_ind", "clim_only", "enso_pca_impute_std_concur", "enso_pca_impute_res_concur", "clim_pca_impute_std_concur", "clim_pca_impute_res_concur")

data_list <- c("observ_annual", "rec_local_m3s")

run_combinations <- expand.grid(predictors=pred_list, data=data_list)	


### Loop through all combations, reading in data and saving
for (c in seq(1,dim(run_combinations)[1])){

pred_c <- as.character(run_combinations$predictors[c])
data_c <- run_combinations$data[c]

if(data_c != "observ_annual"){
	data_name <- as.character(paste0("annual_flow_",data_c))
} else {
	data_name <- as.character(data_c)
}

### Read in CSV file
file_location <- file.path(write_output_path, paste0(site_id,"_",output_name, "_coef_",data_c,"_",pred_c,".csv"))

coef_df <- read.csv(file_location)
rownames(coef_df) <- as.character(coef_df$X)

### Reorganize dataframe, remove the first column, set class
coef_df <- subset(coef_df, select=-c(X))
colnames(coef_df) <- seq(1,12)
coef_df <- data.frame(coef_df)

### Create plot of annual drivers
predictor_levels <- c("annual_norm", "annual_norm_neg_1year", "annual_norm_pos_1year")
predictor_labels <- c("Concurrent Water Year\nStd Normal", "Previous (-1) Water Year\nStd Normal", "Next (+1) Water Year\nStd Normal")
p <- norm_coef_plot(coef_df, predictor_levels=predictor_levels, predictor_labels=predictor_labels)
### Save annual drivers
ggsave(paste0(file.path(write_figures_path,"png/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_annual.png"), p, width=6, height=4, dpi=600)
ggsave(paste0(file.path(write_figures_path,"svg/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_annual.svg"), p, width=6, height=4)
ggsave(paste0(file.path(write_figures_path,"pdf/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_annual.pdf"), p, width=6, height=4)


#### Add a combination of 2 ENSO methods, not needed
#coef_df_enso_comb <- coef_df[rownames(coef_df) == "ENSO",] + coef_df[rownames(coef_df) == "ENSO_short",]
#rownames(coef_df_enso_comb) <- "ENSO_comb"
#coef_df <-  rbind(coef_df, coef_df_enso_comb)

### Create plot of climate drivers
predictor_levels <- c("ENSO", "ENSO_short", "PDO")
predictor_labels <- c("ENSO (NADA)", "ENSO (Pacific Proxy)", "PDO")
p <- norm_coef_plot(coef_df, predictor_levels=predictor_levels, predictor_labels=predictor_labels)
### Reset the color scheme
p <- p + scale_colour_brewer(name="Predictor", type="qual", palette = "Accent")
### Save annual drivers
ggsave(paste0(file.path(write_figures_path,"png/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_clim_ind.png"), p, width=6, height=4, dpi=600)
ggsave(paste0(file.path(write_figures_path,"svg/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_clim_ind.svg"), p, width=6, height=4)
ggsave(paste0(file.path(write_figures_path,"pdf/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_clim_ind.pdf"), p, width=6, height=4)

### If it includes pca, create plots
if (grepl("pca",pred_c)) {
### Create a plot of PC coefficients
predictor_levels <- paste0("PC",seq(1,8))
predictor_labels <- predictor_levels
p <- norm_coef_plot(coef_df, predictor_levels=predictor_levels, predictor_labels=predictor_labels)
### Reset color scale
p <- p + scale_colour_brewer(name="Predictor", type="qual", palette = "Set2")
### Save annual drivers
ggsave(paste0(file.path(write_figures_path,"png/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_pcs.png"), p, width=6, height=4, dpi=600)
ggsave(paste0(file.path(write_figures_path,"svg/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_pcs.svg"), p, width=6, height=4)
ggsave(paste0(file.path(write_figures_path,"pdf/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_pcs.pdf"), p, width=6, height=4)
}

}

}










###########################################################################
###  Set up a loop to run through all site_ides and Plot results
###########################################################################
for (n in seq(1,length(site_id_list))) {

site_id <- site_id_list[n]
site_name <- site_name_list[n]
recons_file_name <- recons_file_name_list[n]

### Create table of all combinations
### Chose not to plot all the lags (too much) 
pred_list <- c("clim_only_lag","clim_pca_impute_lag")

if (site_id == "10109001") {
	data_list <- c("observ_annual", "rec_local_m3s", "rec_region_m3s")
} else {
	data_list <- c("observ_annual", "rec_m3s")
}

run_combinations <- expand.grid(predictors=pred_list, data=data_list)	

### Loop through all combations, reading in data and saving
for (c in seq(1,dim(run_combinations)[1])){

pred_c <- as.character(run_combinations$predictors[c])
data_c <- run_combinations$data[c]

if(data_c != "observ_annual"){
	data_name <- as.character(paste0("annual_flow_",data_c))
} else {
	data_name <- as.character(data_c)
}

### Read in CSV file
file_location <- file.path(write_output_path, paste0(site_id,"_",output_name, "_coef_",data_c,"_",pred_c,".csv"))
coef_df <- read.csv(file_location)
rownames(coef_df) <- as.character(coef_df$X)

### Reorganize dataframe, remove the first column, set class
coef_df <- subset(coef_df, select=-c(X))
colnames(coef_df) <- seq(1,12)
coef_df <- data.frame(coef_df)

### Create plot of annual drivers
predictor_levels <- c("annual_norm", "annual_norm_neg_1year", "annual_norm_pos_1year")
predictor_labels <- c("Concurrent Water Year\nStd Normal", "Previous (-1) Water Year\nStd Normal", "Next (+1) Water Year\nStd Normal")
p <- norm_coef_plot(coef_df, predictor_levels=predictor_levels, predictor_labels=predictor_labels)
### Save annual drivers
ggsave(paste0(file.path(write_figures_path,"png/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_annual.png"), p, width=6, height=4, dpi=600)
ggsave(paste0(file.path(write_figures_path,"svg/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_annual.svg"), p, width=6, height=4)
ggsave(paste0(file.path(write_figures_path,"pdf/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_annual.pdf"), p, width=6, height=4)


#### Add a combination of 2 ENSO methods, not needed
#coef_df_enso_comb <- coef_df[rownames(coef_df) == "ENSO",] + coef_df[rownames(coef_df) == "ENSO_short",]
#rownames(coef_df_enso_comb) <- "ENSO_comb"
#coef_df <-  rbind(coef_df, coef_df_enso_comb)

### Create plot of climate drivers
predictor_levels <- c("ENSO_0", "ENSO_short_0", "PDO_0")
predictor_labels <- c("ENSO (NADA)", "ENSO (Pacific Proxy)", "PDO")
p <- norm_coef_plot(coef_df, predictor_levels=predictor_levels, predictor_labels=predictor_labels)
### Reset the color scheme
p <- p + scale_colour_brewer(name="Predictor", type="qual", palette = "Accent")
### Save annual drivers
ggsave(paste0(file.path(write_figures_path,"png/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_clim_ind.png"), p, width=6, height=4, dpi=600)
ggsave(paste0(file.path(write_figures_path,"svg/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_clim_ind.svg"), p, width=6, height=4)
ggsave(paste0(file.path(write_figures_path,"pdf/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_clim_ind.pdf"), p, width=6, height=4)

### Create plots of PCs at different lags
lag_names <- c("-2", "-1", "0", "1", "2")
lag_numbers <- c(".2", ".1", "0", "1", "2")

### If it includes pca, create plots
if (grepl("pca",pred_c)) {

for (k in seq(1,length(lag_numbers))) {
### Create a plot of PC coefficients
predictor_labels <- paste0("PC",seq(1,8))
predictor_levels <- paste0(predictor_labels, "_",lag_numbers[k])
predictor_labels <- paste0(predictor_labels, " Lag ", lag_names[k])
p <- norm_coef_plot(coef_df, predictor_levels=predictor_levels, predictor_labels=predictor_labels)
### Reset color scale
p <- p + scale_colour_brewer(name="Predictor", type="qual", palette = "Set2")
p <- p + coord_cartesian(ylim=c(-0.2,0.2))
#
### Save annual drivers
ggsave(paste0(file.path(write_figures_path,"png/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_pcs_lag_",lag_names[k],".png"), p, width=6, height=4, dpi=600)
ggsave(paste0(file.path(write_figures_path,"svg/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_pcs_lag_",lag_names[k],".svg"), p, width=6, height=4)
ggsave(paste0(file.path(write_figures_path,"pdf/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_pcs_lag_",lag_names[k],".pdf"), p, width=6, height=4)
}


### Extract PC names
coef_df_names <- rownames(coef_df)
coef_df_pc <- substr(coef_df_names,1,2) == "PC"
coef_df_pc[[1]] <- TRUE
coef_df_pc <- coef_df[coef_df_pc,]
coef_df_names <- rownames(coef_df_pc)

for (k in seq(1,8)) {
### Subset to each PC
pc_test <- substr(coef_df_names,1,3) == paste0("PC", k)
pc_test[[1]] <- TRUE
pc_subset <- coef_df_pc[pc_test,]

### Create plot of climate drivers
predictor_levels <- rownames(pc_subset)
predictor_levels <- predictor_levels[seq(2,length(predictor_levels))]
predictor_labels <- c("Lag -2", "Lag -1", "Concur", "Lag +1", "Lag +2")
p <- norm_coef_plot(pc_subset, predictor_levels=predictor_levels, predictor_labels=predictor_labels)

### Save plot
ggsave(paste0(file.path(write_figures_path,"png/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_pc",k,"_lag",".png"), p, width=6, height=4, dpi=600)
ggsave(paste0(file.path(write_figures_path,"svg/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_pc",k,"_lag",".svg"), p, width=6, height=4)
ggsave(paste0(file.path(write_figures_path,"pdf/"), site_id, "_percent_pred_coef_", data_name,"_",pred_c,"_pc",k,"_lag",".pdf"), p, width=6, height=4)

}

}
}
}


