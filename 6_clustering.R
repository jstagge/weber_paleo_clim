
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
require(staggefuncs)
require(tidyverse)

### Load these functions for this unique project
require(ggrepel)
require(cluster)
require(ggdendro)
require(ade4)
library(dendextend)
require(ggsci)

### Load project specific functions
file.sources = list.files(function_path, pattern="*.R", recursive=TRUE)
sapply(file.path(function_path, file.sources),source)

### Load global functions
file.sources = list.files(global_path, pattern="*.R", recursive=TRUE)
sapply(file.path(global_path, file.sources),source)


### Create theme for dendrogram
theme_dendro <- theme_classic_correct(10)+  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())



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

cluster.list <- 2:18

###########################################################################
###  Read in Data
###########################################################################
### Read in observed flow
#read_location <- file.path(write_output_base_path, paste0(site_id,"_obs_perc_ts.csv"))
#flow_obs <- read.csv(file = read_location)

### Read in climate change
read_location <- file.path(write_output_base_path, paste0(site_id,"_drought_details.csv"))

drought_event_summary <- read.csv(file = read_location)




###########################################################################
###  Prepare to cluster
###########################################################################

### Cut to only paleoclimate and observed
drought_event_allcols <- subset(drought_event_summary, data=="paleo" | data=="observed")
### Extract data for clustering
row.names(drought_event_allcols) <- paste0(substr(drought_event_allcols$begin,1,7), "_",drought_event_allcols$data)
 
data_clustering <- drought_event_allcols[,seq(3,16)]



###########################################################################
## Calculate Distance Matrix
###########################################################################
### Handle Circular variables
data_circ <- data_clustering[,c(2,3)]
data_circ <- prep.circular(data_circ, c(1,1), c(12,12))

### Handle Ordinal variables
data_ord <- data.frame(data_clustering[,c(1)])

### Handle Quantitative variables
data_quant <- data_clustering[,seq(4,dim(data_clustering)[2])]  #dim(data_clustering)[2] #[,seq(1,5)]

### Create distance matrix
ktab1 <- ktab.list.df(list(data_ord, data_quant))  #data_circ, 
### Apply Gowers distance to a mix of circular, ordinal, and quantitative columns
### Scale each column by range before calculating distances
### Apply the method of Pavoine S., Vallet, J., Dufour, A.-B., Gachet, S. and Daniel, H. (2009) On the challenge of treating various types of variables: Application for improving the measurement of functional diversity. Oikos, 118, 391–402.
distMatrix <- dist.ktab(ktab1, c("O", "Q"), option="scaledBYrange")  #"C", 


### Handle Quantitative variables
data_quant <- data_clustering[,seq(4,5)]  #dim(data_clustering)[2] #[,seq(1,5)]

### Create distance matrix
ktab1 <- ktab.list.df(list(data_ord, data_quant))  #data_circ, 
### Apply Gowers distance to a mix of circular, ordinal, and quantitative columns
### Scale each column by range before calculating distances
### Apply the method of Pavoine S., Vallet, J., Dufour, A.-B., Gachet, S. and Daniel, H. (2009) On the challenge of treating various types of variables: Application for improving the measurement of functional diversity. Oikos, 118, 391–402.
distMatrix <- dist.ktab(ktab1, c("O", "Q"), option="scaledBYrange")  #"C", 




###########################################################################
## Run Agglomorative Heirarchical Clustering
###########################################################################
### Name cluster method and create folder for results
cluster_method <- "heir_agglom"

### Create folder to save results in
save_folder <- file.path(write_output_base_path, cluster_method)
dir.create(save_folder, showWarnings=FALSE)

##  Calculate agglomerative clustering
hclust_result <- hclust (d=distMatrix, method='ward.D2')
	
###########################################################################
## Determine clusters
###########################################################################	
#clusters_df <- calc_clusters(clustering_tree=hclust_result, cluster.list = cluster.list, drought_df=drought_event_allcols)
clusters_df <- calc_clusters(clustering_tree=hclust_result, cluster.list = cluster.list)

### Separate cluster results
clusters <- clusters_df$clust
cluster_goodness <- clusters_df$goodness

### Plot goodness of fit
dir.create(file.path(save_folder, "clusters"), showWarnings=FALSE)
p <- plot_cluster_goodness(cluster_goodness=cluster_goodness, save_directory=save_folder) 

### Create master dataframe with cluster information
cluster_event_df <- cbind(drought_event_allcols,clusters)


### Save cluster information
write.csv(cluster_event_df, file.path(save_folder, "clusters/cluster_df.csv"), row.names=FALSE)
write.csv(cluster_goodness, file.path(save_folder, "clusters/cluster_goodness.csv"), row.names=FALSE)

###########################################################################
## Plot Clusters
###########################################################################
for (k in cluster.list) {
	p <- plot_clusters(plot_df=cluster_event_df, k=k, save_directory=save_folder)
}

###########################################################################
## Plot Dendrogram
###########################################################################
### Loop through each of the cluster numbers
for (k in cluster.list) {
	p <- plot_dendrogram(clustering_tree=hclust_result, k=k, save_directory=save_folder, yscale=c(-0.2,2, 0.2, 2.05), label_offsets=c(0,-0.2,-0.4))
}


