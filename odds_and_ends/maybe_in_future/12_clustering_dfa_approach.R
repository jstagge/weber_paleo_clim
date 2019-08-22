
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
require(tidyverse)

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
## Open Output File so I can add demand shortages
###########################################################################
#e1 <- new.env() 
#load(file.path(weber_output_path, "weber_storage_output.RData"), e1) 
#drought_event_system <- get('drought_event_system', e1) 
#rm(e1)
### This does a weird thing where data is wrong

load(file.path(weber_output_path, "weber_storage_output.RData"))
load(file.path(weber_output_path, "weber_delivery_output.RData"))


###########################################################################
###  New Approach to  Drought Events
###########################################################################
data_and_responses <- expand.grid(data=data_levels, response=response_levels)

data_counter <- 0
	
for (i in seq(1,dim(data_and_responses)[1])){

	data_i <- data_and_responses$data[i]
	response_i <- data_and_responses$response[i]
	
	stor_subset <- stor_all %>%
		filter(data==data_i, response == response_i) %>% 
		arrange(date)
		
	shortage_subset <- demand_deliv_df %>%
		filter(data==data_i, response == response_i) %>%
		select(date, data, response, node, demand_shortage) %>%
		spread(node, demand_shortage) %>% 	
		arrange(date)	
	
	all_subset <- stor_subset %>% 
		full_join(shortage_subset, by="date", suffix=c("_stor","_short")) %>%
		arrange(date)	

	for (j in seq(1, dim(drought_event_summary)[1])){
		event_j <- drought_event_summary[j,]
	
		event_data <- all_subset %>%
			filter(date >= as.Date(event_j$begin), date <= as.Date(event_j$end))
		
		if (dim(event_data)[1] > 0) {
			to_min_stor <- event_data %>% select(starts_with("res_"), upper_weber_stor, upper_ogden_stor, lower_stor, total_res, current_res)
			to_min_stor <- apply(to_min_stor, 2, min, na.rm=TRUE)
			names(to_min_stor) <- paste0(names(to_min_stor), "_min")
		
			to_max_stor <- event_data %>% select(moderate, severe, extreme)
			to_max_stor <- apply(to_max_stor, 2, max, na.rm=TRUE)
			names(to_max_stor) <- paste0(names(to_max_stor), "_max")
			
			to_max_short <- event_data %>% select(starts_with("SA"), upper_weber_short, upper_ogden_short, lower_short, system)
			to_max_short <- apply(to_max_short, 2, max, na.rm=TRUE)
			names(to_max_short) <- paste0(names(to_max_short), "_max")	
		
			event_j <- cbind(event_j, response=response_i, t(to_min_stor), t(to_max_stor), t(to_max_short))
		
			if (data_counter == 1){
				drought_event_system <- rbind(drought_event_system, event_j)
			} else {
				drought_event_system <- event_j
				data_counter <- 1
			}
		}
	}
}
		

### Need to refactor data column
drought_event_system$data <- factor(drought_event_system$data, levels=c("paleo", "observed", "base", "WWN5", "HWN5", "CTN5", "WDN5", "HDN5"), labels=c("paleo", "observed", "base", "ww", "hw", "ctn", "wd", "hd"))








ggplot(subset(drought_event_system, response=="Base"), aes(x=dura_months/12, y=min_perc, colour=current_res_min)) + geom_point(size=4) + scale_colour_distiller(type = "seq", palette = "YlOrRd") + theme_classic_new()

ggplot(subset(drought_event_system, response=="Base"), aes(x=dura_months/12, y=min_perc, colour=current_res_min)) + geom_point(size=4) + scale_color_viridis(option="magma") + theme_classic_new()


ggplot(subset(drought_event_system, response=="Base"), aes(x=dura_months/12, y=min_perc, colour=current_res_min)) + geom_point(size=4) + scale_color_viridis(option="cividis") + theme_classic_new()

ggplot(subset(drought_event_system, response=="Base"), aes(x=dura_months/12, y=min_perc, colour=current_res_min)) + geom_point(size=4) + scale_color_cividis(direction=-1) + theme_classic_new()
 

ggplot(subset(drought_event_system, response=="Base"), aes(x=dura_months/12, y=min_perc, colour=upper_ogden_stor_min)) + geom_point(size=6) + scale_color_cividis(direction=-1) + theme_classic_new()
		

ggplot(subset(drought_event_system, response=="Base" & data=="paleo"), aes(x=dura_months/12, y=min_perc, colour=log10(system_max))) + geom_point(size=4) + scale_colour_viridis(option="magma")


ggplot(subset(drought_event_system, response=="Base" & data=="paleo"), aes(x=dura_months/12, y=min_perc, colour=log10(system_max))) + geom_point(size=4) + scale_color_cividis(direction=1)


ggplot(subset(drought_event_system, response=="Base"), aes(x=dura_months/12, y=min_perc, colour=log10(SA11_max))) + geom_point(size=4) + scale_colour_distiller(type = "seq", palette = "YlOrRd", direction = 1) + theme_classic_new()



###########################################################################
###  Prepare to cluster
###########################################################################
### Select columns to retain (everything to the right of response column)
response_column <- which(names(drought_event_system) == "response")
select_columns <- seq(response_column + 1, dim(drought_event_system)[2])

### Extract data for clustering
data_clustering <- drought_event_system %>% 
	filter(data %in% c("observed", "paleo")) %>%
	filter(response == "Base") %>%
	#arrange(begin) # %>%
	select(select_columns)

head(data_clustering)
tail(data_clustering)

###########################################################################
## Calculate Distance Matrix
###########################################################################
### Create distance matrix
scaled_data <- scale(data_clustering)
rownames(scaled_data) <- seq(1,dim(data_clustering)[1])

distMatrix <- dist(scaled_data, method="euclidean")

###########################################################################
## Run Agglomorative Heirarchical Clustering
###########################################################################
### Name cluster method and create folder for results
cluster_method <- "heir_agglom_storage_delivery"
	
### Create folder to save results in
save_folder <- file.path(write_output_base_path, "clustering")
save_folder <- file.path(save_folder, cluster_method)
dir.create(save_folder, recursive=TRUE, showWarnings=FALSE)

##  Calculate agglomerative clustering
hclust_result <- hclust (d=distMatrix, method='ward.D2')
	
###########################################################################
## Determine clusters
###########################################################################	
clusters_df <- calc_clusters(clustering_tree=hclust_result, cluster.list = cluster.list)

### Separate cluster results
clusters <- clusters_df$clust
cluster_goodness <- clusters_df$goodness

### Plot goodness of fit
dir.create(file.path(save_folder, "clusters"), showWarnings=FALSE)
p <- plot_cluster_goodness(cluster_goodness=cluster_goodness, save_directory=save_folder) 

### Create master dataframe with cluster information
cluster_event_df <- cbind(data_clustering,clusters)

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







####### A different approach
## generate 25 objects, divided into 2 clusters.
require(factoextra)
fviz_nbclust(scaled_data, pam, method = "wss", diss=distMatrix) +
  geom_vline(xintercept = 3, linetype = 2)
  
  
 fviz_nbclust(scaled_data, pam, method = "silhouette", diss=distMatrix) +
  geom_vline(xintercept = 3, linetype = 2) 
  
  
 fviz_nbclust(scaled_data, pam, method = "gap_stat", diss=distMatrix) +
  geom_vline(xintercept = 3, linetype = 2) 
  
  
  
#yup <- pam(data_clustering, k=6, diss=FALSE, stand=TRUE)
yup <- pam(distMatrix, k=5, diss=TRUE, stand=FALSE)

sil = silhouette (yup$clustering, distMatrix)
plot(sil)


summary(yup)

plot(yup$data, col = yup$clustering)
points(yup$medoids, col = 1:2, pch = 4)

plot (yup$data, col = yup$clustering)
#add the medoids to the plot
points(yup$medoids, col = 1:3, pch = 4)



cluster_results <- data_clustering
cluster_results$cluster <- yup$clustering

toplot <- cluster_results %>%
   gather(variable, value, -cluster)

toplot$variable <- factor(toplot$variable, levels=colnames(cluster_results))
toplot$cluster <- factor(toplot$cluster)
  
ggplot(toplot, aes(x=cluster, y=value)) + geom_boxplot() + geom_jitter(width = 0.2, aes(colour=cluster), alpha=0.5) + theme_classic_new() + facet_wrap(~ variable, scales="free")
  


,dis)






 data_clustering <- drought_event_system %>% 
					filter(response == "Base") %>%
					select(starts_with("res_"), upper_weber_min, upper_ogden_min, lower_min, total_res_min, current_res_min, moderate_max, severe_max, extreme_max) %>%
					select(-contains("res_9"), -contains("res_10"))
			

### Come back and think about this
mydata <- scale(data_clustering)
distMatrix <- dist(mydata)

##  Calculate agglomerative clustering
hclust_result <- hclust (d=distMatrix, method='ward.D2')
	

require(cluster)
hc2 <- agnes(data_clustering, method = "complete")
hc2$ac

# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(data_clustering, method = x)$ac
}

map_dbl(m, ac)

hc3 <- agnes(data_clustering, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 

sub_grp <- cutree(hc3, k = 6)

cluster_results <- data_clustering %>%
  mutate(cluster = sub_grp) 
  
require(factoextra)
fviz_cluster(list(data = data_clustering, cluster = sub_grp))

fviz_nbclust(data_clustering, FUN = hcut, method = "wss", k.max=20)

fviz_nbclust(data_clustering, FUN = hcut, method = "silhouette", k.max=20)

gap_stat <- clusGap(data_clustering, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)




