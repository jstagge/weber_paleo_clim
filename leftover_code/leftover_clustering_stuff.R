

###########################################################################
## Run Divisive Heirarchical Clustering
###########################################################################
### Name cluster method and create folder for results
cluster_method <- "heir_divis"

### Create folder to save results in
save_folder <- file.path(write_output_path, cluster_method)
dir.create(save_folder, showWarnings=FALSE)

##  Calculate divisive clustering
dv <- diana(x=distMatrix, diss=TRUE, stand=FALSE)
hclust_result <- as.hclust(dv)
	
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
cluster_event_df <- cbind(drought_event_df,clusters)

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
	p <- plot_dendrogram(clustering_tree=hclust_result, k=k, save_directory=save_folder, yscale=c(-0.05,0.55, 0.05, 0.57), label_offsets=c(0,-0.06,-0.12))
}

###########################################################################
## Create Parallel Coordinates Plot
###########################################################################
### Peak at 7 clusters, 9, 13

### Prepare data for plot
rownames(cluster_event_df) <- cluster_event_df$begin

### Only works up to 10 clusters 
plot_par_coord(full_cluster_df=plot_data, k=8)

### 9 is also good
plot_par_coord(full_cluster_df=plot_data, k=9)

### 10 is also good
plot_par_coord(full_cluster_df=plot_data, k=10)

### Run 15
plot_data <- data.frame(cluster=cluster_event_df$k_15, cluster_event_df[,c(seq(3,13), 15, 14, 16)])
plot_data <- plot_data[complete.cases(plot_data),]

### Create plot
parcoords(plot_data, reorderable = TRUE, brushMode = "1D-axes-multi",rownames = F, width=2000,alphaOnBrushed=0.2, color="#000000"
)


### Run 15
plot_data <- data.frame(cluster=cluster_event_df$k_15, cluster_event_df[,c(seq(3,13), 15, 14, 16)])
plot_data <- plot_data[complete.cases(plot_data),]

### Create plot
parcoords(plot_data, reorderable = TRUE, brushMode = "1D-axes-multi",rownames = F, width=2000,alphaOnBrushed=0.2, color="#000000"
)












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
clusters_df <- calc_clusters(clustering_tree=hclust_result, cluster.list = cluster.list, drought_df=drought_event_allcols)

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







####### A different approach






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



