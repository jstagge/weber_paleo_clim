
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

#load(file.path(weber_output_path, "weber_storage_output.RData"))
#load(file.path(weber_output_path, "weber_delivery_output.RData"))



###########################################################################
## Calculate Distance Matrix
###########################################################################
### Extract data for clustering
data_clustering <- drought_event_summary %>% 
	filter(data %in% c("observed", "paleo")) #%>%
	#arrange(begin) # %>%
	#select(-begin, -end, -data)

head(data_clustering)
tail(data_clustering)

#row.names(data_clustering) <- substr(drought_event_df$begin,1,7)
 
### Handle Circular variables
data_circ <- data_clustering %>% select(begin_month, end_month)
head(data_circ)
plot(data_circ$begin_month, data_circ$end_month)
data_circ <- prep.circular(data_circ, c(1,1), c(12,12))
head(data_circ)
### This changes the numbers - strange, actually just makes them go from 0 to 11

### Handle Ordinal variables
data_ord <- data.frame(data_clustering$dura_months)

### Handle Quantitative variables
data_quant <- data_clustering %>% select(-begin, -end, -begin_month, -end_month, -dura_months, -data)

### Create distance matrix
ktab1 <- ktab.list.df(list(data_ord, data_quant)) #data_circ, 
### Apply Gowers distance to a mix of circular, ordinal, and quantitative columns
### Scale each column by range before calculating distances
### Apply the method of Pavoine S., Vallet, J., Dufour, A.-B., Gachet, S. and Daniel, H. (2009) On the challenge of treating various types of variables: Application for improving the measurement of functional diversity. Oikos, 118, 391–402.

distMatrix <- dist.ktab(ktab1, c("O", "Q"), option="scaledBYrange")

#ktab1 <- ktab.list.df(list(data_circ, data_ord, data_quant)) 
#distMatrix <- dist.ktab(ktab1, c("C", "O", "Q"), option="scaledBYrange")

###########################################################################
## Run Agglomorative Heirarchical Clustering
###########################################################################
### Name cluster method and create folder for results
cluster_method <- "heir_agglom"

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
dir.create(file.path(save_folder, "goodness"), showWarnings=FALSE)
p <- plot_cluster_goodness(cluster_goodness=cluster_goodness, save_directory=save_folder) 

### Create master dataframe with cluster information
cluster_event_df <- cbind(data_clustering,clusters)

### Save cluster information
dir.create(file.path(save_folder, "clusters"), showWarnings=FALSE)
write.csv(cluster_event_df, file.path(save_folder, "clusters/cluster_df.csv"), row.names=FALSE)
write.csv(cluster_goodness, file.path(save_folder, "goodness/cluster_goodness.csv"), row.names=FALSE)


### Silhouette plot
for (k in cluster.list) {
	p <- plot_sil(cluster_df=clusters, dist=distMatrix, k=k, save_directory=save_folder)
}


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

cluster_event_df$k_11 <- as.character(cluster_event_df$k_11)
cluster_event_df %>%
  group_by(k_11) %>%
  summarize(avg_min_perc = mean(min_perc))
 
ggplot(cluster_event_df, aes(x=k_11, y=min_perc)) + geom_boxplot() + geom_jitter(width = 0.2, aes(colour=k_11)) + theme_classic_new() + scale_color_d3(name="Cluster", palette = "category20")
	

yup <- cluster_event_df %>%
	mutate(cluster = k_11) %>%
   select(-begin, -end, -data, -starts_with("k_")) %>%
   gather(variable, value, -cluster)

yup$variable <- factor(yup$variable, levels=colnames(cluster_event_df))
  
yup$cluster <- factor(yup$cluster, levels=as.character(seq(1,11)))
ggplot(yup, aes(x=cluster, y=value)) + geom_boxplot() + geom_jitter(width = 0.2, aes(colour=cluster), alpha=0.5) + theme_classic_new() + facet_wrap(~ variable, scales="free") + scale_color_d3(name="Cluster", palette = "category20") 

  
  
###########################################################################
## Create Parallel Coordinates Plot
###########################################################################
### Prepare data for plot
rownames(cluster_event_df) <- cluster_event_df$begin

### Only works up to 10 clusters 
plot_par_coord(full_cluster_df=plot_data, k=6)

### 9 is also good
plot_par_coord(full_cluster_df=plot_data, k=9)

### Run 12
plot_data <- data.frame(cluster=cluster_event_df$k_12, cluster_event_df[,c(seq(3,13), 15, 14, 16)])
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
## Find Centroid
###########################################################################
year_only <- as.numeric(substr(rownames(clusters),1,4))
### Choose 1920 to capture the dust bowl years
obs_test <- year_only >= 1920
pre_test <- year_only < 1920

### 6 Cluster solution
sapply(unique(clusters[,5]), clust.medoid, as.matrix(distMatrix), clusters[,5])

### Post 1920 and Pre 1920 clusters
sapply(unique(clusters[obs_test,5]), clust.medoid, as.matrix(distMatrix)[obs_test,obs_test], clusters[obs_test,5])
sapply(unique(clusters[pre_test,5]), clust.medoid, as.matrix(distMatrix)[pre_test,pre_test], clusters[pre_test,5])




 
###########################################################################
## Manually create final plot of 6
###########################################################################
require(RColorBrewer)
#manual_pal <- c(pal_d3(palette="category20")(5)[c(1,2,4,5)], brewer.pal(7, "YlGn")[seq(3,7)])
#manual_pal <- manual_pal[c(1, 5, 6, 2, 7, 3, 8, 9, 4)]

#manual_pal <- c(pal_d3(palette="category20")(5)[c(1,2,4,5)], brewer.pal(7, "YlGn")[seq(3,7)])
#manual_pal <- manual_pal[c(1, 5, 6, 2, 7, 3, 8, 9, 4)]
#cb_check_plot(pal=manual_pal)


manual_pal <- c(pal_d3(palette="category20")(4), brewer.pal(9, "Purples")[seq(3,9)])
manual_pal <- manual_pal[c(1, 5, 6, 7, 8, 2, 9, 10, 4, 11, 3)]

manual_shapes <- c(25, 20, 20, 20, 20, 23, 20, 20, 24, 20, 22)

cb_check_plot(pal_d3(palette="category20")(10))

cb_check_plot(pal=manual_pal)

cb_check_plot(pal=pal_d3(palette="category20")(6))

#1 - blue, moderate
#6 - brown, extreme
#9 - poop green, moderate weird
#11 - corn blue, small

#4 - red, shortish
#2 - orange, short

#1, 6 (most extreme), 9, 11

#4? 2?
k <- 11
cluster_event_df$k_11 <- factor(cluster_event_df$k_11, levels=as.character(seq(1,11)))

bold_test <- cluster_event_df$k_11 %in% c("1", "6", "9", "11")
plot_bold <- cluster_event_df[bold_test, ]
plot_alpha <- cluster_event_df[!bold_test, ]

p <- ggplot(cluster_event_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4), colour=k_11, fill=k_11, shape=k_11))
	#p <- p + geom_point(aes(colour=timeperiod))
	p <- p + geom_point(size=2.5, alpha=0)
	p <- p + geom_point(data=plot_alpha, size=2.5, alpha=0.4)
	p <- p + geom_point(data=plot_bold, size=2.5, alpha=0.85)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,16,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	p <- p + scale_color_manual(name="Cluster", values=manual_pal)
	p <- p + scale_fill_manual(name="Cluster", values=manual_pal)
	p <- p + scale_shape_manual(name="Cluster", values=manual_shapes)
	p <- p + theme_classic(10)
	p <- p + theme(legend.position = c(0.85, 0.85))
p <- p + guides(colour=guide_legend(ncol=3))
	p



### Save results
ggsave(file.path(save_folder, paste0("clusters/",cluster_method,"_clust_k_",k,"_points_new_colors.png")),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(save_folder, paste0("clusters/",cluster_method,"_clust_k_",k,"_points_new_colors.pdf")),  p, width=6.5, height=5)
ggsave(file.path(save_folder, paste0("clusters/",cluster_method,"_clust_k_",k,"_points_new_colors.svg")),  p, width=6.5, height=5)


k <- 11
### Cut to every other year
dend_labels <- hclust_result$labels[hclust_result$order]
#dend_labels <- substr(cluster_event_df$begin[as.numeric(hclust_result$labels)], 1, 7)

### Add to reorganize cluster numbers
clusters_k <- clusters[,(11-1)]
cluster_number <- clusters_k[match(dend_labels, names(clusters_k))]

### Continue processing labels
dend_labels <- substr(cluster_event_df$begin[as.numeric(dend_labels)], 1, 4)
length_labels <- length(dend_labels)
dend_labels_everynth <- rep(NA, length_labels)
dend_labels_everynth[seq(1, length_labels, 3)] <- dend_labels[seq(1, length_labels, 3)]

#color_list <- pal_d3("category20")(k)[cluster_number]
#color_list <- unique(color_list)
color_list <- manual_pal[cluster_number]
color_list <- unique(color_list)

length_clusters <- length(hclust_result$height)

cut_line <- (hclust_result$height[length_clusters-k+1] + hclust_result$height[length_clusters-k+2])/2


### Create Dendrogram of all years with a cut at 6 clusters
dend <- hclust_result %>% as.dendrogram %>%
   set("branches_k_color", k=k, value = color_list ) %>%
   set("labels_colors", k=k, value = color_list ) %>%
   set("labels", dend_labels_everynth) %>% set("labels_to_char") %>%
   set("labels_cex", c(0.6)) %>%
    set("branches_lwd", c(0.5))

# Convert to ggdend object
ggd1 <- as.ggdend(dend)

# Change the theme to minimal
p <- ggplot(ggd1, horiz = TRUE, theme = theme_dendro)
p <- p + geom_hline(yintercept=cut_line, linetype="dashed")
p <- p + labs(y = "Sum of Within-Group Distance Squared")
#p <- p + coord_flip(ylim=yscale[3:4])
p

### Save clusters
ggsave(file.path(save_folder, paste0("dendro/",cluster_method,"_dendro_k_",k,"_everythird_new_colors.png")),  width=7.5, height=6, dpi=300)
ggsave(file.path(save_folder, paste0("dendro/",cluster_method,"_dendro_k_",k,"_everythird_new_colors.pdf")),  width=7.5, height=6)
ggsave(file.path(save_folder, paste0("dendro/",cluster_method,"_dendro_k_",k,"_everythird_new_colors.svg")),  width=7.5, height=6)





k <- 2
### Cut to every other year
dend_labels <- hclust_result$labels[hclust_result$order]
#dend_labels <- substr(cluster_event_df$begin[as.numeric(hclust_result$labels)], 1, 7)

### Add to reorganize cluster numbers
clusters_k <- clusters[,(k-1)]
cluster_number <- clusters_k[match(dend_labels, names(clusters_k))]

### Continue processing labels
dend_labels <- substr(cluster_event_df$begin[as.numeric(dend_labels)], 1, 4)
length_labels <- length(dend_labels)
dend_labels_everynth <- rep(NA, length_labels)
dend_labels_everynth[seq(1, length_labels, 3)] <- dend_labels[seq(1, length_labels, 3)]

#color_list <- pal_d3("category20")(k)[cluster_number]
#color_list <- unique(color_list)
#color_list <- manual_pal[cluster_number]
#color_list <- unique(color_list)
color_list <- c("#FF7F0EFF", "#6A51A3")
#54278F
#6A51A3

length_clusters <- length(hclust_result$height)

cut_line <- (hclust_result$height[length_clusters-k+1] + hclust_result$height[length_clusters-k+2])/2


### Create Dendrogram of all years with a cut at 6 clusters
dend <- hclust_result %>% as.dendrogram %>%
   set("branches_k_color", k=k, value = color_list ) %>%
   set("labels_colors", k=k, value = color_list ) %>%
   set("labels", dend_labels_everynth) %>% set("labels_to_char") %>%
   set("labels_cex", c(0.6)) %>%
    set("branches_lwd", c(0.5))

# Convert to ggdend object
ggd1 <- as.ggdend(dend)

# Change the theme to minimal
p <- ggplot(ggd1, horiz = TRUE, theme = theme_dendro)
p <- p + geom_hline(yintercept=cut_line, linetype="dashed")
p <- p + labs(y = "Sum of Within-Group Distance Squared")
#p <- p + coord_flip(ylim=yscale[3:4])
p

### Save clusters
ggsave(file.path(save_folder, paste0("dendro/",cluster_method,"_dendro_k_",k,"_everythird_new_colors.png")),  width=7.5, height=6, dpi=300)
ggsave(file.path(save_folder, paste0("dendro/",cluster_method,"_dendro_k_",k,"_everythird_new_colors.pdf")),  width=7.5, height=6)
ggsave(file.path(save_folder, paste0("dendro/",cluster_method,"_dendro_k_",k,"_everythird_new_colors.svg")),  width=7.5, height=6)



###########################################################################
## Save Workspace
###########################################################################
save.image(file.path(weber_output_path, "weber_clustering.RData"))

#load(file.path(weber_output_path, "weber_clustering.RData"))


###########################################################################
## Calculate failures for each drought event
###########################################################################
load(file.path(weber_output_path, "weber_storage_output.RData"))
load(file.path(weber_output_path, "weber_delivery_output.RData"))

cluster_stor_df <- inner_join(x=cluster_event_df, y=scenario_df[,c(1,2)], by="data")

cluster_stor_df$upper_ogden_min <- NA
cluster_stor_df$upper_weber_min <- NA
cluster_stor_df$lower_min <- NA
cluster_stor_df$system_min <- NA
cluster_stor_df$trigger <- "None"

for (i in seq(1, dim(cluster_stor_df)[1])){

	event_i <- cluster_stor_df[i,]
	begin_i <- as.Date(event_i$begin)
	end_i <- as.Date(event_i$end)
	
	data_i <- as.character(event_i$data)
	response_i <- as.character(event_i$response)
	
	stor_i <- subset(stor_all, date >= begin_i-60 & date <=end_i+60 & data==data_i & response == response_i)
	demand_i <- subset(demand_deliv_df, date >= begin_i-60 & date <=end_i+60 & data==data_i & response == response_i)
	
	if (length(stor_i$total_res)>0){
	cluster_stor_df$upper_ogden_min[i] <- min(stor_i$upper_ogden, na.rm=TRUE)
	cluster_stor_df$upper_weber_min[i] <- min(stor_i$upper_weber, na.rm=TRUE)
	cluster_stor_df$lower_min[i] <- min(stor_i$lower, na.rm=TRUE)
	cluster_stor_df$system_min[i] <- min(stor_i$total_res, na.rm=TRUE)
	
		if(sum(stor_i$moderate > 0) > 0) {cluster_stor_df$trigger[i] <- "Moderate"}
		if(sum(stor_i$severe > 0) > 0) {cluster_stor_df$trigger[i] <- "Severe"}
		if(sum(stor_i$extreme > 0) > 0) {cluster_stor_df$trigger[i] <- "Extreme"}
	
	cluster_stor_df$upper_ogden_short[i] <- max(subset(demand_i, node=="upper_ogden")$demand_shortage, na.rm=TRUE)
	cluster_stor_df$upper_weber_short[i] <- max(subset(demand_i, node=="upper_weber")$demand_shortage, na.rm=TRUE)
	cluster_stor_df$lower_short[i] <- max(subset(demand_i, node=="lower")$demand_shortage, na.rm=TRUE)
	cluster_stor_df$system_short[i] <- max(subset(demand_i, node=="system")$demand_shortage, na.rm=TRUE)
	
	}
		
}


cluster_stor_df$k_11 <- as.character(cluster_stor_df$k_11)
cluster_stor_df$trigger <- factor(cluster_stor_df$trigger, levels=c("None", "Extreme", "Severe", "Moderate")) 
 
ggplot(subset(cluster_stor_df, response=="Base"), aes(x=k_11, fill=trigger)) + geom_bar() + scale_fill_manual(values=c(NA, "red", "orange", "yellow"))

### 21 of 30 have storage failures

plot_df <- subset(cluster_stor_df, response=="Base")
plot_df$k_11 <- factor(plot_df$k_11, as.character(seq(1,11))) 
ggplot(plot_df, aes(x=k_11, y=system_min, fill=k_11)) + geom_boxplot()  + scale_fill_manual(values=manual_pal)


ggplot(plot_df, aes(x=k_11, y=system_min*0.02831685, fill=k_11)) +  geom_hline(yintercept = c(trigger_df$moderate[9], trigger_df$severe[9], trigger_df$extreme[9])*0.02831685, linetype="dotted") + geom_boxplot()  + geom_jitter(width = 0.1, colour="grey20", shape=1, alpha=0.4) + scale_fill_manual(values=manual_pal) + theme_classic_new() + scale_x_discrete(name="Cluster") + scale_y_continuous(name="Minimum System Storage (m3)", breaks=seq(0,16000,2000)) + theme(legend.position="none")

dir.create(file.path(save_folder, "storage"))

ggsave(file.path(save_folder, paste0("storage/","cluster_11_total_storage.png")),  width=7.5, height=5.5, dpi=300)

ggplot(plot_df, aes(x=k_11, y=system_short*0.02831685, fill=k_11)) + geom_boxplot()  + geom_jitter(width = 0.1, colour="grey20", shape=1, alpha=0.4) + scale_fill_manual(values=manual_pal) + theme_classic_new() + scale_y_continuous(name="Maximum System Delivery Failure (m3)") + scale_x_discrete(name="Cluster") +  theme(legend.position="none")

ggsave(file.path(save_folder, paste0("storage/","cluster_11_total_shortage.png")),  width=7.5, height=5.5, dpi=300)


### Lower Weber delivery only shorts for orange cluster
### Upper Weber for orange, blue, and just barely for green.  Interesting, not red even though it shows higher shortage and storage failures
### Upper Ogden even for purple


 + coord_cartesian(ylim=c(0,50))


scale_y_continuous(name="Maximum System Shortage (m3)", breaks=seq(0,16000,2000)) +



ggplot(plot_df, aes(x=k_11, y=lower_min, fill=k_11)) + geom_boxplot()  + geom_jitter(width = 0.1, colour="grey20", shape=1, alpha=0.4) + scale_fill_manual(values=manual_pal) + theme_classic_new() + scale_x_discrete(name="Cluster") + scale_y_continuous(name="Lower Weber Storage") + theme(legend.position="none")


ggplot(plot_df, aes(x=k_11, y=upper_ogden_min, fill=k_11)) + geom_boxplot()  + scale_fill_manual(values=manual_pal)


ggplot(plot_df, aes(x=k_11, y=upper_weber_min, fill=k_11)) + geom_boxplot()  + scale_fill_manual(values=manual_pal)

### 11 made most sense without begin/end






head(cluster_event_df)

load(file.path(weber_output_path, "weber_storage_output.RData"))

plot_df <- drought_event_summary[!is.na(drought_event_summary$system_min), ]
plot_df_subset <- plot_df[plot_df$data %in% c("paleo", "observed", "hd"),]


events_base <- subset(drought_event_summary, )
head(cluster_event_df)


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








	
### Plot time series	
	
char_years_list <- c("1932-08", "1999-11", "1976-06", "2009-07",  "1987-05", "1957-10", "2012-05", "1993-11")

#manual_pal <- c(pal_d3(palette="category20")(4), brewer.pal(7, "Purples")[seq(3,7)])
#manual_pal <- manual_pal[c(1, 5, 6, 2, 7, 4, 8, 9, 3)]
#manual_pal <- c(manual_pal, manual_pal[1], manual_pal[4])

cluster_event_df <- drought_event_df
rownames(cluster_event_df) <- substr(cluster_event_df$begin,1,7)

for (n in seq(1,length(char_years_list))){

row_name <- char_years_list[n]

row_n <- rownames(cluster_event_df) %in% row_name

yup <- cluster_event_df[row_n,]

begin_test <- which(recon_df$date %in% yup$begin) -1
end_test <- which(recon_df$date %in% yup$end) + 1

subset_ts <- recon_df[seq(begin_test,end_test),]

base_date <- as.Date(paste0(subset_ts$water_year[1]-1, "-10-01"))

yup <- as.Date(subset_ts$date) - base_date
subset_ts$relative_date <-  yup + as.Date("0000-10-01")
subset_ts$date <- as.Date(subset_ts$date)

first_date <- min(subset_ts$date)
last_date <- first_date + years(14)

breaks_qtr <- seq(as.Date("1400-1-1"), by = "6 months", length.out = 6000)
labels_year = format(seq(from = min(breaks_qtr), to = max(breaks_qtr), by = "1 years"), "%Y")
labs = c(sapply(labels_year, function(x) {
    c(x, rep("", 1))
    }))    
#labs <- labs[seq(1,30)]
    
theme_ts <- theme_classic_correct()+ theme(axis.text.x = element_text(angle = 30, hjust = 1))

p <- ggplot(subset_ts, aes(x=date, y=percentile))
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="grey40")
p <- p + geom_line(colour="blue")#colour=manual_pal[n])
p <- p  + theme_classic()
#p <- p + scale_x_date(name="Date", date_labels = "%Y", breaks = seq(as.Date("1400-1-1"), by = "24 months", length.out = 1000),  date_minor_breaks = "12 months")
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent)
p <- p + coord_cartesian(ylim=c(0,0.62),  xlim=c(first_date, last_date))
p
#
save_folder <- file.path(write_output_path, "heir_agglom/ts")
dir.create(save_folder, showWarnings=FALSE)

### Save results
ggsave(paste0(row_name, "_obs_ts.png"),  p, width=6, height=3, dpi=300)
ggsave(paste0(row_name, "_obs_ts.pdf"),  p, width=6, height=3)
#ggsave(paste0(row_name, "_obs_ts.svg"),  p, width=6, height=3)

}




weber_recon_df <- read.csv("/run/media/jhstagge/Data/Documents/work/output/paleo_weber/apr_model/10128500_apr_model_enso_pca_impute_std_concur_reconst_ts_rec_local.csv")
weber_recon_df$percentile <- pnorm(weber_recon_df$monthly_norm)
recon_orig_date <- as.Date(weber_recon_df$date) 

weber_recon_df$date <- recon_orig_date + years(2012-1706)
weber_recon_df_plot <- subset(weber_recon_df, date < as.Date("2030-01-01"))
#recon_adjust_df <- recon_df 
#recon_adjust_df$date <- recon_adjust_df$date - years(306)

p <- ggplot(weber_recon_df_plot, aes(x=date, y=percentile))
p <- p + geom_hline(yintercept=0.5, linetype="longdash", color="black")
p <- p + geom_line(data=recon_df, colour="grey20", alpha=0.5)
p <- p + geom_line(colour="#1F77B4FF", size=0.7)#colour=manual_pal[n])
p <- p  + theme_classic_new(14)
p <- p + scale_x_date(labels = labs, breaks = breaks_qtr, name = "Year") 
p <- p + scale_y_continuous(name="Flow Percentile",  labels = scales::percent, expand=c(0,0))
p <- p + coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2012-01-01"), as.Date("2025-01-01")))
p
#

ggsave("1706_drought_vs_2012.png",  p, width=8, height=4, dpi=600)

