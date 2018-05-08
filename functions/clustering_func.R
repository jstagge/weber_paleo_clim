###########################################################################
## Determine clusters
###########################################################################

calc_clusters <- function(clustering_tree, cluster.list) {

### Calculate length of sites
site_length <- length(clustering_tree$labels)

### Create dataframe to hold cluster results
clusters <- matrix(NA,site_length,length(cluster.list))
colnames(clusters) <- paste("k_",cluster.list, sep="")
rownames(clusters) <- seq(1,site_length)
#This causes failure if data is not in order
#rownames(clusters) <- substr(drought_event_df$begin, 1, 7)

cluster_goodness <- matrix(NA,length(cluster.list), 2)
rownames(cluster_goodness) <- paste("k_",cluster.list, sep="")
colnames(cluster_goodness) <- c("k", "sil_width")

### Loop through the cluster numbers in cluster.list and save results
for (k in cluster.list) {
### Name clusters in order of dendrogram
cluster_number <- cutree(clustering_tree, k, order_clusters_as_data=FALSE )
### Reorder back to original order
clusters[,k-1] <- as.factor(cluster_number[order(as.numeric(names(cluster_number)))])

cluster_goodness[k-1,] <- c(k, summary(silhouette(cutree(clustering_tree, k), distMatrix))$avg.width)
}

### Special re-ordering step to ensure that each step uses the same clustering numbers as
### before. Newest cluster number is assigned to the smaller of the split clusters
for (k in cluster.list[seq(2,length(cluster.list))]) {

### Extract clusters for k and k-1
clusters_k <- data.frame(clusters[,c(k-2,k-1)])
clusters_k$new <- NA

 for (j in seq(1,k-1)) {
 	cluster_test <- clusters_k[,1]==j 
 	insert_cluster <- clusters_k[cluster_test,]
 	
 	### if all of new column is the same, use same name as previous
	if (sum(cluster_test)==1) {
		clusters_k$new[cluster_test] <- insert_cluster[1,1]
	} else if (var(insert_cluster[,2]) == 0 ){
 		clusters_k$new[cluster_test] <- insert_cluster[1,1]
 	} else {
 		freq_table <- table(insert_cluster[,2])
 		freq_table_sorted <- sort(freq_table)
 		
 		former_names <- names(freq_table_sorted)
		
		### Replace less frequent with k and more frequent with j
 		clusters_k$new[clusters_k[,2] == former_names[1]]  <- k
 		clusters_k$new[clusters_k[,2] == former_names[2]] <- j
 	}	
 }
clusters[,(k-1)] <- clusters_k$new
}

return(list(clust=clusters, goodness=cluster_goodness))

}



###########################################################################
## Goodness of fit plot
###########################################################################
plot_cluster_goodness <- function (cluster_goodness, save_directory) {
### Plot goodness of fit
cluster_goodness <- as.data.frame(cluster_goodness)
p <- ggplot(cluster_goodness, aes(k,sil_width)) 
p <- p + theme_classic_correct()
p <- p + geom_line()
p <- p + geom_point()
p <- p + xlab("Number of Clusters")  #insert the x-axis title
p <- p + ylab("Avg. Silhouette Width")  #insert the y-axis title
p <- p + scale_x_continuous(breaks=seq(0,19,2))

### Save results
ggsave(file.path(save_directory, "goodness/silWidth_bycluster.png"),  p, width=3.5, height=2.7, dpi=300)
ggsave(file.path(save_directory, "goodness/silWidth_bycluster.pdf"),  p, width=3.5, height=2.7)
ggsave(file.path(save_directory, "goodness/silWidth_bycluster.svg"),  p, width=3.5, height=2.7)

return(p)
}


###########################################################################
## Silhouette plot
###########################################################################
plot_sil <- function(cluster_df, dist, k, save_directory) {
	
	cluster <- cluster_df[,(k-1)]
	
	## Generate colors
	plot_col <- pal_d3("category20")(max(cluster))
	
	### Set save location
	save_location <- file.path(save_directory, paste0("goodness/",cluster_method,"_sil_k_",k))
	
	### Save png
	png(filename=paste0(save_location, ".png"), 
    	type="cairo",
    	units="in", 
    	width=5, 
    	height=5, 
    	pointsize=12, 
    	res=300)
    plot(silhouette(cluster, distMatrix), col = plot_col, border=NA)
	dev.off()
	
	### Save pdf
	pdf(file=paste0(save_location, ".pdf"), 
    	width=5, 
    	height=5, 
    	pointsize=12)
    plot(silhouette(cluster, distMatrix), col = plot_col, border=NA)
	dev.off()	
		
	### Save svg
	svg(filename=paste0(save_location, ".svg"), 
    	width=5, 
    	height=5, 
    	pointsize=12)
    plot(silhouette(cluster, distMatrix), col = plot_col, border=NA)
	dev.off()
}


###########################################################################
## Cluster plot
###########################################################################
plot_clusters <- function(plot_df, k, save_directory) {
	p <- ggplot(cluster_event_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	#p <- p + geom_point(aes(colour=timeperiod))
	p <- p + geom_text(aes(colour=as.factor(get(paste0("k_",k)))), size=2.7)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,16,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	p <- p + scale_color_d3(name="Cluster", palette = "category20")
	p <- p + theme_classic(10)
	p <- p + theme(legend.position = c(0.85, 0.85))
	
	if (k > 4) {
	p <- p + guides(colour=guide_legend(ncol=2))
	}
	if (k > 8) {
	p <- p + guides(colour=guide_legend(ncol=3))
	}	
	if (k > 12) {
	p <- p + guides(colour=guide_legend(ncol=4))
	}		

### Save results
ggsave(file.path(save_directory, paste0("clusters/",cluster_method,"_clust_k_",k,".png")),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(save_directory, paste0("clusters/",cluster_method,"_clust_k_",k,".pdf")),  p, width=6.5, height=5)
ggsave(file.path(save_directory, paste0("clusters/",cluster_method,"_clust_k_",k,".svg")),  p, width=6.5, height=5)

### Plot with only points
	p <- ggplot(cluster_event_df, aes(x=dura_months/12, y=min_perc*100, label=substr(begin,1,4)))
	p <- p + geom_point(aes(colour=as.factor(get(paste0("k_",k)))), size=2.5)
	#p <- p + geom_text(aes(colour=as.factor(get(paste0("k_",k)))), size=2.7)
	p <- p + scale_x_continuous(name="Drought Duration (Years)", breaks=seq(0,16,2))
	p <- p + scale_y_continuous(name="Min Flow Percentile")
	p <- p + scale_color_d3(name="Cluster", palette = "category20")
	#p <- p + scale_shape_discrete(solid=T)
	p <- p + theme_classic_new(12)
	p <- p + theme(legend.position = c(0.85, 0.85))
		
	if (k > 4) {
	p <- p + guides(colour=guide_legend(ncol=2))
	}
	if (k > 8) {
	p <- p + guides(colour=guide_legend(ncol=3))
	}	
	if (k > 12) {
	p <- p + guides(colour=guide_legend(ncol=4))
	}		

### Save results
ggsave(file.path(save_directory, paste0("clusters/",cluster_method,"_clust_k_",k,"_points.png")),  p, width=6.5, height=5, dpi=300)
ggsave(file.path(save_directory, paste0("clusters/",cluster_method,"_clust_k_",k,"_points.pdf")),  p, width=6.5, height=5)
ggsave(file.path(save_directory, paste0("clusters/",cluster_method,"_clust_k_",k,"_points.svg")),  p, width=6.5, height=5)


	return(p)
}
	
	

###########################################################################
## Dendrogram  plot
###########################################################################
plot_dendrogram <- function(clustering_tree, k, save_directory, yscale, label_offsets) {

### Cut to every other year
dend_labels <- clustering_tree$labels[clustering_tree$order]

### Add to reorganize cluster numbers
clusters_k <- clusters[,(k-1)]
cluster_number <- clusters_k[match(dend_labels, names(clusters_k))]

### Continue processing labels
dend_labels <- substr(dend_labels,1,4)
length_labels <- length(dend_labels)
dend_labels_everynth <- rep(NA, length_labels)
dend_labels_everynth[seq(1, length_labels, 3)] <- dend_labels[seq(1, length_labels, 3)]

color_list <- pal_d3("category20")(k)[cluster_number]
color_list <- unique(color_list)

### Create Dendrogram of all years with a cut at 6 clusters
dend <- clustering_tree %>% as.dendrogram %>%
  set("branches_k_color", k=k, value = color_list ) %>%
   set("labels_colors", k=k, value = color_list ) %>%
   set("labels", dend_labels) %>% set("labels_to_char") %>%
   set("labels_cex", c(0.6)) %>%
    set("branches_lwd", c(0.5))
# plot the dend in usual "base" plotting engine:
#plot(dend)
# Convert to ggdend object
ggd1 <- as.ggdend(dend)
 
length_clusters <- length(clustering_tree$height)

cut_line <- (clustering_tree$height[length_clusters-k+1] + clustering_tree$height[length_clusters-k+2])/2

# Change the theme to minimal
p <- ggplot(ggd1, horiz = TRUE, theme = theme_dendro, offset_labels=label_offsets)
p <- p + geom_hline(yintercept=cut_line, linetype="dashed")
p <- p + labs(y = "Sum of Within-Group Distance Squared")
p <- p + coord_flip(ylim=yscale[1:2])
#p

### Save clusters
dir.create(file.path(save_directory, "dendro"), showWarnings=FALSE)
ggsave(file.path(save_directory, paste0("dendro/",cluster_method,"_dendro_k_",k,"_allyears.png")),  width=7.5, height=6, dpi=300)
ggsave(file.path(save_directory, paste0("dendro/",cluster_method,"_dendro_k_",k,"_allyears.pdf")),  width=7.5, height=6)
ggsave(file.path(save_directory, paste0("dendro/",cluster_method,"_dendro_k_",k,"_allyears.svg")),  width=7.5, height=6)


### Create Dendrogram of all years with a cut at 6 clusters
dend <- clustering_tree %>% as.dendrogram %>%
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
p <- p + coord_flip(ylim=yscale[3:4])
#p

### Save clusters
ggsave(file.path(save_directory, paste0("dendro/",cluster_method,"_dendro_k_",k,"_everythird.png")),  width=7.5, height=6, dpi=300)
ggsave(file.path(save_directory, paste0("dendro/",cluster_method,"_dendro_k_",k,"_everythird.pdf")),  width=7.5, height=6)
ggsave(file.path(save_directory, paste0("dendro/",cluster_method,"_dendro_k_",k,"_everythird.svg")),  width=7.5, height=6)

return(p)

}


	

###########################################################################
## Parallel Coordinates  plot
###########################################################################
plot_par_coord <- function(full_cluster_df, k) {

### Determine column index for cluster
column_index <- which(names(cluster_event_df) == paste0("k_", k))

### Create plotting dataframe
plot_data <- data.frame(cluster=cluster_event_df[, column_index], cluster_event_df[,c(seq(3,13), 15, 14, 16)])
plot_data <- plot_data[complete.cases(plot_data),]


### Create plot
parcoords(plot_data, reorderable = TRUE, brushMode = "1D-axes-multi",rownames = F, width=2000,alphaOnBrushed=0.1, 
color = list(
      colorBy = "cluster"
      ,colorScale = htmlwidgets::JS("d3.scale.category10()")
      )
)

}



# function to find medoid in cluster i
### Copied from here: https://www.biostars.org/p/11987/
clust.medoid = function(i, distmat, clusters) {
    ind = (clusters == i)

	if (sum(ind) == 0 ){
		return(NA)
	} else if (sum(ind) ==1 ){
		names(clusters[ind])
	} else {
    	names(which.min(rowSums( distmat[ind, ind] )))
    # c(min(rowMeans( distmat[ind, ind] )))
    }
}



