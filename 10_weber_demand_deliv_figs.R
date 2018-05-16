
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
## Open Output File
###########################################################################
load(file.path(weber_output_path, "weber_delivery_output.RData"))

load(file.path(weber_output_path, "weber_clustering.RData"))



###########################################################################
## Set initial values
###########################################################################
res_colors <- cb_pal(pal="wong", 3, sort=FALSE)


cc_colors <- c("#0072B2", "#56B4E9", "#E69F00" , "#D55E00")

#############################################################
###  Quick Check Plots
#############################################################
demand_deliv_plot <- demand_deliv_df[demand_deliv_df$data %in% c("paleo", "observed", "hd"),]

### To plot
p <- ggplot(demand_deliv_plot, aes(x=date, y=shortage))
p <- p + geom_line(aes(group=data), size=0.12)
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (ac-ft)")
p <- p + theme(legend.position="bottom")
p <- p + theme(axis.text.x = element_text(angle = 40, hjust = 1))
p <- p + facet_wrap(~variable, nrow = 4)
p


### Save figures
ggsave(file.path(write_output_base_path,"Shortage_all_facets.png"),  p, width=14, height=6.5, dpi=300)
ggsave(file.path(write_output_base_path,"Shortage_all_facets.pdf"),  p, width=14, height=6.5)
ggsave(file.path(write_output_base_path,"Shortage_all_facets.svg"),  p, width=14, height=6.5)


### To plot
p <- ggplot(demand_deliv_plot, aes(x=date, y=shortage_perc))
p <- p + geom_line(aes(group=data), size=0.12)
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shorage (% of Demand)", breaks=seq(0,1,0.1), labels=percent)
p <- p + theme(legend.position="bottom")
p <- p + facet_wrap(~variable, nrow = 4)
p


### Save figures
ggsave(file.path(write_output_base_path,"Shortage_perc_all_facets.png"),  p, width=14, height=6.5, dpi=300)
ggsave(file.path(write_output_base_path,"Shortage_perc_all_facets.pdf"),  p, width=14, height=6.5)
ggsave(file.path(write_output_base_path,"Shortage_perc_all_facets.svg"),  p, width=14, height=6.5)


### To plot
p <- ggplot(demand_deliv_plot, aes(x=date))
p <- p + geom_line(aes(group=data, y=demand), size=0.12, colour="grey20")
p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
#p <- p + geom_line(aes(group=data, y=request), size=0.12, colour="red")
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
p <- p + facet_wrap(~variable, nrow = 4)
p



#############################################################
###  Check SA1
#############################################################
sa1_plot <- subset(demand_deliv_plot, variable=="SA1")

### To plot
p <- ggplot(sa1_plot, aes(x=date))
p <- p + geom_line(aes(group=data, y=demand), size=0.12, colour="grey20")
p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p




#############################################################
###  Check SA12
#############################################################
sa1_plot <- subset(demand_deliv_plot, variable=="SA12")

### To plot
p <- ggplot(sa1_plot, aes(x=date))
p <- p + geom_line(aes(group=data, y=demand), size=0.12, colour="grey20")
p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p


### To plot
p <- ggplot(sa1_plot, aes(x=month, group=year))
p <- p + geom_line(aes(y=shortage), size=0.12, colour="grey20")
p <- p + theme_classic_new()
#p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p



#############################################################
###  Check SA14
#############################################################
sa1_plot <- subset(demand_deliv_plot, variable=="SA14")

### To plot
p <- ggplot(sa1_plot, aes(x=date))
p <- p + geom_line(aes(group=data, y=demand), size=0.12, colour="grey20")
p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p


### To plot
p <- ggplot(sa1_plot, aes(x=month, group=year))
p <- p + geom_line(aes(y=shortage), size=0.12, colour="grey20")
p <- p + theme_classic_new()
#p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (Demand - Delivery)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p




#############################################################
###  Check SA10
#############################################################
sa10_plot <- subset(demand_deliv_plot, variable=="SA10")
#sa10_plot <- sa10_plot[,c(4, 5, 8, 10)] #9, 
#sa10_plot <- melt(sa10_plot, c("date", "data"))

### To plot
p <- ggplot(sa10_plot, aes(x=date))
p <- p + geom_line(aes(group=data, y=demand), size=0.12, colour="grey20")
p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
p <- p + theme_classic_new()
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Delivery (Ac-ft)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p

### Save figures
ggsave(file.path(write_output_base_path,"SA10_deliv_demand.png"),  p, width=8, height=4.5, dpi=300)
ggsave(file.path(write_output_base_path,"SA10_deliv_demand.pdf"),  p, width=8, height=4.5)
ggsave(file.path(write_output_base_path,"SA10_deliv_demand.svg"),  p, width=8, height=4.5)



sa10_plot <- subset(demand_deliv_plot, variable=="SA10")
sa10_plot <- sa10_plot[,c(4, 5, 8, 10)] #9, 
sa10_plot <- melt(sa10_plot, c("date", "data"))
### Sort by node and then id
sa10_plot <- sa10_plot[order(sa10_plot$variable, sa10_plot$data ,sa10_plot$date),] 

### To plot
p <- ggplot(sa10_plot, aes(x=date))
p <- p + geom_line(aes(group=data, colour=variable, y=value), size=0.12)
#p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
p <- p + theme_classic_new()
p <- p + scale_colour_manual(values=c("black", "red"))
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Delivery (ac-ft)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p




#############################################################
###  Plot against clusters
#############################################################













ggplot(subset(hash_full, response=="Base" & data=="paleo"), aes(x=node, y=reliability)) + geom_col(position="identity") + theme_classic() + coord_cartesian(ylim=c(0.75, 1))








ggplot(subset(hash_full, response=="Base" & data=="paleo"& short_type=="demand"), aes(x=node, y=reliability)) + geom_col(position="identity") + theme_classic() + coord_cartesian(ylim=c(0.75, 1))

ggplot(subset(hash_full, node %in% c("upper_weber", "upper_ogden", "lower") & data=="paleo"& short_type=="demand"), aes(x=node, y=reliability, fill=response)) + geom_bar(position = "dodge", stat="identity") + theme_classic() + coord_cartesian(ylim=c(0.75, 1))

p <- ggplot(subset(hash_full, node %in% c("upper_weber", "upper_ogden", "lower") & response=="Base"& short_type=="demand"), aes(x=node, y=reliability, fill=data)) + geom_bar(position = "dodge", stat="identity") + theme_classic() + coord_cartesian(ylim=c(0.5, 1))
p + scale_fill_manual(values=c("#7fc97f", "grey40", "grey80", cc_colors))




p <- ggplot(subset(hash_full, node %in% c("upper_weber", "upper_ogden", "lower") & response=="Base" & short_type=="demand"), aes(x=node, y=vulnerability, fill=data)) + geom_bar(position = "dodge", stat="identity") + theme_classic()
p + scale_fill_manual(values=c("#7fc97f", "grey40", "grey80", cc_colors))



p <- ggplot(subset(hash_full, node %in% c("upper_weber", "upper_ogden", "lower") & response=="Base" & short_type=="demand"), aes(x=node, y=resilience, fill=data)) + geom_bar(position = "dodge", stat="identity") + theme_classic()
p + scale_fill_manual(values=c("#7fc97f", "grey40", "grey80", cc_colors))



p <- ggplot(subset(hash_full, node %in% c("upper_weber", "upper_ogden", "lower") & response=="Base"), aes(x=node, y=reliability, fill=data)) + geom_bar(position = "dodge", stat="identity") + theme_classic()
p + scale_fill_manual(values=c("#7fc97f", "grey40", "grey80", cc_colors)) + facet_wrap(~short_type, nrow=2)




p <- ggplot(subset(hash_full, response=="Base"), aes(x=node, y=reliability, fill=data)) + geom_bar(position = "dodge", stat="identity") + theme_classic()
p + scale_fill_manual(values=c("#7fc97f", "grey40", "grey80", cc_colors)) + facet_wrap(~short_type, nrow=2)






p <- ggplot(subset(demand_deliv_df, response=="Base"), aes(x=date, y=demand_shortage_perc, group=data)) + facet_wrap(~node) + geom_line() + theme_classic_new()
p

p + geom_line(data=subset(demand_deliv_df, response=="ChalkCreekRes"), aes(y=demand_shortage_perc), colour="red")

p <- ggplot(subset(demand_deliv_df, response=="Base"), aes(x=date, y=demand_shortage, group=data)) + facet_wrap(~node) + geom_line()
p + geom_line(data=subset(demand_deliv_df, response=="ChalkCreekRes"), aes(y=demand_shortage), colour="red")

p <- ggplot(subset(demand_deliv_df, response=="Base"), aes(x=date, y=request_shortage_perc, group=data)) + facet_wrap(~node) + geom_line()
p
p + geom_line(data=subset(demand_deliv_df, response=="ChalkCreekRes"), aes(y=request_shortage_perc), colour="red")









p <- ggplot(subset(shortage_system, data=="hd" | data=="paleo" | data=="observed"), aes(x=date, colour=data, y=./1000))
p <- p + geom_line()
#p <- p + geom_line(aes(group=data, y=deliv), size=0.12, colour="blue")
p <- p + theme_classic_new()
#p <- p + scale_colour_manual(name="Scenario", values= c("grey30", "#CC79A7", "#E69F00", "#0072B2"), labels=c("Observed", "Reconstr", "WD", "WW"), guide = guide_legend())
p <- p + scale_colour_manual(name="Scenario", values= c("#D55E00", "grey30", "#CC79A7"), breaks=c("hd", "observed", "paleo"), labels=c("HD", "Observed", "Reconstr"), guide = guide_legend())
#p <- p + scale_colour_manual(values=c("black", "red"))
p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Shortage (1,000 ac-ft)")
p <- p + theme(legend.position="bottom")
#p <- p + facet_wrap(~variable, nrow = 4)
p


### Save figures
ggsave(file.path(write_output_base_path,"paleo_future_system_shortage_acft_hot.png"),  p, width=8, height=4.5, dpi=300)
ggsave(file.path(write_output_base_path,"paleo_future_system_shortage_acft_hot.pdf"),  p, width=8, height=3.5)
ggsave(file.path(write_output_base_path,"paleo_future_system_shortage_acft_hot.svg"),  p, width=8, height=3.5)



#############################################################
###  Plot requests
#############################################################
### It looks like there is a significant difference between paleo and climate change requests at SA 1 and SA10

p <- ggplot(request_all, aes(x=month, group=wy))
p <- p + geom_line(aes(y=SA1, colour=data), size=0.12)
p <- p + theme_classic_new()
#p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Request")
p <- p + theme(legend.position="bottom")
p <- p + facet_wrap(~data, nrow = 4)
p


### Make months a factor
demand_df$month_adj_forplot <- c(seq(4, 12), seq(1,3))

### Make months a factor
demand_df$month_adj_forplot <- c(seq(4, 12), seq(1,3))

p <- ggplot(request_all, aes(x=month))
p <- p + geom_boxplot(aes(y=SA1, colour=data))
#p <- p + geom_line(data=demand_df, aes(x=month_adj_forplot, y=SA1), colour="black")
p <- p + theme_classic_new()
#p <- p + coord_cartesian(xlim=c(as.Date("1428-01-01"), as.Date("2070-01-01")), expand=FALSE) #ylim=c(0,1), 
#p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="Request")
p <- p + theme(legend.position="bottom")
p

### Funky - April requests are always exactly the same for entire Paleo run.
### Climate change runs vary across scenarios, but each climate change run has the same number, but different for each scenario.


#############################################################
###  Do calculations regionally
#############################################################


#############################################################
###  Hashimoto Figures
#############################################################
p <- ggplot(subset(hash_df, data=="SA2"), aes(x=site, y=resilience))
p <- p + geom_bar()
p <- p + theme_classic()
p


p <- ggplot(hash_df, aes(x=site, y=resilience, fill=data))
p <- p +  geom_bar(position="dodge", stat="identity")
p <- p + theme_classic()
p


p <- ggplot(hash_df, aes(x=site, y=reliability, fill=data))
p <- p +  geom_bar(position="dodge", stat="identity")
p <- p + scale_fill_manual(name="Scenario", values= c("green", "black", "grey30", "orange", "red"), limits=c("paleo", "observed", "base", "hw", "hd"), guide = guide_legend())
p <- p + theme_classic()
p

p <- ggplot(subset(hash_df, data= c("paleo", "observed", "base", "hw", "hd")), aes(x=site, y=resilience, fill=data))
p <- p +  geom_bar(position="dodge", stat="identity")
p <- p + scale_fill_manual(name="Scenario", values= c("green", "black", "grey30", "orange", "red"), limits=c("paleo", "observed", "base", "hw", "hd"), guide = guide_legend())
p <- p + theme_classic()
p



p <- ggplot(hash_df[hash_df$data %in% c("paleo", "observed", "base", "hw", "hd"), ], aes(x=site, y=resilience, fill=data))
p <- p +  geom_bar(position="dodge", stat="identity")
p <- p + scale_fill_manual(name="Scenario", values= c("green", "black", "grey30", "orange", "red"), limits=c("paleo", "observed", "base", "hw", "hd"), guide = guide_legend())
p <- p + theme_classic()
p


p <- ggplot(hash_df[hash_df$data %in% c("paleo", "observed", "base", "hw", "hd"), ], aes(x=site, y=reliability, fill=data))
p <- p +  geom_bar(position="dodge", stat="identity")
p <- p + scale_fill_manual(name="Scenario", values= c("green", "black", "grey30", "orange", "red"), limits=c("paleo", "observed", "base", "hw", "hd"), guide = guide_legend())
p <- p + theme_classic()
p


p <- ggplot(hash_df[hash_df$data %in% c("paleo", "observed", "base", "hw", "hd"), ], aes(x=site, y=vulnerability, fill=data))
p <- p +  geom_bar(position="dodge", stat="identity")
p <- p + scale_fill_manual(name="Scenario", values= c("green", "black", "grey30", "orange", "red"), limits=c("paleo", "observed", "base", "hw", "hd"), guide = guide_legend())
p <- p + theme_classic()
p











