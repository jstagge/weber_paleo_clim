
yuppers <- area_df[, names(area_df) %in% c("date", "data", "upper_ogden", "upper_weber", "lower")]
yuppers <- melt(yuppers, id.vars=c("date", "data"))#, measure.vars=c("upper_ogden", "upper_weber", "lower"))
yuppers$data <- factor(yuppers$data)
yuppers$variable <- factor(yuppers$variable)

p <- ggplot(subset(yuppers,data=="paleo"), aes(x=date, y=value, fill=variable))
p <- p + geom_area()
p


p <- ggplot(subset(yuppers,data=="observed"), aes(x=date, y=value, colour=variable))
p <- p + geom_line()
p




p <- p + geom_hline(yintercept=0, size=0.2)
p <- p + geom_hline(yintercept=-max_stor/1000, size=0.4, colour="black", linetype="longdash")
p <- p + geom_line(data=base_df, colour="grey30")
#p <- p + geom_area( data=base_df, position = "identity", alpha=0.5)
p <- p + theme_classic_new()
p <- p + scale_fill_manual(name="Scenario", values= c("grey30", "#CC79A7", "#E69F00", "#0072B2"), labels=c("Observed", "Reconstr", "WD", "WW"), guide = guide_legend())
#p <- p + coord_cartesian(xlim=c(as.Date("1920-01-01"), as.Date("2018-01-01")))
p <- p + scale_x_date(name="Date", breaks=seq(as.Date("1200-01-01"), as.Date("2100-01-01"), by="50 years"), date_labels = "%Y")
p <- p + scale_y_continuous(name="System Storage Deficit (1,000 ac-ft)", breaks=seq(-10000,500,50))
p <- p + theme(legend.position="bottom")
p




