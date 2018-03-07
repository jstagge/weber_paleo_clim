readNWISdv

readNWISdv




 site_id <- '04085427'
     startDate <- '2012-01-01'
     endDate <- '2012-06-30'
     pCode <- '00060'
     ## Not run:
     
     rawDailyQ <- readNWISdv(site_id,pCode, startDate, endDate)




site_id_list <- c("10128500")
site_name_list <- c("Weber River")
param_cd <- "00060"

     startDate <- '1800-01-01'
     endDate <- '2018-01-01'


n <- 1
site_id <- site_id_list[n]

rawDailyQ <- readNWISdv(site_id,param_cd, startDate, endDate)


plot(rawDailyQ$Date, rawDailyQ$X_00060_00003, type="l", ylim=c(0,5000))
