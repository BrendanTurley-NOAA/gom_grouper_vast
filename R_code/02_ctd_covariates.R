library(lubridate)

# setwd("C:/Users/brendan.turley/Desktop/FL_habs/data/groomed")
setwd("~/Desktop/professional/projects/Postdoc_FL/data/ctd/groomed")

### version 3 of combined_groomed data standardized subsetting and removal of bad data
combined <- read.csv('ctd_combined_groomed_v6.csv')
combined$date <- as.Date(combined$date)
combined_stations <- read.csv('ctd_combined_stations_groomed_v6.csv')
combined_stations$date <- as.Date(combined_stations$date)
combined <- combined[which(year(combined$date)>=2004),]
combined_stations <- combined_stations[which(year(combined_stations$date)>=2004),]
ind <- which(combined$z_crm<=50)
combined <- combined[ind,]
combined_stations <- combined_stations[ind,]
names(combined)

covariate_data <- combined[which(year(combined$date)>=2010 & month(combined$date)==6 |
                                   year(combined$date)>=2010 & month(combined$date)==7),]
ind <- c(14:15,20:22)
names(covariate_data)[ind]
covariate_data <- covariate_data[,ind]
covariate_data <- covariate_data[,c(4,5,3,1,2)]
covariate_data$date <- year(covariate_data$date)
names(covariate_data) <- c('Lat','Lon','Year','Bot_DO','Bot_Temp')

setwd('~/Desktop/professional/projects/Postdoc_FL/data/grouper/')
write.csv(covariate_data,'covariate_data.csv')
