rm(list=ls())
gc()

library(DHARMa)
library(fields)
library(lubridate)
library(ncdf4)
library(raster)
library(rgdal)
library(sp)
library(sf)
library(tidyr)
library(VAST)

### set model run version and create new folder to save results
run <- 'vmod6'
# pthwy <- paste0('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/',run)
# if (file.exists(pthwy)) {
#   cat("The folder already exists")
# } else {
#   dir.create(pthwy)
# }

setwd('~/Desktop/professional/projects/Postdoc_FL/data/grouper/')
# data <- read.csv('grp_snp_2019.csv')
data <- read.csv('grp_snp_2022.csv')
data$year <- year(data$date)
sort(unique(data$year))
data <- data[data$year!=2021,]
data$VESSEL <- as.factor(data$VESSEL)
data$rg_wt[which(is.na(data$rg_wt))] <- 0
data$rs_wt[which(is.na(data$rs_wt))] <- 0

# par(mfrow=c(2,5))
# for(i in sort(unique(data$year))){
#   tmp <- data[data$year==i,]
#   plot(tmp$DECSLON,tmp$DECSLAT,
#        xlim=range(data$DECSLON,na.rm=T),ylim=range(data$DECSLAT,na.rm=T),
#        cex=(tmp$rg_wt+1))
#   mtext(i)
# }

region <- read.csv('region.csv')


### strata limits
# https://github.com/James-Thorson-NOAA/VAST/issues/176
strata_lon <- aggregate(data$DECSLON,by=list(data$STAT_ZONE),range,na.rm=T)
strata_lon$x <- round(strata_lon$x)
strata_lon$x[5:7,2] <- -82
strata_lon$x[1,2] <- -81
strata_lat <- aggregate(data$DECSLAT,by=list(data$STAT_ZONE),range,na.rm=T)
strata_lat$x <- round(strata_lat$x)
strata_lat$x[1,1] <- 24
strata_lat$x[8,1] <- 28.5

plot(strata_lon$x[,1],strata_lat$x[,1],asp=1,
     xlim=range(strata_lon$x),ylim=range(strata_lat$x))
for(i in 1:nrow(strata_lat)){
  rect(strata_lon$x[i,1],strata_lat$x[i,1],
       strata_lon$x[i,2],strata_lat$x[i,2])
  text(strata_lon$x[i,1],strata_lat$x[i,1],strata_lat$Group.1[i],adj=1,cex=2)
}

strata.limits <- data.frame( 'STRATA' = strata_lon$Group.1,
                             'north_border' = strata_lat$x[,2], 'south_border' = strata_lat$x[,1],
                             'east_border' = strata_lon$x[,2], 'west_border' = strata_lon$x[,1] )

Extrapolation_List = make_extrapolation_info( Region='User',
                                              input_grid = region[,c(2:3,5)],
                                              strata.limits=strata.limits )
print(strata.limits)
colSums( Extrapolation_List$a_el )
which(colSums(Extrapolation_List$a_el)==0)

### ----------------- VAST modeling -----------------
# Make settings (turning off bias.correct to save time for example)
settings = make_settings(n_x = 300, # set higher to avoid logKappa2 boundary issues
                         Region = 'User', 
                         purpose = "index2", 
                         bias.correct = FALSE,
                         knot_method = 'grid',
                         ObsModel = c(1,0), # ?make_data for details
                         strata.limits = strata.limits) # can be used to sum indices of abundance over strata

# setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
# Run model
fit <- fit_model(settings = settings, 
                Lat_i = data$DECSLAT, 
                Lon_i = data$DECSLON, 
                t_i = data$year, 
                b_i = as_units(data$rg_wt,'kg'), 
                a_i = as_units(data$AreaSwept_km2,'km^2'),
                v_i = as.numeric(data$VESSEL)-1,
                # covariate_data #https://github.com/James-Thorson-NOAA/VAST/issues/262
                input_grid=region,
                run_model = T, # make FALSE to see upper (fit$tmb_list$Upper) and lower (fit$tmb_list$Lpper) starting bounds
                working_dir = pthwy
                ) 

results <- plot(fit)

file_save <- paste0('~/Desktop/professional/projects/Postdoc_FL/data/grouper/',as.Date(Sys.time()),'_',run,'_results.RData') 
save.image(file_save)
