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


### ----------------- VAST modeling -----------------
# Make settings (turning off bias.correct to save time for example)
settings = make_settings(n_x = 300, # set higher to avoid logKappa2 boundary issues
                         Region = 'User', 
                         purpose = "index2", 
                         bias.correct = FALSE,
                         knot_method = 'grid',
                         ObsModel = c(1,0)) # ?make_data for details
                         #strata.limits = data$STAT_ZONE) # can be used to sum indices of abundance over strata

setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
# Run model
fit = fit_model(settings = settings, 
                Lat_i = data$DECSLAT, 
                Lon_i = data$DECSLON, 
                t_i = data$year, 
                b_i = as_units(data$rg_wt,'kg'), 
                a_i = as_units(data$AreaSwept_km2,'km^2'),
                # v_i = as.numeric(data$VESSEL)-1,
                input_grid=region,
                run_model = T # make FALSE to see upper (fit$tmb_list$Upper) and lower (fit$tmb_list$Lpper) starting bounds
                ) 

vmod3 <- plot(fit)

file_save <- paste0('~/Desktop/professional/projects/Postdoc_FL/data/grouper/',as.Date(Sys.time()),'_vmod3_results.RData') 
save.image(file_save)
