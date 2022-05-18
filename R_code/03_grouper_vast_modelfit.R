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
run <- 'vmod14'
# pthwy <- paste0('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/',run)
# if (file.exists(pthwy)) {
#   cat("The folder already exists")
# } else {
#   dir.create(pthwy)
# }

setwd('~/Desktop/professional/projects/Postdoc_FL/data/grouper/')
# data <- read.csv('grp_snp_2019.csv')
# data <- read.csv('grp_snp_2022.csv')
data <- read.csv('grp_snp_2022_v2.csv')
data$year <- year(data$date)
data$VESSEL <- as.factor(data$VESSEL)
# sort(unique(data$STAT_ZONE)) ### should not be stat zone 1
# sort(unique(data$year))
### should I remove 2021? Yes, spatial distribution of effort different
data <- data[-which(data$year==2021),]

### adding covariate data and formula
data <- data[-which(is.na(data$bot_do)),]
# https://github.com/James-Thorson-NOAA/VAST/wiki/Specify-covariates-and-visualize-responses
# see also: https://github.com/James-Thorson-NOAA/VAST/issues/262
covariate_data <- data.frame(Lat=data$DECSLAT,
                             Lon=data$DECSLON,
                             Year=year(data$date),
                             bot_do=data$bot_do)
                             # bot_temp=data$TEMP_BOT)
# ### vast doesn't like NAs, but these data don't have to exactly match the input biomass data
# covariate_data <- covariate_data[-which(is.na(covariate_data$bot_do)),]
X1_formula <- ~ bot_do #+ bot_temp
X2_formula <- ~ 1

### alternate covariate
# covariate_data <- read.csv('covariate_data.csv')
### bottom DO
# if(length(which(is.na(covariate_data$Bot_DO)))>0){
#   covariate_data <- covariate_data[-which(is.na(covariate_data$Bot_DO)),]
# }
# X1_formula <- ~ Bot_DO #+ bot_temp
### bottom temp
# if(length(which(is.na(covariate_data$Bot_Temp)))>0){
#   covariate_data <- covariate_data[-which(is.na(covariate_data$Bot_Temp)),]
# }
# X1_formula <- ~ Bot_Temp
# X2_formula <- ~ 1

# par(mfrow=c(2,5))
# for(i in sort(unique(data$year))){
#   tmp <- data[data$year==i,]
#   plot(tmp$DECSLON,tmp$DECSLAT,
#        xlim=range(data$DECSLON,na.rm=T),ylim=range(data$DECSLAT,na.rm=T),
#        cex=(tmp$rg_wt+1))
#   mtext(i)
# }

### custom defined region
region <- read.csv('region.csv')

### strata limits
# https://github.com/James-Thorson-NOAA/VAST/issues/176
strata_lon <- aggregate(data$DECSLON,by=list(data$STAT_ZONE),range,na.rm=T)
strata_lon$x <- round(strata_lon$x)
strata_lon$x[2:3,2] <- -81.5
strata_lon$x[4:6,2] <- -82
strata_lat <- aggregate(data$DECSLAT,by=list(data$STAT_ZONE),range,na.rm=T)
strata_lat$x <- round(strata_lat$x)
strata_lat$x[7,1] <- 28.5

# plot(strata_lon$x[,1],strata_lat$x[,1],asp=1,
#      xlim=range(strata_lon$x),ylim=range(strata_lat$x))
# for(i in 1:nrow(strata_lat)){
#   rect(strata_lon$x[i,1],strata_lat$x[i,1],
#        strata_lon$x[i,2],strata_lat$x[i,2])
#   text(strata_lon$x[i,1],strata_lat$x[i,1],strata_lat$Group.1[i],adj=1,cex=2)
# }

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
settings = make_settings(n_x = 250, # set higher to avoid logKappa2 boundary issues; 250 worked with Delta-Lognormal; 400 for Delta-Generalized Gamma
                         Region = 'User', 
                         purpose = "index2", 
                         bias.correct = FALSE,
                         knot_method = 'grid',
                         use_anisotropy = T,
                         FieldConfig = matrix(c('IID','IID','IID','IID','IID','IID'), ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2"))),
                         RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=4, "Epsilon2"=4),
                         ObsModel = c(4,0) # try Delta-Gamma as alternative to Delta-Lognormal
                         # strata.limits = strata.limits # can be used to sum indices of abundance over strata
                         )
### change the Beta and Epsilon; fixed or random with AR1
# settings$RhoConfig
## L_epsilon2_z changed to run Delta-Generalized Gamma model
# settings$FieldConfig[2,2]=0

### Percent deviance explained
# https://github.com/James-Thorson-NOAA/VAST/wiki/Percent-deviance-explained

pthwy <- '~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/'
# Run model
fit <- fit_model(settings = settings, 
                Lat_i = data$DECSLAT, 
                Lon_i = data$DECSLON, 
                t_i = data$year, 
                b_i = as_units(data$rg_wt,'kg'), 
                a_i = as_units(data$AreaSwept_km2,'km^2'),
                v_i = as.numeric(data$VESSEL)-1,
                X1_formula = X1_formula,
                X2_formula = X2_formula,
                covariate_data = covariate_data,
                # X1config_cp = array(3,dim=c(1,1)), ### added vmod14; not working
                input_grid = region,
                run_model = T, # make FALSE to see upper (fit$tmb_list$Upper) and lower (fit$tmb_list$Lpper) starting bounds
                working_dir = pthwy
                ) 

setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
results <- plot(fit)

file_save <- paste0('~/Desktop/professional/projects/Postdoc_FL/data/grouper/',as.Date(Sys.time()),'_',run,'_results.RData') 
save.image(file_save)
