rm(list=ls())
gc()

library(fields)
library(lubridate)
library(ncdf4)
library(raster)
library(rgdal)
library(sp)
library(sf)
library(tidyr)
library(VAST)

# # load data set
# # see `?load_example` for list of stocks with example data 
# # that are installed automatically with `FishStatsUtils`. 
# example = load_example( data_set="EBS_pollock" )
# 
# for(i in sort(unique(EBS_pollock_data$sampling_data$year))){
#   temp <- EBS_pollock_data$sampling_data[which(EBS_pollock_data$sampling_data$year==i),]
#   plot(temp$long,temp$lat,cex=log10(temp$catch),asp=1)
#   mtext(i)
# }
# # Make settings (turning off bias.correct to save time for example)
# settings = make_settings( n_x = 100, 
#                           Region = example$Region, 
#                           purpose = "index2", 
#                           bias.correct = FALSE )
# 
# # Run model
# fit = fit_model( settings = settings, 
#                  Lat_i = example$sampling_data[,'Lat'], 
#                  Lon_i = example$sampling_data[,'Lon'], 
#                  t_i = example$sampling_data[,'Year'], 
#                  b_i = example$sampling_data[,'Catch_KG'], 
#                  a_i = example$sampling_data[,'AreaSwept_km2'] )
# 
# # Plot results
# plot( fit )

# ?make_data
# ?make_mesh
# ?make_settings
# ?fit_model

example = load_example( data_set="multimodal_red_snapper" )


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
settings = make_settings(n_x = 100, 
                         Region = 'User', 
                         purpose = "index2", 
                         bias.correct = FALSE,
                         knot_method = 'grid',
                         ObsModel = c(1,0))
                         #strata.limits = data$STAT_ZONE) # can be used to sum indices of abundance over strata
# settings$FieldConfig[2,2]=0

setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
# Run model
fit = fit_model(settings = settings, 
                Lat_i = data$DECSLAT, 
                Lon_i = data$DECSLON, 
                t_i = data$year, 
                b_i = as_units(data$rg_wt,'kg'), 
                a_i = as_units(data$AreaSwept_km2,'km^2'),
                # v_i = as.numeric(data$VESSEL)-1,
                input_grid=region)

# Check bounds for the following parameters:
    # Param starting_value     Lower       MLE     Upper final_gradient
# 28 logkappa2     -0.1053605 -5.476512 -2.108411 -2.108411     -0.2129862

vmod4 <- plot(fit)

### custom plots: https://github.com/James-Thorson-NOAA/VAST/wiki/Plots-using-ggplot
mdl <- make_map_info(Region = settings$Region,
                     spatial_list = fit$spatial_list,
                     Extrapolation_List = fit$extrapolation_list)

D_gt <- fit$Report$D_gct[,1,] # drop the category
yrs <- sort(unique(data$year))
dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=yrs)
## tidy way of doing this, reshape2::melt() does
## it cleanly but is deprecated
D_gt <- D_gt %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Year", values_to='D')
D <- merge(D_gt, mdl$PlotDF, by.x='cell', by.y='x2i')


plot(D$Lon[D$Year==2014],D$Lat[D$Year==2014],asp=1,cex=log10(strip_units(D$D[D$Year==2014])+1))
plot(D$Lon[D$Year==2015],D$Lat[D$Year==2015],asp=1,cex=log10(strip_units(D$D[D$Year==2015])+1))
plot(D$Lon[D$Year==2018],D$Lat[D$Year==2018],asp=1,cex=log10(strip_units(D$D[D$Year==2018])+1))
plot(D$Lon[D$Year==2019],D$Lat[D$Year==2019],asp=1,cex=log10(strip_units(D$D[D$Year==2019])+1))


# save.image("~/Desktop/professional/projects/Postdoc_FL/data/grouper/20220330_vmod3_results.RData")
load("~/Desktop/professional/projects/Postdoc_FL/data/grouper/20220330_vmod2_results.RData")


index_yr <- as.numeric(paste(vmod2$Index$Table$Time))
index <- vmod2$Index$Table$Estimate
index_se <- vmod2$Index$Table$`Std. Error for Estimate`

plot(index_yr,index,typ='b',
     ylim=range(c(index+index_se,index-index_se)))
polygon(c(index_yr,rev(index_yr)),
        c(index+index_se,rev(index-index_se)))

barplot(c(NA,diff(index)),names.arg = index_yr)
mtext('Change in index')

cog <- (vmod2$Range$SD_mean_Z_ctm)

coord_utm <- SpatialPoints(cbind(cog[,,1,1], cog[,,2,1]), proj4string = CRS("+proj=utm +zone=17 +units=km"))
coord_ll <- spTransform(coord_utm, CRS("+proj=longlat"))

coord_utm_ucl <- SpatialPoints(cbind(cog[,,1,1]+cog[,,1,2], cog[,,2,1]+cog[,,2,2]), proj4string = CRS("+proj=utm +zone=17 +units=km"))
coord_ll_ucl <- spTransform(coord_utm_ucl, CRS("+proj=longlat"))

coord_utm_lcl <- SpatialPoints(cbind(cog[,,1,1]-cog[,,1,2], cog[,,2,1]-cog[,,2,2]), proj4string = CRS("+proj=utm +zone=17 +units=km"))
coord_ll_lcl <- spTransform(coord_utm_lcl, CRS("+proj=longlat"))


plot(index_yr,coord_ll@coords[,1],typ='b',
     ylim=range(c(coord_ll_ucl@coords[,1],coord_ll_lcl@coords[,1])))
polygon(c(index_yr,rev(index_yr)),
        c(coord_ll_ucl@coords[,1],rev(coord_ll_lcl@coords[,1])))

plot(index_yr,coord_ll@coords[,2],typ='b',
     ylim=range(c(coord_ll_ucl@coords[,2],coord_ll_lcl@coords[,2])))
polygon(c(index_yr,rev(index_yr)),
        c(coord_ll_ucl@coords[,2],rev(coord_ll_lcl@coords[,2])))

plot(coord_ll@coords,typ='b')
text(coord_ll@coords,labels=index_yr)


area_occ <- (vmod2$Range$SD_effective_area_ctl[1,,1,])

plot(index_yr,area_occ[,1],typ='b',
     ylim=range(c(area_occ[,1]+area_occ[,2],area_occ[,1]-area_occ[,2])))
polygon(c(index_yr,rev(index_yr)),
        c(area_occ[,1]+area_occ[,2],rev(area_occ[,1]-area_occ[,2])))

barplot(c(NA,diff(area_occ[,1])),names.arg = index_yr)
mtext('Change in area occupied')



par(mfrow=c(2,5))
for(i in 1:10){
  # quilt.plot(vmod2$map_list$PlotDF[,2:1],vmod2$Dens_xt[,i],asp=1,breaks=seq(-3,3,.1),col = tim.colors(60))
  bubblePlot(vmod2$map_list$PlotDF[,2:1],vmod2$Dens_xt[,i],asp=1,size=.5,highlight=F,zlim=c(-3,3))
  mtext(2009+i)
}



r <- raster(ncols=300,nrow=300,
            xmn=min(vmod2$map_list$PlotDF[,2]),xmx=max(vmod2$map_list$PlotDF[,2]),
            ymn=min(vmod2$map_list$PlotDF[,1]),ymx=max(vmod2$map_list$PlotDF[,1]),
            res=.11)
lon <- seq(r@extent@xmin,r@extent@xmax,len=r@ncols)
lat <- seq(r@extent@ymin,r@extent@ymax,len=r@nrows)

col_pal <- colorRampPalette(c('gray20','purple','darkorange','gold'))
col_pal <- colorRampPalette(c('purple4','dodgerblue4','seagreen3','khaki1'))
col_pal <- colorRampPalette(c('honeydew2','darkseagreen3','forestgreen','darkslategrey'))


dens <- exp(vmod2$Dens_xt)
cutoff <- quantile((dens),.95)
dens[which(dens>=cutoff)] <- cutoff
breaks <- pretty(dens,n=20)
par(mfrow=c(2,5))
for(i in 1:10){
  r1 <- raster::rasterize(vmod2$map_list$PlotDF[,2:1],r,dens[,i],fun=mean)
  plot(r1,col=col_pal(length(breaks)-1),breaks=breaks,asp=1)
  mtext(2009+i)
}

r2 <- matrix(r1@data@values,r1@ncols,r1@nrows)
image(lon,lat,r2[,ncol(r2):1],col=col_pal(length(breaks)-1),breaks=breaks,asp=1)
