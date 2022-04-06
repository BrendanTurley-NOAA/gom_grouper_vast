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

load("~/Desktop/professional/projects/Postdoc_FL/data/grouper/2022-04-06_vmod3_results.RData")

### index
index_yr <- as.numeric(paste(vmod2$Index$Table$Time))
index <- vmod2$Index$Table$Estimate
index_st <- vmod2$Index$Table$Estimate/mean(vmod2$Index$Table$Estimate)
index_se <- vmod2$Index$Table$`Std. Error for Estimate`
### how to standardize the standard errors?

par(mfrow=c(2,1))
plot(index_yr,index,typ='b',
     ylim=range(c(index+index_se,index-index_se)))
polygon(c(index_yr,rev(index_yr)),
        c(index+index_se,rev(index-index_se)))

barplot(c(NA,diff(index)),names.arg = index_yr)
mtext('Change in index')

plot(index_yr,index_st,typ='b')


### center of gravity
cog <- (vmod2$Range$SD_mean_Z_ctm)

coord_utm <- SpatialPoints(cbind(cog[,,1,1], cog[,,2,1]), proj4string = CRS("+proj=utm +zone=17 +units=km"))
coord_ll <- spTransform(coord_utm, CRS("+proj=longlat"))

coord_utm_ucl <- SpatialPoints(cbind(cog[,,1,1]+cog[,,1,2], cog[,,2,1]+cog[,,2,2]), proj4string = CRS("+proj=utm +zone=17 +units=km"))
coord_ll_ucl <- spTransform(coord_utm_ucl, CRS("+proj=longlat"))

coord_utm_lcl <- SpatialPoints(cbind(cog[,,1,1]-cog[,,1,2], cog[,,2,1]-cog[,,2,2]), proj4string = CRS("+proj=utm +zone=17 +units=km"))
coord_ll_lcl <- spTransform(coord_utm_lcl, CRS("+proj=longlat"))

par(mfrow=c(2,1))
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


### area occupied
area_occ <- (vmod2$Range$SD_effective_area_ctl[1,,1,])

par(mfrow=c(2,1))
plot(index_yr,area_occ[,1],typ='b',
     ylim=range(c(area_occ[,1]+area_occ[,2],area_occ[,1]-area_occ[,2])))
polygon(c(index_yr,rev(index_yr)),
        c(area_occ[,1]+area_occ[,2],rev(area_occ[,1]-area_occ[,2])))

barplot(c(NA,diff(area_occ[,1])),names.arg = index_yr)
mtext('Change in area occupied')


### density plots

# par(mfrow=c(2,5))
# for(i in 1:10){
#   # quilt.plot(vmod2$map_list$PlotDF[,2:1],vmod2$Dens_xt[,i],asp=1,breaks=seq(-3,3,.1),col = tim.colors(60))
#   bubblePlot(vmod2$map_list$PlotDF[,2:1],vmod2$Dens_xt[,i],asp=1,size=.5,highlight=F,zlim=c(-3,3))
#   mtext(2009+i)
# }


r <- raster(ncols=300,nrow=300,
            xmn=min(vmod2$map_list$PlotDF[,2]),xmx=max(vmod2$map_list$PlotDF[,2]),
            ymn=min(vmod2$map_list$PlotDF[,1]),ymx=max(vmod2$map_list$PlotDF[,1]),
            res=.11)
lon <- seq(r@extent@xmin,r@extent@xmax,len=r@ncols)
lat <- seq(r@extent@ymin,r@extent@ymax,len=r@nrows)

col_pal <- colorRampPalette(c('gray20','purple','darkorange','gold'))
# col_pal <- colorRampPalette(c('purple4','dodgerblue4','seagreen3','khaki1'))
# col_pal <- colorRampPalette(c('honeydew2','darkseagreen3','forestgreen','darkslategrey'))


dens <- exp(vmod2$Dens_xt)
cutoff <- quantile((dens),.95)
dens[which(dens>=cutoff)] <- cutoff
breaks <- pretty(dens,n=20)

par(mfrow=c(2,5),mar=c(4,4,1.5,1))
for(i in 1:10){
  r1 <- raster::rasterize(vmod2$map_list$PlotDF[,2:1],r,dens[,i],fun=mean)
  # plot(r1,col=col_pal(length(breaks)-1),breaks=breaks,asp=1)
  r2 <- matrix(r1@data@values,r1@ncols,r1@nrows)
  image(lon,lat,r2[,ncol(r2):1],col=col_pal(length(breaks)-1),breaks=breaks,asp=1)
  points(coord_ll@coords[i,1],coord_ll@coords[i,2],col=3,pch=16)
  # if(i>1)(
  #   arrows(coord_ll@coords[i-1,1],coord_ll@coords[i-1,2],
  #          coord_ll@coords[i,1],coord_ll@coords[i,2],col='gray90',length=.05)
  # )
  mtext(2009+i)
}

r2 <- matrix(r1@data@values,r1@ncols,r1@nrows)
image(lon,lat,r2[,ncol(r2):1],col=col_pal(length(breaks)-1),breaks=breaks,asp=1)

