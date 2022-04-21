rm(list=ls())
gc()

library(fields)
library(lubridate)
library(ncdf4)
library(raster)
library(rgdal)
library(scales)
library(sp)
library(sf)
library(tidyr)
library(VAST)

### load map
# setwd("C:/Users/brendan.turley/Desktop/FL_habs/ne_10m_admin_0_countries")
# setwd("~/Desktop/professional/biblioteca/data/shapefiles/ne_10m_admin_0_countries")
# world <- readOGR('ne_10m_admin_0_countries.shp')
setwd("~/Desktop/professional/biblioteca/data/shapefiles/gshhg-shp-2.3.7/GSHHS_shp/h/")
world <- readOGR('GSHHS_h_L1.shp')
world <- crop(world, extent(-86.5, -80, 24.5, 31))

### load VAST results
load("~/Desktop/professional/projects/Postdoc_FL/data/grouper/2022-04-21_vmod8_results.RData")
model <- 'vmod8'

### index
index_yr <- as.numeric(paste(results$Index$Table$Time))
index <- results$Index$Table$Estimate
index_st <- results$Index$Table$Estimate/mean(results$Index$Table$Estimate)
index_se <- results$Index$Table$`Std. Error for Estimate`
### how to standardize the standard errors?

setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
png("rg_index_vast.png", height = 8, width = 8, units = 'in', res=300)
par(mfrow=c(2,1),mar=c(5,5,1,1),oma=c(1,0,0,0))
plot(index_yr,index,
     typ='o',las=1,
     xlab='Year',ylab='',
     col=4,pch=16,xaxt='n',
     ylim=range(c(index+index_se,index-index_se)))
axis(1,2010:20119)
mtext('Red grouper index (kg)',2,4)
polygon(c(index_yr,rev(index_yr)),
        c(index+index_se,rev(index-index_se)),
        col=alpha(4,.3))

barplot(c(NA,diff(index)),names.arg = index_yr,las=1)
abline(h=0,lty=5)
mtext(expression(paste(Delta,' index')),2,4)
mtext(paste('Run:',model),1,adj=1,outer=T,col='red',line=-1)
dev.off()

par(mfrow=c(1,1))
plot(index_yr,index_st,typ='b')
points(index_yr,index_st,typ='b',col=2,pch=16)

### center of gravity
cog <- (results$Range$SD_mean_Z_ctm)

coord_utm <- SpatialPoints(cbind(cog[,,1,1], cog[,,2,1]), proj4string = CRS("+proj=utm +zone=17 +units=km"))
coord_ll <- spTransform(coord_utm, CRS("+proj=longlat"))

coord_utm_ucl <- SpatialPoints(cbind(cog[,,1,1]+cog[,,1,2], cog[,,2,1]+cog[,,2,2]), proj4string = CRS("+proj=utm +zone=17 +units=km"))
coord_ll_ucl <- spTransform(coord_utm_ucl, CRS("+proj=longlat"))

coord_utm_lcl <- SpatialPoints(cbind(cog[,,1,1]-cog[,,1,2], cog[,,2,1]-cog[,,2,2]), proj4string = CRS("+proj=utm +zone=17 +units=km"))
coord_ll_lcl <- spTransform(coord_utm_lcl, CRS("+proj=longlat"))


setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
png("rg_cog_vast.png", height = 8, width = 8, units = 'in', res=300)
par(mfrow=c(2,1),mar=c(5,5,1,1),oma=c(1,0,0,0))
plot(index_yr,coord_ll@coords[,1],
     typ='o',col=2,pch=16,las=1,xaxt='n',
     xlab='Year',ylab='Longitude',
     ylim=range(c(coord_ll_ucl@coords[,1],coord_ll_lcl@coords[,1])))
axis(1,2010:2019)
polygon(c(index_yr,rev(index_yr)),
        col=alpha(2,.3),
        c(coord_ll_ucl@coords[,1],rev(coord_ll_lcl@coords[,1])))
mtext('Center of gravity')

plot(index_yr,coord_ll@coords[,2],
     typ='o',col=4,pch=16,las=1,xaxt='n',
     xlab='Year',ylab='Latitude',
     ylim=range(c(coord_ll_ucl@coords[,2],coord_ll_lcl@coords[,2])))
axis(1,2010:2019)
polygon(c(index_yr,rev(index_yr)),
        col=alpha(4,.3),
        c(coord_ll_ucl@coords[,2],rev(coord_ll_lcl@coords[,2])))
mtext(paste('Run:',model),1,adj=1,outer=T,col='red',line=-1)
dev.off()

par(mfrow=c(1,1))
plot(coord_ll@coords,typ='b')
text(coord_ll@coords,labels=index_yr)


par(mfrow=c(2,1))
plot(index_yr,cog[,,1,1]-min(cog[,,1,1]),
# plot(index_yr,cog[,,1,1]-cog[,1,1,1],
     las=1,typ='o',col=2,pch=16,
     xlab='Year',ylab='Shift Eastward (km)')
plot(index_yr,cog[,,2,1]-min(cog[,,2,1]),
# plot(index_yr,cog[,,2,1]-cog[,1,2,1],
     las=1,typ='o',col=4,pch=16,
     xlab='Year',ylab='Shift Northward (km)')


### area occupied
area_occ <- (results$Range$SD_effective_area_ctl[1,,1,])

setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
png("rg_area_vast.png", height = 8, width = 8, units = 'in', res=300)
par(mfrow=c(2,1),mar=c(5,5,1,1),oma=c(1,0,0,0))
plot(index_yr,area_occ[,1],
     typ='o',col=4,pch=16,las=1,
     xlab='Year',ylab='',xaxt='n',
     ylim=range(c(area_occ[,1]+area_occ[,2],area_occ[,1]-area_occ[,2])))
axis(1,2010:2019)
mtext(expression(paste('Area occupied (km'^2,')')),2,3.5)
polygon(c(index_yr,rev(index_yr)),
        col=alpha(4,.3),
        c(area_occ[,1]+area_occ[,2],rev(area_occ[,1]-area_occ[,2])))

barplot(c(NA,diff(area_occ[,1])),names.arg = index_yr,las=1)
mtext(expression(paste(Delta,'area occupied (km'^2,')')),2,3.5)
mtext(paste('Run:',model),1,adj=1,outer=T,col='red',line=-1)
dev.off()


### density plots

plot(index_yr,index/area_occ[,1],
     typ='o',pch=16,col=2,
     xlab='Year',ylab='Density (kg/km^2)')

barplot(c(NA,diff(index/area_occ[,1])),names.arg = index_yr)
mtext('Change in mean density')

# par(mfrow=c(2,5))
# for(i in 1:10){
#   # quilt.plot(results$map_list$PlotDF[,2:1],results$Dens_xt[,i],asp=1,breaks=seq(-3,3,.1),col = tim.colors(60))
#   bubblePlot(results$map_list$PlotDF[,2:1],results$Dens_xt[,i],asp=1,size=.5,highlight=F,zlim=c(-3,3))
#   mtext(2009+i)
# }


r <- raster(ncols=300,nrow=300,
            xmn=min(results$map_list$PlotDF[,2]),xmx=max(results$map_list$PlotDF[,2]),
            ymn=min(results$map_list$PlotDF[,1]),ymx=max(results$map_list$PlotDF[,1]),
            res=.11)
lon <- seq(r@extent@xmin,r@extent@xmax,len=r@ncols)
lat <- seq(r@extent@ymin,r@extent@ymax,len=r@nrows)


# col_pal <- colorRampPalette(c('gray20','dodgerblue4','indianred3','gold1'))
col_pal <- colorRampPalette(rev(c('khaki1','cadetblue2','dodgerblue3','slateblue4')))
lm_neg <- colorRampPalette(c('dodgerblue4','deepskyblue3','lightskyblue1','gray95'))
lm_pos <- colorRampPalette(c('gray95','rosybrown1','tomato2','red4'))


dens <- exp(results$Dens_xt)
# dens <- (results$Dens_xt) # log space
cutoff <- quantile((dens),.99)
dens[which(dens>=cutoff)] <- cutoff
breaks <- pretty(dens,n=20)


setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
png("rg_den_vast.png", height = 5, width = 10, units = 'in', res=300)
par(mfrow=c(2,5),mar=c(4,4,1.5,1),oma=c(0,0,0,5))
for(i in 1:10){
  r1 <- raster::rasterize(results$map_list$PlotDF[,2:1],r,dens[,i],fun=mean)
  # plot(r1,col=col_pal(length(breaks)-1),breaks=breaks,asp=1)
  r2 <- matrix(r1@data@values,r1@ncols,r1@nrows)
  image(lon,lat,r2[,ncol(r2):1],
        col=col_pal(length(breaks)-1),breaks=breaks,asp=1,
        xlab='',ylab='')
  points(coord_ll@coords[i,1],coord_ll@coords[i,2],col=2,pch=16)
  plot(world,add=T,col='gray80')
  # if(i>1)(
  #   arrows(coord_ll@coords[i-1,1],coord_ll@coords[i-1,2],
  #          coord_ll@coords[i,1],coord_ll@coords[i,2],col='gray90',length=.05)
  # )
  mtext(2009+i)
}
par(oma=c(1,0,0,3))
imagePlot(legend.only=TRUE, zlim=range(dens),col=col_pal(20),
          legend.lab=expression(paste('Red grouper density (kg km'^-2,')')),
          legend.line=2.5,legend.cex = .7,legend.width = 1.75)
mtext(paste('Run:',model),1,adj=1,outer=T,col='red',line=-1)
dev.off()



breaks <- seq(-5.5,5.5,.5)
# breaks <- seq(-2,2,.1) # if log space
cols <- c(lm_neg(length(which(breaks<0))),lm_pos(length(which(breaks>0))))

setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
png("rg_diffvast.png", height = 5, width = 10, units = 'in', res=300)
par(mfrow=c(2,5),mar=c(4,4,1.5,1),oma=c(0,0,0,5))
for(i in 1:10){
  r1 <- raster::rasterize(results$map_list$PlotDF[,2:1],r,dens[,i],fun=mean)
  # plot(r1,col=col_pal(length(breaks)-1),breaks=breaks,asp=1)
  r2 <- matrix(r1@data@values,r1@ncols,r1@nrows)
  r2 <- r2[,ncol(r2):1]
  if(i>=2){
    r_diff <- r2-r3
    if(max(r_diff,na.rm=T)>max(breaks)){
      r_diff[which(r_diff>max(breaks))] <- max(breaks)
    }
    if(min(r_diff,na.rm=T)<min(breaks)){
      r_diff[which(r_diff<min(breaks))] <- min(breaks)
    }
    image(lon,lat,r_diff,
          asp=1,col=cols,breaks=breaks,xlab='',ylab='')
    plot(world,add=T,col='gray80')
    mtext(paste(2009+i,'-',2009+i-1))
    print(quantile(r_diff,na.rm=T))
  }
  r3 <- r2
}
par(oma=c(1,0,0,3))
imagePlot(legend.only=TRUE, zlim=range(breaks),col=cols,
          legend.lab=expression(paste(Delta,'Red grouper density (kg km'^-2,')')),
          legend.line=2.5,legend.cex = .7,legend.width = 1.75)
mtext(paste('Run:',model),1,adj=1,outer=T,col='red',line=-1)
dev.off()


##### log space
dens <- (results$Dens_xt) # log space
cutoff <- quantile((dens),.99)
dens[which(dens>=cutoff)] <- cutoff
breaks <- pretty(dens,n=20)

setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
png("rg_den_vast_log.png", height = 5, width = 10, units = 'in', res=300)
par(mfrow=c(2,5),mar=c(4,4,1.5,1),oma=c(0,0,0,5))
for(i in 1:10){
  r1 <- raster::rasterize(results$map_list$PlotDF[,2:1],r,dens[,i],fun=mean)
  # plot(r1,col=col_pal(length(breaks)-1),breaks=breaks,asp=1)
  r2 <- matrix(r1@data@values,r1@ncols,r1@nrows)
  image(lon,lat,r2[,ncol(r2):1],
        col=col_pal(length(breaks)-1),breaks=breaks,asp=1,
        xlab='',ylab='')
  points(coord_ll@coords[i,1],coord_ll@coords[i,2],col=2,pch=16)
  plot(world,add=T,col='gray80')
  # if(i>1)(
  #   arrows(coord_ll@coords[i-1,1],coord_ll@coords[i-1,2],
  #          coord_ll@coords[i,1],coord_ll@coords[i,2],col='gray90',length=.05)
  # )
  mtext(2009+i)
}
par(oma=c(1,0,0,3))
imagePlot(legend.only=TRUE, zlim=range(dens),col=col_pal(20),
          legend.lab=expression(paste('Red grouper density (kg km'^-2,')')),
          legend.line=2.5,legend.cex = .7,legend.width = 1.75)
mtext(paste('Run:',model),1,adj=1,outer=T,col='red',line=-1)
dev.off()


breaks <- seq(-2,2,.1) # if log space
cols <- c(lm_neg(length(which(breaks<0))),lm_pos(length(which(breaks>0))))

setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
png("rg_diffvast_log.png", height = 5, width = 10, units = 'in', res=300)
par(mfrow=c(2,5),mar=c(4,4,1.5,1),oma=c(0,0,0,5))
for(i in 1:10){
  r1 <- raster::rasterize(results$map_list$PlotDF[,2:1],r,dens[,i],fun=mean)
  # plot(r1,col=col_pal(length(breaks)-1),breaks=breaks,asp=1)
  r2 <- matrix(r1@data@values,r1@ncols,r1@nrows)
  r2 <- r2[,ncol(r2):1]
  if(i>=2){
    r_diff <- r2-r3
    if(max(r_diff,na.rm=T)>max(breaks)){
      r_diff[which(r_diff>max(breaks))] <- max(breaks)
    }
    if(min(r_diff,na.rm=T)<min(breaks)){
      r_diff[which(r_diff<min(breaks))] <- min(breaks)
    }
    image(lon,lat,r_diff,
          asp=1,col=cols,breaks=breaks,xlab='',ylab='')
    plot(world,add=T,col='gray80')
    mtext(paste(2009+i,'-',2009+i-1))
    print(quantile(r_diff,na.rm=T))
  }
  r3 <- r2
}
par(oma=c(1,0,0,3))
imagePlot(legend.only=TRUE, zlim=range(breaks),col=cols,
          legend.lab=expression(paste(Delta,'Red grouper density (kg km'^-2,')')),
          legend.line=2.5,legend.cex = .7,legend.width = 1.75)
mtext(paste('Run:',model),1,adj=1,outer=T,col='red',line=-1)
dev.off()
