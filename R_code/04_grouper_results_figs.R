rm(list=ls())
gc()

library(fields)
library(lubridate)
library(MASS)
library(ncdf4)
library(raster)
library(rgdal)
library(scales)
library(sp)
library(sf)
library(tidyr)
library(VAST)

### takes a linear index of a matrix and returns row (r) and column (c) indices
### requires some linear index (ind)
### e.g., ind <- which(a==0) or ind <- which(is.na(a)), if a is a matrix m x n
ind2sub <- function(ind,a){
  m <- nrow(a)
  r <- ((ind-1) %% m) + 1
  c <- floor((ind-1) / m) + 1
  return(cbind(r,c))
}

### load bathymetry
setwd("~/Desktop/professional/biblioteca/data")
bathy <- nc_open('etopo1.nc')
topo <- ncvar_get(bathy, 'Band1')
topo_lat <- ncvar_get(bathy, 'lat')
topo_lon <- ncvar_get(bathy, 'lon')
nc_close(bathy)

### load map
# setwd("C:/Users/brendan.turley/Desktop/FL_habs/ne_10m_admin_0_countries")
# setwd("~/Desktop/professional/biblioteca/data/shapefiles/ne_10m_admin_0_countries")
# world <- readOGR('ne_10m_admin_0_countries.shp')
setwd("~/Desktop/professional/biblioteca/data/shapefiles/gshhg-shp-2.3.7/GSHHS_shp/h/")
world <- readOGR('GSHHS_h_L1.shp')
world <- crop(world, extent(-86.5, -80, 24.5, 31))

### load VAST results
load("~/Desktop/professional/projects/Postdoc_FL/data/grouper/2022-05-17_vmod14_results.RData")

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
mtext(paste('Run:',run),1,adj=1,outer=T,col='red',line=-1)
dev.off()


### ----------------- SEDAR 61 indices -----------------
add_years <- function(df){
  merge(data.frame(Year=seq(min(df$Year),max(df$Year))),
        df,by='Year',all=T)
}

setwd('~/Desktop/professional/projects/Postdoc_FL/data/grouper')
bbl <- read.table('rg_bbl_index.txt',header=T)
trawl <- read.table('rg_trawl_index.txt',header=T)
hb <- read.table('rg_headboat_index.txt',header=T)
video <- read.table('rg_video_index.txt',header=T)

bbl <- add_years(bbl)
trawl <- add_years(trawl)
hb <- add_years(hb)
video <- add_years(video)
vast_indx <- data.frame('Year'=index_yr,'Index'=index_st)
vast_indx <- add_years(vast_indx)
indx <- merge(bbl,trawl,by=c('Year'),all=T)
indx <- merge(indx,hb,by=c('Year'),all=T)
indx <- merge(indx,video,by=c('Year'),all=T)
indx <- merge(indx,vast_indx,by=c('Year'),all=T)
indx_en <- apply(indx[,c(5,12,20,27)],1,median,na.rm=T)
indx_en2 <- apply(indx[,c(5,12,20,27,29)],1,median,na.rm=T)

setwd('~/Desktop/professional/projects/Postdoc_FL/figures/grouper/vast/')
png("rg_indices_vast.png", height = 5, width = 8, units = 'in', res=300)
par(mfrow=c(1,1),mar=c(5,5,1,1),oma=c(1,0,0,0))
plot(bbl$Year,bbl$Scaled_Index,
     typ='o',lwd=2,pch=16,lty=1,
     xlim=c(2000,2019),ylim=c(0,2.5),
     xaxt='n',xlab='',ylab='Scaled index')
axis(1,2000:2020)
points(trawl$Year,trawl$Scaled_Index,
       typ='o',lwd=2,pch=17,lty=2,col='gray20')
points(hb$Year,hb$relative_index,
       typ='o',lwd=2,pch=15,lty=3,col='gray40')
points(video$Year,video$std_indx,
       typ='o',lwd=2,pch=18,lty=4,col='gray60')
points(vast_indx$Year,vast_indx$Index,
       typ='o',lwd=2,pch=20,lty=1,col=4)
legend('topright',c('LL','Trawl','HB','Video','Ensemble','VAST'),
       lty=c(1,2,3,4,1,1),pch=c(16,17,15,18,NA,20),
       col=c(1,'gray20','gray40','gray60',2,4),bty='n',lwd=2)
points(indx$Year,indx_en,typ='l',col=2,lwd=2,pch=20)
# points(indx$Year,indx_en2,typ='l',col=3,lwd=2,pch=20)
mtext('RG indices - SEDAR 61')
abline(v=seq(2000,2020,5),lty=5,col='gray60')
dev.off()

par(mfrow=c(1,1))
plot(index_yr,index_st,typ='o',col=2,pch=16)


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
mtext(paste('Run:',run),1,adj=1,outer=T,col='red',line=-1)
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


### center of depth

lon_i <- unlist(lapply(results$map_list$PlotDF[,2], function(x) which.min(abs(x-topo_lon))))
lat_i <- unlist(lapply(results$map_list$PlotDF[,1], function(x) which.min(abs(x-topo_lat))))

bottom_z <- rep(NA,length(lon_i))
for(i in 1:length(lon_i)){
  bottom_z[i] <- topo[lon_i[i],lat_i[i]]
}

# plot(results$map_list$PlotDF[,2],results$map_list$PlotDF[,1],
#      cex=-bottom_z/100,asp=1)

dens <- exp(results$Dens_xt)

cod <- matrix(NA,10,3)
par(mfrow=c(1,2))
for(i in 1:10){
  # center of depth
  cod[i,1] <- sum(dens[,i]*bottom_z)/sum(dens[,i])
  # range
  # kde1 <- kde2d(dens[,i],bottom_z,n=50)
  # image(kde1)
  # ci_95 <- quantile(kde1$z,c(.95))
  # kde3 <- kde1
  # kde3$z[which(kde3$z<ci_95)] <- NA
  # image(kde3)
  # mtext(2009+i)
  # ind <- which(!is.na(kde3$z))
  # rc <- ind2sub(ind,kde3$z)
  # cod[i,2:3] <- range(kde3$y[rc[,2]])
 
  dens_sort_sum <- cumsum(dens[,i][order(bottom_z)])/sum(dens[,i])
  z_sort <- bottom_z[order(bottom_z)]
  
  bind <- which(dens_sort_sum<=.05)
  bot <- z_sort[max(bind)]
  
  tind <- which(dens_sort_sum>=.95)
  top <- z_sort[min(tind)]
  
  cod[i,2:3] <- c(bot,top)
}

par(mfrow=c(1,1))
plot(2010:2019,cod[,1],ylim=range(cod[,2:3]))
polygon(c(2010:2019,rev(2010:2019)),
        c(cod[,2],rev(cod[,3])))


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
mtext(paste('Run:',run),1,adj=1,outer=T,col='red',line=-1)
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
mtext(paste('Run:',run),1,adj=1,outer=T,col='red',line=-1)
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
mtext(paste('Run:',run),1,adj=1,outer=T,col='red',line=-1)
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
mtext(paste('Run:',run),1,adj=1,outer=T,col='red',line=-1)
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
mtext(paste('Run:',run),1,adj=1,outer=T,col='red',line=-1)
dev.off()
