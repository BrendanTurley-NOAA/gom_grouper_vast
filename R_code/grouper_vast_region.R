library(lubridate)
library(ncdf4)
library(sp)
library(sf)


setwd("~/Desktop/professional/biblioteca/data")
# bathy <- nc_open('crm_vol3.nc')
# topo <- ncvar_get(bathy,'z')
# topo_lat <- ncvar_get(bathy, 'y')
# topo_lon <- ncvar_get(bathy, 'x')
# nc_close(bathy)
# bathy <- nc_open('etopo1.nc')
# topo <- ncvar_get(bathy, 'Band1')
# topo_lat <- ncvar_get(bathy, 'lat')
# topo_lon <- ncvar_get(bathy, 'lon')
# nc_close(bathy)
# rm(bathy)

### subset for region of interest
lonbox_e <- -80.6 ### Florida Bay
lonbox_w <- -88 ### mouth of Mississippi River
latbox_n <- 30.5 ### northern coast
latbox_s <- 24.3 ### southern edge of Ket West
# ind_lat <- which(topo_lat<=latbox_n & topo_lat>=latbox_s)
# ind_lon <- which(topo_lon<=lonbox_e & topo_lon>=lonbox_w)
# topo_lat <- topo_lat[ind_lat]
# topo_lon <- topo_lon[ind_lon]
# topo <- topo[ind_lon,ind_lat]
# topo[which(topo>=0)] <- 0
# setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
# png("test.png", height = 10, width = 10, units = 'in', res=300)
# image(topo_lon,topo_lat,topo,asp=1)
# dev.off()

### bathymetry
b_url <- 'https://tds.hycom.org/thredds/dodsC/datasets/GLBy0.08/expt_93.0/topo/depth_GLBy0.08_09m11.nc'
data <- nc_open(b_url)
topo_lat <- ncvar_get(data,'Latitude')
topo_lon <- ncvar_get(data,'Longitude')
ind_lat <- which(topo_lat<=latbox_n & topo_lat>=latbox_s)
ind_lon <- which(topo_lon<=(lonbox_e+360) & topo_lon>=(lonbox_w+360))
topo_lat <- topo_lat[ind_lat]
topo_lon <- -(360-topo_lon[ind_lon])
topo <- ncvar_get(data,'bathymetry',
                  start=c(ind_lon[1],ind_lat[1],1),
                  count=c(length(ind_lon),length(ind_lat),1))
topo <- -topo
topo[which(is.na(topo))] <- 0
nc_close(data)


setwd('~/Desktop/professional/projects/Postdoc_FL/data/grouper/')
data <- read.csv('grp_snp.csv')
data$year <- year(data$date)
data$rg_wt[which(is.na(data$rg_wt))] <- 0
data$rs_wt[which(is.na(data$rs_wt))] <- 0

# for(i in sort(unique(data$year))){
#   temp <- data[which(data$year==i),]
#   plot(temp$DECSLON,temp$DECSLAT,cex=(temp$rg_wt),asp=1)
#   mtext(i)
# }

### ----------------- user defined grid -----------------
#### Turn it into a spatial polygon object
region_extent <- data.frame(long=data$DECSLON, lat=data$DECSLAT)
region_extent <- region_extent[chull(region_extent),]
str(region_extent)
## Need to duplicate a point so that it is connected
region_extent <- rbind(region_extent, region_extent[1,])
## https://www.maths.lancs.ac.uk/~rowlings/Teaching/Sheffield2013/cheatsheet.html
poly <- Polygon(region_extent)
polys <- Polygons(list(poly), ID='all')
sps <- SpatialPolygons(list(polys))
## I think the F_AREA could be dropped here
sps <- SpatialPolygonsDataFrame(sps, data.frame(Id=factor('all'), F_AREA=1, row.names='all'))
proj4string(sps)<- CRS("+proj=longlat +datum=WGS84")
sps <- spTransform(sps, CRS("+proj=longlat +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
### Get UTM zone for conversion to UTM projection
## retrieves spatial bounding box from spatial data [,1] is
## longitude
lon <- sum(bbox(sps)[1,])/2
## convert decimal degrees to utm zone for average longitude, use
## for new CRS
utmzone <- floor((lon + 180)/6)+1
crs_LL <- CRS('+proj=longlat +ellps=WGS84 +no_defs')
sps@proj4string <- crs_LL

### Create the VAST extroplation grid for method 1 and 2
## Convert the final in polygon to UTM
crs_UTM <- CRS(paste0("+proj=utm +zone=",utmzone," +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
region_polygon <- spTransform(sps, crs_UTM)

### Construct the extroplation grid for VAST using sf package
## Size of grid **in meters** (since working in UTM). Controls
## the resolution of the grid.
cell_size <- 5000 # 10000 also works well
## This step is slow at high resolutions
region_grid <- st_make_grid(region_polygon, cellsize = cell_size, what = "centers")
## Convert region_grid to Spatial Points to SpatialPointsDataFrame
region_grid <- as(region_grid, "Spatial")
region_grid_sp <- as(region_grid, "SpatialPointsDataFrame")
## combine shapefile data (region_polygon) with Spatial Points
## (region_grid_spatial) & place in SpatialPointsDataFrame data
## (this provides you with your strata identifier (here called
## Id) in your data frame))
region_grid_sp@data <- over(region_grid, region_polygon)

## Convert back to lon/lat coordinates as that is what VAST uses
region_grid_LL <- as.data.frame(spTransform(region_grid_sp, crs_LL))
region_df <- with(region_grid_LL,
                  data.frame(Lon=coords.x1,
                             Lat=coords.x2, Id,
                             Area_km2=( (cell_size/1000)^2),
                             row=1:nrow(region_grid_LL)))
## Filter out the grid that does not overlap (outside extent)
region <- subset(region_df, !is.na(Id))
## This is the final file needed.
str(region)

# par(mfrow=c(2,2))
# with(region_extent, plot(long, lat, main='Extent in points in LL'))
# plot(region_polygon, main='Polygon in UTM', axes=TRUE)
# plot(region_grid, col=ifelse(is.na(region_df$Id), 'red', 'black'),
# axes=TRUE, main='Extrapolation area UTM')
# with(region, plot(Lon, Lat, main='Extrapolation region in LL', pch='.',asp=1))
# rm(crs_LL,crs_UTM,region_grid_LL,region_df,region_grid,region_grid_sp,region_polygon,sps,poly,polys,region_extent)

### find bottom depth from CRM bathymetry
lon_i <- unlist(lapply(region$Lon, function(x) which.min(abs(x-topo_lon))))
lat_i <- unlist(lapply(region$Lat, function(x) which.min(abs(x-topo_lat))))

bottom_z <- rep(NA,length(lon_i))
for(i in 1:length(lon_i)){
  bottom_z[i] <- topo[lon_i[i],lat_i[i]]
}
rm(topo,topo_lon,topo_lat,lon_i,lat_i)

ind <- which(bottom_z<=(-110) | bottom_z>=(-10))
region2 <- region[-ind,]
with(region2, plot(Lon, Lat, main='Extrapolation region in LL', pch='.',asp=1))
region <- region2
rm(region2)
### locations to remove
# rm_ll <- locator()

# setwd('~/Desktop/professional/projects/Postdoc_FL/figures/')
# png('bottom5.png',width=10,height=10, units = 'in', res=300)
plot(region$Lon,region$Lat,cex=-bottom_z[-ind]/100,asp=1)
points(region$Lon,region$Lat,pch='.',col=2)
# dev.off()


setwd('~/Desktop/professional/projects/Postdoc_FL/data/grouper/')
write.csv(region,'region.csv')
