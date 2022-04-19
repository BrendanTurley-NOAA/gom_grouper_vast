# from David Hansiko
# Under SEAMAP, experimental trawling on the West Florida Shelf was initiated in 2008 by the state of FL.   In 2010, the trawl surveys were expanded in that area to include NMFS and FL sampling effort.

### SEDAR 61 - bottom trawl index
# https://sedarweb.org/docs/wpapers/S61_WP_12_SEAMAP_groundfish.pdf

### other things to consider
# 1 sample and select weight are separate and mutual exclusive; if one happens then the other does not? (used the is_sample column)
# 2 there is vessel speed and time towed to get distance; need to find net dimensions
# 3 there is good/bad trawl codes; use them!
# 4 BGSREC has weights; use this instead of GLFREC
# 5 INVREC has bottom type
# 6 there is also gear type in STAREC; use them!

### to estimate area swept by gear
## https://doi.org/10.1093/icesjms/fsv099
## https://doi.org/10.1093/icesjms/fsu235
## https://doi.org/10.1126/science.aat6713
## https://sustainablefisheries-uw.org/global-footprint-of-fishing/

library(lubridate)
library(MASS)
library(nlstools)
library(raster)
library(rgdal)
library(ROCR)
library(scales)

### standard error
std.er <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))
}

### subset for region of interest
lonbox_e <- -80.6 ### Florida Bay
lonbox_w <- -86 ### mouth of Mississippi River
latbox_n <- 30.5 ### northern coast
latbox_s <- 24.3 ### southern edge of Ket West

# setwd("~/Desktop/professional/biblioteca/data/shapefiles/gshhg-shp-2.3.7/GSHHS_shp/h/")
# world <- readOGR('GSHHS_h_L1.shp')
# world <- crop(world, extent(-87, -79.5, 24.3, 31))

# setwd('~/Desktop/professional/projects/Postdoc_FL/data/ctd/seamap/public_seamap_csvs/')
setwd('~/Desktop/professional/projects/Postdoc_FL/data/ctd/seamap/public_seamap_csvs/public_seamap_csvs_2022/')

catch <- read.csv('BGSREC.csv')
cruises <- read.csv('CRUISES.csv')
ctd <- read.csv('CTDREC.csv')
ctd_data <- read.csv('CTDCASTREC.csv')
# station <- read.csv('STAREC_3.csv') ### I had to manually modify the file by deleting the comment column because the comments had end of file (EOF) character in the string
station <- read.csv('STAREC.csv')
vessels <- read.csv('VESSELS.csv')
envrec <- read.csv('ENVREC.csv')
# lth_freq <- read.csv('GLFREC.csv')

### removed 'BAD' trawls
ind <- which(station$HAULVALUE=='B')
station <- station[-ind,]
### shrimp trawls used
ind <- grep('ST',station$GEARS)
station <- station[ind,]

### subset for location of interest
ind <- which(station$DECSLON<=lonbox_e & station$DECSLON>=lonbox_w & station$DECSLAT<=latbox_n & station$DECSLAT>=latbox_s)
station <- station[ind,]
### keep only columns of interest
ind <- c(1:3,5:7,16:17,25:33,36:38,42:43,45:47)
names(station)[ind]
station_redux <- station[,ind]
### time trawled
start_time <- mdy_hm(station_redux$START_DATE)
end_time <- mdy_hm(station_redux$END_DATE)
station_redux$trawl_hrs <- as.numeric(end_time-start_time)/3600
trawl_km <- station_redux$VESSEL_SPD*station_redux$trawl_hrs*1.852 # km
trawl_wdth <- 50*.6*2.54/1000 # 60% footrope spread in inches converted to cm then km
station_redux$AreaSwept_km2 <- trawl_km*trawl_wdth # km^2 area swept
### subset CTD data for corresponding station data already subsetted
ctd_data <- ctd_data[is.element(ctd_data$STATIONID,unique(station$STATIONID)),]


### ----------------- pull out bottom data from CTD -----------------
ids <- sort(unique(ctd_data$CTDID))
bot_do <- bot_temp <- rep(NA,length(ids))
station_id <- cruise_id <- rep(NA,length(ids))
n <- 1
for(i in ids){
  temp <- ctd_data[which(ctd_data$CTDID==i),]
  bot_do[n] <- temp$OXY_MG[which.max(temp$DEPTH)]
  bot_temp[n] <- temp$TEMP[which.max(temp$DEPTH)]
  station_id[n] <- unique(temp$STATIONID)
  cruise_id[n] <- unique(temp$CRUISEID)
  n <- n + 1
}
bot_do[which(bot_do>14)] <- NA
bot_temp[which(bot_temp>40)] <- NA


### ----------------- subset species -----------------
spcde <- 170021211 #red grouper
# spcde <- 170022104 # gag; not many instances
ind <- which(catch$BIO_BGS==spcde)
grouper <- catch[ind,]

spcde <- 170151107 # red snapper
ind <- which(catch$BIO_BGS==spcde)
snapper <- catch[ind,]

spcde1 <- 613000000 # sponge/porifera
spcde2 <- 613000010
spcde3 <- 613000020
spcde4 <- 613000030
spcde5 <- 613000040
spcde6 <- 613000050
ind <- which(catch$BIO_BGS==spcde1 | 
               catch$BIO_BGS==spcde2 | 
               catch$BIO_BGS==spcde3 | 
               catch$BIO_BGS==spcde4 | 
               catch$BIO_BGS==spcde5 |
               catch$BIO_BGS==spcde6)
sponge <- catch[ind,]
sponge_agg <- aggregate(sponge$SELECT_BGS,by=list(sponge$CRUISEID,sponge$STATIONID),sum,na.rm=T)
names(sponge_agg) <- c('CRUISEID','STATIONID','sp_wt')

### rename count and count extrapolation and total weight
grouper <- grouper[,c(2:3,11:12,14)]
names(grouper)[c(3:5)] <- c('rg_cnt','rg_cntexp','rg_wt')
snapper <- snapper[,c(2:3,11:12,14)]
names(snapper)[c(3:5)] <- c('rs_cnt','rs_cntexp','rs_wt')

### merge species
grp_snp <- merge(grouper,snapper,by=c('STATIONID','CRUISEID'),all=T)
grp_snp <- merge(grp_snp,sponge_agg,by=c('STATIONID','CRUISEID'),all=T)

### subset for only FL and NOAA cruises and 2010 onward and trawl survey data
grp_cruises <- cruises[is.element(cruises$CRUISEID,unique(grp_snp$CRUISEID)),]
FL_cruises <- grp_cruises[which(grp_cruises$SOURCE=='FL' | grp_cruises$SOURCE=='US'),]
FL_cruises <- FL_cruises[which(FL_cruises$YR>=2010 & grepl('summer',FL_cruises$TITLE,ignore.case = T) &
                                 grepl('groundfish',FL_cruises$TITLE,ignore.case = T)),]
### only summer due to fall converage issues following Pollock et al. 2018; SEDAR61-WP-12
### 2008 and 2009 excluded due to poor coverage
# FL_cruises <- FL_cruises[which(FL_cruises$YR>=2008 & grepl('groundfish',FL_cruises$TITLE,ignore.case = T)),]


### ----------------- merge data for logistic regression -----------------
grp_do <- data.frame(STATIONID=station_id,
                     CRUISEID=cruise_id,
                     bot_do=bot_do)

grp <- merge(grp_do,grp_snp,by=c('STATIONID','CRUISEID'),all=T)
new <- merge(grp,station_redux,by=c('STATIONID','CRUISEID'),all.x=T)
new <- new[is.element(new$CRUISEID,FL_cruises$CRUISEID),]
new$date <- mdy_hm(paste(new$MO_DAY_YR,
                         paste(substr(sprintf('%04d',as.numeric(new$TIME_MIL)),1,2),
                               substr(sprintf('%04d',as.numeric(new$TIME_MIL)),3,4),sep=':')))
# new$date <- mdy_hm(paste(new$MO_DAY_YR,sprintf('%04d',as.numeric(new$TIME_MIL))))
ind <- which(is.na(new$MO_DAY_YR))
new <- new[-ind,]
### make presence/absence
new$rg_pr_ab <- ifelse(!is.na(new$rg_cnt),1,0)
new$rs_pr_ab <- ifelse(!is.na(new$rs_cnt),1,0)
### cpue
new$rg_cpue <- new$rg_cntexp/new$trawl_hrs
new$rs_cpue <- new$rs_cntexp/new$trawl_hrs
### water depth
new <- merge(new,envrec,by=c('STATIONID','CRUISEID'),all.x=T)
### which depth
par(mfrow=c(1,2))
plot(new$DECSLON,new$DECSLAT,cex=new$DEPTH_EMAX/100,asp=1)
plot(new$DECSLON,new$DECSLAT,cex=new$DEPTH_EWTR/100,asp=1)
par(mfrow=c(1,1))
plot(new$DEPTH_EMAX,new$DEPTH_EWTR)
abline(0,1,col=2)
### EMAX seems more reasonable
new$DEPTH_EMAX <- (-new$DEPTH_EMAX)


ind <- c(1:10,18,27:29,33:39,52)
ind[10:22] <- ind[10:22]+1 # add for sponge data
ind <- c(ind,10) # add for sponge data
names(new)[ind]
# [1] "STATIONID"     "CRUISEID"      "bot_do"        "rg_cnt"        "rg_cntexp"    
# [6] "rg_wt"         "rs_cnt"        "rs_cntexp"     "rs_wt"         "VESSEL.x"     
# [11] "TEMP_BOT"      "STAT_ZONE"     "DECSLAT"       "DECSLON"       "trawl_hrs"    
# [16] "AreaSwept_km2" "date"          "rg_pr_ab"      "rs_pr_ab"      "rg_cpue"      
# [21] "rs_cpue"       "DEPTH_EMAX" 
grp_snp_cov <- new[,ind]
names(grp_snp_cov)
ind <- c(2,1,10,17,12:16,22,4:6,18,20,7:9,19,21,3,11,23) # 23 is sponge wt
names(grp_snp_cov)[ind]
grp_snp_cov <- grp_snp_cov[,ind]
vessel_m <- match(grp_snp_cov$VESSEL.x,vessels$VESSELID)
grp_snp_cov$VESSEL.x <- vessels$NAME[vessel_m]
names(grp_snp_cov)[3] <- 'VESSEL'

grp_snp_cov <- grp_snp_cov[which(!is.na(grp_snp_cov$AreaSwept_km2)),]

aggregate(grp_snp_cov$STAT_ZONE,by=list(year(grp_snp_cov$date)),unique)
plot(grp_snp_cov$DECSLON[which(year(grp_snp_cov$date)==2021)],
     grp_snp_cov$DECSLAT[which(year(grp_snp_cov$date)==2021)],
     xlim=range(grp_snp_cov$DECSLON),ylim=range(grp_snp_cov$DECSLAT),
     cex=grp_snp_cov$rg_pr_ab[which(year(grp_snp_cov$date)==2021)]+1,asp=1)

sort(unique(grp_snp_cov$STAT_ZONE))
strata_lon <- aggregate(grp_snp_cov$DECSLON,by=list(grp_snp_cov$STAT_ZONE),range,na.rm=T)
# strata_lon$x <- round(strata_lon$x)
strata_lat <- aggregate(grp_snp_cov$DECSLAT,by=list(grp_snp_cov$STAT_ZONE),range,na.rm=T)
# strata_lat$x <- round(strata_lat$x)

plot(strata_lon$x[,1],strata_lat$x[,1],asp=1,
     xlim=range(strata_lon$x),ylim=range(strata_lat$x))
for(i in 1:9){
  rect(strata_lon$x[i,1],strata_lat$x[i,1],
       strata_lon$x[i,2],strata_lat$x[i,2])
  text(strata_lon$x[i,1],strata_lat$x[i,1],strata_lat$Group.1[i],adj=1,cex=2)
}
### assign those in statzone 0 to 4
grp_snp_cov$STAT_ZONE[which(grp_snp_cov$STAT_ZONE==0)] <- 4
### remove stat zone 1
grp_snp_cov <- grp_snp_cov[-which(grp_snp_cov$STAT_ZONE==1),]
### make zeros
grp_snp_cov$rg_wt[which(is.na(grp_snp_cov$rg_wt))] <- 0
grp_snp_cov$rs_wt[which(is.na(grp_snp_cov$rs_wt))] <- 0
grp_snp_cov$sp_wt[which(is.na(grp_snp_cov$sp_wt))] <- 0

setwd('~/Desktop/professional/projects/Postdoc_FL/data/grouper/')
write.csv(grp_snp_cov,'grp_snp_2022_v2.csv',row.names = F) # v2 has sponge weight

yr <- sort(unique(year(grp_snp_cov$date)))
for(i in yr){
  tmp <- grp_snp_cov[which(year(grp_snp_cov$date)==i),]
  bubblePlot(tmp$DECSLON,tmp$DECSLAT,tmp$TEMP_BOT,
             xlim=range(grp_snp_cov$DECSLON),ylim=range(grp_snp_cov$DECSLAT),
             zlim=range(grp_snp_cov$TEMP_BOT,na.rm=T),asp=1)
  mtext(i)
}

yr <- sort(unique(year(grp_snp_cov$date)))
for(i in yr){
  tmp <- grp_snp_cov[which(year(grp_snp_cov$date)==i),]
  bubblePlot(tmp$DECSLON,tmp$DECSLAT,tmp$bot_do,
             xlim=range(grp_snp_cov$DECSLON),ylim=range(grp_snp_cov$DECSLAT),
             zlim=range(grp_snp_cov$bot_do,na.rm=T),asp=1)
  mtext(i)
}

yr <- sort(unique(year(grp_snp_cov$date)))
grp_snp_cov$lsp_wt <- log10(grp_snp_cov$sp_wt+1)
for(i in yr){
  tmp <- grp_snp_cov[which(year(grp_snp_cov$date)==i),]
  bubblePlot(tmp$DECSLON,tmp$DECSLAT,tmp$lsp_wt,
             xlim=range(grp_snp_cov$DECSLON),ylim=range(grp_snp_cov$DECSLAT),
             zlim=range(grp_snp_cov$lsp_wt,na.rm=T),asp=1)
  mtext(i)
}

sp_yr <- aggregate(grp_snp_cov$sp_wt,by=list(year(grp_snp_cov$date)),mean,na.rm=T)
barplot(sp_yr$x)
