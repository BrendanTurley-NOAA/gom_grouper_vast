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

setwd("~/Desktop/professional/biblioteca/data/shapefiles/gshhg-shp-2.3.7/GSHHS_shp/h/")
world <- readOGR('GSHHS_h_L1.shp')
world <- crop(world, extent(-87, -79.5, 24.3, 31))

setwd('~/Desktop/professional/projects/Postdoc_FL/data/ctd/seamap/public_seamap_csvs/')
# setwd('~/Desktop/professional/projects/Postdoc_FL/data/ctd/seamap/public_seamap_csvs/public_seamap_csvs_2022/')

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
bot_do <- bot_dos <- bot_temp <- rep(NA,length(ids))
station_id <- cruise_id <- rep(NA,length(ids))
n <- 1
for(i in ids){
  temp <- ctd_data[which(ctd_data$CTDID==i),]
  bot_do[n] <- temp$OXY_MG[which.max(temp$DEPTH)]
  bot_dos[n] <- temp$OXSAT[which.max(temp$DEPTH)]
  bot_temp[n] <- temp$TEMP[which.max(temp$DEPTH)]
  station_id[n] <- unique(temp$STATIONID)
  cruise_id[n] <- unique(temp$CRUISEID)
  n <- n + 1
}
bot_do[which(bot_do>14)] <- NA
bot_temp[which(bot_temp>40)] <- NA
bot_dos[which(bot_do!=0 & bot_dos==0)] <- NA


### ----------------- subset species -----------------
spcde <- 170021211 #red grouper
# spcde <- 170022104 # gag; not many instances
ind <- which(catch$BIO_BGS==spcde)
grouper <- catch[ind,]

spcde <-170151107 # red snapper
ind <- which(catch$BIO_BGS==spcde)
snapper <- catch[ind,]

### rename count and count extrapolation and total weight
grouper <- grouper[,c(2:3,11:12,14)]
names(grouper)[c(3:5)] <- c('rg_cnt','rg_cntexp','rg_wt')
snapper <- snapper[,c(2:3,11:12,14)]
names(snapper)[c(3:5)] <- c('rs_cnt','rs_cntexp','rs_wt')

### merge species
grp_snp <- merge(grouper,snapper,by=c('STATIONID','CRUISEID'),all=T)

### subset for only FL and NOAA cruises and 2010 onward and trawl survey data
grp_cruises <- cruises[is.element(cruises$CRUISEID,unique(grp_snp$CRUISEID)),]
FL_cruises <- grp_cruises[which(grp_cruises$SOURCE=='FL' | grp_cruises$SOURCE=='US'),]
FL_cruises <- FL_cruises[which(FL_cruises$YR>=2010 & grepl('groundfish',FL_cruises$TITLE,ignore.case = T)),]
### 2008 and 2009 excluded due to poor coverage
# FL_cruises <- FL_cruises[which(FL_cruises$YR>=2008 & grepl('groundfish',FL_cruises$TITLE,ignore.case = T)),]


### ----------------- merge data for logistic regression -----------------
grp_do <- data.frame(STATIONID=station_id,
                     CRUISEID=cruise_id,
                     bot_do=bot_do,
                     bot_dos=bot_dos)

grp <- merge(grp_do,grp_snp,by=c('STATIONID','CRUISEID'),all=T)
new <- merge(grp,station_redux,by=c('STATIONID','CRUISEID'),all.x=T)
new <- new[is.element(new$CRUISEID,FL_cruises$CRUISEID),]
new$date <- mdy_hm(paste(new$MO_DAY_YR,sprintf('%04d',as.numeric(new$TIME_MIL))))
ind <- which(is.na(new$MO_DAY_YR))
new <- new[-ind,]
### make presence/absence
new$rg_pr_ab <- ifelse(!is.na(new$rg_cnt),1,0)
new$rs_pr_ab <- ifelse(!is.na(new$rs_cnt),1,0)
### cpue by number
new$rg_cpue <- new$rg_cntexp/new$trawl_hrs
new$rs_cpue <- new$rs_cntexp/new$trawl_hrs
### cpue by weight
new$rg_cpue_wt <- new$rg_wt/new$trawl_hrs
new$rs_cpue_wt <- new$rs_wt/new$trawl_hrs
### water depth
new <- merge(new,envrec,by=c('STATIONID','CRUISEID'),all.x=T)
### which depth
par(mfrow=c(1,2))
plot(new$DECSLON,new$DECSLAT,cex=new$DEPTH_EMAX/100,asp=1)
plot(new$DECSLON,new$DECSLAT,cex=-new$DEPTH_EWTR/100,asp=1)
### EMAX seems more reasonable
new$DEPTH_EMAX <- (-new$DEPTH_EMAX)


# ind <- c(1:10,18,27:29,33:39,52)
ind <- c(1:11,19,28:30,34:40,53) # with bot_dos
names(new)[ind]
# [1] "STATIONID"     "CRUISEID"      "bot_do"        "rg_cnt"        "rg_cntexp"     "rg_wt"        
# [7] "rs_cnt"        "rs_cntexp"     "rs_wt"         "VESSEL.x"      "TEMP_BOT"      "STAT_ZONE"    
# [13] "DECSLAT"       "DECSLON"       "trawl_hrs"     "AreaSwept_km2" "date"          "rg_pr_ab"     
# [19] "rs_pr_ab"      "rg_cpue"       "rs_cpue"       "DEPTH_ESRF"  
grp_snp_cov <- new[,ind]
# ind <- c(2,1,10,17,12:16,22,4:6,18,20,7:9,19,21)
ind <- c(2,1,11,18,13:17,23,5:7,19,21,8:10,20,22) # with bot_dos
names(grp_snp_cov)[ind]
# [1] "CRUISEID"      "STATIONID"     "VESSEL.x"      "date"          "STAT_ZONE"     "DECSLAT"      
# [7] "DECSLON"       "trawl_hrs"     "AreaSwept_km2" "DEPTH_ESRF"    "rg_cnt"        "rg_cntexp"    
# [13] "rg_wt"         "rg_pr_ab"      "rg_cpue"       "rs_cnt"        "rs_cntexp"     "rs_wt"        
# [19] "rs_pr_ab"      "rs_cpue"   
grp_snp_cov <- grp_snp_cov[,ind]
vessel_m <- match(grp_snp_cov$VESSEL.x,vessels$VESSELID)
grp_snp_cov$VESSEL.x <- vessels$NAME[vessel_m]
names(grp_snp_cov)[3] <- 'VESSEL'

grp_snp_cov <- grp_snp_cov[which(!is.na(grp_snp_cov$AreaSwept_km2)),]

# aggregate(grp_snp_cov$STAT_ZONE,by=list(year(grp_snp_cov$date)),unique)
# plot(grp_snp_cov$DECSLON[which(year(grp_snp_cov$date)==2011)],
#      grp_snp_cov$DECSLAT[which(year(grp_snp_cov$date)==2011)],
#      xlim=range(grp_snp_cov$DECSLON),ylim=range(grp_snp_cov$DECSLAT))

setwd('~/Desktop/professional/projects/Postdoc_FL/data/grouper/')
write.csv(grp_snp_cov,'grp_snp.csv',row.names = F)

### ----------------- plot some distributions -----------------
### aggregate mean and SE
cnts_yr <- aggregate(new$rg_cntexp,by=list(year(new$date)),mean,na.rm=T)
se_yr <- aggregate(new$rg_cntexp,by=list(year(new$date)),std.er,na.rm=T)

# setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
# png("rg_cnts_stations.png", height = 5, width = 7, units = 'in', res=300)
b <- barplot(cnts_yr$x,names=cnts_yr$Group.1,
             las=1,ylim=c(0,5),
             ylab='Mean count per station',main='Red grouper - SEAMAP trawl')
arrows(b,cnts_yr$x-se_yr$x,b,cnts_yr$x+se_yr$x,angle=90,length=.05,code=3)
# dev.off()

### proportion of positive stations per year
pos <- aggregate(new$rg_pr_ab,by=list(year(new$date)),sum,na.rm=T)
tot <- aggregate(new$rg_pr_ab,by=list(year(new$date)),length)
barplot(pos$x/tot$x,names=pos$Group.1,ylab='Proportion postive',las=1)

rg_pos <- new[which(new$rg_pr_ab>0),]
pos_cnt_yr <- aggregate(rg_pos$rg_cpue_wt,by=list(year(rg_pos$date)),mean,na.rm=T)
pos_se_yr <- aggregate(rg_pos$rg_cpue_wt,by=list(year(rg_pos$date)),std.er,na.rm=T)

b <- barplot(pos_cnt_yr$x,names=pos_cnt_yr$Group.1,
             las=1,ylim=c(0,6),
             ylab='Mean cpue (kg/hr) per station',main='Red grouper - SEAMAP trawl')
arrows(b,pos_cnt_yr$x-pos_se_yr$x,b,pos_cnt_yr$x+pos_se_yr$x,angle=90,length=.05,code=3)



plot(new$DECSLON,new$DECSLAT,asp=1)
plot(new$bot_do,new$rg_pr_ab)
plot(new$bot_dos,new$rg_pr_ab)
boxplot(new$bot_do[which(new$rg_pr_ab>0)],new$bot_do[which(new$rg_pr_ab==0)])
boxplot(new$bot_dos[which(new$rg_pr_ab>0)],new$bot_dos[which(new$rg_pr_ab==0)])
boxplot(new$TEMP_BOT[which(new$rg_pr_ab>0)],new$TEMP_BOT[which(new$rg_pr_ab==0)])
plot(new$bot_do,new$rg_cntexp)
plot(new$bot_do,new$rg_cpue)
plot(new$bot_dos,new$rg_cntexp)
plot(new$bot_dos,new$rg_cpue_wt)
plot(new$bot_dos,new$rg_cpue)

plot(new$bot_do,new$rs_pr_ab)
boxplot(new$bot_do[which(new$rs_pr_ab>0)],new$bot_do[which(new$rs_pr_ab==0)])
boxplot(new$TEMP_BOT[which(new$rs_pr_ab>0)],new$TEMP_BOT[which(new$rs_pr_ab==0)])
plot(new$bot_do,new$rs_cntexp)
plot(new$bot_do,new$rs_cpue)


boxplot(new$bot_do[which(new$rg_pr_ab>0)]~
          year(new$date[which(new$rg_pr_ab>0)]),
        at=seq(1,20,2),xlim=c(0,20),ylim=c(0,10),lty=1)
boxplot(new$bot_do[which(new$rg_pr_ab==0)]~
          year(new$date[which(new$rg_pr_ab==0)]),
        at=seq(2,20,2),add=T,xaxt='n',col=2,lty=1,lwd=2,yaxt='n')

boxplot(new$bot_dos[which(new$rg_pr_ab>0)]~
          year(new$date[which(new$rg_pr_ab>0)]),
        at=seq(1,20,2),xlim=c(0,20),ylim=c(0,150),lty=1)
boxplot(new$bot_dos[which(new$rg_pr_ab==0)]~
          year(new$date[which(new$rg_pr_ab==0)]),
        at=seq(2,20,2),add=T,xaxt='n',col=2,lty=1,lwd=2,yaxt='n')


### bottom DO
d1 <- density(new$bot_do[which(new$rg_pr_ab>0)],na.rm=T)
d2 <- density(new$bot_do[which(new$rg_pr_ab==0)],na.rm=T)

# setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
# png("rg_botdo.png", height = 5, width = 8, units = 'in', res=300)
plot(d1$x,d1$y,
     xlim=c(0,9.5),ylim=c(0,.6),las=1,
     typ='l',lwd=2,
     xlab='Bottom dissolved oxygen',ylab='Density')
points(d2$x,d2$y,typ='l',col=2,lwd=2)
legend('topleft',
       c(paste0('Present (n = ',d1$n,')'),paste0('Absent (n = ',d2$n,')')),
       col=c(1,2),lty=1,lwd=2,bty='n')
abline(v=c(median(new$bot_do[which(new$rg_pr_ab>0)],na.rm=T),
           median(new$bot_do[which(new$rg_pr_ab==0)],na.rm=T)),
       lty=5,col=c(1,2))
mtext('Red grouper - SEAMAP trawl')
# dev.off()

### bottom DO sat
d1 <- density(new$bot_dos[which(new$rg_pr_ab>0)],na.rm=T)
d2 <- density(new$bot_dos[which(new$rg_pr_ab==0)],na.rm=T)

# setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
# png("rg_botdo.png", height = 5, width = 8, units = 'in', res=300)
plot(d1$x,d1$y,
     xlim=c(0,150),ylim=c(0,.05),las=1,
     typ='l',lwd=2,
     xlab='Bottom dissolved oxygen',ylab='Density')
points(d2$x,d2$y,typ='l',col=2,lwd=2)
legend('topleft',
       c(paste0('Present (n = ',d1$n,')'),paste0('Absent (n = ',d2$n,')')),
       col=c(1,2),lty=1,lwd=2,bty='n')
abline(v=c(median(new$bot_dos[which(new$rg_pr_ab>0)],na.rm=T),
           median(new$bot_dos[which(new$rg_pr_ab==0)],na.rm=T)),
       lty=5,col=c(1,2))
mtext('Red grouper - SEAMAP trawl')
# dev.off()

### bottom temperature
d1 <- density(new$TEMP_BOT[which(new$rg_pr_ab>0)],na.rm=T)
d2 <- density(new$TEMP_BOT[which(new$rg_pr_ab==0)],na.rm=T)

# setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
# png("rg_bottemp.png", height = 5, width = 8, units = 'in', res=300)
plot(d1$x,d1$y,
     xlim=c(12,34),ylim=c(0,.13),las=1,
     typ='l',lwd=2,
     xlab='Bottom temperature',ylab='Density')
points(d2$x,d2$y,typ='l',col=2,lwd=2)
legend('topleft',
       c(paste0('Present (n = ',d1$n,')'),paste0('Absent (n = ',d2$n,')')),
       col=c(1,2),lty=1,lwd=2,bty='n')
abline(v=c(median(new$TEMP_BOT[which(new$rg_pr_ab>0)],na.rm=T),
           median(new$TEMP_BOT[which(new$rg_pr_ab==0)],na.rm=T)),
       lty=5,col=c(1,2))
mtext('Red grouper - SEAMAP trawl')
# dev.off()

plot(density(new$TEMP_BOT[which(new$rg_pr_ab>0)],na.rm=T))
lines(density(new$TEMP_BOT[which(new$rg_pr_ab==0)],na.rm=T),col=2)

plot(density(new$bot_do[which(new$rs_pr_ab>0)],na.rm=T))
lines(density(new$bot_do[which(new$rs_pr_ab==0)],na.rm=T),col=2)

plot(density(new$bot_dos[which(new$rs_pr_ab>0)],na.rm=T))
lines(density(new$bot_dos[which(new$rs_pr_ab==0)],na.rm=T),col=2)

plot(density(new$TEMP_BOT[which(new$rs_pr_ab>0)],na.rm=T))
lines(density(new$TEMP_BOT[which(new$rs_pr_ab==0)],na.rm=T),col=2)

plot(density(new$rg_cnt[which(new$rs_pr_ab>0)],na.rm=T))
lines(density(new$rg_cnt[which(new$rs_pr_ab==0)],na.rm=T),col=2)

plot(density(new$rs_cnt[which(new$rg_pr_ab>0)],na.rm=T))
lines(density(new$rs_cnt[which(new$rg_pr_ab==0)],na.rm=T),col=2)

plot(density(new$DEPTH_EMAX[which(new$rg_pr_ab>0)],na.rm=T))
lines(density(new$DEPTH_EMAX[which(new$rg_pr_ab==0)],na.rm=T),col=2)

which(!is.na(new$rg_cntexp) &!is.na(new$rs_cntexp))
plot(new$rg_cntexp,new$rs_cntexp)

boxplot((new$rg_wt/new$rg_cntexp)~year(new$date))

### ----------------- logistic regression -----------------
### could model prescence/absence in response to DO
### data.frame: station ID, date, lon, lat, grouper i/o, bottom DO
pos <- aggregate(new$rg_pr_ab,by=list(year(new$date)),sum,na.rm=T)
tot <- aggregate(new$rg_pr_ab,by=list(year(new$date)),length)
barplot(pos$x/tot$x,names=pos$Group.1,ylab='Proportion postive',las=1)
### set aside
new2 <- new
new$year <- as.factor(year(new$date))
# new <- new[which(new$bot_do<4),]
plot(new$bot_do,new$rg_pr_ab)
model <- glm(rg_pr_ab ~ bot_do + TEMP_BOT + DEPTH_EMAX + year, data=new, family=binomial(link='logit'))
# model2 <- glm(rg_pr_ab ~ bot_do + TEMP_BOT,data=new, family=binomial(link='logit'))
# model3 <- glm(rg_pr_ab ~ bot_do,data=new, family=binomial(link='logit'))
m <- summary(model)
m$null.deviance-m$deviance
# m2 <- summary(model2)
# m2$null.deviance-m2$deviance
# m3 <- summary(model3)
# m3$null.deviance-m3$deviance
# AIC(model,model2,model3)
anova(model, test="Chisq")
# plot(model)

newdata <- data.frame(bot_do=seq(min(new$bot_do,na.rm=T),
                                 max(new$bot_do,na.rm=T),
                                 length.out=nrow(new)),
                      TEMP_BOT=seq(min(new$TEMP_BOT,na.rm=T),
                                 max(new$TEMP_BOT,na.rm=T),
                                 length.out=nrow(new)),
                      DEPTH_EMAX=seq(min(new$DEPTH_EWTR,na.rm=T),
                                   max(new$DEPTH_EWTR,na.rm=T),
                                   length.out=nrow(new)),
                      year=levels(new$year))
newdata$rg_pr_ab = predict(model, newdata, type="response")
plot(new$bot_do,new$rg_pr_ab)
for(i in 2010:2019){
  temp <- newdata[which(newdata$year==i),]
  points(temp$bot_do,temp$rg_pr_ab,typ='l',col=2)
}
points(newdata$bot_do,newdata$rg_pr_ab,typ='l')

boxplot(newdata$rg_pr_ab~newdata$year)

cols <- colorRampPalette(c('forestgreen','gray80','purple'))
ncols <- cols(length(2010:2019))
setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png("rg_logistic_do_temp_z_yr.png", height = 5, width = 5, units = 'in', res=300)
plot(new$bot_do,new$rg_pr_ab,col=alpha(1,.3),xlab='Bottom DO',ylab='Probability of RG',pch='|')
for(i in 2010:2019){
  temp <- newdata[which(newdata$year==i),]
  points(temp$bot_do,temp$rg_pr_ab,col=ncols[i-2009],typ='l',lwd=2)
}
# points(newdata$bot_do,newdata$rg_pr_ab,col='darkorange',typ='l',lwd=2)
mtext('Logistic regression results')
legend('topleft',legend=c(2010:2019),col=ncols,lty=1,bty='n',lwd=1.5,cex=.7)
dev.off()


roc <- prediction(newdata$rg_pr_ab,new$rg_pr_ab)
perf <- performance(roc, "tpr", "fpr")
plot(perf, colorize=TRUE)
unlist(slot(performance(roc, "auc"), "y.values"))


model2 <- glm(rg_pr_ab ~ bot_do + TEMP_BOT + DEPTH_EMAX, data=new, family=binomial(link='logit'))
m2 <- summary(model2)
m2$null.deviance-m2$deviance

newdata2 <- data.frame(bot_do=seq(min(new$bot_do,na.rm=T),
                                 max(new$bot_do,na.rm=T),
                                 length.out=nrow(new)),
                      TEMP_BOT=seq(min(new$TEMP_BOT,na.rm=T),
                                   max(new$TEMP_BOT,na.rm=T),
                                   length.out=nrow(new)),
                      DEPTH_EMAX=seq(min(new$DEPTH_EWTR,na.rm=T),
                                     max(new$DEPTH_EWTR,na.rm=T),
                                     length.out=nrow(new)))
newdata2$rg_pr_ab = predict(model2, newdata2, type="response")

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png("rg_logistic_do_temp_z.png", height = 5, width = 5, units = 'in', res=300)
plot(new$bot_do,new$rg_pr_ab,col=alpha(1,.3),xlab='Bottom DO',ylab='Probability of RG',pch='|')
points(newdata2$bot_do,newdata2$rg_pr_ab,col='darkorange',typ='l',lwd=2)
mtext('Logistic regression results')
dev.off()

roc <- prediction(newdata2$rg_pr_ab,new$rg_pr_ab)
perf <- performance(roc, "tpr", "fpr")
plot(perf, colorize=TRUE)
unlist(slot(performance(roc, "auc"), "y.values"))

### ----------------- load CTD data -----------------
### load CTD groomed data
setwd("~/Desktop/professional/projects/Postdoc_FL/data/ctd/groomed")

### version 3 of combined_groomed data standardized subsetting and removal of bad data
combined <- read.csv('ctd_combined_groomed_v6.csv')
combined$date <- as.Date(combined$date)
combined_stations <- read.csv('ctd_combined_stations_groomed_v6.csv')
combined_stations$date <- as.Date(combined_stations$date)
combined <- combined[which(year(combined$date)>=2004),]
combined_stations <- combined_stations[which(year(combined_stations$date)>=2004),]
### remove gross differences
# combined <- combined[-which(abs(combined$z_crm-combined$z_max)>15),]
gr_diff <- which(abs(combined$z_crm-combined$z_station)>15)
combined$bot_do[gr_diff] <- NA
combined$bot_t[gr_diff] <- NA
### CRM depth and CRM depth with hypoxia
ind_crm <- which(combined$bot_code=='crm')
combined$bot_code[-ind_crm] <- 'not'
# combined$bot_code[which(is.na(combined$bot_code))] <- 'not'
# ind_hyp <- which(combined$bot_do<=3.5 & combined$bot_code=='crm')
ind_nhyp <- which(combined$bot_do>2 & combined$bot_do<=3.5)
ind_hyp <- which(combined$bot_do<=2)
ind_lodo <- which(combined$bot_do<=3.5)
lo_do <- combined[ind_lodo,]


### ----------------- plots some years -----------------
# setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
# png("rg_trawl.png", height = 5, width = 10, units = 'in', res=300)
# par(mfrow=c(2,5),mar=c(2.5,2.5,1.5,1),oma=c(0,0,2,0))
# for(i in 2010:2019){
#   temp <- new[which(year(new$date)==i & month(new$date)<9),]
#   temp2 <- lo_do[which(year(lo_do$date)==i & month(lo_do$date)>5 & month(lo_do$date)<11),]
#   plot(temp$DECSLON[which(is.na(temp$rg_cntexp))],
#        temp$DECSLAT[which(is.na(temp$rg_cntexp))],asp=1,pch=18,col='deepskyblue',
#        xlim=c(-86,-81),ylim=c(24.5,30.5),
#        xlab='',ylab='')
#   points(temp$DECSLON[which(!is.na(temp$rg_cntexp))],
#          temp$DECSLAT[which(!is.na(temp$rg_cntexp))],pch=21,bg='deeppink',
#          cex=log10(temp$rg_cntexp[which(!is.na(temp$rg_cntexp))])+1)
#   points(temp2$lon,temp2$lat,bg='gold',pch=24)
#   plot(world,add=T)
#   mtext(i)
#   mtext('Red grouper - SEAMAP summer trawl',outer=T)
# }
# dev.off()

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png("rg_trawl.png", height = 5, width = 10, units = 'in', res=300)
par(mfrow=c(2,5),mar=c(2.5,2.5,1.5,1),oma=c(0,0,2,0))
for(i in 2010:2019){
  temp <- new[which(year(new$date)==i & month(new$date)<9),]
  temp2 <- lo_do[which(year(lo_do$date)==i & month(lo_do$date)>5 & month(lo_do$date)<11),]
  plot(temp$DECSLON[which(is.na(temp$rg_cpue))],
       temp$DECSLAT[which(is.na(temp$rg_cpue))],asp=1,pch=18,col='deepskyblue',
       xlim=c(-86,-81),ylim=c(24.5,30.5),
       xlab='',ylab='')
  points(temp$DECSLON[which(!is.na(temp$rg_cpue))],
         temp$DECSLAT[which(!is.na(temp$rg_cpue))],pch=21,bg='deeppink',
         cex=log10(temp$rg_cpue[which(!is.na(temp$rg_cpue))])+1)
  points(temp2$lon,temp2$lat,bg='gold',pch=24)
  plot(world,add=T)
  mtext(i)
  mtext('Red grouper - SEAMAP summer trawl',outer=T)
}
dev.off()

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png("rg_trawl_wt.png", height = 5, width = 10, units = 'in', res=300)
par(mfrow=c(2,5),mar=c(2.5,2.5,1.5,1),oma=c(0,0,2,0))
for(i in 2010:2019){
  temp <- new[which(year(new$date)==i & month(new$date)<9),]
  temp2 <- lo_do[which(year(lo_do$date)==i & month(lo_do$date)>5 & month(lo_do$date)<11),]
  plot(temp$DECSLON[which(is.na(temp$rg_wt))],
       temp$DECSLAT[which(is.na(temp$rg_wt))],asp=1,pch=18,col='deepskyblue',
       xlim=c(-86,-81),ylim=c(24.5,30.5),
       xlab='',ylab='')
  points(temp$DECSLON[which(!is.na(temp$rg_wt))],
         temp$DECSLAT[which(!is.na(temp$rg_wt))],pch=21,bg='deeppink',
         cex=log10(temp$rg_wt[which(!is.na(temp$rg_wt))])+1)
  points(temp2$lon,temp2$lat,bg='gold',pch=24)
  plot(world,add=T)
  mtext(i)
  mtext('Red grouper - SEAMAP summer trawl',outer=T)
}
dev.off()

# setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
# png("rs_trawl.png", height = 5, width = 10, units = 'in', res=300)
# par(mfrow=c(2,5),mar=c(2.5,2.5,1.5,1),oma=c(0,0,2,0))
# for(i in 2010:2019){
#   temp <- new[which(year(new$date)==i & month(new$date)<9),]
#   temp2 <- lo_do[which(year(lo_do$date)==i & month(lo_do$date)>5 & month(lo_do$date)<11),]
#   plot(temp$DECSLON[which(is.na(temp$rs_cntexp))],
#        temp$DECSLAT[which(is.na(temp$rs_cntexp))],asp=1,pch=18,col='deepskyblue',
#        xlim=c(-86,-81),ylim=c(24.5,30.5),
#        xlab='',ylab='')
#   points(temp$DECSLON[which(!is.na(temp$rs_cntexp))],
#          temp$DECSLAT[which(!is.na(temp$rs_cntexp))],pch=21,bg='deeppink',
#          cex=log10(temp$rs_cntexp[which(!is.na(temp$rs_cntexp))])+1)
#   points(temp2$lon,temp2$lat,bg='gold',pch=24)
#   plot(world,add=T)
#   mtext(i)
#   mtext('Red snapper - SEAMAP summer trawl',outer=T)
# }
# dev.off()

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png("rs_trawl.png", height = 5, width = 10, units = 'in', res=300)
par(mfrow=c(2,5),mar=c(2.5,2.5,1.5,1),oma=c(0,0,2,0))
for(i in 2010:2019){
  temp <- new[which(year(new$date)==i & month(new$date)<9),]
  temp2 <- lo_do[which(year(lo_do$date)==i & month(lo_do$date)>5 & month(lo_do$date)<11),]
  plot(temp$DECSLON[which(is.na(temp$rs_cpue))],
       temp$DECSLAT[which(is.na(temp$rs_cpue))],asp=1,pch=18,col='deepskyblue',
       xlim=c(-86,-81),ylim=c(24.5,30.5),
       xlab='',ylab='')
  points(temp$DECSLON[which(!is.na(temp$rs_cpue))],
         temp$DECSLAT[which(!is.na(temp$rs_cpue))],pch=21,bg='deeppink',
         cex=log10(temp$rs_cpue[which(!is.na(temp$rs_cpue))])+1)
  points(temp2$lon,temp2$lat,bg='gold',pch=24)
  plot(world,add=T)
  mtext(i)
  mtext('Red snapper - SEAMAP summer trawl',outer=T)
}
dev.off()


### ----------------- center of gravity -----------------
### by CPUE
com <- new[which(new$rg_pr_ab>0),]
# w <- com$rg_cntexp/sum(com$rg_cntexp,na.rm=T)
w <- com$rg_cpue/sum(com$rg_cpue,na.rm=T)
lon_wm <- weighted.mean(com$DECSLON,w)
lat_wm <- weighted.mean(com$DECSLAT,w)

plot(com$DECSLON,com$DECSLAT,asp=1)
points(lon_wm,lat_wm,col=2,pch=16,cex=2)

lat_com <- rep(NA,length(2010:2019))
par(mfrow=c(2,5))
for(i in 2010:2019){
  com <- new[which(new$rg_pr_ab>0 & year(new$date)==i),]
  # w <- com$rg_cntexp/sum(com$rg_cntexp,na.rm=T)
  w <- com$rg_cpue/sum(com$rg_cpue,na.rm=T)
  lon_wm <- weighted.mean(com$DECSLON,w)
  lat_wm <- weighted.mean(com$DECSLAT,w)
  # lon_wm <- mean(com$DECSLON,na.rm=T)
  # lat_wm <- mean(com$DECSLAT,na.rm=T)
  lat_com[i-2009] <- lat_wm
  
  plot(com$DECSLON,com$DECSLAT,asp=1,xlim=range(new$DECSLON),ylim=range(new$DECSLAT))
  points(lon_wm,lat_wm,col=2,pch=16,cex=2)
  mtext(i)
}

com <- new[which(new$rg_pr_ab>0),]
lat_range <- aggregate(com$DECSLAT,by=list(year(com$date)),range,na.rm=T)
names(lat_range) <- c('year','range')

# lat_quants <- aggregate(com$DECSLAT,by=list(year(com$date)),quantile,c(.25,.75),na.rm=T)

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png("rg_cog_range.png", height = 5, width = 7, units = 'in', res=300)
par(mfrow=c(1,1))
plot(2010:2019,lat_com,
     ylim=range(lat_range[,2]),
     typ='b',las=1,xlab='',ylab='Latitude',xaxt='n',lwd=2)
polygon(c(2010:2019,rev(2010:2019)),
        c(lat_range$range[,1],rev(lat_range$range[,2])),col=alpha(3,.2))
# polygon(c(2010:2019,rev(2010:2019)),
#         c(lat_quants$x[,1],rev(lat_quants$x[,2])),col=alpha(2,.3))
# arrows(2010:2019,lat_range$range[,1],2010:2019,lat_range$range[,2],code=3,angle=90,length=.05)
axis(1,2010:2019)
mtext('Latitude center of gravity (by CPUE)')
dev.off()

### by weight
com <- new[which(new$rg_pr_ab>0),]
w <- com$rg_wt/sum(com$rg_wt,na.rm=T)
lon_wm <- weighted.mean(com$DECSLON,w)
lat_wm <- weighted.mean(com$DECSLAT,w)

plot(com$DECSLON,com$DECSLAT,asp=1)
points(lon_wm,lat_wm,col=2,pch=16,cex=2)

lat_com <- rep(NA,length(2010:2019))
par(mfrow=c(2,5))
for(i in 2010:2019){
  com <- new[which(new$rg_pr_ab>0 & year(new$date)==i),]
  w <- com$rg_wt/sum(com$rg_wt,na.rm=T)
  lon_wm <- weighted.mean(com$DECSLON,w)
  lat_wm <- weighted.mean(com$DECSLAT,w)
  # lon_wm <- mean(com$DECSLON,na.rm=T)
  # lat_wm <- mean(com$DECSLAT,na.rm=T)
  lat_com[i-2009] <- lat_wm
  
  plot(com$DECSLON,com$DECSLAT,asp=1,xlim=range(new$DECSLON),ylim=range(new$DECSLAT))
  points(lon_wm,lat_wm,col=2,pch=16,cex=2)
  mtext(i)
}

com <- new[which(new$rg_pr_ab>0),]
lat_range <- aggregate(com$DECSLAT,by=list(year(com$date)),range,na.rm=T)
names(lat_range) <- c('year','range')

# lat_quants <- aggregate(com$DECSLAT,by=list(year(com$date)),quantile,c(.25,.75),na.rm=T)

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png("rg_cog_range_wt.png", height = 5, width = 7, units = 'in', res=300)
par(mfrow=c(1,1))
plot(2010:2019,lat_com,
     ylim=range(lat_range[,2]),
     typ='b',las=1,xlab='',ylab='Latitude',xaxt='n',lwd=2)
polygon(c(2010:2019,rev(2010:2019)),
        c(lat_range$range[,1],rev(lat_range$range[,2])),col=alpha(3,.2))
# polygon(c(2010:2019,rev(2010:2019)),
#         c(lat_quants$x[,1],rev(lat_quants$x[,2])),col=alpha(2,.3))
# arrows(2010:2019,lat_range$range[,1],2010:2019,lat_range$range[,2],code=3,angle=90,length=.05)
axis(1,2010:2019)
mtext('Latitude center of gravity (by weight)')
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

plot(bbl$Year,bbl$Scaled_Index,
     ylim=c(min(bbl$LCL),max(bbl$UCL)),
     typ='b',las=1)
polygon(c(bbl$Year,rev(bbl$Year)),
        c(bbl$LCL,rev(bbl$UCL)),
        col=alpha(4,.2))
points(bbl$Year,bbl$Scaled_Index,
       typ='o',lwd=2,pch=16,lty=1)

polygon(c(trawl$Year,rev(trawl$Year)),
        c(trawl$LCL,rev(trawl$UCL)),
        col=alpha(2,.2))
points(trawl$Year,trawl$Scaled_Index,
       typ='o',lwd=2,pch=17,lty=2)

polygon(c(hb$Year,rev(hb$Year)),
        c(hb$LCL,rev(hb$UCL)),
        col=alpha(3,.2))
points(hb$Year,hb$relative_index,
       typ='o',lwd=2,pch=15,lty=3)

points(video$Year,video$std_indx,
       typ='o',lwd=2,pch=18,lty=4)

bbl <- add_years(bbl)
trawl <- add_years(trawl)
hb <- add_years(hb)
video <- add_years(video)
indx <- merge(bbl,trawl,by=c('Year'),all=T)
indx <- merge(indx,hb,by=c('Year'),all=T)
indx <- merge(indx,video,by=c('Year'),all=T)
indx_en <- apply(indx[,c(5,12,20,27)],1,median,na.rm=T)

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png("rg_indices.png", height = 5, width = 8, units = 'in', res=300)
plot(bbl$Year,bbl$Scaled_Index,
     typ='o',lwd=2,pch=16,lty=1,
     ylim=c(0,2.5),xaxt='n',xlab='',ylab='Scaled index')
axis(1,2000:2017)
points(trawl$Year,trawl$Scaled_Index,
       typ='o',lwd=2,pch=17,lty=2,col='gray20')
points(hb$Year,hb$relative_index,
       typ='o',lwd=2,pch=15,lty=3,col='gray40')
points(video$Year,video$std_indx,
       typ='o',lwd=2,pch=18,lty=4,col='gray60')
legend('topright',c('LL','Trawl','HB','Video','Ensemble'),
       lty=c(1,2,3,4,1),pch=c(16,17,15,18,NA),
       col=c(1,'gray20','gray40','gray60',2),bty='n',lwd=2)
points(indx$Year,indx_en,typ='l',col=2,lwd=2,pch=20)
mtext('RG indices - SEDAR 61')
abline(v=seq(2000,2017,5),lty=5,col='gray60')
dev.off()


### ----------------- length frequency data -----------------
setwd('~/Desktop/professional/projects/Postdoc_FL/data/ctd/seamap/public_seamap_csvs')
lth_freq <- read.csv('GLFREC.csv')
lth_freq <- lth_freq[is.element(lth_freq$CRUISEID,new$CRUISEID),]
spcde <- 170021211 #red grouper
# spcde <- 170022104 # gag; not many instances
ind <- which(lth_freq$BIO_GLF==spcde)
rg_lth_freq <- lth_freq[ind,]
### subset for red grouper
spcde <-170151107 # red snapper
ind <- which(lth_freq$BIO_GLF==spcde)
rs_lth_freq <- lth_freq[ind,]
rm(lth_freq)
### find total weight per station
rg_wts <- aggregate(rg_lth_freq$INDVL_WT,by=list(rg_lth_freq$CRUISEID,rg_lth_freq$STATIONID),sum,na.rm=T)
names(rg_wts) <- c('CRUISEID','STATIONID','tot_wts')
new3 <- merge(new,rg_wts,by=c('CRUISEID','STATIONID'),all.x=T)
plot(new3$rg_cnt,new3$tot_wts)

hist(rg_lth_freq$LEN_GLF,breaks=seq(0,900,50),xlab='Length (mm)',main='Red Grouper',las=1)
abline(v=18*25.4,col=2,lwd=2,lty=5)

hist(rg_lth_freq$INDVL_WT,xlab='Length (mm)',main='Red Grouper',las=1)

rg_lth_freq_n <- merge(rg_lth_freq,station_redux,by=c('STATIONID','CRUISEID'),all.x=T)
rg_lth_freq_n$date <- mdy_hm(paste(rg_lth_freq_n$MO_DAY_YR,sprintf('%04d',as.numeric(rg_lth_freq_n$TIME_MIL))))

rg_lth_yr <- aggregate(rg_lth_freq_n$LEN_GLF,list(year(rg_lth_freq_n$date)),median,na.rm=T)
rg_lth_yrl <- aggregate(rg_lth_freq_n$LEN_GLF,list(year(rg_lth_freq_n$date)),length)

b <- barplot(rg_lth_yr$x,names=rg_lth_yr$Group.1,ylim=c(0,375),
             ylab='Median Length (mm)',las=1,main='Red grouper - SEAMAP trawl')
text(b,rg_lth_yr$x+10,paste('n=',rg_lth_yrl$x))

# setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
# png("rg_lth_freq.png", height = 5, width = 8, units = 'in', res=300)
boxplot(rg_lth_freq_n$LEN_GLF~year(rg_lth_freq_n$date),
        xlab='Year',ylab='Length (mm)',main='Red grouper - SEAMAP trawl',
        varwidth=T,las=1,lty=1)
abline(h=c(median(rg_lth_freq_n$LEN_GLF,na.rm=T),18*25.4),
       col=c(1,2),lwd=2,lty=5)
# dev.off()

y_d <- matrix(NA,10,512)
setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png("rg_lth_freq2.png", height = 5, width = 8, units = 'in', res=300)
plot(1,1,col='white',xlim=c(2009,2019),ylim=c(0,1000),
     xlab='Year',ylab='Red grouper length (mm)')
abline(v=2010:2019,col='gray50')
for(i in 2010:2019){
  tmp <- rg_lth_freq_n[which(year(rg_lth_freq_n$date)==i),]
  dt <- density(tmp$LEN_GLF,from=0,to=1000)
  y_d[i-2009,] <- dt$y
  print(range(dt$x))
  points(i-dt$y*200,dt$x,typ='l')
  # polygon(c(i-dt$y*200,i-dt$y[1]*200),c(dt$x,dt$x[1]),col=2)
  polygon(c(i,i-dt$y*200,i,i),c(dt$x[1],dt$x,dt$x[length(dt$x)],dt$x[1]),col=2)
  # points(i-dt$y[which.max(dt$y)]*200,median(tmp$LEN_GLF,na.rm=T),pch='-',cex=2)
  xxx <- which.min(abs(dt$x-median(tmp$LEN_GLF,na.rm=T)))
  points(i-dt$y[xxx]*200,median(tmp$LEN_GLF,na.rm=T),pch='_',cex=2)
}
abline(h=c(median(rg_lth_freq_n$LEN_GLF,na.rm=T),18*25.4),
       col=c('gray80',1),lwd=2,lty=5)
dev.off()

image(2010:2019,dt$x,y_d)


### mean wt/yr
rg_lth_freq_n$sz_class <- cut(rg_lth_freq_n$LEN_GLF,breaks=seq(0,1000,200))
aggregate(rg_lth_freq_n$LEN_GLF,by=list(year(rg_lth_freq_n$date)),mean,na.rm=T)

mwt_sz <- aggregate(rg_lth_freq_n$LEN_GLF,
                    by=list(year(rg_lth_freq_n$date),rg_lth_freq_n$sz_class),
                    length)

sz_class <- levels(mwt_sz$Group.2)
plot(mwt_sz$Group.1,mwt_sz$x,col='white')
for(i in 1:length(sz_class)){
  points(mwt_sz$Group.1[which(mwt_sz$Group.2==sz_class[i])],
         mwt_sz$x[which(mwt_sz$Group.2==sz_class[i])],
         col=i,typ='b')  
}

par(mfrow=c(2,1))
boxplot(rg_lth_freq_n$INDVL_WT~year(rg_lth_freq_n$date))
boxplot(rg_lth_freq_n$LEN_GLF~year(rg_lth_freq_n$date))

# rg_lth_freq_n$lth_wt <- rg_lth_freq_n$LEN_GLF/rg_lth_freq_n$INDVL_WT/1000
# aggregate(rg_lth_freq_n$lth_wt,
#           by=list(year(rg_lth_freq_n$date)),mean,na.rm=T)
# 
# mwt_sz <- aggregate(rg_lth_freq_n$lth_wt,
#                     by=list(year(rg_lth_freq_n$date),rg_lth_freq_n$sz_class),
#                     mean,na.rm=T)
# 
# sz_class <- levels(mwt_sz$Group.2)
# plot(mwt_sz$Group.1,mwt_sz$x,col='white')
# for(i in 1:length(sz_class)){
#   points(mwt_sz$Group.1[which(mwt_sz$Group.2==sz_class[i])],
#          mwt_sz$x[which(mwt_sz$Group.2==sz_class[i])],
#          col=i,typ='b')  
# }


plot(rg_lth_freq$LEN_GLF,rg_lth_freq$INDVL_WT,col=year(rg_lth_freq_n$date),log='xy')

lth_wt <- data.frame(wt=rg_lth_freq$INDVL_WT,
                     len=rg_lth_freq$LEN_GLF,
                     l_wt=log10(rg_lth_freq$INDVL_WT),
                     l_len=log10(rg_lth_freq$LEN_GLF),
                     year=year(rg_lth_freq_n$date),
                     month=month(rg_lth_freq_n$date))

mod1 <- lm(l_wt~l_len+year,data=lth_wt)
summary(mod1)
confint(mod1)

mod2 <- nls(wt~a*len^b,
            data=lth_wt,
            start=list(a=.01,b=3),
            control=nls.control(maxiter = 100),
            model=T)
residuals(mod2)
lth_wt <- na.omit(lth_wt)
boxplot(residuals(mod2)~lth_wt$year)
abline(h=0,lty=5,col=2)


preds <- data.frame(len=seq(0,900,5))
preds$wts <- coef(mod2)[1]*preds$len^coef(mod2)[2]

predict(mod2)
confint2(mod2)

plot(lth_wt$len,lth_wt$wt,col=lth_wt$year)
points(preds,lty=1,col=4,typ='l')


# setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
# png("rg_lth_freq.png", height = 5, width = 8, units = 'in', res=300)
boxplot(resid(mod2)~lth_wt$year,
        xlab='Year',ylab='Length (mm)',main='Red grouper - SEAMAP trawl',
        varwidth=T,las=1,lty=1)
# dev.off()

a_i <- data.frame(a=rep(NA,length(2010:2019)),
                  a_lcl=rep(NA,length(2010:2019)),
                  a_ucl=rep(NA,length(2010:2019)))
b_i <- data.frame(b=rep(NA,length(2010:2019)),
                  b_lcl=rep(NA,length(2010:2019)),
                  b_ucl=rep(NA,length(2010:2019)))
resids <- list()
setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
# png("rg_len_wt.png", height = 5, width = 10, units = 'in', res=300)
par(mfrow=c(2,5),mar=c(4,4,1.5,1),oma=c(0,0,2,0))
for(i in 2010:2019){
  temp <- lth_wt[which(lth_wt$year==i),]
  tmod <- try(nls(wt~a*len^b,
              data=temp,
              start=list(a=.01,b=3),
              control=nls.control(maxiter = 100),
              algorithm = "port"))
  if(class(tmod)!='try-error'){
    a_i[i-2009,1] <- coef(tmod)[1]
    a_i[i-2009,2:3] <- confint2(tmod)[c(1,3)]
    b_i[i-2009,1] <- coef(tmod)[2]
    b_i[i-2009,2:3] <- confint2(tmod)[c(2,4)] 
    
    preds <- data.frame(len=seq(0,900,5))
    preds$wts <- coef(tmod)[1]*preds$len^coef(tmod)[2]
    
    resids[[i-2009]] <- resid(tmod)
    
    plot(lth_wt$len,lth_wt$wt,col='gray80',pch=20,
         xlab='',ylab='')
    points(temp$len,temp$wt,col=2,pch=16)
    points(preds,lty=1,col=4,typ='l')
    mtext(i)
    if(i==2010){
      legend('topleft',c('All data','Year i','Model'),
             pch=c(16,16,NA),col=c('gray80',2,4),lty=c(NA,NA,1),
             bty='n',cex=.75)
      mtext('Red grouper - SEAMAP summer trawl',outer=T)
    }
    if(i==2010 | i==2015){
      mtext('Weight (kg)',2,line=2.5,cex=.75)
    }
    if(i>2014){
      mtext('Length (mm)',1,line=2.5,cex=.75)
    }
  }
}
# dev.off()

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png("rg_len_wt_parms.png", height = 6, width = 7, units = 'in', res=300)
par(mfrow=c(2,1),mar=c(4,5,1,1),oma=c(0,0,2,0))
plot(2010:2019,a_i$a,
     ylim=range(a_i),pch=16,
     xlab='',ylab='',
     las=1,typ='b',xaxt='n')
mtext('a',2,line=4)
polygon(c(2010:2019,rev(2010:2019)),
        c(a_i$a_lcl,rev(a_i$a_ucl)),col=alpha('darkorange',.2))
axis(1,2010:2019)
# arrows(2010:2019,a_i$a_lcl,2010:2019,a_i$a_ucl,code=3,angle=90,length=.05)
mtext(expression(paste('Red grouper: L ~ aW'^'b')),outer=T)

plot(2010:2019,b_i$b,
     ylim=range(b_i),pch=16,
     xlab='Year',ylab='',
     las=1,typ='b',xaxt='n')
mtext('b',2,line=4)
polygon(c(2010:2019,rev(2010:2019)),
        c(b_i$b_lcl,rev(b_i$b_ucl)),col=alpha(4,.2))
axis(1,2010:2019)
# arrows(2010:2019,b_i$b_lcl,2010:2019,b_i$b_ucl,code=3,angle=90,length=.05)
dev.off()

b <- boxplot(resids[[1]],xlim=c(1,10))
for(i in 2:10){
  boxplot(resids[[i]],add=T,at=i)
}

### ----------------- catch curve -----------------
## VBGM
# L_t ~ L_inf (1-exp(-k(t-t0)))
# t ~ (log(1-(L_t/L_inf))/-k)+t0

### from SEDAR12: Red Grouper (Epinephelus morio) age - length structure and description of growth from the eastern Gulf of Mexico: 1992-2001 
### Table 4
L_inf <- 923
k <- 0.11
t0 <- -3.21

age <- (log(1-(lth_wt$len/L_inf))/-k)+t0
age[which(age<0)] <- 0
age_r <- floor(age)
age_r2 <- round(age)

plot(age,age_r)
points(age,age_r2,col=2)
abline(v=1:30,lty=5,col='gray80')

age_agg <- aggregate(age_r,by=list(lth_wt$year,age_r),length)
names(age_agg) <- c('year','age','count')
age_agg$logcount <- log(age_agg$count)
# age_agg <- age_agg[-which(age_agg$age==0),]

plot(age_agg$age,age_agg$logcount)

z_i <- data.frame(z=rep(NA,length(2010:2019)),
                  z_lcl=rep(NA,length(2010:2019)),
                  z_ucl=rep(NA,length(2010:2019)))
for(i in 2010:2019){
  temp <- age_agg[which(age_agg$year==i),]
  mod <- lm(logcount~age,data=temp)
  z_i[i-2009,1] <- coef(mod)[2]
  z_i[i-2009,2:3] <- confint(mod)[c(2,4)]
  
  plot(temp$age,temp$logcount)
  mtext(i)
}


age_cl <- unique(age_agg$age)
plot(age_agg$year,age_agg$logcount,col='white')
for(i in 1:length(age_cl)){
  points(age_agg$year[which(age_agg$age==age_cl[i])],
         age_agg$logcount[which(age_agg$age==age_cl[i])],
         col=i,typ='b')  
}

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png("rg_mortality.png", height = 6, width = 7, units = 'in', res=300)
plot(2010:2019,-z_i$z,
     ylim=range(-z_i),pch=16,
     xlab='Year',ylab='Mortality (z)',
     las=1,typ='o',xaxt='n')
polygon(c(2010:2019,rev(2010:2019)),
        c(-z_i$z_lcl,rev(-z_i$z_ucl)),col=alpha(2,.2))
axis(1,2010:2019)
# arrows(2010:2019,-z_i$z_lcl,2010:2019,-z_i$z_ucl,code=3,angle=90,length=.05)
mtext(expression(paste('Red grouper: N'['t+1']*' ~ N'['t']*'e'^'-z')))
dev.off()


