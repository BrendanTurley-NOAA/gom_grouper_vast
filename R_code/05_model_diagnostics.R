### help with diagnostics
### https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

library(ape)
library(DHARMa)

load("~/Desktop/professional/projects/Postdoc_FL/data/grouper/2022-04-21_vmod8_results.RData")

plot(results$dharmaRes)
tq <- testQuantiles(results$dharmaRes)

testZeroInflation(results$dharmaRes)

testDispersion(results$dharmaRes)
testDispersion(results$dharmaRes,alternative='less')
### under-dispersion may be indicative of autocorrelation ot model misspecification

data[outliers(results$dharmaRes),]

residuals(results$dharmaRes)

plotResiduals(results$dharmaRes,form=results$dharmaRes$fittedPredictedResponse)

plotResiduals(results$dharmaRes,data$year)
plotResiduals(results$dharmaRes,data$VESSEL)
plotResiduals(results$dharmaRes,data$DEPTH_EMAX)
plotResiduals(results$dharmaRes,data$STAT_ZONE)
plotResiduals(results$dharmaRes,data$AreaSwept_km2)


ta_res <- recalculateResiduals(results$dharmaRes,fit$data_frame$t_i)
testTemporalAutocorrelation(ta_res$scaledResiduals,2010:2019)

# testSpatialAutocorrelation(results$dharmaRes,fit$data_frame$t_i)

## grouping by LAT/LON
## https://github.com/florianhartig/DHARMa/issues/297
data$g <- paste0(data$DECSLON,data$DECSLAT)
n <- recalculateResiduals(results$dharmaRes, group= data$g)
x <- aggregate(data$DECSLON, list(data$g), mean)$x
y <- aggregate(data$DECSLAT, list(data$g), mean)$x
testSpatialAutocorrelation(n, x = x, y = y,alternative = 'l')
### interpretation
# A negative z-value: data is clustered in a competitive way. 
# The spatial distribution of high values and low values in the dataset is more spatially dispersed than would be expected if underlying spatial processes were random. 
# For example, high values may be repelling high values or negative values may be repelling negative values.


par(mfrow=c(2,5))
for(i in 2010:2019){
tmp <- which(fit$data_frame$t_i==i)

dM = 1/as.matrix(dist(cbind(fit$data_frame$Lon_i[tmp], fit$data_frame$Lat_i[tmp])))
diag(dM) <- 0
dM[which(is.infinite(dM))] <- 0
m.i <- Moran.I(residuals(results$dharmaRes)[tmp],dM,na.rm=T)

plot(fit$data_frame$Lon_i[tmp],
     fit$data_frame$Lat_i[tmp],
     asp=1,col=ifelse(residuals(results$dharmaRes)[tmp]>=.5,2,4))
mtext(i,adj=0)
mtext(paste0('pval = ',signif(m.i$p.value,2)),adj=1)
}
