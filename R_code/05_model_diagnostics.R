library(DHARMa)

load("~/Desktop/professional/projects/Postdoc_FL/data/grouper/2022-04-06_vmod3_results.RData")

plot(vmod3$dharmaRes)

testDispersion(vmod3$dharmaRes)
testDispersion(vmod3$dharmaRes,alternative='less')
### underdispersion may be indicative of autocorrelation

data[outliers(vmod3$dharmaRes),]

residuals(vmod3$dharmaRes)

plotResiduals(vmod3$dharmaRes,form=vmod3$dharmaRes$fittedPredictedResponse)


ta_res <- recalculateResiduals(vmod3$dharmaRes,fit$data_frame$t_i)
testTemporalAutocorrelation(ta_res$scaledResiduals,2010:2019)

testSpatialAutocorrelation(vmod3$dharmaRes,fit$data_frame$t_i)

## grouping by LAT/LON
## https://github.com/florianhartig/DHARMa/issues/297
data$g <- paste0(data$DECSLON,data$DECSLAT)
n <- recalculateResiduals(vmod3$dharmaRes, group= data$g)
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
m.i <- Moran.I(residuals(vmod3$dharmaRes)[tmp],dM,na.rm=T)

plot(fit$data_frame$Lon_i[tmp],
     fit$data_frame$Lat_i[tmp],
     asp=1,col=ifelse(residuals(vmod3$dharmaRes)[tmp]>=.5,2,4))
mtext(i,adj=0)
mtext(paste0('pval = ',signif(m.i$p.value,2)),adj=1)
}
