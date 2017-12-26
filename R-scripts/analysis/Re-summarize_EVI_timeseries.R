### This script: 
## 1) loads the EVI data matrix created by "extract_and_summarize_EVI.R" script. The matrix contains EVI values for filtered pixels from 2000-2013
## 2) Summarizes these time series to get features to use to predict mortality in 2015-16. 
## 3) loads the climate (and x,y coordinate) data matrix created by "extract_and_summarize_EVI.R"


library(sp)
library(raster)
library(rgdal)
library(lattice)
library(fields)
library(VGAM)
library(viridis)
library(car)
library(lubridate)
library(INLA)
library(parallel)


#### Load data ####

# Enter the locations of files to work with 
# geotifs are the MODIS EVI data that Mike K exported from Google Earth Engine, they are stored in a single folder (geotif_folder) and are all names consistently with the prefix geotif_filename and a date code. 
geotif_folder = "./features/ee-sn_jep_modis_ts_quality_mask_epsg3310/"
geotif_filename =  "sn_jep_modis_ts_quality_mask_epsg3310_"
filenames = dir(geotif_folder, pattern=geotif_filename)

# Projection info
albers.proj <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
wgs.proj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Get date information from the geotif filenames
date_codes = sapply(filenames, substr, start=nchar(geotif_filename)+5, stop=nchar(geotif_filename)+12) 
dates = strptime(date_codes, "%Y%m%d")

# Get the locations of the target pixels 
# EVI template raster 
evi_template = raster("features/sierra-nevada-250m-evi-template.tif")
# Get target pixels
load("features/target_pixels.Rdata")

# load one mortality layer as a template
mort_template = raster("features/ADS-rasterized/Y2015_sp122.tif")
mort_2015_2016 = raster("features/ADS-rasterized/Y2015_spALL.tif") + raster("features/ADS-rasterized/Y2016_spALL.tif")
# reproject to same projection as EVI data
mort_albers <- projectRaster(mort_2015_2016, evi_template)

# Load EVI matrix
load("features/working-files/evi_data_matrix_jepson_PPN+SMC_central+south.Rdata")


#### Summarize the EVI time series ####

# how many "good" values are there per month? 
obs_by_mon = rep(NA, 12)
for (i in 0:11) obs_by_mon[i+1] = sum(!is.na(evi_mat[,dates$mon==i]))
barplot(obs_by_mon, names.arg=as.character(1:12))
# all values present June-Sept, good number in May. April and Oct missing ~40% of values

# In previous summary version, we masked out all EVI values from months other than May-Sept. This time, keep all the information in. 

# Seasonality -- characterize seasonal fluctuation and its amplitude (and whether it changes through time?)

# Try fitting an ARIMA model 
# plot some example time series
n_pixels <- nrow(evi_mat)
n_times <- ncol(evi_mat)
par(mfrow=c(4, 4), mar=rep(2, 4))
for (i in 1:16) plot(evi_mat[i*1000,], pch=16, cex=0.5, col="slateblue")

# how many observations per year? 
table(year(dates))

time_subset <- year(dates)>1999 & year(dates) < 2013

m <- arima(evi_mat[3000,time_subset], order = c(1, 0, 0), seasonal = list(order=c(0,1,1), period=23))
coef(m)
tsdiag(m)
AIC(m)

# Fit a sine wave to data, given that the period is 23 obs/year
sin_fit <- function(x) {
  yhat <- x[1] + x[2]*sin(seas.time+x[3]) 
  return(sum((y-yhat)^2, na.rm=T))
}
wave_fit <- function(x) {
  yhat <- x[1] + x[2]*sin(seas.time) + x[3]*cos(seas.time)
  return(sum((y-yhat)^2, na.rm=T))
}

y <- evi_mat[1000,time_subset]
seas.time <- yday(dates[time_subset])/365 * 2*pi #((1:n) %% 23)/23 * (2*pi)
n <- length(y)
fit1 <- optim(c(0.3, 0.5, 1), sin_fit, method="BFGS")
fit1
# note that it tends to slide the sine wave back by 1/4 year, making jan 1 the trough and jul 1 the peak, which is what we'd expect

# plot result
plot(1:n, y, cex=0.7)
lines(1:n, fit1$par[1] + fit1$par[2]*sin(seas.time+fit1$par[3]), col="cyan4", lwd=3)



# INLA version
y <- evi_mat[3000, time_subset]
n <- length(y)
d = data.frame(y=y, time=1:n)
d$trend = scale((1:n)/230, center=T, scale=F) # make trend into a rate per decade, centered at 0
d$month = month(dates[time_subset]) 
d$doy = yday(dates[time_subset])
d$sinwave = sin((d$doy-91)/365 * 2*pi) # slide sine wave back 1/4 year

# trend as random effect
#formula  = y ~  f(trend, model="ar1", cyclic=FALSE, param=c(1,0.0001)) + f(month, model="seasonal", season.length=12, param=c(1,0.001))
# trend as fixed effect
formula  = y ~  trend + sinwave + f(month, model="seasonal", season.length=12, param=c(1,0.0001)) #+ f(time, model="rw1", param=c(100, 0.0001))

mod = inla(formula, family="gaussian", data=d)
summary(mod)

plot(mod, plot.fixed.effects=F, plot.random.effects=T, plot.lincomb = FALSE, plot.hyperparameters = FALSE, plot.predictor = FALSE, plot.q = FALSE, plot.cpo = FALSE)

#plot(y~month, d)
#lines(d$month, d$sinwave*mod$summary.fixed$mean[3] + mod$summary.fixed$mean[1], col="blue" )

# amplitude as difference of min and max REs
#max(mod$summary.random$month$mean) - min(mod$summary.random$month$mean)

# amplitude as effect of sine wave * 2 (effect at max minus min value of the sine wave)
mod$summary.fixed$mean[3]*2

# level (intercept, which is mean EVI at middle of the time series and at the mean seasonality value)
mod$summary.fixed[1,]

# trend (change in EVI per decade)
mod$summary.fixed[2,]



#### Apply INLA model to all pixels and extract these coef values 
n.cells <- nrow(evi_mat)
seas.amp <- seas.amp.025 <- seas.amp.975  <- mean.evi <- mean.evi.025 <- mean.evi.975 <- trend.evi <- trend.evi.025 <- trend.evi.975 <- rep(NA, n.cells)
time_subset <- year(dates) < 2013

# specify model
formula  = y ~  trend + sinwave + f(month, model="seasonal", season.length=12, param=c(1,0.0001)) #+ f(time, model="rw1", param=c(100, 0.0001))
  
# function to fit inla model in to one cell 
inla.ts.fit <- function(y, ts_dates, time_scalar, time_shift, formula) {
  require(lubridate)
  require(INLA)
  n <- length(y)
  d = data.frame(y=y, time=1:n)
  d$trend = scale((1:n)/time_scalar, center=T, scale=F) # make trend into a rate per decade, centered at 0
  d$month = month(ts_dates) 
  d$doy = yday(ts_dates)
  d$sinwave = sin((d$doy-time_shift)/365 * 2*pi) # slide sine wave back 1/4 year
  mod = inla(formula, family="gaussian", data=d)
  seas.amp <- mod$summary.fixed$mean[3]
  seas.amp.025 <- mod$summary.fixed$`0.025quant`[3]
  seas.amp.975 <- mod$summary.fixed$`0.975quant`[3]
  mean.evi <- mod$summary.fixed$mean[1] 
  mean.evi.025 <- mod$summary.fixed$`0.025quant`[1]
  mean.evi.975 <- mod$summary.fixed$`0.975quant`[1]
  trend.evi <- mod$summary.fixed$mean[2]
  trend.evi.025 <- mod$summary.fixed$`0.025quant`[2]
  trend.evi.975 <- mod$summary.fixed$`0.975quant`[2]
  return(c(seas.amp, seas.amp.025, seas.amp.975, mean.evi, mean.evi.025, mean.evi.975, trend.evi, trend.evi.025, trend.evi.975))
}

# Test on a few cells 
system.time(apply(evi_mat[1:14,time_subset], 1, inla.ts.fit, ts_dates=dates[time_subset], time_scalar=230, time_shift=91, formula=formula))

# Test in parallel 
no.cores <- detectCores() - 1
cl <- makeCluster(no.cores) # for Windows
cl <- makeCluster(no.cores, type="FORK")
system.time(parRapply(cl=cl, evi_mat[1:14, time_subset], inla.ts.fit, ts_dates=dates[time_subset], time_scalar=230, time_shift=91, formul=formula))
# runs 3x as fast, whole analysis should take ~3.75 hours  
# note it returns all the values in a long vector that's n.cells * 9 long

# Run in parallel
fit.vec <- parRapply(cl=cl, evi_mat[, time_subset], inla.ts.fit, ts_dates=dates[time_subset], time_scalar=230, time_shift=91, formula=formula)

# Reformat output 
evi_summary <- matrix(fit.vec, nrow=n.cells, byrow=TRUE)
evi_summary <- as.data.frame(evi_summary)
rownames(evi_summary) = rownames(evi_mat)
names(evi_summary) = c("seas.amp", "seas.amp.025", "seas.amp.975", "mean.evi", "mean.evi.025", "mean.evi.975", "trend.evi", "trend.evi.025", "trend.evi.975")


## Add the mortality data 
mort_masked = mask(mort_albers, target_pixels, maskvalue=0)
evi_summary$mort = getValues(mort_masked)[as.integer(rownames(evi_summary))]

dim(evi_summary)

# Store intermediate file 
save(evi_summary, file="features/working-files/evi_summary_new_PPN+SMC_jepson_central+south.Rdata")


### AML STOPPED HERE 12/22

### DY Resumed here 12/26 ###

# load climate data (clim_mat), which also contains cell x,y coordinates (Albers I believe)
load("features/working-files/climate_data_matrix_jepson_PPN+SMC_central+south.Rdata")

# the rows of climate data (which also contains cell x,y coordinates) should match the rows of the evi_summary matrix
# you could bind them by row number (simply cbind them) or by row names (which should also match)
# I (Derek) did not do this yet because I could not get INLA to finish running for all cells (and thus generate evi_summary)




############################################
# Do some simple regressions

load("features/working-files/evi_summary_new_PPN+SMC_jepson_central+south.Rdata")
evi_template = raster("features/sierra-nevada-250m-evi-template.tif")


### Run a simple model to check associations -- use tobit model in vgam library
hist(evi_summary$mort)
cols_to_standardize = c("evi_mean", "seas_change", "within_year_sd", "among_year_sd", "wet_dry_diff", "linear_trend")
for (i in 1:length(cols_to_standardize)) evi_summary[,cols_to_standardize[i]] = scale(evi_summary[,cols_to_standardize[i]])

# check for correlation in explanatory variables
vif(lm(mort~evi_mean+seas_change+within_year_sd + among_year_sd+wet_dry_diff+linear_trend, data=evi_summary))
cor(evi_summary[,cols_to_standardize], use="pairwise.complete")

# fit model
m = vglm(mort~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, tobit, data=evi_summary, trace=TRUE)
summary(m)

barplot(coef(m)[3:8], horiz=T, las=2, main="tobit model coefficients", col=ifelse(coef(m)[3:8]<0, "red", "blue"), cex.names=0.5)

plot(evi_summary$mort[!is.na(evi_summary$mort)]~predict(m, type="response"))
abline(0,1)


#####################
# Make output plots

plot_to_region <- function(cell.values, cell.index, crop_layer) { # index is the row numbers of the cells to plot, and indexes grid cells in the original evi_template
  r_tmp = evi_template
  r_tmp = setValues(r_tmp, rep(NA, length(r_tmp)))
  r_tmp[cell.index] = cell.values
  r_plot = crop(r_tmp, crop_layer)
  plot(r_plot,col=tim.colors(16))#viridis(10))
}


###########################
# Visualization

# plot some EVI summaries
# make a template for plotting 
load("features/target_pixels.Rdata")
target_pixels_na = target_pixels
target_pixels_na[target_pixels_na==0] = NA
plot_to_region(evi_summary$evi_mean, evi_summary$cell_number,subset_layer_albers); title("mean EVI")

plot_to_region(evi_summary$seas_change, evi_summary$cell_number, subset_layer_albers); title("Early-season EVI minus late-season EVI", cex.main=0.7)
plot_to_region(sqrt(evi_summary$among_year_var-min(evi_summary$among_year_var)), evi_summary$cell_number, subset_layer_albers); title("among-year sd")
plot_to_region(evi_summary$wet_dry_diff, evi_summary$cell_number, subset_layer_albers); title("Wet-year mean EVI minus dry-year mean EVI", cex.main=0.7)

# observed and predicted mortality 
mort_pred = fitted(m)
# As a quick fix for visualization, truncate mortality prediction / fit at 0
mort_pred[mort_pred<0] = 0
par(mfrow=c(1,2))
plot_to_region(sqrt(evi_summary$mort[!is.na(evi_summary$mort)]), evi_summary$cell_number[!is.na(evi_summary$mort)], subset_layer_albers)
title("Sqrt observed mortality")
plot_to_region(sqrt(mort_pred)[!is.na(evi_summary$mort)], evi_summary$cell_number[!is.na(evi_summary$mort)], subset_layer_albers)
title("Sqrt model fit")


# look at random individual pixels
par(mfrow=c(4,4), mar=rep(2, 4))
for (i in 1:16) plot(evi_mat[sample(1:nrow(evi_mat), 1),dates$mon %in% c(5,6,7,8,9)]~linear_time[dates$mon %in% c(5,6,7,8,9)], pch=16, cex=0.4, ylim=c(0, 0.9))

# all pixels averaged
evi_mean_all = apply(evi_mat, 2, mean, na.rm=T)
# long time series
plot(evi_mean_all, ylim=c(0.2, 0.5), type="l",lwd=2, col="cyan4")
# I'd say this shows that 2013 was low, clearly a drought year, but not out of the normal range for the rest of the years. So for model fitting, seems ok to go through 2013. The later years are drastically low, esp 2016. Will be interesting to see the rebound in 2017, if any. 
# plotting the spatial average for each year. 
plot(evi_mean_all[dates$year==100]~dates$yday[dates$year==100], type="l", lwd=2, col="cyan4", ylim=c(0.2, 0.5))
for (i in 101:112) lines(evi_mean_all[dates$year==i]~dates$yday[dates$year==i], type="l", lwd=2, col="cyan4")
for (i in 113:116) lines(evi_mean_all[dates$year==i]~dates$yday[dates$year==i ], type="l", lwd=2, col="orange3")

m_lin = lm(sqrt(mort)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, data=evi_summary[!is.na(evi_summary$mort),])
summary(m_lin)

mort_pred = predict(m_lin)
mort_pred[mort_pred<0]= 0
par(mfrow=c(1,2))
plot_to_region(sqrt(evi_summary$mort), evi_summary$cell_number, subset_layer_albers)
title("Sqrt observed mortality")
# stupid hack to get the same color scale in both plots
mort_pred[1] = sqrt(max(evi_summary$mort))
plot_to_region(mort_pred[!is.na(evi_summary$mort)], evi_summary$cell_number[!is.na(evi_summary$mort)],subset_layer_albers)
title("Sqrt model fit")

#### Cross-validation of simple linear models ##### 

# split data into fitting (90%) and test (10%) data sets
n <- nrow(evi_summary)
holdout.ind <- rep(TRUE, n)
holdout.ind[sample(1:nrow(evi_summary), size=round(0.1*n))] = FALSE
evi_fit <- evi_summary[holdout.ind,]
evi_test <- evi_summary[!holdout.ind,]

# simple cv of tobit model
m <- vglm(mort~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, tobit, data=evi_fit, trace=TRUE)
summary(m)

m.pred <- predict(m, newdata=evi_test)[,1]
plot(evi_test$mort ~ m.pred)
cor(evi_test$mort, m.pred) # predictive R2=0.11

# simple cv of linear model
m_lin <- lm(sqrt(mort)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, data=evi_fit)
summary(m_lin)
m_lin.pred <- predict(m_lin, newdata=evi_test)
plot(m.pred ~ I(m_lin.pred^2))
plot(evi_test$mort ~ I(m_lin.pred^2))
cor(evi_test$mort, m_lin.pred) # predictive R2=0.11

# Same for gam 
gam.formula <- sqrt(mort) ~ s(seas_change)+s(wet_dry_diff)+s(within_year_sd) + s(among_year_sd)+s(linear_trend) #s(evi_mean)
m_gam <- gam(gam.formula, data=evi_fit)
summary(m_gam)

qqnorm(resid(m_gam))

m_gam.pred <- predict(m_gam, newdata=evi_test)
plot(evi_test$mort ~ m_gam.pred)
cor(evi_test$mort, m_gam.pred)^2 # predictive R2 = 0.13, but without mean evi, it's only 0.07




