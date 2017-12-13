### This script: 
## 1) loads the data matrix creates by "extract_and_summarize_EVI.R" script. The matrix contains EVI values for filtered pixels from 2000-2013
## 2) Summarizes these time series to get features to use to predict mortality in 2015-16. 


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

# load one mortality layer as a template
mort_template = raster("features/ADS-rasterized/Y2015_sp122.tif")
mort_2015_2016 = raster("features/ADS-rasterized/Y2015_spALL.tif") + raster("features/ADS-rasterized/Y2016_spALL.tif")

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

time_subset <- year(dates)>2000 & year(dates) < 2014

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
  yhat <- x[1] + x[2]*sin(seas.time+x[3]) + x[4]*cos(seas.time+x[5])
  return(sum((y-yhat)^2, na.rm=T))
}

y <- evi_mat[4000,time_subset]
seas.time <- ((1:n) %% 23)/23 * (2*pi)
n <- length(y)
fit1 <- optim(c(2, 1, 0.1), sin_fit)
fit1

# plot result
plot(1:n, y, cex=0.7)
lines(1:n, fit1$par[1] + fit1$par[2]*sin(seas.time+fit1$par[3]), col="cyan4", lwd=3)

# INLA version
y <- evi_mat[2000,time_subset]
n <- length(y)
d = data.frame(y=y, trend = 1:n, seasonal=1:n)
formula  = y ~  f(trend, model="rw2", cyclic=FALSE, param=c(1,0.0001)) + f(seasonal,model="seasonal",season.length=23,param=c(1,0.1))
mod = inla(formula, family="gaussian", data=d, control.family=list(param=c(4,4)), control.predictor = list(link = 1))
summary(mod)

plot(mod$summary.fitted.values$mean, pch=16, cex=0.7, col="blue")
points(y,  pch=16, cex=0.7, col="red")
lines(mod$summary.random$seasonal$mean + mod$summary.fixed$mean)
lines(mod$summary.random$trend$mean + mod$summary.fixed$mean, col="darkgray", lwd=2)




plot(mod, plot.fixed.effects=F, plot.random.effects=T, plot.lincomb = FALSE, plot.hyperparameters = FALSE, plot.predictor = FALSE, plot.q = FALSE, plot.cpo = FALSE)
### Make single-number summaries of EVI time series and store in data frame

# Summarize extracted values into data frame
evi_summary = data.frame(cell_number=as.integer(rownames(evi_mat)))
evi_summary$evi_mean = apply(evi_mat[,time_index], 1, mean, na.rm=T)
evi_summary$evi_mayjun = apply(evi_mat[,dates$mon %in% c(4, 5) & dates$year <= (endyear-1900) & dates$year >= (startyear-1900) ], 1, mean, na.rm=T)
evi_summary$evi_augsept = apply(evi_mat[,dates$mon %in% c(7,8) & dates$year <= (endyear-1900) & dates$year >= (startyear-1900) ], 1, mean, na.rm=T)
evi_summary$seas_change = evi_summary$evi_mayjun - evi_summary$evi_augsept
evi_summary$seas_change_prop = (evi_summary$evi_augsept/evi_summary$evi_mayjun)-1
evi_summary$total_var = apply(evi_mat[,time_index], 1, var, na.rm=T)

# wet-year vs dry-year difference in late-season EVI
# define wet years as 2000, 2005, 2006, 2010, 2011
# define dry years as 2002, 2007 ( could also include 2013 if that year's in the training data)
wetmean = apply(evi_mat[,dates$year %in% c(100,105, 106, 2010, 2011) & dates$mon %in% c(7,8)], 1, mean, na.rm=T)
drymean = apply(evi_mat[,dates$year %in% c(102,107, 113) & dates$mon %in% c(7,8)], 1, mean, na.rm=T) # Q whether to include 2013 here
evi_summary$wet_dry_diff = wetmean-drymean
evi_summary$wet_dry_propdiff = drymean/wetmean-1

# add trends
linear_time = scale(as.integer(dates[1:length(dates)]-dates[1]))
# Note this is slow with many pixels and should be parallelized! 
system.time(for (i in 1:nrow(evi_mat)) evi_summary$linear_trend[i] = coef(lm(evi_mat[i,time_index]~linear_time[time_index]))[2])

# partition variance
years = (startyear-1900):(endyear-1900)
ncells = nrow(evi_mat)
annual.mean = matrix(NA, nrow=ncells, ncol=length(years))
annual.var = matrix(NA, nrow=ncells, ncol=length(years))
for (i in 1:length(years)) {
    timeind = dates$year == years[i] & dates$mon %in% c(5,6,7,8,9)
  for (j in 1:ncells) {
    annual.mean[j,i] = mean(evi_mat[j,timeind], na.rm=T)
    annual.var[j,i] = var(evi_mat[j,timeind], na.rm=T)
  }
}

# within-year variance
evi_summary$within_year_var = apply(annual.var, 1, mean, na.rm=T)
evi_summary$within_year_sd = sqrt(evi_summary$within_year_var)

# among-year variance
evi_summary$among_year_var = apply(annual.mean, 1, var, na.rm=T)
evi_summary$among_year_sd = sqrt(evi_summary$among_year_var)

# ratio of among-to-within-year variance 
evi_summary$among_to_within_ratio = evi_summary$among_year_var / evi_summary$within_year_var




### Quick look 
sub = sample(1:nrow(evi_mat), size=500, replace=FALSE)
#pairs(evi_summary[sub,])

## Check variance outliers
# Note removing the disturbed pixels seems to have removed the very high variance outliers
# look at cells with high within-year variance 
outliers = which(evi_summary$within_year_var>0.01)
par(mfrow=c(5,5), mar=rep(2,4))
for (i in 1:25) plot(evi_mat[outliers[i],dates$mon %in% c(5,6,7,8,9)]~linear_time[dates$mon %in% c(5,6,7,8,9)], pch=16, cex=0.5, ylim=c(0, 1))
title("timeseries of EVI 2000-2016")
# same for high among-year variance 
outliers = which(evi_summary$among_year_var>0.002)
par(mfrow=c(5,5), mar=rep(2,4))
for (i in 1:25) plot(evi_mat[outliers[i],dates$mon %in% c(5,6,7,8,9)]~linear_time[dates$mon %in% c(5,6,7,8,9)], pch=16, cex=0.5, ylim=c(0, 1))
title("timeseries of EVI 2000-2016")
# these seem now to represent mainly cell with strong time-trends


## Add the mortality data 
mort_masked = mask(mort_albers, target_pixels, maskvalue=0)
evi_summary$mort = getValues(mort_masked)[as.integer(rownames(evi_mat))]
pairs(evi_summary[sub,])

x = as.matrix(cor(evi_summary[evi_summary$among_year_var<0.002 & evi_summary$within_year_var<0.005,], use="pairwise.complete"))
heatmap(x, col=viridis(12))
x

# Clean out NAs
evi_summary = evi_summary[complete.cases(evi_summary),]
dim(evi_summary)

# Store intermediate file 
save(evi_summary, file="features/working-files/evi_summary_PPN+SMC_jepson_central+south.Rdata")


############################################
# Do some simple regressions

load("features/working-files/evi_summary_PPN+SMC_jepson_central+south.Rdata")
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




