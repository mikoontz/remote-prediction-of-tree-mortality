### This script: 
## 1) creates a raster stack from the 16-day EVI files exported by Mike
## 2) pulls the data values out of these to create a data matrix.
## 3) generates summary statistics for the time series of EVI values for each pixel from 2000-2013.
## 4) Does simple correlation and regression analysis to test which of these summary stats predicts mortality in 2015-16.
## 5) displays summary stats and model predictions as a plotted raster.


library(sp)
library(raster)
library(rgdal)
library(lattice)
library(fields)
library(VGAM)
library(viridis)
library(car)
library(spdep)
library(mgcv)

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

# Load selected mortality response variables
mort_2011 = raster("features/ADS-rasterized/Y2011_sp122.tif") 
mort_2012 = raster("features/ADS-rasterized/Y2012_sp122.tif") 
mort_2013 = raster("features/ADS-rasterized/Y2013_sp122.tif") 
mort_2014 = raster("features/ADS-rasterized/Y2014_sp122.tif") 
mort_2015 = raster("features/ADS-rasterized/Y2015_sp122.tif") 
mort_2016 = raster("features/ADS-rasterized/Y2016_sp122.tif") 
sum(getValues(mort_2011)>0, na.rm=T); sum(getValues(mort_2012)>0, na.rm=T);sum(getValues(mort_2013)>0, na.rm=T);sum(getValues(mort_2014)>0, na.rm=T);sum(getValues(mort_2015)>0, na.rm=T);sum(getValues(mort_2016)>0, na.rm=T)

plot(mort_2013)
plot(mort_2014)
plot(mort_2015)
plot(mort_2016)

# Create a target cover layer for pixels with specified percent forest type cover
# Start by focusing only on PPN for test run
cover_min = 80 # Minumum cover of WHR type
target_cover = raster("features/calveg-pct-cover-rasters/sierra_nevada_250m_calveg_cover_whr_type_PPN.tif") + raster("features/calveg-pct-cover-rasters/sierra_nevada_250m_calveg_cover_whr_type_SMC.tif")

target_pixels = target_cover
target_pixels[target_cover<80] = 0
target_pixels[target_cover>=80] = 1
target_pixels[is.na(target_pixels)] = 0 

# Optionally further geographically subset the target pixels
subset_layer = shapefile("features/jepson-central+southern-outline.shp")
subset_layer_albers = spTransform(subset_layer, albers.proj)
target_pixels = mask(target_pixels, subset_layer_albers, updatevalue=0)

# Reproject mortality raster to match the EVI geotiffs (in Albers projection). 
# Note target veg raster is already in this projection. 
mort_2011_albers = projectRaster(mort_2011, evi_template)
mort_2012_albers = projectRaster(mort_2012, evi_template)
mort_2013_albers = projectRaster(mort_2013, evi_template)
mort_2014_albers = projectRaster(mort_2014, evi_template)
mort_2015_albers = projectRaster(mort_2015, evi_template)
mort_2016_albers = projectRaster(mort_2016, evi_template)

# Remover numerical artefacts from projection
mort_2011_albers[which(getValues(mort_2011_albers)<0)] = 0
mort_2012_albers[which(getValues(mort_2012_albers)<0)] = 0
mort_2013_albers[which(getValues(mort_2013_albers)<0)] = 0
mort_2014_albers[which(getValues(mort_2014_albers)<0)] = 0
mort_2015_albers[which(getValues(mort_2015_albers)<0)] = 0
mort_2016_albers[which(getValues(mort_2016_albers)<0)] = 0



######################################################
# Explore summaries of EVI time series

load("features/working-files/evi_summary_PPN+SMC_jepson_central+south.Rdata")

## Add the mortality data 
mort_2011_masked = mask(mort_2011_albers, target_pixels, maskvalue=0)
evi_summary$mort_2011 = getValues(mort_2011_masked)[evi_summary$cell_number]
mort_2012_masked = mask(mort_2012_albers, target_pixels, maskvalue=0)
evi_summary$mort_2012 = getValues(mort_2012_masked)[evi_summary$cell_number]
mort_2013_masked = mask(mort_2013_albers, target_pixels, maskvalue=0)
evi_summary$mort_2013 = getValues(mort_2013_masked)[evi_summary$cell_number]
mort_2014_masked = mask(mort_2014_albers, target_pixels, maskvalue=0)
evi_summary$mort_2014 = getValues(mort_2014_masked)[evi_summary$cell_number]
mort_2015_masked = mask(mort_2015_albers, target_pixels, maskvalue=0)
evi_summary$mort_2015 = getValues(mort_2015_masked)[evi_summary$cell_number]
mort_2016_masked = mask(mort_2016_albers, target_pixels, maskvalue=0)
evi_summary$mort_2016 = getValues(mort_2016_masked)[evi_summary$cell_number]

# Clean out NAs
evi_summary = evi_summary[complete.cases(evi_summary),]
dim(evi_summary)

# create a variable that describes average mortality in a grid cell's neighbors
maxdist = 740 # set size of neighborhood -- started with 740 for small neighborhood
              # results below don't seem different when this is increased to 1480
evi_points = rasterToPoints(evi_template)
summary_points = evi_points[evi_summary$cell_number,1:2]
evi_neigh = dnearneigh(summary_points, d1=100, d2=740, longlat=FALSE)
mort_neigh = lapply(evi_neigh, f<- function(x){return(mean(evi_summary$mort_2011[x], na.rm=T))})
evi_summary$mort_neigh_2011 = unlist(mort_neigh)
mort_neigh = lapply(evi_neigh, f<- function(x){return(mean(evi_summary$mort_2012[x], na.rm=T))})
evi_summary$mort_neigh_2012 = unlist(mort_neigh)
mort_neigh = lapply(evi_neigh, f<- function(x){return(mean(evi_summary$mort_2013[x], na.rm=T))})
evi_summary$mort_neigh_2013 = unlist(mort_neigh)
mort_neigh = lapply(evi_neigh, f<- function(x){return(mean(evi_summary$mort_2014[x], na.rm=T))})
evi_summary$mort_neigh_2014 = unlist(mort_neigh)
mort_neigh = lapply(evi_neigh, f<- function(x){return(mean(evi_summary$mort_2015[x], na.rm=T))})
evi_summary$mort_neigh_2015 = unlist(mort_neigh)

# Store intermediate file 
save(evi_summary, file="features/working-files/evi_summary_with_mort_PPN+SMC_jepson_central+south.Rdata")


############################################
# Do some simple regressions for different years

### standardize x matrix
cols_to_standardize = c("evi_mean", "seas_change", "within_year_sd", "among_year_sd", "wet_dry_diff", "linear_trend")
for (i in 1:length(cols_to_standardize)) evi_summary[,cols_to_standardize[i]] = scale(evi_summary[,cols_to_standardize[i]])

# Weird outliers for wet_dry_diff -- toss these
evi_summary$wet_dry_diff[evi_summary$wet_dry_diff < -5] = NA
evi_summary = evi_summary[complete.cases(evi_summary),]


### Check out shapes of responses
# Combine the 2 major mortality years
evi_summary$mort_2015_16 = evi_summary$mort_2015 + evi_summary$mort_2016
# Combine 2014-16
evi_summary$mort_2014_16 = evi_summary$mort_2015+evi_summary$mort_2015 + evi_summary$mort_2016


# fit gam model to look at response shapes
m1 = gam(mort_2015_16~s(evi_mean)+s(seas_change)+wet_dry_diff+within_year_sd + s(among_year_sd)+s(linear_trend), data=evi_summary, trace=TRUE)
summary(m1)
plot(m)
 
# does adding lag help?
m2 = gam(mort_2015_16~s(evi_mean)+s(seas_change)+wet_dry_diff+within_year_sd + s(among_year_sd)+s(linear_trend) + s(mort_2014), data=evi_summary, trace=TRUE)
summary(m2)
# not much 

# does adding lagged local neighborhood help?
m3 = gam(mort_2015_16~s(evi_mean)+s(seas_change)+wet_dry_diff+within_year_sd + s(among_year_sd)+s(linear_trend) +s(mort_neigh_2014), data=evi_summary, trace=TRUE)
summary(m3) 
plot(m3) # a little, at very high mortality levels 

### Convert this over to a tobit model for 2014
m = vglm(mort_2014~evi_mean+I(evi_mean^2) + seas_change + I(seas_change^2)+wet_dry_diff+I(wet_dry_diff^2)+within_year_sd+linear_trend + I(linear_trend^2)+mort_neigh_2013, tobit, data=evi_summary, trace=TRUE)
#m = vglm(mort_2014~evi_mean+I(evi_mean^2) + seas_change +wet_dry_diff+within_year_sd +mort_neigh_2013, tobit, data=evi_summary, trace=TRUE)
BIC(m) # full model plus neighborhood mortality has best BIC 
summary(m)
m0 = vglm(mort_2014~1, tobit, data=evi_summary, trace=TRUE)
1 - (-2*logLik(m)) / (-2*logLik(m0))


# PLot model coefficients
barplot(coef(m)[3:length(coef(m))], horiz=T, las=2, main="Tobit model coefficients", col=ifelse(coef(m)[3:length(coef(m))]<0, "red", "blue"), cex.names=0.5)


### Convert this over to a tobit model  for 2015-16
m = vglm(mort_2015_16~evi_mean+I(evi_mean^2) + seas_change + I(seas_change^2)+wet_dry_diff+I(wet_dry_diff^2)+within_year_sd +linear_trend + I(linear_trend^2)+mort_neigh_2014  + mort_2014, tobit, data=evi_summary, trace=TRUE)
BIC(m) # full model plus neighborhood mortality has best BIC 
summary(m)
m0 = vglm(mort_2015_16~1, tobit, data=evi_summary, trace=TRUE)
1 - (-2*logLik(m)) / (-2*logLik(m0))
barplot(coef(m)[3:length(coef(m))], horiz=T, las=2, main="Tobit model coefficients", col=ifelse(coef(m)[3:length(coef(m))]<0, "red", "blue"), cex.names=0.5)


# Same for all 3 years
m = vglm(mort_2014_16~evi_mean+I(evi_mean^2) + seas_change + I(seas_change^2)+wet_dry_diff+I(wet_dry_diff^2)+within_year_sd +linear_trend + I(linear_trend^2)+mort_neigh_2013, tobit, data=evi_summary, trace=TRUE)
BIC(m) # full model plus neighborhood mortality has best BIC 
summary(m)
m0 = vglm(mort_2014_16~1, tobit, data=evi_summary, trace=TRUE)
1 - (-2*logLik(m)) / (-2*logLik(m0))
barplot(coef(m)[3:length(coef(m))], horiz=T, las=2, main="Tobit model coefficients", col=ifelse(coef(m)[3:length(coef(m))]<0, "red", "blue"), cex.names=0.5)


# Note wet_dry_diff matters in 2014, 2015 but not 2016

# look at predictions 
# Scatterplot
mort_pred = predict(m, type="response")
mort_pred[mort_pred<0] = 0
plot(evi_summary$mort_2014~mort_pred)
abline(0,1)

# map
plot_to_region <- function(cell.values, cell.index, crop_layer) { # index is the row numbers of the cells to plot, and indexes grid cells in the original evi_template
  r_tmp = evi_template
  r_tmp = setValues(r_tmp, rep(NA, length(r_tmp)))
  r_tmp[cell.index] = cell.values
  r_plot = crop(r_tmp, crop_layer)
  plot(r_plot,col= tim.colors(16))#rev(viridis(16)))# #
}
par(mfrow=c(1,2))
plot_to_region(sqrt(evi_summary$mort_2014), evi_summary$cell_number, subset_layer_albers)
title("Sqrt observed mortality")
# stupid hack to get the same color scale in both plots
mort_pred[1] = max(evi_summary$mort_2014)
plot_to_region(sqrt(mort_pred), evi_summary$cell_number,subset_layer_albers)
title("Sqrt model fit")

plot_to_region(mort_pred-evi_summary$mort_2014, evi_summary$cell_number, subset_layer_albers)

###########
# For presentation perhaps stop here with modeling. 
# Then talk about implications by showing maps of the "important" explanatory variables and what they might indicate. 
# Compare "early" drought response in 2014 (pretty predictable, relates to past sensitivity of forest patches to dry years, and to seasonal dry-down). 
# To the later "extreme" drought mortality of 2015-16 (not very predictable except in broad terms, possibly due to beetle population dynamics, but also possibly to finer-scale factors that aren't in this analysis)
# About beetle dynamics: ask if adding in a forest heterogeneity layer (from Mike) improves prediction. The expectationis that it does, because heterogeneity at some intermediate scale should moderate impacts of beetles as well as fire. Though possibly there is an exception in extreme conditions whereby the fire generates its own weather/wind and where the beetles switch to different attack behavior. 




### Convert this over to a hurdle model

# Define mortality year(s) and threshold
evi_summary$mort = evi_summary$mort_2015+evi_summary$mort_2016
evi_summary$mort_pa = as.integer((evi_summary$mort)>0.1) # set threshold: 1 is 1 tree per year per acre
sum(evi_summary$mort_pa)/nrow(evi_summary)

# Presence or absence of mortality
m.pa = glm(mort_pa~evi_mean+I(evi_mean^2) + seas_change + I(seas_change^2)+wet_dry_diff+I(wet_dry_diff^2)+within_year_sd + among_year_sd + I(among_year_sd^2)+linear_trend + I(linear_trend^2), data=evi_summary, family="binomial")
summary(m.pa)
# pct deviance explained
m0.pa = glm(mort_pa~1, data=evi_summary, family="binomial")
1 - (-2*logLik(m.pa)[1]) / (-2*logLik(m0.pa)[1]) # a lot -- 24%
# ROC AUC
colAUC(fitted(m.pa), evi_summary$mort_pa)

# plot coefs
barplot(coef(m.pa), horiz=T, las=2, main="binomial model coefficients", col=ifelse(coef(m.pa)<0, "red", "blue"), cex.names=0.5)

# Amount of mortality, given mortality occurred
m.dens = lm(sqrt(mort)~evi_mean+I(evi_mean^2) + seas_change + I(seas_change^2)+wet_dry_diff+I(wet_dry_diff^2)+within_year_sd + among_year_sd + I(among_year_sd^2)+linear_trend + I(linear_trend^2), data=evi_summary, subset=evi_summary$mort_pa==1)
summary(m.dens)
# R2 is pretty weak -- 0.11

# plot coefs
barplot(coef(m.dens)[2:16], horiz=T, las=2, main="linear model coefficients", col=ifelse(coef(m.dens)[2:16]<0, "red", "blue"), cex.names=0.5)


# map of presence absence
par(mfrow=c(1,2))
plot_to_region(evi_summary$mort_pa, evi_summary$cell_number, subset_layer_albers)
title("observed mortality (pres/abs)")
plot_to_region(fitted(m.pa), evi_summary$cell_number,subset_layer_albers)
title("model fit (prob. mortality)")

# map of amount of mortality, given present
par(mfrow=c(1,2)) 
plot_to_region(sqrt(evi_summary$mort_2015_16), evi_summary$cell_number, subset_layer_albers)
title("sqrt observed mortality (TPA)")
mort_pred  = fitted(m.dens)
#mort_pred[1] = max(sqrt(evi_summary$mort_2015_16))
plot_to_region(mort_pred, evi_summary$cell_number,subset_layer_albers)
title("model fit")


plot(fitted(m.pa),jitter(evi_summary$mort_pa), pch=".") # not too bad at predicting where mortality occurred    

plot(mort_pred[evi_summary$mort_2015_16>0]~sqrt(evi_summary$mort_2015_16[evi_summary$mort_2015_16>0]))
# Terrible at predicting how much, given it occurred







# fit model
m = vglm(mort_2011~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, tobit, data=evi_summary, trace=TRUE)
summary(m)
m = vglm(mort_2012~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, tobit, data=evi_summary, trace=TRUE)
summary(m)
m = vglm(mort_2013~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, tobit, data=evi_summary, trace=TRUE)
summary(m)
m = vglm(mort_2014~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, tobit, data=evi_summary, trace=TRUE)
summary(m)
m = vglm(mort_2015~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, tobit, data=evi_summary, trace=TRUE)
summary(m)
m = vglm(mort_2016~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, tobit, data=evi_summary, trace=TRUE)
summary(m)

# Test for neighborhood effect from 2013
m = vglm(mort_2014~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, tobit, data=evi_summary, trace=TRUE)
summary(m)
m = vglm(mort_2014~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend+sqrt(mort_neigh), tobit, data=evi_summary, trace=TRUE)
summary(m) # vastly better with neighborhood effect
m = vglm(mort_2015~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend+sqrt(mort_neigh), tobit, data=evi_summary, trace=TRUE)
summary(m) # weak neighborhood effect
m = vglm(mort_2016~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend+sqrt(mort_neigh), tobit, data=evi_summary, trace=TRUE)
summary(m) # weirdly neighborhood effect is strong again

barplot(coef(m)[3:15], horiz=T, las=2, main="tobit model coefficients", col=ifelse(coef(m)[3:8]<0, "red", "blue"), cex.names=0.5)

plot(evi_summary$mort[!is.na(evi_summary$mort)]~predict(m, type="response"))
abline(0,1)

# dumb linear version
m_lin = lm(sqrt(mort_2016)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend , data=evi_summary)
summary(m_lin)
m_lin = lm(sqrt(mort_2015)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend , data=evi_summary)
summary(m_lin)
m_lin = lm(sqrt(mort_2014)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend , data=evi_summary)
summary(m_lin)
m_lin = lm(sqrt(mort_2013)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend , data=evi_summary)
summary(m_lin)
m_lin = lm(sqrt(mort_2012)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend , data=evi_summary)
summary(m_lin)
m_lin = lm(sqrt(mort_2011)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend , data=evi_summary)
summary(m_lin)

# Linear version with neighborhood effect
m_lin = lm(sqrt(mort_2014)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend+sqrt(mort_neigh), data=evi_summary)
summary(m_lin) # neighborhood effect is large, but doesn't explain much variation
m_lin = lm(sqrt(mort_2015)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend +sqrt(mort_neigh), data=evi_summary)
summary(m_lin) # no neighborhood effect
m_lin = lm(sqrt(mort_2016)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend +sqrt(mort_neigh), data=evi_summary)
summary(m_lin) # neighborhood effect signif but not very explanatory, mostly it's climate





###############################################################
# Hurdle model for presence/absence of mortality, then for amount where it occurs

# Define mortality year(s) and threshold
evi_summary$mort = evi_summary$mort_2015+evi_summary$mort_2016
evi_summary$mort_pa = as.integer((evi_summary$mort)>1) # set threshold: 1 is 1 tree per year per acre
sum(evi_summary$mort_pa)/nrow(evi_summary)

# standardize variables if not yet standardized
if (!exists("cols_to_standardize")) {
  cols_to_standardize = c("evi_mean", "seas_change", "within_year_sd", "among_year_sd", "wet_dry_diff", "linear_trend")
  for (i in 1:length(cols_to_standardize)) evi_summary[,cols_to_standardize[i]] = scale(evi_summary[,cols_to_standardize[i]])
}

# Presence or absence of mortality
m.pa = glm(mort_pa~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend+mort_neigh+mort_2014, data=evi_summary, family="binomial", link="probit")
summary(m.pa)
# probit version
m.dens = glm(mort_pa~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend+mort_2014, data=evi_summary, family=binomial(link="logit"))
summary(m.dens)
colAUC(fitted(m.dens), evi_summary$mort_pa)

# Amount of mortality, given mortality occurred
m.dens = lm(sqrt(mort)~evi_mean+seas_change+wet_dry_diff+within_year_sd + among_year_sd+linear_trend, data=evi_summary, subset=evi_summary$mort_pa==1)
summary(m.dens)







#####################
# Make output plots

plot_to_region <- function(cell.values, cell.index, crop_layer) { # index is the row numbers of the cells to plot, and indexes grid cells in the original evi_template
  r_tmp = evi_template
  plot(subset_layer_albers, col=gray(0.9), lty=0)
  r_tmp = setValues(r_tmp, rep(NA, length(r_tmp)))
  r_tmp[cell.index] = cell.values
  r_plot = crop(r_tmp, crop_layer)
  plot(r_plot,col=tim.colors(16), add=T)#col=viridis(16))#
}


###########################
# Visualization

# plot some EVI summaries
# make a template for plotting 
target_pixels_na = target_pixels
target_pixels_na[target_pixels_na==0] = NA
plot_to_region(evi_summary$evi_mean, evi_summary$cell_number,subset_layer_albers); title("mean EVI")


plot_to_region(evi_summary$seas_change, evi_summary$cell_number, subset_layer_albers); title("Early-season EVI minus late-season EVI", cex.main=0.7)
plot_to_region(sqrt(evi_summary$among_year_var-min(evi_summary$among_year_var)), evi_summary$cell_number, subset_layer_albers); title("among-year sd")
plot_to_region(evi_summary$linear_trend, evi_summary$cell_number, subset_layer_albers); title("linear trend in EVI", cex.main=0.7)
plot_to_region(evi_summary$wet_dry_diff, evi_summary$cell_number, subset_layer_albers); title("Wet-year mean EVI minus dry-year mean EVI", cex.main=0.7)
plot_to_region(sqrt(evi_summary$mort_neigh_2014), evi_summary$cell_number, subset_layer_albers); title("Neighboring mortality in 2014", cex.main=0.7)

# Here the mean EVI, seasonal diff, and wet-dry diff are most interesting, esp if shown in extreme colors that highlight negative vs positive values 
# A lot of the prediction amounts to the sensitivity plus the mean layers. 
# What still needs to be done: some kind of interpretation of what these are reflecting -- e.g. type of veg, growing season patterns / phenology at different elevations and possibly different parts of the region. Possibly different local drought intensity differences? Forest density differences? 
# Maybe check simple correlations with density, elevation, drought intensity, and also visually compare with a forest-type map. 

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
plot(evi_mean_all, ylim=c(0.2, 0.4), type="l",lwd=2, col="cyan4", ylab="Regional mean EVI",xlab="Time step", cex.axis=1.1, cex.lab=1.5)
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





