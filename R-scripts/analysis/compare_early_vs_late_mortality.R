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


# create a variable that describes average mortality in a grid cell's neighbors
evi_points = rasterToPoints(evi_template)
summary_points = evi_points[evi_summary$cell_number,1:2]
evi_neigh = dnearneigh(summary_points, d1=100, d2=740, longlat=FALSE)
mort_neigh = lapply(evi_neigh, f<- function(x){return(mean(evi_summary$mort_2014[x], na.rm=T))})
evi_summary$mort_neigh = unlist(mort_neigh)
hist(log(evi_summary$mort_neigh+0.1))

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

# Store intermediate file 
save(evi_summary, file="features/working-files/evi_summary_with_mort_PPN+SMC_jepson_central+south.Rdata")


############################################
# Do some simple regressions for different years

### Run a simple model to check associations -- use tobit model in vgam library
cols_to_standardize = c("evi_mean", "seas_change", "within_year_sd", "among_year_sd", "wet_dry_diff", "linear_trend")
for (i in 1:length(cols_to_standardize)) evi_summary[,cols_to_standardize[i]] = scale(evi_summary[,cols_to_standardize[i]])

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

barplot(coef(m)[3:8], horiz=T, las=2, main="tobit model coefficients", col=ifelse(coef(m)[3:8]<0, "red", "blue"), cex.names=0.5)

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

#####################
# Make output plots

plot_to_region <- function(cell.values, cell.index, crop_layer) { # index is the row numbers of the cells to plot, and indexes grid cells in the original evi_template
  r_tmp = evi_template
  r_tmp = setValues(r_tmp, rep(NA, length(r_tmp)))
  r_tmp[cell.index] = cell.values
  r_plot = crop(r_tmp, crop_layer)
  plot(r_plot,col=tim.colors(16))#col=viridis(10))
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



