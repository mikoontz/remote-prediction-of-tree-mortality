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
mort_2015_2016 = raster("features/ADS-rasterized/Y2015_sp122.tif") + raster("features/ADS-rasterized/Y2016_sp122.tif")

# Create a target cover layer for pixels with specified percent forest type cover
# Start by focusing only on PPN for test run
cover_min = 80 # Minumum cover of WHR type
target_cover = raster("features/calveg-pct-cover-rasters/sierra_nevada_250m_calveg_cover_whr_type_PPN.tif")
target_pixels = target_cover >= cover_min

# Reproject rasters to match the EVI geotiffs (in Albers projection)
target_albers <- projectRaster(target_pixels, evi_template)
mort_albers = projectRaster(mort_2015_2016, evi_template)

# check 
extent(target_albers) == extent(evi_template)
extent(mort_albers) == extent(evi_template)
length(getValues(target_albers))
length(getValues(mort_albers))
length(getValues(evi_template))

testpoint = SpatialPoints(coords=matrix(c(-120.75, 38.9), ncol=2), proj4string = wgs.proj)
testpoint_alb = spTransform(testpoint, albers.proj)
plot(evi_template); points(testpoint_alb)
extract(evi_template, testpoint_alb, cellnumbers=T)
plot(mort_albers); points(testpoint_alb)
extract(mort_albers, testpoint_alb, cellnumbers=T)
plot(target_albers); points(testpoint_alb)
extract(target_albers, testpoint_alb, cellnumbers=T)
# OK, returns the same cell index for all rasters. 


# Function to extract time series of EVI values for cells identified in the target_pixels layer
# extracts from the set of geotifs in geotif_folder with filename listed in geotif_filenames

stack_evi_layers <- function(target_pixels, geotif_folder, geotif_filenames) {
  r = raster(paste(geotif_folder, geotif_filenames[1], sep=""))
  evi_stack = stack(r)
  for (i in 2:length(geotif_filenames)) {
    r = raster(paste(geotif_folder, geotif_filenames[i], sep=""))
    evi_stack = stack(evi_stack, r)
  }
  return(evi_stack)
}



evi_stack = stack_evi_layers(target_pixels=target_albers, geotif_folder=geotif_folder, geotif_filenames=filenames)
n_times = nlayers(evi_stack)

# turn the EVI values into a matrix with pixels on the rows and times on the columns
evi_mat = getValues(evi_stack)
evi_target_index = which(!is.na(getValues(target_albers)))
evi_mat = evi_mat[evi_target_index,]
evi_mat = evi_mat/10000 # rescale to standard EVI values
colnames(evi_mat) = date_codes # columns indicate date of observation
rownames(evi_mat) = evi_target_index # rows indicate location of pixel within source raster

######################################################
# Summarize the EVI time series

# how many "good" values are there per month? 
#obs_by_mon = rep(NA, 12)
#for (i in 0:11) obs_by_mon[i+1] = sum(!is.na(evi_mat[,dates$mon==i]))
#barplot(obs_by_mon, names.arg=as.character(1:12))
# all values present June-Sept, almost all in May too

# temporally mask out all months but May-Sept
# and the the years after 2012 (i.e. the drought years)
# Q should we include the "early" drought years of 2013-14?
startyear = 2000
endyear = 2012
evi_months = c(4,5,6,7,8) # which months -- note month numbers arew 0-11
time_index = as.integer(which(dates$mon %in% evi_months & dates$year <= (endyear-1900) & dates$year >= (startyear-1900)))
plot(evi_mat[10,time_index])

# Check how many pixels have EVI data at all 
n_pixels = nrow(evi_mat)
missing_index = rep(0, n_pixels)
for (i in 1:n_pixels) missing_index[i] = sum(!is.na(evi_mat[i,]))
sum(missing_index>0) # 1203 pixels -- about half -- contain EVI data

# Further subset the data to areas that have EVI information
evi_mat = evi_mat[missing_index>0,]

### Make single-number summaries of EVI time series and store in data frame


# Summarize extracted values into data frame
evi_summary = data.frame(cell_number=as.integer(rownames(evi_mat)))
evi_summary$evi_mean = apply(evi_mat[,time_index], 1, mean, na.rm=T)
evi_summary$evi_mayjun = apply(evi_mat[,dates$mon %in% c(4, 5) & dates$year <= (endyear-1900) & dates$year >= (startyear-1900) ], 1, mean, na.rm=T)
evi_summary$evi_sept = apply(evi_mat[,dates$mon == 9 & dates$year <= (endyear-1900) & dates$year >= (startyear-1900) ], 1, mean, na.rm=T)
evi_summary$seas_change = evi_summary$evi_sept - evi_summary$evi_mayjun
evi_summary$seas_change_prop = (evi_summary$evi_sept/evi_summary$evi_may)-1
evi_summary$total_var = apply(evi_mat[,time_index], 1, var, na.rm=T)

# wet-year vs dry-year difference in late-season EVI
# define wet years as 2000, 2005, 2006, 2010, 2011
# define dry years as 2002, 2007 ( could also include 2013 if that year's in the training data)
wetmean = apply(evi_mat[,dates$year %in% c(100,105, 106, 2010, 2011) & dates$mon %in% c(7,8)], 1, mean, na.rm=T)
drymean = apply(evi_mat[,dates$year %in% c(102,107, 113) & dates$mon %in% c(7,8)], 1, mean, na.rm=T)
evi_summary$wet_dry_diff = drymean-wetmean
evi_summary$wet_dry_propdiff = drymean/wetmean-1

# add trends
linear_time = scale(as.integer(dates[1:length(dates)]-dates[1]))
for (i in 1:nrow(evi_mat)) evi_summary$linear_trend[i] = coef(lm(evi_mat[i,time_index]~linear_time[time_index]))[2]


# within-year variance
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

# among-year variance
evi_summary$among_year_var = apply(annual.mean, 1, var, na.rm=T)

# ratio of among-to-within-year variance 
evi_summary$among_to_within_ratio = evi_summary$among_year_var / evi_summary$within_year_var

### Quick look 
#pairs(evi_summary[evi_summary$among_year_var<0.002 & evi_summary$within_year_var<0.005,])
# higher within-year variance associated with lower EVI. Maybe reflecting herbaceous cover? Higher among-year variance associated with linear trends (both positive and negative.)


# a few cells have outlier-high variance -- look at these
outliers = which(evi_summary$among_year_var>0.002)
par(mfrow=c(3,3))
for (i in 1:length(outliers)) plot(evi_mat[outliers[i],dates$mon %in% c(5,6,7,8,9)]~linear_time[dates$mon %in% c(5,6,7,8,9)], pch=16, cex=0.5, ylim=c(0.4, 0.9))
title("timeseries of EVI 2000-2016")
# interesting -- two of the high-variance outliers have sharp break in the middle suggesting fire or logging. Others show strong trends, maybe representing areas recovering from disturbance? 


# add the mortality data 
### NEED TO CHECK WHETHER THIS IS GETTING THE RIGHT PIXELS FROM MORTALITY DATA! 
mort_masked = mask(mort_albers, target_pixels)
evi_summary$mort = getValues(mort_masked)[as.integer(rownames(evi_mat))]
pairs(evi_summary[evi_summary$among_year_var<0.002 & evi_summary$within_year_var<0.005,])

x = as.matrix(cor(evi_summary[evi_summary$among_year_var<0.002 & evi_summary$within_year_var<0.005,], use="pairwise.complete"))
heatmap(x, col=viridis(12))
# Note when we include all years, among-year variance has strongest correlation with mortality (0.29)
# When we include just the pre-drought years 2000-2012, high evi, especially in early season, is positively correlated with mortality. Difference in late-season EVI in dry versus wet years is also strongly associated with mortality (sites that showed a drop earlier tended to have more mortality later). 


### Run a simple model to check associations -- use tobit model in vgam library
hist(evi_summary$mort)
cols_to_standardize = c("evi_mayjun", "seas_change_prop", "within_year_var", "among_year_var", "linear_trend", "wet_dry_diff")
for (i in 1:length(cols_to_standardize)) evi_summary[,cols_to_standardize[i]] = scale(evi_summary[,cols_to_standardize[i]])

# check for correlation in explanatory variables
vif(lm(mort~evi_mayjun+seas_change_prop+within_year_var + among_year_var+ linear_trend+wet_dry_diff, data=evi_summary))
cor(evi_summary[,cols_to_standardize])

# fit model
m = vglm(mort~evi_mayjun+seas_change_prop+within_year_var + among_year_var+ linear_trend+wet_dry_diff, tobit, data=evi_summary, trace=TRUE)
summary(m)

plot(evi_summary$mort[!is.na(evi_summary$mort)]~predict(m, type="response"))
abline(0,1)


#####################
# Make output plots

plot_EVI_to_subregion <- function(values, index, target_pixels, target_cover_sub) { # index is the row numbers of the cells to plot, and indexes grid cells in the original evi_template and target_pixels rasters
  # values is the values to assign to these 
  # it uses target_pixels as the template rasters, and target_cover_sub as the extent and coordinate system to display the plot in
  plotraster = target_pixels
  plotraster[index] = values
  plotraster[plotraster>0.95] = NA # get rid of excess indicator values in the template raster
  plotraster = projectRaster(plotraster, target_cover_sub)
  plot(plotraster, col=viridis(12))
}

plot_to_subregion <- function(values, index, target_pixels, target_cover_sub) { # index is the row numbers of the cells to plot, and indexes grid cells in the original evi_template and target_pixels rasters
  # values is the values to assign to these 
  # it uses target_pixels as the template rasters, and target_cover_sub as the extent and coordinate system to display the plot in
  plotraster = target_pixels
  plotraster[index] = values
  plotraster = projectRaster(plotraster, target_cover_sub)
  plot(plotraster, col=viridis(12))
}


# plot some EVI summaries
#par(mfrow=c(2,2))
plot_EVI_to_subregion(evi_summary$evi_mayjun, evi_summary$cell_number, target_pixels, target_cover_sub); title("May-Jun mean EVI")
plot_EVI_to_subregion(evi_summary$seas_change, evi_summary$cell_number, target_pixels, target_cover_sub); title("early- to late-season change in EVI")
plot_EVI_to_subregion(evi_summary$linear_trend, evi_summary$cell_number, target_pixels, target_cover_sub); title("Linear trend 2000-2012")
plot_EVI_to_subregion(evi_summary$wet_dry_diff, evi_summary$cell_number, target_pixels, target_cover_sub); title("Wet-to-dry-year change in EVI")

# observed and predicted mortality 
plot_to_subregion(evi_summary$mort, evi_summary$cell_number, target_pixels, target_cover_sub)
plot_to_subregion(predict(m, type="response"), evi_summary$cell_number[!is.na(evi_summary$mort)], target_pixels, target_cover_sub)

# look at random individual pixels
par(mfrow=c(4,4), mar=rep(2, 4))
for (i in 1:16) plot(evi_mat[sample(1:1203, 1),dates$mon %in% c(5,6,7,8,9)]~linear_time[dates$mon %in% c(5,6,7,8,9)], pch=16, cex=0.4, ylim=c(0.4, 0.9))

# all pixels averaged
evi_mean_all = apply(evi_mat, 2, mean, na.rm=T)
# long time series
plot(evi_mean_all~linear_time)
# I'd say this shows that 2013 was low, clearly a drought year, but not out of the normal range for the rest of the years. So for model fitting, seems ok to go through 2013. The later years are drastically low, esp 2016. Will be interesting to see the rebound in 2017, if any. 
# plotting the spatial average for each year. 
plot(evi_mean_all[dates$year==100]~dates$yday[dates$year==100], type="l", lwd=2, col="cyan4", ylim=c(0.4, 0.9))
for (i in 101:112) lines(evi_mean_all[dates$year==i]~dates$yday[dates$year==i], type="l", lwd=2, col="cyan4")
for (i in 113:116) lines(evi_mean_all[dates$year==i]~dates$yday[dates$year==i ], type="l", lwd=2, col="orange3")
# For this region, 2016 looks pretty flat -- mortality mainly happened in 2015 it appears.

# linear model just to vaguely assess fit
summary(lm(sqrt(mort)~evi_mayjun+seas_change_prop+within_year_var + among_year_var+ linear_trend+wet_dry_diff, data=evi_summary))
