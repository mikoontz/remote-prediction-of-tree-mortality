library(sp)
library(raster)
library(rgdal)
library(lattice)
library(lme4)
library(MASS)
library(fields)

# Enter the locations of files to work with 
# geotifs are the MODIS EVI data that Mike K exported from Google Earth Engine, they are stored in a single folder (geotif_folder) and are all names consistently with the prefix geotif_filename and a date code. 
geotif_folder = "./features/ee_sierra-nevada-forest-quality-mask-modis-time-series/"
geotif_filename =  "sn-whole-ts-modis-forest-quality-mask-"
tmp = dir(geotif_folder)
filenames = tmp[grep(geotif_filename, tmp)]

# Get date information from the geotif filenames
date_codes = sapply(filenames, substr, start=nchar(geotif_filename)+1, stop=nchar(geotif_filename)+8) 
dates = strptime(date_codes, "%Y%m%d")

# Get the locations of the target pixels 
# target cover raster 
target_cover = raster("features/sierra_nevada_250m_calveg_pipo_forest-cover_whr_type.tif")
# load one mortality layer as a template
mort_template = raster("features/ADS-rasterized/Y2015_sp122.tif")
mort_2015_2016 = raster("features/ADS-rasterized/Y2015_sp122.tif") + raster("features/ADS-rasterized/Y2016_sp122.tif")
# subset the target cover raster to the mortality area 
target_cover_sub = crop(target_cover, mort_template)
# load one EVI layer as a template
evi_template = raster(paste(geotif_folder, geotif_filename, date_codes[1], ".tif", sep=""))

# Create a target cover layer for pixels with specified percent PIPO cover
PIPO_cover_min = 80
target_pixels = target_cover_sub >= PIPO_cover_min
target_pixels[target_pixels==0] = NA # set the non-target values to NA which is default value for masking
#plot(target_pixels)

# Reproject target mask to match the EVI geotiffs (in Albers projection)
target_pixels <- projectRaster(target_pixels, evi_template)
mort
# ISSUE- this changes the extent to the extent of the whole evi raster

# Also reproject mortality data to match the EVI geotifs
mort_albers = projectRaster(mort_2015_2016, evi_template)

# Function to extract time series of EVI values for pixels identified in the target_pixels layer
# extracts from the set of geotifs in geotif_folder with filename starting with geotif_filename and ending in an integer date code, as specified in geotif_date_codes.
# Note this is very slow, since it reads in the whole raster for each time step, but for this reason also requires little memory. 
# Probably should make one that first assembles a rasterbrick, then drills through it to get the time series. 
stack_evi_layers <- function(target_pixels, geotif_folder, geotif_filename, geotif_date_codes) {
  r = raster(paste(geotif_folder, geotif_filename, geotif_date_codes[1], ".tif", sep=""))
  evi_stack = stack(r)
  for (i in 2:length(geotif_date_codes)) {
    r = raster(paste(geotif_folder, geotif_filename, geotif_date_codes[i], ".tif", sep=""))
    evi_stack = stack(evi_stack, r)
  }
  return(evi_stack)
}



evi_stack = stack_evi_layers(target_pixels, geotif_folder, geotif_filename, date_codes)
n_times = nlayers(evi_stack)

# turn the EVI values into a matrix with pixels on the rows and times on the columns
evi_mat = getValues(evi_stack)
evi_target_index = which(!is.na(getValues(target_pixels)))
evi_mat = evi_mat[evi_target_index,]
sum(is.na(evi_mat)); sum(evi_mat==1)
evi_mat[evi_mat==1] = NA # replace the NA values
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
m = vglm(mort~evi_mayjun+seas_change_prop+within_year_var + among_year_var+ linear_trend+wet_dry_diff, tobit, data=evi_summary, trace=TRUE)
summary(m)
plot(evi_summary$mort[!is.na(evi_summary$mort)]~predict(m, type="response"))
abline(0,1)


#####################
# Make output plots

plot_to_subregion <- function(values, index, target_pixels, target_cover_sub) { # index is the row numbers of the cells to plot, and indexes grid cells in the original evi_template and target_pixels rasters
  # values is the values to assign to these 
  # it uses target_pixels as the template rasters, and target_cover_sub as the extent and coordinate system to display the plot in
  plotraster = target_pixels
  plotraster[index] = values
  plotraster[plotraster>0.95] = NA # get rid of excess indicator values in the template raster
  plotraster = projectRaster(plotraster, target_cover_sub)
  plot(plotraster, col=viridis(12))
}

# plot some EVI summaries
#par(mfrow=c(2,2))
plot_to_subregion(evi_summary$evi_mayjun, evi_summary$cell_number, target_pixels, target_cover_sub); title("May-Jun mean EVI")
plot_to_subregion(evi_summary$seas_change, evi_summary$cell_number, target_pixels, target_cover_sub); title("early- to late-season change in EVI")
plot_to_subregion(evi_summary$linear_trend, evi_summary$cell_number, target_pixels, target_cover_sub); title("Linear trend 2000-2012")
plot_to_subregion(evi_summary$wet_dry_diff, evi_summary$cell_number, target_pixels, target_cover_sub); title("Wet-to-dry-year change in EVI")


plot_to_subregion(evi_summary$mort, evi_summary$cell_number, target_pixels, target_cover_sub)

# look at random individual pixels
par(mfrow=c(4,4), mar=rep(2, 4))
for (i in 1:16) plot(evi_mat[sample(1:1203, 1),dates$mon %in% c(5,6,7,8,9)]~linear_time[dates$mon %in% c(5,6,7,8,9)], pch=16, cex=0.4, ylim=c(0.4, 0.9))

# all pixels averaged
evi_mean_all = apply(evi_mat[,dates$mon %in% c(5,6,7,8,9)], 2, mean, na.rm=T)
plot(evi_mean_all~linear_time[dates$mon %in% c(5,6,7,8,9)])

# linear model just to vaguely assess fit
summary(lm(mort~evi_mayjun+seas_change_prop+within_year_var+wet_dry_diff, data=evi_summary))

