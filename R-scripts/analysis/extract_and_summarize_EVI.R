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

# check that template includes the same region as the first evi data raster
#evi_000 = raster("features/ee-sn_jep_modis_ts_quality_mask_epsg3310/sn_jep_modis_ts_quality_mask_epsg3310_000_20000218.tif")
#sn = shapefile("features/SierraEcoregion_Jepson/SierraEcoregion_Jepson.shp")
#sn = spTransform(sn, albers.proj)
#cells_in_region = extract(evi_template, sn, cellnumbers=T)[[1]]
#vals_outside_region = getValues(evi_template)[-cells_in_region[,1]]
#length(vals_outside_region); sum(vals_outside_region)
#unique(vals_outside_region)
#vals_outside_region2 = getValues(evi_000)[-cells_in_region[,1]]
#length(vals_outside_region2); sum(vals_outside_region2>0)
#par(mfrow=c(1, 2)); plot(evi_template); plot(evi_000)
#unique(vals_outside_region2)
# results: almost identical, but there are a few cells that "leak out" of the sierra nevada polygon presumably as a result of reprojection. 
#rm(evi_000)

# load one mortality layer as a template
mort_template = raster("features/ADS-rasterized/Y2015_sp122.tif")
mort_2015_2016 = raster("features/ADS-rasterized/Y2015_spALL.tif") + raster("features/ADS-rasterized/Y2016_spALL.tif")

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
mort_albers = projectRaster(mort_2015_2016, evi_template)

# check 
extent(target_pixels) == extent(evi_template)
extent(mort_albers) == extent(evi_template)
length(getValues(target_pixels))
length(getValues(mort_albers))
length(getValues(evi_template))

testpoint = SpatialPoints(coords=matrix(c(-120.75, 38.9), ncol=2), proj4string = wgs.proj)
testpoint_alb = spTransform(testpoint, albers.proj)
plot(evi_template); points(testpoint_alb)
extract(evi_template, testpoint_alb, cellnumbers=T)
plot(mort_albers); points(testpoint_alb)
extract(mort_albers, testpoint_alb, cellnumbers=T)
plot(target_pixels); points(testpoint_alb)
extract(target_pixels, testpoint_alb, cellnumbers=T)
extract(raster("features/ee-sn_jep_modis_ts_quality_mask_epsg3310/sn_jep_modis_ts_quality_mask_epsg3310_000_20000218.tif"), testpoint_alb, cellnumbers=T)
# OK, returns the same cell index for all rasters. 


# Function to extract time series of EVI values for cells identified in the target_pixels layer
# extracts from the set of geotifs in geotif_folder with filename listed in geotif_filenames

stack_evi_layers <- function(geotif_folder, geotif_filenames) {
  r = raster(paste(geotif_folder, geotif_filenames[1], sep=""))
  evi_stack = stack(r)
  for (i in 2:length(geotif_filenames)) {
    r = raster(paste(geotif_folder, geotif_filenames[i], sep=""))
    evi_stack = stack(evi_stack, r)
  }
  return(evi_stack)
}


# Build the stack
evi_stack = stack_evi_layers(geotif_folder=geotif_folder, geotif_filenames=filenames)

# Extract EVI stack values into a matrix with pixels on the rows and times on the columns
evi_mat = getValues(evi_stack)
colnames(evi_mat) = date_codes # rename columns to indicate date of observation
rownames(evi_mat) = as.character(1:length(evi_template)) # rename rows to indicate location of pixel within source raster. This is key for later for extracting mortality values and for matching analysis / summary results back to cells in the reference rasters. 

# retain only the cells that fall within target veg type
# and also within the EVI template that defines the region
# and also were not disturbed from 2000 onward
evi_target_index = getValues(target_pixels)==1
evi_mask_index = getValues(evi_template)==1
fire_dates = raster("features/sierra_nevada_250m_most_recent_fire.tif")
mgt_dates = raster("features/sierra_nevada_250m_most_recent_management.tif")
fire_dates[is.na(fire_dates)] = 0 # put values in NA cells for ease of indexing
mgt_dates[is.na(mgt_dates)] = 0 
disturb_index = !(getValues(mgt_dates)>=2000 | getValues(fire_dates)>=2000)

evi_mat = evi_mat[evi_mask_index & evi_target_index & disturb_index,] # Subset to the rows that are within the domain as defined in the template (cells outside this domain are 0s because that's how EarthEngine works, and there are a few NAs we can also exclude -- check this)

# Convert the missing values from Earth Engine (zeros) into NA's. 
# also convert negative EVI values into NA's
evi_mat[evi_mat<=0] = NA

# Remove rows that have no EVI data, if any
z = apply(evi_mat, 1, f<-function(x){return(sum(!is.na(x)))})
sum(z==0)
evi_mat = evi_mat[z>0,]
object.size(evi_mat)

# Rescale
evi_mat = evi_mat/10000 # rescale to standard EVI values

save(evi_mat, file="features/working-files/evi_data_matrix_jepson_PPN+SMC.Rdata")



######################################################
# Summarize the EVI time series

# load evi data matrix
load(file="features/working-files/evi_data_matrix_jepson_PPN+SMC.Rdata")

# how many "good" values are there per month? 
obs_by_mon = rep(NA, 12)
for (i in 0:11) obs_by_mon[i+1] = sum(!is.na(evi_mat[,dates$mon==i]))
barplot(obs_by_mon, names.arg=as.character(1:12))
# all values present June-Sept, good number in May. April and Oct missing ~40% of values

# Temporally mask out all months but May-Sept
# and the the years after 2013 (i.e. the later drought years)
# Q should we also include 2014?
startyear = 2000
endyear = 2013
evi_months = c(4,5,6,7,8) # which months -- note month numbers are 0-11
time_index = as.integer(which(dates$mon %in% evi_months & dates$year <= (endyear-1900) & dates$year >= (startyear-1900)))
plot(evi_mat[100,time_index])


### Make single-number summaries of EVI time series and store in data frame

# Summarize extracted values into data frame
evi_summary = data.frame(cell_number=as.integer(rownames(evi_mat)))
evi_summary$evi_mean = apply(evi_mat[,time_index], 1, mean, na.rm=T)
evi_summary$evi_mayjun = apply(evi_mat[,dates$mon %in% c(4, 5) & dates$year <= (endyear-1900) & dates$year >= (startyear-1900) ], 1, mean, na.rm=T)
evi_summary$evi_augsept = apply(evi_mat[,dates$mon %in% c(7,8) & dates$year <= (endyear-1900) & dates$year >= (startyear-1900) ], 1, mean, na.rm=T)
evi_summary$seas_change = evi_summary$evi_augsept - evi_summary$evi_mayjun
evi_summary$seas_change_prop = (evi_summary$evi_augsept/evi_summary$evi_mayjun)-1
evi_summary$total_var = apply(evi_mat[,time_index], 1, var, na.rm=T)

# wet-year vs dry-year difference in late-season EVI
# define wet years as 2000, 2005, 2006, 2010, 2011
# define dry years as 2002, 2007 ( could also include 2013 if that year's in the training data)
wetmean = apply(evi_mat[,dates$year %in% c(100,105, 106, 2010, 2011) & dates$mon %in% c(7,8)], 1, mean, na.rm=T)
drymean = apply(evi_mat[,dates$year %in% c(102,107, 113) & dates$mon %in% c(7,8)], 1, mean, na.rm=T) # Q whether to include 2013 here
evi_summary$wet_dry_diff = drymean-wetmean
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

# among-year variance
evi_summary$among_year_var = apply(annual.mean, 1, var, na.rm=T)

# ratio of among-to-within-year variance 
evi_summary$among_to_within_ratio = evi_summary$among_year_var / evi_summary$within_year_var




### Quick look 
sub = sample(1:nrow(evi_mat), size=500, replace=FALSE)
#pairs(evi_summary[sub,])

# Note removing the disturbed pixels seems to have removed the very high variance outliers
# look at cells with high within-year variance 
outliers = which(evi_summary$within_year_var>0.01)
par(mfrow=c(5,5), mar=rep(2,4))
for (i in 1:25) plot(evi_mat[outliers[i],dates$mon %in% c(5,6,7,8,9)]~linear_time[dates$mon %in% c(5,6,7,8,9)], pch=16, cex=0.5, ylim=c(0, 1))
title("timeseries of EVI 2000-2016")
# same for high among-year variance 
outliers = which(evi_summary$among_year_var>0.005)
par(mfrow=c(5,5), mar=rep(2,4))
for (i in 1:25) plot(evi_mat[outliers[i],dates$mon %in% c(5,6,7,8,9)]~linear_time[dates$mon %in% c(5,6,7,8,9)], pch=16, cex=0.5, ylim=c(0, 1))
title("timeseries of EVI 2000-2016")
# these seem now to represent mainly cell with strong time-trends


# add the mortality data 
mort_masked = mask(mort_albers, target_pixels, maskvalue=0)
evi_summary$mort = getValues(mort_masked)[as.integer(rownames(evi_mat))]
pairs(evi_summary[sub,])

x = as.matrix(cor(evi_summary[evi_summary$among_year_var<0.002 & evi_summary$within_year_var<0.005,], use="pairwise.complete"))
heatmap(x, col=viridis(12))
x


# Store intermediate file 
save(evi_summary, file="features/working-files/evi_summary_PPN+SMC_jepson.Rdata")


############################################
# Do some simple regressions

load("features/working-files/evi_summary_PPN+SMC_jepson.Rdata")

### Run a simple model to check associations -- use tobit model in vgam library
hist(evi_summary$mort)
cols_to_standardize = c("evi_mayjun", "seas_change", "within_year_var", "among_year_var", "wet_dry_diff", "linear_trend")
for (i in 1:length(cols_to_standardize)) evi_summary[,cols_to_standardize[i]] = scale(evi_summary[,cols_to_standardize[i]])

# check for correlation in explanatory variables
vif(lm(mort~evi_mayjun+seas_change+within_year_var + among_year_var+wet_dry_diff+linear_trend, data=evi_summary))
cor(evi_summary[,cols_to_standardize], use="pairwise.complete")

# fit model
m = vglm(mort~evi_mayjun+seas_change+wet_dry_diff+within_year_var + among_year_var+linear_trend, tobit, data=evi_summary, trace=TRUE)
summary(m)

plot(evi_summary$mort[!is.na(evi_summary$mort)]~predict(m, type="response"))
abline(0,1)


#####################
# Make output plots

plot_to_region <- function(values, index, target_pixels) { # index is the row numbers of the cells to plot, and indexes grid cells in the original evi_template and target_pixels rasters
  # values is the values to assign to these, using the cell numbers in index
  # it uses target_pixels as the template rasters
  plotraster = target_pixels
  plotraster[index] = values
  plot(plotraster, col=viridis(12))
}


# plot some EVI summaries
#par(mfrow=c(2,2))
target_pixels_na = target_pixels
target_pixels_na[target_pixels_na==0] = NA
plot_to_region(evi_summary$evi_mayjun, evi_summary$cell_number, target_pixels_na); title("May-Jun mean EVI")
plot_to_region(evi_summary$seas_change, evi_summary$cell_number, target_pixels_na); title("early- to late-season change in EVI")

#
plot_to_region(evi_summary$among_year_var, evi_summary$cell_number, target_pixels); title("among-year variance")
plot_to_region(evi_summary$wet_dry_diff, evi_summary$cell_number, target_pixels_na); title("Wet-to-dry-year change in EVI")

# observed and predicted mortality 
mort_pred = fitted(m)
par(mfrow=c(1,2))
plot_to_region(sqrt(evi_summary$mort), evi_summary$cell_number, target_pixels_na)
plot_to_region(sqrt(mort_pred+121), evi_summary$cell_number[!is.na(evi_summary$mort)], target_pixels_na) # something wrong here! !

plot(evi_summary$mort[!is.na(evi_summary$mort)][1:89506]~sqrt(mort_pred+121)); abline(c(0,1))

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


summary(lm(sqrt(mort)~evi_mayjun+seas_change_prop+wet_dry_diff+within_year_var + among_year_var+linear_trend, data=evi_summary[!is.na(evi_summary$mort),]))



