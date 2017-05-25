library(sp)
library(raster)
library(rgdal)
library(lattice)
library(lme4)

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
extract_target_evi <- function(target_pixels, geotif_folder, geotif_filename, geotif_date_codes) {
  r = raster(paste(geotif_folder, geotif_filename, geotif_date_codes[1], ".tif", sep=""))
  #evi_crop = crop(r, extent(target_pixels))
  #evi_mask = mask(evi_crop, target_pixels)
  evi_mask = mask(r, target_pixels)
  evi_stack = stack(evi_mask)
  for (i in 2:length(geotif_date_codes)) {
    r = raster(paste(geotif_folder, geotif_filename, geotif_date_codes[i], ".tif", sep=""))
    #evi_crop = crop(r, extent(target_pixels))
    evi_mask = mask(r, target_pixels)
    evi_stack = stack(evi_stack, evi_mask)
  }
  evi_stack = evi_stack/10000 # rescale to evi scale
  return(evi_stack)
}



target_evi_stack = extract_target_evi(target_pixels, geotif_folder, geotif_filename, date_codes)

# turn the values into a matrix with pixels on the rows and times on the columns
evi_mat = getValues(target_evi_stack)
evi_target_index = which(!is.na(evi_mat[,1]))
evi_mat = evi_mat[evi_target_index,]
evi_mat[evi_mat==0.0001] = NA # replace the NA values
colnames(evi_mat) = date_codes # columns indicate date of observation
rownames(evi_mat) = evi_target_index # rows indicate location of pixel within source raster

######################################################
# Summarize the EVI time series

# how many pixels have EVI time series? 
n_pixels = nrow(evi_mat)
missing_index = rep(0, n_pixels)
for (i in 1:n_pixels) missing_index[i] = sum(!is.na(evi_mat[i,]))
sum(missing_index>0)
# how many "good" values are there per month? 
obs_by_mon = rep(NA, 12)
for (i in 0:11) obs_by_mon[i+1] = sum(!is.na(evi_mat[,dates$mon==i]))
barplot(obs_by_mon, names.arg=as.character(1:12))
# all values present June-Sept, almost all in May too

# temporally mask out all months but May-Sept
# and the the years after 2012 (i.e. the drought years)
# Q should we include the "early" drought years of 2013-14?
startyear = 2000
endyear = 2012
time_index = which(dates$mon %in% c(4,5,6,7,8) & dates$year <= (endyear-1900) & dates$year >= (startyear-1900))
plot(evi_mat[10,time_index])

# Check how many pixels have EVI data at all 
n_pixels = nrow(evi_mat)
missing_index = rep(0, n_pixels)
for (i in 1:n_pixels) missing_index[i] = sum(!is.na(evi_mat[i,]))
sum(missing_index>0) # 1203 pixels -- about half -- contain EVI data

# Further subset the data to areas that have EVI information
evi_mat = evi_mat[missing_index>0,]

### Make single-number summaries of EVI time series and store in data frame

evi_summary = data.frame(cell_number=as.integer(rownames(evi_mat)))
evi_summary$evi_mean = apply(evi_mat[,time_index], 1, mean, na.rm=T)
evi_summary$evi_mayjun = apply(evi_mat[,dates$mon %in% c(4, 5) & dates$year <= (endyear-1900) & dates$year >= (startyear-1900) ], 1, mean, na.rm=T)
evi_summary$evi_sept = apply(evi_mat[,dates$mon == 9 & dates$year <= (endyear-1900) & dates$year >= (startyear-1900) ], 1, mean, na.rm=T)
evi_summary$seas_change = evi_summary$evi_sept - evi_summary$evi_mayjun
evi_summary$seas_change_prop = (evi_summary$evi_sept/evi_summary$evi_may)-1
evi_summary$total_var = apply(evi_mat[,time_index], 1, var, na.rm=T)

# add trends
linear_time = scale(as.integer(dates[1:length(dates)]-dates[1]))
for (i in 1:nrow(evi_mat)) evi_summary$linear_trend[i] = coef(lm(evi_mat[i,time_index]~linear_time[time_index]))[2]


# add the mortality data 
### NEED TO CHECK WHETHER THIS IS GETTING THE RIGHT PIXELS FROM MORTALITY DATA! 
mort_masked = mask(mort_albers, target_pixels)
evi_summary$mort = getValues(mort_masked)[as.integer(rownames(evi_mat))]




# within-year variance
years = 100:112
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


### Quick look 
pairs(evi_summary)


# ratio of among- to within-year variance 



# look at random individual pixels
par(mfrow=c(4,4), mar=rep(2, 4))
for (i in 1:16) plot(evi_mat[sample(1:1203, 1),dates$mon %in% c(5,6,7,8,9)]~linear_time[dates$mon %in% c(5,6,7,8,9)], type="l", ylim=c(0.4, 0.9))


