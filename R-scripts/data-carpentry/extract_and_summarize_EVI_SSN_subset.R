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
# ISSUE- this changes the extent to the extent of the whole evi raster

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
time_index = which(dates$mon %in% c(4,5,6,7,8) & dates$year<=112)
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
evi_summary$evi_mayjun = apply(evi_mat[,dates$mon %in% c(4, 5) & dates$year <= 112], 1, mean, na.rm=T)
evi_summary$evi_sept = apply(evi_mat[,dates$mon == 9 & dates$year <= 112], 1, mean, na.rm=T)
evi_summary$seas_change = evi_summary$evi_sept - evi_summary$evi_mayjun
evi_summary$seas_change_prop = (evi_summary$evi_sept/evi_summary$evi_may)-1

# add trends
linear_time = scale(as.integer(dates[1:length(dates)]-dates[1]))
for (i in 1:nrow(evi_mat)) evi_summary$linear_trend[i] = coef(lm(evi_mat[i,time_index]~linear_time[time_index]))[2]

# within-year variance
evi_summary_

# within-year CV

# among-year variance

# among-year CV

# ratio of among- to within-year variance 


par(mfrow=c(4,4), mar=rep(2, 4))
for (i in 1:16) plot(evi_mat[i,]~linear_time, type="l", ylim=c(0.4, 0.9))


###################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# code carried over from other script below here. 

# look at annual variation for Shaver Lake points
plot.evi.slopes <- function(z, dates) {
  zsub = z[z$mon %in% 6:9,]
  m=lmer(evi~daystd + (1+daystd|year), data=zsub)
  coefs = coef(m)$year
  plot.new(); curve(coefs[1,1]+coefs[1,2]*x, from=0, to=1.5, ylim=c(0.5, 0.9))
  for(i in 1:length(unique(dates$year))) curve(coefs[i,1]+coefs[i,2]*x, add=T, col=i)
  points(evi~daystd, zsub, col="gray")
}

z = data.frame(evi=apply(evivals_shaver, 1, mean, na.rm=T),  year=dates$year, mon=dates$mon, day=dates$yday, daystd = scale(dates$yday))
plot.evi.slopes(z, dates)

z = data.frame(evi=apply(evivals_ill, 1, mean, na.rm=T), year=dates$year, mon=dates$mon, day=dates$yday, daystd = scale(dates$yday))
plot.evi.slopes(z, dates)

z = data.frame(evi=portlo_mean,  year=dates$year, mon=dates$mon, day=dates$yday, daystd = scale(dates$yday))
plot.evi.slopes(z, dates)

# High vs low mortality near Porterville S Sierras, ~6000 feet
locs.porthi= read_GE_locs("./features/example_locations/Porterville_Hi_Mort.txt", header=F)
locs.porthi = SpatialPoints(project(locs.porthi, proj4string(r)))
evivals_porthi = extract_evi(locs.porthi, geotif_folder, geotif_filename, geotif_numbers = datestr)

locs.portlo= read_GE_locs("./features/example_locations/Porterville_Lomort.txt", header=F)
locs.portlo = SpatialPoints(project(locs.portlo, proj4string(r)))
evivals_portlo = extract_evi(locs.portlo, geotif_folder, geotif_filename, geotif_numbers = datestr)


par(mfrow=c(3, 4))
for (i in 1:10) {
  plot(evivals_portlo[,i], type="l", lwd=2, col="cyan3", ylim=c(0.5, 0.9))
  plot(evivals_porthi[,i], type="l", lwd=2, col="orange2", ylim=c(0.5, 0.9))
}

## Yikes almost all these points are masked out for the entire time series!! 
# Compare smoothed means across sites
porthi_mean = apply(evivals_porthi, 1, mean, na.rm=T)
portlo_mean = apply(evivals_portlo, 1, mean, na.rm=T)
shaver_mean = apply(evivals_shaver, 1, mean, na.rm=T)
ill_mean  = apply(evivals_ill, 1, mean, na.rm=T)
plot(porthi_mean, type="l", col="red")
lines(portlo_mean,  col="blue")
# Not sure there are obvious differences here before 2015 -- possibly late-season browning in 2013 in the high-mortality areas. 

plot(porthi_mean~portlo_mean); abline(0,1) # for these few pixels, "high" mortality pixels have less total variation than "low" mortality pixels


