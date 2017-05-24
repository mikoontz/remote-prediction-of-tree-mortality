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

# Reproject target mask to match the EVI geotiffs
target_pixels <- projectRaster(target_pixels, evi_template)

# Function to extract time series of EVI values for pixels identified in the target_pixels layer
# extracts from the set of geotifs in geotif_folder with filename starting with geotif_filename and ending in an integer date code, as specified in geotif_date_codes.
# Note this is very slow, since it reads in the whole raster for each time step, but for this reason also requires little memory. 
# Probably should make one that first assembles a rasterbrick, then drills through it to get the time series. 
extract_target_evi <- function(target_pixels, geotif_folder, geotif_filename, geotif_date_codes) {
  r = raster(paste(geotif_folder, geotif_filename, geotif_date_codes[1], ".tif", sep=""))
  evi_crop = crop(r, extent(target_pixels))
  evi_mask = mask(evi_crop, target_pixels)
  evi_stack = stack(evi_mask)
  for (i in 2:length(geotif_date_codes)) {
    r = raster(paste(geotif_folder, geotif_filename, geotif_date_codes[i], ".tif", sep=""))
    evi_crop = crop(r, extent(target_pixels))
    evi_mask = mask(evi_crop, target_pixels)
    evi_stack = stack(evi_stack, evi_mask)
  }
  evi_stack = evi_stack/10000 # rescale to evi scale
  return(evi_stack)
}



target_evi_stack = extract_target_evi(target_pixels, geotif_folder, geotif_filename, date_codes)


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


