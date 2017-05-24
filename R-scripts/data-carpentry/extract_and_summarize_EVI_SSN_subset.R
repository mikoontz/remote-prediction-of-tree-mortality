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
mort_template = raster(paste("features/ADS-rasterized/Y2015_sp122.tif", sep=""))
# subset the target cover raster to the mortality area 
target_cover_sub = crop(target_cover, mort_template)

# Create a target cover layer for pixels with specified percent PIPO cover
PIPO_cover_min = 80
target_pixels = target_cover_sub >= PIPO_cover_min
target_pixels[target_pixels==0] = NA # set the non-target values to NA which is default value for masking
#plot(target_pixels)

# Function to extract time series of EVI values for pixels identified in the target_pixels layer
# extracts from the set of geotifs in geotif_folder with filename starting with geotif_filename and ending in an integer date code, as specified in geotif_date_codes.
# Note this is very slow, since it reads in the whole raster for each time step, but for this reason also requires little memory. 
# Probably should make one that first assembles a rasterbrick, then drills through it to get the time series. 
extract_target_evi <- function(target_pixels, geotif_folder, geotif_filename, geotif_date_codes) {
  r = raster(paste(geotif_folder, geotif_filename, geotif_date_codes[1], ".tif", sep=""))
  proj4string(r) == proj4string(target_pixels)
  ## WTF? the coordinates of the geotif seem not to be in lat-lon
  

  evi_crop = crop(r, extent(target_pixels))
  evi_mask = mask(evi_crop, target_pixels)
    
  for (i in 2:length(geotif_date_codes)) {
    r = raster(paste(geotif_folder, geotif_filename, geotif_date_codes[i], ".tif", sep=""))
    evivals = rbind(evivals, extract(r, locs))
  }
  evivals[evivals==1] = NA 
  evivals = evivals/10000 # rescale to evi scale
  return(evivals)
}



## Look at mean EVI response for all PIPO areas in S Sierras
## Then compare to the EVI response for the "high mortality" PIPO areas in S Sierras

evi_average_by_polys <- function(polys, geotif_folder, geotif_filename, geotif_numbers) {
  r = raster(paste(geotif_folder, geotif_filename, geotif_numbers[1], ".tif", sep=""))
  mr = mask(r, polys)
  z = getValues(mr)
  na.index = is.na(z)
  z[z==1] = NA 
  evi_mean = mean(z, na.rm=T)
  for (i in 2:length(geotif_numbers)) {
    r = raster(paste(geotif_folder, geotif_filename, geotif_numbers[i], ".tif", sep=""))
    z = getValues(r)
    z[na.index] = NA
    z[z==1] = NA 
    evi_mean = c(evi_mean, mean(z, na.rm=T))
    if (!i%%10) print(i)
  }
  evi_mean = evi_mean/10000  # rescale to evi scale
  return(evi_mean)
}

evi_cellvalues_by_polys <- function(polys, geotif_folder, geotif_filename, geotif_numbers) {
  r = raster(paste(geotif_folder, geotif_filename, geotif_numbers[1], ".tif", sep=""))
  mr = mask(r, polys)
  z = getValues(mr)
  na.index = is.na(z)
  evi_vals = z[!na.index]
  for (i in 2:length(geotif_numbers)) {
    r = raster(paste(geotif_folder, geotif_filename, geotif_numbers[i], ".tif", sep=""))
    z = getValues(r)
    z[na.index] = NA
    evi_vals = rbind(evi_vals, z)
    if (!i%%10) print(i)
  }
  evi_vals[evi_vals==1] = NA
  evi_vals = evi_vals/10000  # rescale to evi scale
  return(evi_vals)
}

#########################################################
# Extract and look at EVI timeseries for some example locations

# Get locations for Illillouette
locs.ill = read_GE_locs("./features/example_locations/Illillouette_conif_forest_points.txt", header=FALSE)
# Reproject to CRS of the geotifs
r = raster(paste(geotif_folder, geotif_filename, datestr[1], ".tif", sep=""))
locs.ill = SpatialPoints(project(locs.ill, proj4string(r)))
# Check
plot(r)
plot(locs.ill, add=TRUE, col="red")


# Extract EVI for the point locations
evivals_ill = extract_evi(locs.ill, geotif_folder, geotif_filename, geotif_numbers = datestr)


# Display evi for individual pixels
par(mfrow=c(3, 4))
for (i in 1:10) plot(evivals_ill[,i], type="l", lwd=2, col="darkgray", ylim=c(0.5, 0.9))
# note half are blank -- don't fall into a "forested" pixel per GEE mask

# Plot average for the 10 pixels
plot(apply(evivals_ill, 1, mean, na.rm=T), type="l", col="darkgray", lwd=2)
# Interestingly, these Illillouette high-elev forests show increasing green-ness through the summer. Could this be a signal of temperature limitation? 
# But appears to be overall lower EVI in the drought years. 


# Try comparing heuristically to points from lower-elevation confer forest stands

locs.shaver = read_GE_locs("./features/example_locations/Shaver_Lake_conifer_forest_points.txt", header=FALSE)

# Reproject to CRS of the geotifs
r = raster(paste(geotif_folder, geotif_filename, datestr[1], ".tif", sep=""))
locs.shaver = SpatialPoints(project(locs.shaver, proj4string(r)))

# Check
plot(r)
plot(locs.shaver, add=TRUE, col="red")

# Extract EVI
evivals_shaver = extract_evi(locs.shaver, geotif_folder, geotif_filename, geotif_numbers = datestr)

# Plot individual pixels
par(mfrow=c(3, 4))
for (i in 1:10) {
  plot(evivals_shaver[,i], type="l", lwd=2, col="darkgray", ylim=c(0.5, 0.9))
}

# Plot average
plot(apply(evivals_shaver, 1, mean, na.rm=T),  type="l", col="red")


# Compare the means of the Illillouette to the Shaver Lake points

par(mfrow=c(1,1))
datecols = rep("white", length(dates)); datecols[which(dates$mon %in% c(6, 7, 8, 9))] = "cyan4"
plot(apply(evivals_ill, 1, mean, na.rm=T), type="p", col=datecols, lwd=2, ylim = c(0.5, 0.85))
datecols = rep("white", length(dates)); datecols[which(dates$mon %in% c(6, 7, 8, 9))] = "orange3"
points(apply(evivals_shaver, 1, mean, na.rm=T),  col=datecols, lwd=2)

# Both locations show substantial greening up from spring to fall. Clearly the higher elevation sites have lower mean EVI. Hard to tell from this small sample, but it looks like the higher-elevation sites also have lower amount of greening. 
# But during 2015-16, Shaver Lake areas crashed to below Illillouette levels. 

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



#### Not updated from here on! Broken. (Andrew 5/17/17)


## Next look at EVI values for all cells within some subsets of pixels 
evi_PIPO = evi_average_by_polys(polys=v.PIPO, geotif_folder, geotif_filename, datestr)
cols = rep("darkgreen", length(evi_PIPO))
cols[dates$mon %in% c(1,2,3,4,5,11, 12)] = "white"
plot(evi_PIPO, type="p", lwd=2, col=cols)
x = 4:8
plot(evi_PIPO[dates$year==100 & dates$mon %in% x]~dates$yday[dates$year==100 & dates$mon %in% x], type="l", ylim=c(0.6, 0.8), lwd=2, col="darkgray", ylab="EVI", xlab = "julian date", main="Growing season EVI for PIPO areas in S. Sierras")
for (i in 100:112) {
  lines(evi_PIPO[dates$year==i & dates$mon %in% x]~dates$yday[dates$year==i & dates$mon %in% x], lwd=2, col="cyan4")
}
for (i in 113:116) {
  lines(evi_PIPO[dates$year==i & dates$mon %in% x]~dates$yday[dates$year==i & dates$mon %in% x], lwd=2, col="orange3")
}
legend("topleft", c("2000-2012", "2013-2016"), col=c("cyan4", "orange3"), lwd=c(2,2))

evi_late = evi_PIPO[dates$yday==288]
evi_diff = evi_late-evi_PIPO[dates$yday == 144]
plot(evi_late)
plot(evi_diff)

# Q how to get the "less affected" parts of this? Separate by Derek's high mortality polygons (high = within those polygons, lower = elsewhere?)

# Since 2013 was the year where the EVI really dropped from normal to drought levels, check what the variation was like among grid cells that year. 
evi_PIPO_cells = evi_cellvalues_by_polys(polys=v.PIPO, geotif_folder, geotif_filename, datestr[dates$year==113])
z = apply(evi_PIPO_cells, 2, f<-function(x){return(sum(!is.na(x)))})
evi_PIPO_cells = evi_PIPO_cells[,z>12]
dim(evi_PIPO_cells)
x=(8:17)
days = x*16
plot(days, evi_PIPO_cells[x,1], type="l", ylim=c(0.5, 0.85), col="white")
palette(tim.colors(100))
for(i in 1:50) {
  z = sample(1:ncol(evi_PIPO_cells), size=1)
  cols = coef(lm(evi_PIPO_cells[x,z]~days))
  lines(days, evi_PIPO_cells[x,z], col=max(c(30,cols[2]*16000+40)), lwd=2)
  #abline(cols, col="lightgray")
}

