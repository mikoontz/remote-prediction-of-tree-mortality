

# Exploratory look at the MODIS geotifs that Mike extracted using Google Earth Engine

library(sp); library(raster); library(rgdal); library(lattice); library(lme4)


# function to convert strings from Google Earth in format "37°17'23.2"N" to numeric decimal degrees
degminsec2dig <- function(x) {
  tmp = strsplit(x, split="°")
  deg = as.numeric(tmp[[1]][1])
  tmp2 = strsplit(tmp[[1]][2], split="'")
  min = as.numeric(tmp2[[1]][1])
  sec = as.numeric(strsplit(tmp2[[1]][2], "\"")[[1]][1])
  return(deg + min/60 + sec/3600)
}

# Function to read in locations from Google Earth and return coordinates
# Assumes these are in 2-column format with Latitude in the first column, Lon in the second
read_GE_locs <- function(filename, header) {
  rawlocs = read.table(filename, header=header, stringsAsFactors = FALSE)
  lat = sapply(rawlocs[,1], degminsec2dig)
  lon = -sapply(rawlocs[,2], degminsec2dig)
  return(as.matrix(cbind(lon, lat)))
}

# Function to extract time series of EVI values, for SpatialPoints identified in locs, from the set of geotifs in geotif_folder with filename starting with geotif_filename and ending in an integer (set of these integers is in geotif_numbers, per Mike Koontz's file naming system)
# Note this is very slow, since it reads in the whole raster for each time step, but for this reason also requires little memory
extract_evi <- function(locs, geotif_folder, geotif_filename, geotif_numbers) {
  r = raster(paste(geotif_filename, geotif_numbers[1], ".tif", sep=""))
  evivals = extract(r, locs)
  for (i in 2:length(geotif_numbers)) {
    r = raster(paste(geotif_filename, geotif_numbers[i], ".tif", sep=""))
    evivals = rbind(evivals, extract(r, locs))
  }
  evivals[evivals==1] = NA 
  evivals = evivals/10000 # rescale to evi scale
  return(evivals)
}


# Enter the locations of files to work with 

geotif_folder = "./features/ee_sierra-nevada-forest-quality-mask-modis-time-series/"
geotif_filename =  "sn-whole-ts-modis-forest-quality-mask-"
tmp = dir(geotif_folder)
filenames = tmp[grep(geotif_filename, tmp)]
loc_filename = ""


## Look at mean EVI response for all PIPO areas in S Sierras
## Then compare to the EVI response for the "high mortality" PIPO areas in S Sierras

evi_average_by_area <- function(polys, geotif_folder, geotif_filename, geotif_numbers) {
  r = raster(paste(geotif_filename, geotif_numbers[1], ".tif", sep=""))
  mr = mask(r, polys)
  z = getValues(mr)
  na.index = is.na(z)
  z[z==1] = NA 
  evi_mean = mean(z, na.rm=T)
  for (i in 2:length(geotif_numbers)) {
    r = raster(paste(geotif_filename, geotif_numbers[i], ".tif", sep=""))
    z = getValues(r)
    z[na.index] = NA
    z[z==1] = NA 
    evi_mean = c(evi_mean, mean(z, na.rm=T))
    if (!i%%10) print(i)
  }
  evi_mean = evi_mean/10000  # rescale to evi scale
  return(evi_mean)
}

evi_cellvalues_by_area <- function(polys, geotif_folder, geotif_filename, geotif_numbers) {
  r = raster(paste(geotif_filename, geotif_numbers[1], ".tif", sep=""))
  mr = mask(r, polys)
  z = getValues(mr)
  na.index = is.na(z)
  evi_vals = z[!na.index]
  for (i in 2:length(geotif_numbers)) {
    r = raster(paste(geotif_filename, geotif_numbers[i], ".tif", sep=""))
    z = getValues(r)
    z[na.index] = NA
    evi_vals = rbind(evi_vals, z)
    if (!i%%10) print(i)
  }
  evi_vals[evi_vals==1] = NA
  evi_vals = evi_vals/10000  # rescale to evi scale
  return(evi_vals)
}

evi_PIPO = evi_average_by_area(polys=v.PIPO, geotif_folder, geotif_filename, datestr)
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
evi_PIPO_cells = evi_cellvalues_by_area(polys=v.PIPO, geotif_folder, geotif_filename, datestr[dates$year==113])
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
# Get data information
#meta = read.csv("metadata0_49.csv")
#dates = strptime(meta$date, "%Y%m%d")
datestr = sapply(filenames, substr, start=39, stop=46)
dates = strptime(datestr, "%Y%m%d")

# Get locations for Illillouette
locs.ill = read_GE_locs("Illillouette_conif_forest_points.txt", header=FALSE)
# Reproject to CRS of the geotifs
r = raster(paste(geotif_filename, datestr[1], ".tif", sep=""))
locs.ill = SpatialPoints(project(locs.ill, proj4string(r)))
# Check
plot(r)
plot(locs.ill, add=TRUE, col="red")


# Extract EVI for the point locations
evivals_ill = extract_evi(locs.ill, geotif_folder, geotif_filename, geotif_numbers = datestr)


# Display evi for individual pixels
par(mfrow=c(3, 4))
for (i in 1:10) {
  plot(evivals_ill[,i], type="l", lwd=2, col="darkgray", ylim=c(0.5, 0.9))
  #lines(c(15, 15), c(0.5, 0.9), col="black")
  #lines(c(43, 43), c(0.5, 0.9), col="black")
}

# Plot average for the 10 pixels
plot(apply(evivals_ill, 1, mean, na.rm=T), type="l", col="darkgray", lwd=2)
# interestingly, these Illillouette high-elev forests show increasing green-ness through the summer. Could this be a signal of temperature limitation? 


# Try comparing heuristically to points from lower-elevation confer forest stands

locs.shaver = read_GE_locs("Shaver_Lake_conifer_forest_points.txt", header=FALSE)

# Reproject to CRS of the geotifs
r = raster(paste(geotif_filename, datestr[1], ".tif", sep=""))
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
  #lines(c(15, 15), c(0.5, 0.9), col="black")
  #lines(c(43, 43), c(0.5, 0.9), col="black")
}

# Plot average
plot(apply(evivals_shaver, 1, mean, na.rm=T),  type="l", col="red")


# compare the means of the illillouette to the shaver lake points

par(mfrow=c(1,1))
datecols = rep("white", length(dates)); datecols[which(dates$mon %in% c(6, 7, 8, 9))] = "cyan4"
plot(apply(evivals_ill, 1, mean, na.rm=T), type="p", col=datecols, lwd=2, ylim = c(0.5, 0.85))
datecols = rep("white", length(dates)); datecols[which(dates$mon %in% c(6, 7, 8, 9))] = "orange3"
points(apply(evivals_shaver, 1, mean, na.rm=T),  col=datecols, lwd=2)

# Both locations show substantial greening up from spring to fall. Clearly the higher elevation sites have lower mean EVI. Hard to tell from this small sample, but it looks like the higher-elevation sites also have lower amount of greening. 
# Question: how does the pattern look during severe drought? 
# How variable is it across years in general? 

# look at annual variation for Shaver Lake points
plotslopes <- function(z, dates) {
  zsub = z[z$mon %in% 6:9,]
  m=lmer(evi~daystd + (1+daystd|year), data=zsub)
  coefs = coef(m)$year
  plot.new(); curve(coefs[1,1]+coefs[1,2]*x, from=0, to=1.5, ylim=c(0.5, 0.9))
  for(i in 1:length(unique(dates$year))) curve(coefs[i,1]+coefs[i,2]*x, add=T, col=i)
  points(evi~daystd, zsub, col="gray")
}

z = data.frame(evi=apply(evivals_shaver, 1, mean, na.rm=T),  year=dates$year, mon=dates$mon, day=dates$yday, daystd = scale(dates$yday))
plotslopes(z, dates)

z = data.frame(evi=apply(evivals_ill, 1, mean, na.rm=T), year=dates$year, mon=dates$mon, day=dates$yday, daystd = scale(dates$yday))
plotslopes(z, dates)

# Annual precip for Fresno county from PRISM
z = data.frame(evi=apply(evivals_shaver, 1, mean, na.rm=T),  year=dates$year, mon=dates$mon, day=dates$yday, daystd = scale(dates$yday))
zsub = z[z$mon %in% 6:9,]
m=lmer(evi~daystd + (1+daystd|year), data=zsub)
coefs = coef(m)$year
annualppt  = c(16.25, 13.24, 8.42, 10.35, 10.8, 14.1, 15.58, 7.02, 7.53, 8.82, 19.72, 10.67, 9.87, 3.24)
barplot(coefs[,2], names.arg= as.character(2000:2013))
barplot(annualppt, names.arg= as.character(2000:2013))
plot(coefs[,2]~annualppt); cor(coefs[,2],annualppt)
# mild positive correlation


z = data.frame(evi=portlo_mean,  year=dates$year, mon=dates$mon, day=dates$yday, daystd = scale(dates$yday))
plotslopes(z, dates)
zsub = z[z$mon %in% 6:9,]
m=lmer(evi~daystd + (1+daystd|year), data=zsub)
coefs = coef(m)$year
annualppt  = c(16.25, 13.24, 8.42, 10.35, 10.8, 14.1, 15.58, 7.02, 7.53, 8.82, 19.72, 10.67, 9.87, 3.24)
par(mfrow=c(2, 1))
barplot(coefs[,2], names.arg= as.character(2000:2013))
barplot(annualppt, names.arg= as.character(2000:2013))
plot(coefs[],2]~annualppt); cor(coefs[,2],annualppt)
# no correlation

# High vs low mortality near Porterville S Sierras, ~6000 feet
locs.porthi= read_GE_locs("Porterville_Hi_Mort.txt", header=F)
locs.porthi = SpatialPoints(project(locs.porthi, proj4string(r)))
evivals_porthi = extract_evi(locs.porthi, geotif_folder, geotif_filename, geotif_numbers = datestr)

locs.portlo= read_GE_locs("Porterville_Lomort.txt", header=F)
locs.portlo = SpatialPoints(project(locs.portlo, proj4string(r)))
evivals_portlo = extract_evi(locs.portlo, geotif_folder, geotif_filename, geotif_numbers = datestr)


par(mfrow=c(3, 4))
for (i in 1:10) {
  plot(evivals_portlo[,i], type="l", lwd=2, col="cyan3", ylim=c(0.5, 0.9))
  plot(evivals_porthi[,i], type="l", lwd=2, col="orange2", ylim=c(0.5, 0.9))
}

## ISSUE -- almost all these points are masked out for the entire time series!! 



# Compare smoothed means across sites

porthi_mean = apply(evivals_porthi, 1, mean, na.rm=T)
portlo_mean = apply(evivals_portlo, 1, mean, na.rm=T)
shaver_mean = apply(evivals_shaver, 1, mean, na.rm=T)
ill_mean  = apply(evivals_ill, 1, mean, na.rm=T)
plot(porthi_mean, type="l", col="red")
lines(portlo_mean,  col="blue")

x = 1:length(datestr) # choose which dates -- 4-14 is April-Sept
plot(loess(porthi_mean[x]~x), type="l", lwd=2, col="red", ylim=c(0.65, 0.85), main="Smoothed EVI April-Sept 2000")
lines(loess(portlo_mean[x]~x), type="l", lwd=2, col="blue")
lines(loess(shaver_mean[x]~x), type="l", lwd=2, col="darkgreen")
#lines(loess(ill_mean[1:20]~x), type="l", lwd=2, col="darkgray")
legend("topleft", c("Shaver Lake", "Porterville Hi Mort", "Porterville Lo Mort"), lwd=rep(2, 3), col=c("darkgreen", "red", "blue") )

plot(porthi_mean~portlo_mean); abline(0,1) # for these few pixels, "high" mortality pixels have less total variation than "low" mortality pixels


