library(sf)
library(fasterize)
library(raster)
library(viridis)
library(rgdal)

# Read in the target cover raster created by script rasterize_calveg_polygons_percent_cover.R
target_cover = raster("features/sierra_nevada_250m_calveg_pipo_forest-cover_whr_type.tif")

# Read in subset polygon that defines the area where the aerial mortality data has been rasterized
area_subset = readOGR("features/so-sierra-subset-mask/so-sierra-subset-mask.shp")

# Read in areal mortality raster for a specified year

mort2015 = raster(paste("features/ADS-rasterized/Y2015_sp122.tif", sep=""))
mort2016 = raster(paste("features/ADS-rasterized/Y2015_sp122.tif", sep=""))
mort = mort2015 + mort2016

proj4string(area_subset) == proj4string(target_cover) # same projection?
proj4string(mort) == proj4string(target_cover)

# Subset the target cover to this area 
target_cover_sub = crop(target_cover, mort)
#plot(target_cover_sub)

# Create a target cover layer for pixels with specified percent PIPO cover
target_cover_sub$cover80 = target_cover_sub >= 80
#plot(target_cover_sub)

dim(mort)
dim(target_cover_sub) # these are at different resolutions
mort_resamp = resample(mort, target_cover_sub$cover80)

# Evaluate 
#par(mfrow=c(1,2))
#plot(mort_resamp); plot(target_cover_sub$cover80) # doesn't look like much overlap 

mort_resamp$mort_bin = as.integer(mort_resamp[[1]]>50)
target_cover_sub$cover80[target_cover_sub$cover80==0] = NA
z = mask(mort_resamp[[1]], target_cover_sub$cover80)
plot(z)


# plot mortality levels in pixels containing >=80% PIPO WHR
plot(target_cover, main="PIPO WHR pct cover")
plot(area_subset, add=T)
plot(z, main="PIPO mortality in high-PIPO areas", col=rainbow(10))
sum(getValues(z)>0, na.rm=T) # 1435 pixels in this ssn subset that are >80% PIPO and have some mortality 
sum(getValues(z)>30, na.rm=T) # 350, pixels in this ssn subset that are >80% PIPO and have "high" mortality 

# Conclusion: there is some variability in mortality in these forests, but maybe not as much mortality as I would have expected (few cells >50 trees per 250m pixel). Need to check with Derek that this looks reasonable when compared to the actual mortality data set and its range of mortality levels. 
