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
plot(target_cover_sub, col=rev(viridis(10)))

# Create a target cover layer for pixels with specified percent PIPO cover
target_cover_sub$cover80 = target_cover_sub >= 80
plot(target_cover_sub, col=rev(viridis(10)))

# Evaluate 
par(mfrow=c(1,2))
plot(mort); plot(target_cover_sub$cover80) 

mort$mort_bin = as.integer(mort[[1]]>=80)
target_cover_sub$cover80[target_cover_sub$cover80==0] = NA
z = mask(mort[[1]], target_cover_sub$cover80)
plot(z, col=rev(viridis(10)))


# plot mortality levels in pixels containing >=80% PIPO WHR
plot(target_cover, main="PIPO WHR pct cover")
plot(area_subset, add=T)
plot(z, main="PIPO mortality in high-PIPO areas", col=rev(viridis(10)))
mortvals = getValues(z)
sum(!is.na(mortvals))
sum(mortvals==0, na.rm=T) # 917 pixels in this ssn subset that are >80% PIPO and have some mortality 
sum(mortvals>0 & mortvals <15, na.rm=T) # 415 pixels in this ssn subset that are >80% PIPO and have low mortality 
sum(mortvals>=15 & mortvals<40, na.rm=T) # 675, pixels in this ssn subset that are >80% PIPO and have "medium" mortality 
sum(getValues(z)>=40, na.rm=T) # 226, pixels in this ssn subset that are >80% PIPO and have "very high" mortality 

# Conclusion: Seems like good range of mortality levels in this area within the PIPO forest type. Total number of pixels is relatively small 2223, but within this number, there is good stratification across mortality levels. 

# I think we can use this for testing statistical models to relate EVI to mortality and environmental factors. 
