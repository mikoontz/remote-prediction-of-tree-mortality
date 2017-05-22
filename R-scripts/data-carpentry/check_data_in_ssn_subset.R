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
yr = 2015
mort = raster(paste("features/ADS-rasterized/Y", yr, "_sp122.tif", sep=""))

proj4string(area_subset) == proj4string(target_cover) # same projection?
proj4string(mort) == proj4string(target_cover)

# Subset the target cover to this area 
target_cover_sub = crop(target_cover, mort)
plot(target_cover_sub)

# Create a target cover layer for pixels with specified percent PIPO cover
target_cover_sub$cover80 = target_cover_sub >= 80
plot(target_cover_sub)

# Evaluate 
par(mfrow=c(1,2))
plot(mort); plot(target_cover_sub$cover80) # doesn't look like much overlap 

dim(mort)
dim(target_cover_sub)

