library(sf)
library(fasterize)
library(raster)
library(viridis)
library(rgdal)

# Read in the target cover raster created by script rasterize_calveg_polygons_percent_cover.R
target_cover = raster("features/sierra_nevada_250m_calveg_pipo_forest-cover_whr_type.tif")

# Read in subset polygon that defines the area where the aerial mortality data has been rasterized
area_subset = readOGR("features/so-sierra-subset-mask/so-sierra-subset-mask.shp")

proj4string(area_subset) == proj4string(target_cover) # same projection?

# Subset the target cover to this area 
