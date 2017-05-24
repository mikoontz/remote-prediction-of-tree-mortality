library(raster)

ads16 <- raster("features/ADS-rasterized/Y2016_sp122.tif")
evi <- raster("features/ts/sn-whole-ts-modis-forest-quality-mask-20000218.tif")

crs(ads16)
crs(evi)

ads16 <- projectRaster(ads16, evi)

identical(crs(ads16), crs(evi))

plot(evi)
plot(ads16, add = TRUE, legend = FALSE)
