library(raster)

# What should the projections be?
evi_template <- raster("features/ts/sn-whole-ts-modis-forest-quality-mask-20000218.tif")
california_albers <- crs(evi_template)

fileNames <- list.files(path = "features/ADS-rasterized/", full.names = TRUE)
rasterList <- lapply(X = fileNames, FUN = raster)

newFileNames <- gsub(pattern = ".tif", replacement = "_ca.tif", x = fileNames)

for (i in seq_along(rasterList)) {
  tmp <- projectRaster(rasterList[[i]], evi_template)
  writeRaster(x = tmp, filename = newFileNames[i])
}
