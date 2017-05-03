library(raster)
library(viridis)

test <- stack("features/ts/sn-whole-ts-modis-forest-quality-mask-20000524.tif")
cf <- raster("features/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-type-full-cell.tif")
cf[cf == 0] <- NA 

bandNames <- c("NDVI", 
               "EVI", 
               "sur_refl_b01", 
               "sur_refl_b02", 
               "sur_refl_b03", 
               "sur_refl_b07",
               "ViewZenith",
               "SolarZenith",
               "RelativeAzimuth",
               "DayOfYear",
               "SummaryQA",
               "DetailedQA",
               "date")

names(test) <- bandNames
cf <- spTransform(cf, crs(test))


?spTransform
plot(test[["EVI"]], col = viridis(10))

# Each image needs to be remasked (when all bands are == -1, mask it)
str(test)
extract(test, y = cbind(0, -1e5))
?extract
