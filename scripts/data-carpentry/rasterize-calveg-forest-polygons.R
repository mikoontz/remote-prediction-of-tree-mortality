library(raster)
library(rgdal)
library(sf)
library(devtools)
# devtools::install_github("ecohealthalliance/fasterize")
library(fasterize)

nsn <- st_read(dsn = "features/ExistingVegNorSierra2000_2014_v1.gdb/", stringsAsFactors = FALSE)
ssn <- st_read(dsn = "features/ExistingVegSouthSierra2000_2008_v1.gdb/", stringsAsFactors = FALSE)

nsn_con <- subset(nsn, subset = nsn$WHRLIFEFOR == "WHR_CON")
ssn_con <- subset(ssn, subset = ssn$WHRLIFEFOR == "WHR_CON")

sn <- shapefile("features/SierraEcoregion_TNC/SierraEcoregion_TNC.shp")
raster_template <- raster("features/sierra-nevada-250m-evi-template.tif")

raster_template[] <- rep(0, times = ncell(raster_template))
raster_template <- mask(raster_template, sn)

nsn_con <- st_transform(nsn_con, crs = proj4string(raster_template))
ssn_con <- st_transform(ssn_con, crs = proj4string(raster_template))

nsn_forest_raster <- fasterize(sf = nsn_con,
                               raster = raster_template)
ssn_forest_raster <- fasterize(sf = ssn_con,
                               raster = raster_template)

sn_forest_raster <- merge(nsn_forest_raster, ssn_forest_raster)
sn_forest_raster <- mask(sn_forest_raster, sn)
sn_forest_raster[sn_forest_raster == 1] <- 0

plot(sn_forest_raster)
plot(sn, add = TRUE)

writeRaster(sn_forest_raster, filename="features/sierra-nevada-250m-calveg-forested-pixels.tif", format="GTiff", overwrite=TRUE)

test <- raster("features/sierra-nevada-250m-calveg-forested-pixels.tif")
plot(test)
