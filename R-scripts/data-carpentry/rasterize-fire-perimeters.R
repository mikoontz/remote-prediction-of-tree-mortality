library(sf)
library(fasterize)
library(raster)
library(viridis)

# Load project boundary and raster template
sn <- shapefile("features/SierraEcoregion_TNC/SierraEcoregion_TNC.shp")
raster_template <- raster("features/sierra-nevada-250m-evi-template.tif")

# Read in fire perimeters layers
#st_layers(dsn = "features/FRAP-fire-perimeters/fire16_1.gdb") # what layers are aviailable in the geodatabase?
fire.perims <- st_read(dsn = "features/FRAP-fire-perimeters/fire16_1.gdb", layer="firep16_1", stringsAsFactors = FALSE)
rxburn.perims <- st_read(dsn = "features/FRAP-fire-perimeters/fire16_1.gdb", layer="rxburn16_1", stringsAsFactors = FALSE)

#all we need from each is the year; trim to that column and then merge fires and prescribed burns
all.perims <- rbind(fire.perims[,"YEAR_"],rxburn.perims[,"YEAR_"])
all.perims$YEAR_ <- as.numeric(all.perims$YEAR_)

# Reproject to projection of raster template
all.perims <- st_transform(all.perims, crs = proj4string(raster_template))

# Disaggregate the raster_template to get a finer resolution
raster_template_fine <- disaggregate(raster_template, fact = c(10, 10))

# Rasterize the fire polygons, taking the most recent year
perims_target_fine <- fasterize(sf = all.perims, field = "YEAR_", raster = raster_template_fine, fun = "max")

# Get the year of the most recent fire within each 250m cell by taking the max of the 100 mini-cells within it
perims_target_coarse <- aggregate(perims_target_fine, fact = c(10, 10), fun = max)

plot(perims_target_coarse, col = viridis(10))

# Export a GeoTiff for use with R
writeRaster(perims_target_coarse, filename = "features/sierra_nevada_250m_most_recent_fire.tif")