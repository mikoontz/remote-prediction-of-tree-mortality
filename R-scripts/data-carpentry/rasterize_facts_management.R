library(devtools)
#devtools::install_github("ecohealthalliance/fasterize") 
library(sf)
library(fasterize)
library(raster)
library(viridis)

# Load project boundary and raster template
raster_template <- raster("features/sierra-nevada-250m-evi-template.tif")

# Read in facts management polygons
#st_layers(dsn = "features/FRAP-fire-perimeters/fire16_1.gdb") # what layers are aviailable in the geodatabase?
facts <- st_read(dsn = "features/FACTS/CA_Activity_merged.shp",stringsAsFactors = FALSE)

# there are a few features that cover the entire state but contain no info. they have blanks in the SUID column. remove them
facts <- facts[facts$SUID != "",]
facts <- facts[!is.na(facts$SUID),]

# Reproject to projection of raster template
facts <- st_transform(facts, crs = proj4string(raster_template))

# get the year out of the date
facts$year_compl <- as.numeric(substr(facts$DATE_COMPL,1,4))

# Disaggregate the raster_template to get a finer resolution
raster_template_fine <- disaggregate(raster_template, fact = c(10, 10))

# Rasterize the facts polygons, taking the most recent year if overlap
facts_target_fine <- fasterize(sf = facts, field = "year_compl", raster = raster_template_fine, fun = "max")

# Get the year of the most recent management within each 250m cell by taking the max of the 100 mini-cells within it
facts_target_coarse <- aggregate(facts_target_fine, fact = c(10, 10), fun = max)

plot(facts_target_coarse, col = viridis(10))

# Export a GeoTiff for use with R
writeRaster(facts_target_coarse, filename = "features/sierra_nevada_250m_most_recent_management.tif",overwrite=TRUE)
