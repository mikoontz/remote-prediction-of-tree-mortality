# rm(list = ls())
library(raster)
library(rgdal)
library(sf)
library(devtools)
# devtools::install_github("ecohealthalliance/fasterize")
library(fasterize)

nsn <- st_read(dsn = "features/ExistingVegNorSierra2000_2014_v1.gdb", stringsAsFactors = FALSE) # From DY: I had to remove the slash after the filename for this to work
ssn <- st_read(dsn = "features/ExistingVegSouthSierra2000_2008_v1.gdb", stringsAsFactors = FALSE) # From DY: I had to remove the slash after the filename for this to work

# Subsets by Wildlife Habitat Relationship lifeform (https://www.fs.fed.us/r5/rsl/projects/classification/cv-cwhr-xwalk.html)
# nsn_con_forest <- subset(nsn,
#                   subset =
#                     nsn$WHRLIFEFORM == "WHR_CON" | # Conifer forest/woodland
#                     nsn$WHRLIFEFORM == "WHR_MIX")  # Mixed conifer and hardwood forest/woodland
# 
# ssn_con_forest <- subset(ssn,
#                   subset =
#                     ssn$WHRLIFEFORM == "WHR_CON" | # Conifer forest/woodland
#                     ssn$WHRLIFEFORM == "WHR_MIX")  # Mixed conifer and hardwood forest/woodland
# 
# filename <- "features/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-lifeform.tif"

# Subsets by Wildlife Habitat Relationship type (http://frap.fire.ca.gov/projects/frap_veg/classification)
# The end product will be a raster with each forested cell's center having at least some overlap with a conifer forest polygon

forest.whr.types <- c("SMC",  # Sierra mixed conifer
                      "MCN",  # Mixed conifer
                      "MHC",  # Mixed hardwood-conifer
                      "SCN",  # Subalpine conifer
                      "JPN",  # Jeffrey pine
                      "PPN",  # Ponderosa pine
                      "WFR",  # White fir
                      "RFR",  # Red fir
                      "DFR")  # Douglas-fir

nsn_con_forest <- subset(nsn, subset = nsn$WHRTYPE %in% forest.whr.types)
ssn_con_forest <- subset(ssn, subset = ssn$WHRTYPE %in% forest.whr.types)

sn <- shapefile("features/SierraEcoregion_TNC/SierraEcoregion_TNC.shp")
raster_template <- raster("features/sierra-nevada-250m-evi-template.tif")

nsn_con_forest <- st_transform(nsn_con_forest, crs = proj4string(raster_template))
ssn_con_forest <- st_transform(ssn_con_forest, crs = proj4string(raster_template))

# All pixels whose centers overlap a forested polygon will get a value of 1
nsn_con_forest_r <- fasterize(sf = nsn_con_forest,
                               raster = raster_template)
ssn_con_forest_r <- fasterize(sf = ssn_con_forest,
                               raster = raster_template)

# Combines the North Sierra Nevada and the South Sierra Nevada
sn_con_forest_r <- merge(nsn_con_forest_r, ssn_con_forest_r)

# There are some pixels from the South Sierra Nevada region that are not in the
# SierraEcoregion_TNC polygon, so we can mask those out.
sn_con_forest_r <- mask(sn_con_forest_r, sn)

filename <- "features/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-type.tif"
writeRaster(sn_con_forest_r, filename = filename, format="GTiff", overwrite=TRUE)

# We want the "good" (i.e. forested) pixels to have a value of 1
# and "bad" (i.e. non-forested) pixels to have a value of 0. We can't mask out
# non-forested pixels because Earth Engine asset uploads don't seem to like
# that approach.
sn_con_forest_r[is.na(sn_con_forest_r)] <- 0 # Turn all masked pixels to 0

plot(sn_con_forest_r)
plot(sn, add = TRUE)

filename <- "features/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-type-no-mask.tif"
writeRaster(sn_con_forest_r, filename = filename, format="GTiff", overwrite=TRUE)
