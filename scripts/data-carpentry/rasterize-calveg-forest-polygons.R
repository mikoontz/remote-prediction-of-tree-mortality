# rm(list = ls())
library(raster)
library(rgdal)
library(sf)
library(devtools)
# devtools::install_github("ecohealthalliance/fasterize")
library(fasterize)

nsn <- st_read(dsn = "features/ExistingVegNorSierra2000_2014_v1.gdb/", stringsAsFactors = FALSE)
ssn <- st_read(dsn = "features/ExistingVegSouthSierra2000_2008_v1.gdb/", stringsAsFactors = FALSE)

# Subsets by Wildlife Habitat Relationship lifeform (https://www.fs.fed.us/r5/rsl/projects/classification/cv-cwhr-xwalk.html)
nsn_con_forest <- subset(nsn,
                  subset =
                    nsn$WHRLIFEFORM == "WHR_CON" | # Conifer forest/woodland
                    nsn$WHRLIFEFORM == "WHR_MIX")  # Mixed conifer and hardwood forest/woodland

ssn_con_forest <- subset(ssn,
                  subset =
                    ssn$WHRLIFEFORM == "WHR_CON" | # Conifer forest/woodland
                    ssn$WHRLIFEFORM == "WHR_MIX")  # Mixed conifer and hardwood forest/woodland

# Subsets by Wildlife Habitat Relationship type (http://frap.fire.ca.gov/projects/frap_veg/classification)
# nsn_con_forest <- subset(nsn, 
#                   subset = 
#                     nsn$WHRTYPE == "SMC" |  # Sierra mixed conifer
#                     nsn$WHRTYPE == "MCN" |  # Mixed conifer
#                     nsn$WHRTYPE == "MHC" |  # Mixed hardwood-conifer
#                     nsn$WHRTYPE == "SCN" |  # Subalpine conifer
#                     nsn$WHRTYPE == "JPN" |  # Jeffrey pine
#                     nsn$WHRTYPE == "PPN" |  # Ponderosa pine
#                     nsn$WHRTYPE == "WFR" |  # White fir
#                     nsn$WHRTYPE == "RFR" |  # Red fir
#                     nsn$WHRTYPE == "DFR")   # Douglas fir
# 
# ssn_con_forest <- subset(ssn, 
#                   subset = 
#                     ssn$WHRTYPE == "SMC" |  # Sierra mixed conifer
#                     ssn$WHRTYPE == "MCN" |  # Mixed conifer
#                     ssn$WHRTYPE == "MHC" |  # Mixed hardwood-conifer
#                     ssn$WHRTYPE == "SCN" |  # Subalpine conifer
#                     ssn$WHRTYPE == "JPN" |  # Jeffrey pine
#                     ssn$WHRTYPE == "PPN" |  # Ponderosa pine
#                     ssn$WHRTYPE == "WFR" |  # White fir
#                     ssn$WHRTYPE == "RFR" |  # Red fir
#                     ssn$WHRTYPE == "DFR")   # Douglas fir

sn <- shapefile("features/SierraEcoregion_TNC/SierraEcoregion_TNC.shp")
raster_template <- raster("features/sierra-nevada-250m-evi-template.tif")

raster_template[] <- rep(0, times = ncell(raster_template))
raster_template <- mask(raster_template, sn)

nsn_con_forest <- st_transform(nsn_con_forest, crs = proj4string(raster_template))
ssn_con_forest <- st_transform(ssn_con_forest, crs = proj4string(raster_template))

nsn_con_forest_r <- fasterize(sf = nsn_con_forest,
                               raster = raster_template)
ssn_con_forest_r <- fasterize(sf = ssn_con_forest,
                               raster = raster_template)

sn_con_forest_r <- merge(nsn_con_forest_r, ssn_con_forest_r)
sn_con_forest_r <- mask(sn_con_forest_r, sn)
sn_con_forest_r[sn_con_forest_r == 1] <- 0

plot(sn_con_forest_r)
plot(sn, add = TRUE)

writeRaster(sn_con_forest_r, filename="features/sierra-nevada-250m-calveg-conifer-forested-pixels-by-lifeform.tif", format="GTiff", overwrite=TRUE)

test <- raster("features/sierra-nevada-250m-calveg-conifer-forested-pixels-by-lifeform.tif")
plot(test)
