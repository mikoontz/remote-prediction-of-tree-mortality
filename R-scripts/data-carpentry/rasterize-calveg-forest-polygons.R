# rm(list = ls())
library(raster)
library(rgdal)
library(sf)
library(devtools)
# devtools::install_github("ecohealthalliance/fasterize")
library(fasterize)

library(gdalUtils)

#Need these directories to output large intermediate data files, but they're ignored in .gitignore so they aren't populated to start, so create them if they don't exist
dir.create("features/intermediate products")
dir.create("features/intermediate products/CALVEG shapefiles")

# Define conifer WHR types (http://frap.fire.ca.gov/projects/frap_veg/classification)
con.whr.types <- c("SMC",  # Sierra mixed conifer
                   "MCN",  # Mixed conifer
                   "MHC",  # Mixed hardwood-conifer
                   "SCN",  # Subalpine conifer
                   "JPN",  # Jeffrey pine
                   "PPN",  # Ponderosa pine
                   "WFR",  # White fir
                   "RFR",  # Red fir
                   "DFR")  # Douglas-fir

# Read in CALVEG layers
nsn <- st_read(dsn = "features/ExistingVegNorSierra2000_2014_v1.gdb", stringsAsFactors = FALSE) # From DY: I had to remove the slash after the filename for this to work
ssn <- st_read(dsn = "features/ExistingVegSouthSierra2000_2008_v1.gdb", stringsAsFactors = FALSE) # From DY: I had to remove the slash after the filename for this to work

# Keep only the relevant attributes
nsn <- nsn[,c("WHRTYPE","WHRLIFEFORM")]
ssn <- ssn[,c("WHRTYPE","WHRLIFEFORM")]

# Load project boundary and raster template
sn <- shapefile("features/SierraEcoregion_TNC/SierraEcoregion_TNC.shp")
raster_template <- raster("features/sierra-nevada-250m-evi-template.tif")

# Reproject CALVEG layers to projection of raster template
nsn <- st_transform(nsn, crs = proj4string(raster_template))
ssn <- st_transform(ssn, crs = proj4string(raster_template))

# Create new column that indicates whether a polygon contains conifer forest (1 if yes, 2 if no)
nsn$con_forest <- ifelse(nsn$WHRTYPE %in% con.whr.types,1,2)
ssn$con_forest <- ifelse(ssn$WHRTYPE %in% con.whr.types,1,2)

# Create another shapefile that contains only non-conifer polygons
nsn.nonconifer <- nsn[nsn$con_forest != 1,]
ssn.nonconifer <- ssn[ssn$con_forest != 1,]

# Write to disk in order to rasterize using gdal_rasterize in next step
# But first need to delete the files if they already exist because st_write has some problems with overwriting
do.call(file.remove, list(list.files("features/intermediate products/CALVEG shapefiles", full.names = TRUE)))

st_write(nsn,"features/intermediate products/CALVEG shapefiles/CALVEG_nsn.shp") #annoyingly this puts a second ".shp" after the filename
st_write(nsn.nonconifer,"features/intermediate products/CALVEG shapefiles/CALVEG_nsn_nonconifer.shp")
st_write(ssn,"features/intermediate products/CALVEG shapefiles/CALVEG_ssn.shp")
st_write(ssn.nonconifer,"features/intermediate products/CALVEG shapefiles/CALVEG_ssn_nonconifer.shp")

# Get resolution and extent of template raster (needed for rasterization)
template.res <- res(raster_template)
template.extent <- extent(raster_template)[c(1,3,2,4)]

# Make raster indicating whether the center of each cell overlaps a polygon that is a conifer type (1 if yes; 2 if no; nodata if no polygon overlap)
nsn.raster <- gdal_rasterize("features/intermediate products/CALVEG shapefiles/CALVEG_nsn.shp","features/intermediate products/calveg_conifer_nsn.tif",
                              a="con_forest", tr=template.res, te=template.extent,
                               l="CALVEG_nsn",a_nodata=NA,verbose=TRUE,output_Raster=TRUE)
ssn.raster <- gdal_rasterize("features/intermediate products/CALVEG shapefiles/CALVEG_ssn.shp","features/intermediate products/calveg_conifer_ssn.tif",
                             a="con_forest", tr=template.res, te=template.extent,
                             l="CALVEG_ssn",a_nodata=NA,verbose=TRUE,output_Raster=TRUE)

# Make a raster indicating whether ANY PART of the cell overlaps a polygon that is a non-forest type (2 if this is the case; nodata if it is not)
nsn.raster.nonconifer <- gdal_rasterize("features/intermediate products/CALVEG shapefiles/CALVEG_nsn_nonconifer.shp","features/intermediate products/calveg_nonconifer_nsn.tif",
                             a="con_forest", tr=template.res, te=template.extent, at=TRUE,
                             l="CALVEG_nsn_nonconifer",a_nodata=NA,verbose=TRUE,output_Raster=TRUE)
ssn.raster.nonconifer <- gdal_rasterize("features/intermediate products/CALVEG shapefiles/CALVEG_ssn_nonconifer.shp","features/intermediate products/calveg_nonconifer_ssn.tif",
                                        a="con_forest", tr=template.res, te=template.extent, at=TRUE,
                                        l="CALVEG_ssn_nonconifer",a_nodata=NA,verbose=TRUE,output_Raster=TRUE)

# Merge both conifer rasters (north and south Sierra)
sn.raster <- merge(nsn.raster,ssn.raster)

# Merge both nonconifer rasters (north and south Sierra); where they overlap, take the maximum (this means set a cell to "nonconifer" if either the north or south raster said nonconifer)
sn.raster.nonconifer <- mosaic(nsn.raster.nonconifer,ssn.raster.nonconifer,fun=max)

# Identify cells where the center of the cell overlapped a conifer polygon, and no part of the cell overlapped a non-conifer polygon
sn_con_forest_r <- (sn.raster == 1) & (is.na(sn.raster.nonconifer))

# Mask out the non-conifer pixels (they are at this point assigned 0; set them to NA)
sn_con_forest_r[sn_con_forest_r == 0] <- NA


# There are some pixels from the South Sierra Nevada region that are not in the
# SierraEcoregion_TNC polygon, so we can mask those out.
sn_con_forest_r <- mask(sn_con_forest_r, sn)

filename <- "features/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-type_full-cell.tif"
writeRaster(sn_con_forest_r, filename = filename, format="GTiff", overwrite=TRUE)

# We want the "good" (i.e. forested) pixels to have a value of 1
# and "bad" (i.e. non-forested) pixels to have a value of 0. We can't mask out
# non-forested pixels because Earth Engine asset uploads don't seem to like
# that approach.
sn_con_forest_r[is.na(sn_con_forest_r)] <- 0 # Turn all masked pixels to 0

plot(sn_con_forest_r)
plot(sn, add = TRUE)

filename <- "features/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-type-no-mask_full-cell.tif"
writeRaster(sn_con_forest_r, filename = filename, format="GTiff", overwrite=TRUE)
