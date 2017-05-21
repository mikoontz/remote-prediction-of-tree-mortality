# rm(list = ls())
library(raster)
library(rgdal)
library(sf)
library(devtools)
# devtools::install_github("ecohealthalliance/fasterize")
library(fasterize)

library(gdalUtils)

#Need these directories to output large intermediate data files, but they're ignored in .gitignore so they aren't populated to start, so create them if they don't exist
intermProd <- "features/intermediate-products"
intermProdSubDir <- paste(intermProd, "CALVEG-shapefiles", sep = "/")

if (!dir.exists(intermProd)) {
  dir.create(intermProd)
  if (!dir.exists(intermProdSubDir)) {
    dir.create(intermProdSubDir)
  }
}

# Define conifer WHR types (http://frap.fire.ca.gov/projects/frap_veg/classification)
# con.whr.types <- c("SMC",  # Sierra mixed conifer
#                    "MCN",  # Mixed conifer
#                    "MHC",  # Mixed hardwood-conifer
#                    "SCN",  # Subalpine conifer
#                    "JPN",  # Jeffrey pine
#                    "PPN",  # Ponderosa pine
#                    "WFR",  # White fir
#                    "RFR",  # Red fir
#                    "DFR")  # Douglas-fir

# Use just Ponderosa pine WHR Type
con.whr.types <- "PPN"  

# Filenames to change depending on how forest subset works
north_SN_targetWHR_fileName <- "CALVEG_nsn_pipo"
north_SN_nonTargetWHR_fileName <- "CALVEG_nsn_non_pipo"

south_SN_targetWHR_fileName <- "CALVEG_ssn_pipo"
south_SN_nonTargetWHR_fileName <- "CALVEG_ssn_non_pipo"

# File directories to write for the several files that comprise each shapefile
filesToWrite <- c(paste(intermProdSubDir, north_SN_targetWHR_fileName, sep = "/"),
                  paste(intermProdSubDir, north_SN_nonTargetWHR_fileName, sep = "/"),
                  paste(intermProdSubDir, south_SN_targetWHR_fileName, sep = "/"),
                  paste(intermProdSubDir, south_SN_nonTargetWHR_fileName, sep = "/"))

# Target spp shape file path and (to write) .tif file path
nsn_target_shp <- paste0(intermProdSubDir, "/", north_SN_targetWHR_fileName, "/", north_SN_targetWHR_fileName, ".shp")
nsn_target_tif <- paste0(intermProd, "/", north_SN_targetWHR_fileName, ".tif")

ssn_target_shp <- paste0(intermProdSubDir, "/", south_SN_targetWHR_fileName, "/", south_SN_targetWHR_fileName, ".shp")
ssn_target_tif <- paste0(intermProd, "/", south_SN_targetWHR_fileName, ".tif")

# Non-target spp shape file path and (to write) .tif file path
nsn_nonTarget_shp <- paste0(intermProdSubDir, "/", north_SN_nonTargetWHR_fileName, "/", north_SN_nonTargetWHR_fileName, ".shp")
nsn_nonTarget_tif <- paste0(intermProd, "/", north_SN_nonTargetWHR_fileName, ".tif")

ssn_nonTarget_shp <- paste0(intermProdSubDir, "/", south_SN_nonTargetWHR_fileName, "/", south_SN_nonTargetWHR_fileName, ".shp")
ssn_nonTarget_tif <- paste0(intermProd, "/", south_SN_nonTargetWHR_fileName, ".tif")

maskedFilename <- "features/sierra-nevada-250m-calveg-pipo-forested-pixels-by-whr-type-full-cell.tif"
nonMaskedFilename <- "features/sierra-nevada-250m-calveg-pipo-forested-pixels-by-whr-type-no-mask-full-cell.tif"

# Should new intermediate files be written at all? Or will we just read the
# intermediate features from the script directly because the long time step
# of creating them has already been done. Note this is a different question than
# whether previously created files should be overwritten
newFiles <- FALSE

# Load project boundary and raster template
sn <- shapefile("features/SierraEcoregion_TNC/SierraEcoregion_TNC.shp")
raster_template <- raster("features/sierra-nevada-250m-evi-template.tif")

if (newFiles) {
  # Should intermediate files be overwritten? Will take lots more time, but critical to
  # set to TRUE if any changes will be made to "update" a file with the same name
  overwrite <- FALSE
  
  # Read in CALVEG layers
  nsn <- st_read(dsn = "features/ExistingVegNorSierra2000_2014_v1.gdb", stringsAsFactors = FALSE) # From DY: I had to remove the slash after the filename for this to work
  ssn <- st_read(dsn = "features/ExistingVegSouthSierra2000_2008_v1.gdb", stringsAsFactors = FALSE) # From DY: I had to remove the slash after the filename for this to work
  
  # Keep only the relevant attributes
  nsn <- nsn[,c("WHRTYPE","WHRLIFEFORM")]
  ssn <- ssn[,c("WHRTYPE","WHRLIFEFORM")]
  
  # Reproject CALVEG layers to projection of raster template
  nsn <- st_transform(nsn, crs = proj4string(raster_template))
  ssn <- st_transform(ssn, crs = proj4string(raster_template))
  
  # Create new column that indicates whether a polygon contains conifer forest (1 if yes, 2 if no)
  nsn$con_forest <- ifelse(test = nsn$WHRTYPE %in% con.whr.types, yes = 1, no = 2)
  ssn$con_forest <- ifelse(test = ssn$WHRTYPE %in% con.whr.types, yes = 1, no = 2)
  
  # Create another shapefile that contains only non-conifer polygons
  nsn.nonconifer <- nsn[nsn$con_forest != 1, ]
  ssn.nonconifer <- ssn[ssn$con_forest != 1, ]
  
  # Write to disk in order to rasterize using gdal_rasterize in next step
  objsToWrite <- list(nsn, nsn.nonconifer, ssn, ssn.nonconifer)
  names(objsToWrite) <- filesToWrite
  
  # If the file exists already, and you want to make a change, delete them and
  # then rewrite to avoid trouble with st_write() on already existing files
  
  if (overwrite)
    lapply(filesToWrite, FUN = file.remove)
  
  if (any(sapply(filesToWrite, FUN = file.exists)) & !overwrite) {
    print(filesToWrite)
    warning("Some of these files already exist, but overwrite was set to FALSE. 
      If you've made changes to the files and want to overwrite them, change the
      overwrite variable to TRUE so that the st_write() function works. 
      If you don't need to make any changes to the files, and can just read
      them in from your local disk, set newFiles to FALSE.")
    }
  
  # Which objects have files that do not already exist? Only write new files 
  # when they don't already exist
  idx <- which(!sapply(filesToWrite, FUN = file.exists))
  
  if (length(idx > 0)) {
    for (i in idx) {
      if (!file.exists(names(objsToWrite[i]))) {
        st_write(obj = objsToWrite[[i]], 
                 dsn = names(objsToWrite[i]),
                 driver = "ESRI Shapefile")
      }
    }
  }

} # End newFiles

# Get resolution and extent of template raster (needed for rasterization)
template.res <- res(raster_template)
template.extent <- extent(raster_template)[c(1,3,2,4)]

# Make raster indicating whether the center of each cell overlaps a polygon that is a conifer type (1 if yes; 2 if no; nodata if no polygon overlap)
nsn.raster <- gdal_rasterize(nsn_target_shp, nsn_target_tif,
                             a = "con_forest", tr = template.res, te = template.extent,
                             l = north_SN_targetWHR_fileName, a_nodata = NA, verbose = TRUE, output_Raster = TRUE)

ssn.raster <- gdal_rasterize(ssn_target_shp, ssn_target_tif,
                             a = "con_forest", tr = template.res, te = template.extent,
                             l = south_SN_targetWHR_fileName, a_nodata = NA, verbose = TRUE, output_Raster = TRUE)

# Make a raster indicating whether ANY PART of the cell overlaps a polygon that is a non-forest type (2 if this is the case; nodata if it is not)
nsn.raster.nonconifer <- gdal_rasterize(nsn_nonTarget_shp, nsn_nonTarget_tif,
                                        a = "con_forest", tr = template.res, te = template.extent, at = TRUE,
                                        l = north_SN_nonTargetWHR_fileName, a_nodata = NA, verbose = TRUE, output_Raster = TRUE)

ssn.raster.nonconifer <- gdal_rasterize(ssn_nonTarget_shp, ssn_nonTarget_tif,
                                        a = "con_forest", tr = template.res, te = template.extent, at = TRUE,
                                        l = south_SN_nonTargetWHR_fileName, a_nodata = NA, verbose = TRUE, output_Raster = TRUE)

# Merge both conifer rasters (north and south Sierra)
sn.raster <- merge(nsn.raster, ssn.raster)

# Merge both nonconifer rasters (north and south Sierra); where they overlap, take the maximum (this means set a cell to "nonconifer" if either the north or south raster said nonconifer)
sn.raster.nonconifer <- mosaic(nsn.raster.nonconifer, ssn.raster.nonconifer,fun=max)

# Identify cells where the center of the cell overlapped a conifer polygon, and no part of the cell overlapped a non-conifer polygon
sn_con_forest_r <- (sn.raster == 1) & (is.na(sn.raster.nonconifer))

# Mask out the non-conifer pixels (they are at this point assigned 0; set them to NA)
sn_con_forest_r[sn_con_forest_r == 0] <- NA


# There are some pixels from the South Sierra Nevada region that are not in the
# SierraEcoregion_TNC polygon, so we can mask those out.
sn_con_forest_r <- mask(sn_con_forest_r, sn)

writeRaster(sn_con_forest_r, filename = maskedFilename, format="GTiff", overwrite=TRUE)

# We want the "good" (e.g., conifer forested, pipo forested) pixels to have a value of 1
# and "bad" (e.g., non-conifer forested, non-pipo forested) pixels to have a value of 0. We can't mask out
# non-forested pixels because Earth Engine asset uploads don't seem to like
# that approach.
sn_con_forest_r[is.na(sn_con_forest_r)] <- 0 # Turn all masked pixels to 0

plot(sn_con_forest_r)
plot(sn, add = TRUE)

writeRaster(sn_con_forest_r, filename = nonMaskedFilename, format="GTiff", overwrite=TRUE)
