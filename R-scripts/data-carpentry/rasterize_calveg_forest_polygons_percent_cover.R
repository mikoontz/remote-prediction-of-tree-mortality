#### This script:
## 1) defines a set of target vegetation types for analysis, using the CalVeg WHR types
## 2) for each WHR type, selects polygons from the Sierra CalVeg tiles that are type
## 3) rasterizes it to a resolution 100 times as fine as the EVI template raster, then aggregate back to template resolution to produce layer that has percent cover of the target WHR type. 


library(sf)
library(fasterize)
library(raster)
library(viridis)

target_whr_types <- c("SMC",  # Sierra mixed conifer
                    "MCN",  # Mixed conifer
                    "MHC",  # Mixed hardwood-conifer
                    "SCN",  # Subalpine conifer
                    "JPN",  # Jeffrey pine
                    "PPN",  # Ponderosa pine
                    "WFR",  # White fir
                    "RFR",  # Red fir
                    "DFR")  # Douglas-fir

# Use just Ponderosa pine WHR Type
#target_whr_types <- "PPN"  

# Load raster template
# Raster template is an EVI layer exported from Earth Engine 
# It is already masked to the target subregion (Jepson ecoregions)
# CURRENTLY WAITING FOR FINAL VERSION FROM MIKE
raster_template <- raster("features/sierra-nevada-250m-evi-template.tif")
raster_template[] <- 0

# Read in CALVEG layers
nsn <- st_read(dsn = "features/ExistingVegNorSierra2000_2014_v1.gdb", stringsAsFactors = FALSE)
ssn <- st_read(dsn = "features/ExistingVegSouthSierra2000_2008_v1.gdb", stringsAsFactors = FALSE)

# Keep only the relevant attributes
nsn <- nsn[,c("WHRTYPE","WHRLIFEFORM")]
ssn <- ssn[,c("WHRTYPE","WHRLIFEFORM")]

# Reproject CALVEG layers to projection of raster template
nsn <- st_transform(nsn, crs = proj4string(raster_template))
ssn <- st_transform(ssn, crs = proj4string(raster_template))

# Create a new layer for each target forest type that indicates what percentage of each pixel contains that forest type.
n_whr_types = length(target_whr_types)

for (i in 1:n_whr_types) {
  # Create new column that indicates whether a polygon contains target WHR type i (1 if yes, 2 if no)
  nsn_tmp = nsn
  ssn_tmp = ssn
  nsn_tmp$target_forest <- ifelse(test = nsn_tmp$WHRTYPE == target_whr_types[i], yes = 1, no = 0)
  ssn_tmp$target_forest <- ifelse(test = ssn_tmp$WHRTYPE %in% target_whr_types, yes = 1, no = 0)

  # Disaggregate the raster_template to get a finer resolution
  raster_template_fine <- disaggregate(raster_template, fact = c(10, 10))

  # Rasterize the north and south Sierra Nevada polygons representing target forest types separately
  nsn_target_raster <- fasterize(sf = nsn_tmp, field = "target_forest", raster = raster_template_fine, fun = "sum")
  ssn_target_raster <- fasterize(sf = ssn_tmp, field = "target_forest", raster = raster_template_fine, fun = "sum")

  # Get cover estimate for each 250m cell by summing the 100 mini-cells within it
  nsn_target_cover <- aggregate(nsn_target_raster, fact = c(10, 10), fun = sum)
  ssn_target_cover <- aggregate(ssn_target_raster, fact = c(10, 10), fun = sum)

  # Merge the North and South rasters back together
  target_cover <- merge(nsn_target_cover, ssn_target_cover)
  #plot(target_cover, col = viridis(10))

  # Export a GeoTiff to use to construct vegetation type masks within target region
  writeRaster(target_cover, filename = paste("features/sierra_nevada_250m_calveg_cover_whr_type_", target_whr_types[i], ".tif", sep=""))

}
