library(devtools)
#devtools::install_github("ecohealthalliance/fasterize") 
library(sf)
library(fasterize)
library(raster)
library(viridis)

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

# Load project boundary and raster template
sn <- shapefile("features/SierraEcoregion_TNC/SierraEcoregion_TNC.shp")
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

# Create new column that indicates whether a polygon contains conifer forest (1 if yes, 2 if no)
nsn$target_forest <- ifelse(test = nsn$WHRTYPE %in% con.whr.types, yes = 1, no = 0)
ssn$target_forest <- ifelse(test = ssn$WHRTYPE %in% con.whr.types, yes = 1, no = 0)

# Disaggregate the raster_template to get a finer resolution
raster_template_fine <- disaggregate(raster_template, fact = c(10, 10))

# Rasterize the north and south Sierra Nevada polygons representing target forest
# types separately
nsn_target_raster <- fasterize(sf = nsn, field = "target_forest", raster = raster_template_fine, fun = "sum")
ssn_target_raster <- fasterize(sf = ssn, field = "target_forest", raster = raster_template_fine, fun = "sum")

# Get cover estimate for each 250m cell by summing the 100 mini-cells within it
nsn_target_cover <- aggregate(nsn_target_raster, fact = c(10, 10), fun = sum)
ssn_target_cover <- aggregate(ssn_target_raster, fact = c(10, 10), fun = sum)

# Merge the North and South rasters back together
target_cover <- merge(nsn_target_cover, ssn_target_cover)
plot(target_cover, col = viridis(10))

# Export a GeoTiff for use with R (including mask around Sierra Nevada)
writeRaster(target_cover, filename = "features/sierra_nevada_250m_calveg_pipo_forest-cover_whr_type.tif")

# Export a GeoTiff for use with Earth Engine (all masked pixels get a value of 0)
target_cover_ee <- target_cover
target_cover_ee[is.na(target_cover_ee)] <- 0

writeRaster(target_cover_ee, filename = "features/sierra_nevada_250m_calveg_pipo_forest-cover_whr_type_no_mask.tif")
