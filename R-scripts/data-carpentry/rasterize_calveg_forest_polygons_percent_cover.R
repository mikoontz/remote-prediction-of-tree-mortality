library(sf)
library(fasterize)

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

# Create another shapefile that contains only non-conifer polygons
nsn.nonconifer <- nsn[nsn$target_forest != 1, ]
ssn.nonconifer <- ssn[ssn$target_forest != 1, ]