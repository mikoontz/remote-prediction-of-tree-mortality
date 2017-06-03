library(viridis)
library(raster)
library(sf)

jp <- st_read("features/jepcodes-v7.kml")

# The KML attributes to extract
kml_attributes <- c("FID", "JEPCODE", "JEP_REG", "AREA_ACRES", "AREA_SQ_MI")

# Helper function to extract the first list element
first <- function(x) {
  x[[1]]
}

# Use regular expressions to find the values for the different kml attributes
# that are buried in the HTML within the $Description column of the simple
# feature spatial object

for (i in seq_along(kml_attributes)) {
  for (j in seq_along(row.names(jp))) {
    
    # Starting character within the $Description column for a particular
    # kml attribute
    start_char <-
      kml_attributes[i] %>%
      paste0("</td> <td>") %>%
      gregexpr(text = jp$Description[j]) %>%
      first() %>%
      sum(attr(., which = "match.length"))

    # Stopping character within the $Description column for a particular
    # kml attribute
    stop_char <- 
      kml_attributes[i] %>%
      paste0("</td> <td>[A-Z|a-z| |0-9|.]+</td>") %>%
      gregexpr(text = jp$Description[j]) %>%
      first() %>%
      sum(attr(., which = "match.length")) %>%
      '-'(6) 
    # Note we subtract 6 characters because we searched with "</td>"
    # at the end
    
    # Assign a new column with the extracted data
    jp[j, kml_attributes[i]] <- substr(x = jp$Description[j], start = start_char, stop = stop_char)
  }
}

# Delete the old columns to reduce clutter
jp$Name <- NULL
jp$Description <- NULL

# Convert some jp columns to numeric
jp$AREA_ACRES <- as.numeric(jp$AREA_ACRES)
jp$AREA_SQ_MI <- as.numeric(jp$AREA_SQ_MI)

# These are the Jepson Regions that we will call "The Sierra Nevada"
sn_districts <- c("Central High Sierra Nevada District", 
                  "Central Sierra Nevada Foothills District", 
                  "Northern Sierra Nevada Foothills District", 
                  "Southern Sierra Nevada Foothills District", 
                  "Northern High Sierra Nevada District", 
                  "Southern High Sierra Nevada District", 
                  "Tehachapi Mountain Area Subregion")

sn <- jp[jp$JEP_REG %in% sn_districts, ]

sn %<>%
  st_zm(drop = TRUE) %>%
  as("Spatial") %>%
  aggregate() %>%
  st_as_sf()

st_write(obj = sn, dsn = "features/SierraEcoregion_Jepson/SierraEcoregion_Jepson.shp")
st_write(obj = sn, dsn = "features/SierraEcoregion_Jepson.kml", driver = "KML")

# Next step is to upload the sn as a kml to a Fusion Table, use it to clip a 
# MODIS EVI image, then export that clipped EVI image back here for proper masking
