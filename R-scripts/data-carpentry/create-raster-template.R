library(viridis)
library(raster)
library(sf)

sn <- shapefile("features/SierraEcoregion_TNC/SierraEcoregion_TNC.shp")
sn <- st_read("features/SierraEcoregion_Jepson/")

kml_attributes <- c("FID", "JEPCODE", "JEP_REG", "AREA_ACRES", "AREA_SQ_MI")


first <- function(x) {
  x[[1]]
}

for (i in seq_along(kml_attributes)) {
  for (j in seq_along(row.names(sn))) {
    
    start_char <-
      kml_attributes[i] %>%
      paste0("</td> <td>") %>%
      gregexpr(text = sn$Description[j]) %>%
      first() %>%
      sum(attr(., which = "match.length"))
    
    stop_char <- 
      kml_attributes[i] %>%
      paste0("</td> <td>[A-Z|a-z| |0-9|.]+</td>") %>%
      gregexpr(text = sn$Description[j]) %>%
      first() %>%
      sum(attr(., which = "match.length")) %>%
      '-'(6)
    
    sn[j, kml_attributes[i]] <- substr(x = sn$Description[j], start = start_char, stop = stop_char)
  }
}

tail(sn, 5)
plot(sn$geometry)
start_char <- gregexpr(pattern = paste0(kml_attributes[i], "</td> <td>"), text = sn$Description[j])[[1]]
attr_length <- attr(start_char, which = "match.length")
start_char <- start_char + attr_length

gregexpr(pattern = paste0(kml_attributes[i], "</td> <td>[A-Z|a-z| |0-9]+</td>"), text = sn$Description[j])[[1]]

gregexpr(pattern = "JEPCODE</td> <td>", text = sn$Description[1])
gregexpr(pattern = "JEPCODE</td> <td>[A-Z|a-z| |0-9]+</td>", text = sn$Description[1])

gregexpr(pattern = "JEP_REG</td> <td>", text = sn$Description[1])
gregexpr(pattern = "JEP_REG</td> <td>[A-Z|a-z| |0-9]+</td>", text = sn$Description[1])

gregexpr(pattern = "FID</td> <td>", text = sn$Description[1])
gregexpr(pattern = "FID</td> <td>[A-Z|a-z| |0-9]+</td>", text = sn$Description[1])

gregexpr(pattern = "AREA_ACRES</td> <td>", text = sn$Description[1])
gregexpr(pattern = "AREA_ACRES</td> <td>[A-Z|a-z| |0-9|.]+</td>", text = sn$Description[1])

gregexpr(pattern = "AREA_SQ_MI</td> <td>", text = sn$Description[1])
gregexpr(pattern = "AREA_SQ_MI</td> <td>[A-Z|a-z| |0-9|.]+</td>", text = sn$Description[1])

as.character(sn$Description[1])
?gregexpr
str(sn)
sn$Description[1]
evi <- raster("features/sierra-nevada-250m-evi-template.tif")
  
plot(sn$geometry)
plot(evi)
plot(sn, add = TRUE)

m_evi <- mask(evi, sn)
plot(m_evi)

writeRaster(m_evi, filename="features/sierra-nevada-250m-evi-template.tif", format="GTiff", overwrite=TRUE)

test <- raster("features/sierra-nevada-250m-evi-template.tif")
plot(test)
