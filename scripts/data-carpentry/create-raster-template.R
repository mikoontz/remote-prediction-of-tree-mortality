library(viridis)
library(raster)

sn <- shapefile("features/SierraEcoregion_TNC/SierraEcoregion_TNC.shp")
evi <- raster("features/sierra-nevada-250m-evi-template.tif")
  
plot(sn)
plot(evi)
plot(sn, add = TRUE)

m_evi <- mask(evi, sn)
plot(m_evi)

writeRaster(m_evi, filename="features/sierra-nevada-250m-evi-template.tif", format="GTiff", overwrite=TRUE)

test <- raster("features/sierra-nevada-250m-evi-template.tif")
plot(test)