library(rgdal)


facts.fueltrt <- readOGR(dsn = "features/FACTS/S_USA.Activity_HazFuelTrt_PL",layer="S_USA.Activity_HazFuelTrt_PL", stringsAsFactors = FALSE)

writeOGR(facts.fueltrt,getwd(),"features/FACTS/cleaned/fueltrt",driver="ESRI Shapefile",overwrite=TRUE)



facts.reforest <- readOGR(dsn = "features/FACTS/S_USA.Activity_SilvReforestation",layer="S_USA.Activity_SilvReforestation", stringsAsFactors = FALSE)
