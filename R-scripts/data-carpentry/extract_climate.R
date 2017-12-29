### This script: 
## 1) converts monthly temperature and ppt data from raw PRISM and TopoWx layers to multilayer rasters using
##    same grid (resolution and origin) as the EVI layers, with one layer for each month of each year from 2000 to 2016

library(sp)
library(raster)
library(stringr)

# Projection info
albers.proj <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
wgs.proj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Get the locations (grid resolution, origin, ectent) of the target pixels 
# EVI template raster 
evi_template = raster("features/sierra-nevada-250m-evi-template.tif")

## dates (layer names) need to be in format as YYYYMMDD


#### Process PRISM ####

## open and stack all PRISM layers


years <- 2000:2016
months <- 1:12
months.chr <- str_pad(months,width=2,side="left",pad="0")

date <- paste0(rep(years,each=12),months.chr)

folder <- "features/raw-climate-layers/prism-ppt/"
prefix <- "PRISM_ppt_stable_4kmM3_"
suffix <- "_bil.bil"

prism.ppt.files <- paste0(folder,prefix,date,suffix)

prism.ppt.stack <- stack(prism.ppt.files)

##resample PRISM ppt to EVI template
extent.template.albers <- extent(evi_template)
extent.template.albers <- as(extent.template.albers,"SpatialPolygons")
proj4string(extent.template.albers) <- proj4string(evi_template)
extent.template.geo <- spTransform(extent.template.albers,proj4string(prism.ppt.stack))

prism.ppt.stack <- raster::crop(prism.ppt.stack,extent.template.geo)
layer.names <- paste0("ppt_",date)
names(prism.ppt.stack) <- layer.names

prism.ppt.template <- raster::projectRaster(prism.ppt.stack,evi_template)

## write projected prism data as a grd (which preserves layer names)
writeRaster(prism.ppt.template,"features/climate-layers-stacked-projected/ppt-monthly.grd")








#### Process TopoWx tmax (tmx) ####

## open and stack all TopoWx layers


years <- 2000:2016
months <- 1:12
months.chr <- str_pad(months,width=2,side="left",pad="0")
date <- paste0(rep(years,each=12),months.chr)


folder <- "features/climate-layers-raw/topowx-tmax/"
prefix <- "tmax_"
suffix <- ".nc"

topowx.tmax.files <- paste0(folder,prefix,years,suffix)

topowx.tmax.stack <- stack(topowx.tmax.files)

##resample PRISM ppt to EVI template
extent.template.albers <- extent(evi_template)
extent.template.albers <- as(extent.template.albers,"SpatialPolygons")
proj4string(extent.template.albers) <- proj4string(evi_template)
extent.template.geo <- spTransform(extent.template.albers,proj4string(topowx.tmax.stack))

topowx.tmax.stack <- raster::crop(topowx.tmax.stack,extent.template.geo)
layer.names <- paste0("tmx_",date)
names(topowx.tmax.stack) <- layer.names

topowx.tmax.template <- raster::projectRaster(topowx.tmax.stack,evi_template)

## write projected prism data as a grd (which preserves layer names)
# writeRaster(topowx.tmax.template,"features/climate-layers-stacked-projected/tmx-monthly.grd")





#### Process TopoWx tmin (tmn) ####

## open and stack all TopoWx layers


years <- 2000:2016
months <- 1:12
months.chr <- str_pad(months,width=2,side="left",pad="0")
date <- paste0(rep(years,each=12),months.chr)


folder <- "features/climate-layers-raw/topowx-tmin/"
prefix <- "tmin_"
suffix <- ".nc"

topowx.tmn.files <- paste0(folder,prefix,years,suffix)

topowx.tmn.stack <- stack(topowx.tmn.files)

##resample PRISM ppt to EVI template
extent.template.albers <- extent(evi_template)
extent.template.albers <- as(extent.template.albers,"SpatialPolygons")
proj4string(extent.template.albers) <- proj4string(evi_template)
extent.template.geo <- spTransform(extent.template.albers,proj4string(topowx.tmn.stack))

topowx.tmn.stack <- raster::crop(topowx.tmn.stack,extent.template.geo)
layer.names <- paste0("tmn_",date)
names(topowx.tmn.stack) <- layer.names

topowx.tmn.template <- raster::projectRaster(topowx.tmn.stack,evi_template)

## write projected prism data as a grd (which preserves layer names)
# writeRaster(topowx.tmn.template,"features/climate-layers-stacked-projected/tmn-monthly.grd")




#### Compute and save TopoWx tmean (tmp) ####

## open and stack all TopoWx layers

topowx.tmp.stack <- (topowx.tmn.stack + topowx.tmax.stack)/2
layer.names <- paste0("tmp_",date)
names(topowx.tmp.stack) <- layer.names

topowx.tmp.template <- raster::projectRaster(topowx.tmp.stack,evi_template)

## write projected prism data as a grd (which preserves layer names)
writeRaster(topowx.tmp.template,"features/climate-layers-stacked-projected/tmp-monthly.grd")

