#-----------------------------------------------------#
#### 1. Load packages and define utility functions ####
#-----------------------------------------------------#

# Specify packages
p <- c("sp","rgdal","raster","ncdf4","gdalUtils","spatial.tools","plyr","reshape2","lme4","ggplot2","rgeos","data.table")

# If necessary: install/update packages
# install.packages(p) # To install `ncdf4` it may be necessary to follow instructions here: http://cirrus.ucsd.edu/~pierce/ncdf/

# Load packages
lapply(p,require,character.only=TRUE)

# Utility functions
fn.bin <- function(x){ifelse(x >0,1,0)} # is a value over 0 or not?



#-----------------------------------------------------#
#### 2. Grid-level variables (not dependent on year) ####
#-----------------------------------------------------#

### Create master rasters (the grid to hold all variables) ###

albers.proj <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

## Load project boundary and raster template
#for restricted area
project.area <- shapefile("features/so-sierra-subset-mask/so-sierra-subset-mask.shp")

#for entire area
project.area <- shapefile("features/SierraEcoregion_Jepson/sierra_ecoregion.shp")


raster_template <- raster("features/sierra-nevada-250m-evi-template.tif")

project.area <- spTransform(project.area,albers.proj)

# Define resolutions
agg.factor <- 10 #how much finer should the fine-resolution raster be, in each dimension. must be integer

# create coarese and fine res raster templates
master.coarse <- mask(crop(raster_template,project.area),project.area)
master.coarse <- master.coarse*0 # set all cells to 0 because we don't need the values that came with the template
master.fine <- disaggregate(master.coarse,fact=agg.factor)

## Extract the parameters needed to re-create the same grids from scratch
master.fine.res <- res(master.fine)
master.fine.extent <- extent(master.fine)[c(1,3,2,4)]

master.coarse.res <- res(master.coarse)
master.coarse.extent <- extent(master.coarse)[c(1,3,2,4)]

master.proj <- projection(master.coarse)


### Create coarse grid that represents the proportion of each coarse grid cell that is land
# Open waterbody layer
waterbody <- readOGR("features/waterbodies","CA_Lakes") # Polygon shapefile layer of major waterbodies. https://www.wildlife.ca.gov/Data/GIS/Clearinghouse
waterbody$water <- 1 # Everywhere that has a polygon is water. Coded as 1.
waterbody <- spTransform(waterbody,projection(master.fine)) # Project to the same projection as the fine-resolution raster template
writeOGR(waterbody, getwd(),"features/waterbodies/waterbody_proj",driver="ESRI Shapefile",overwrite=TRUE) # Write to shapefile for rasterization

rm(waterbody)

# Rasterize it at fine resolution
water.raster <- gdal_rasterize("features/waterbodies/waterbody_proj.shp","features/waterbodies/waterbodyraster.tif",
                               a="water", tr=master.fine.res, te=master.fine.extent,
                               l="waterbody_proj",a_nodata=NA,verbose=TRUE,output_Raster=TRUE)

#Aggregate water raster to the scale for statistical analysis, so that the values in coarse cells reflect the proportion of the coarse cell that is water
waterraster.coarse <- aggregate(water.raster,agg.factor,fun=sum,na.rm=TRUE)
propwater <- waterraster.coarse/(agg.factor^2) # calculate proportion as number of fine cells that were water divided by total number of cells
propland <- 1-propwater # calculate the proportion that is land
propland[is.na(propland)] <- 1 # if the value was NA, that means there was no water there, so it was completely land
propland[propland < .5] <- NA # if it's more than half water, set it to NA so that we do not use those cells as data points.



#-----------------------------------------------------#
####3. Rasterize mortality survey polygons ############
#-----------------------------------------------------#



species.set <- c(122,108,15,20,202,807,631) #which species to rasterize
species.set <- c(122)

### Rasterize mortality data for each year
year.set <- c(2009:2016)
for(i in 1:length(year.set)){
  
  #define year of interest; and aerial survey file for that year
  year <- year.set[i]
  mort.file <- paste("ADS_",year,sep="")
  flight.file <- paste("Flown_area_",year,sep="")
  
  ## Open corresponding shapefiles
  # Survey shapefiles should be in ESRI shapefile format, in the folder "ADS-shapefiles" which is a subfolder of the working directory, with naming "ADS_2009.shp" etc.
  # There should also be flown area shapefiles (same format) in the same subfolder, one file for each year, with naming "Flown_area_2009.shp" etc.
  # Annual survey spatial data accessible here: http://www.fs.usda.gov/detail/r5/forest-grasslandhealth/?cid=fsbdev3_046696 (needs to be converted from geodatabases into shapefiles before opening in this script)
  mort.polygon <- readOGR("features/ADS-shapefiles",mort.file,stringsAsFactors=FALSE)
  flight.polygon <- readOGR("features/ADS-shapefiles",flight.file,stringsAsFactors=FALSE)
  
  # Reproject to projection of master raster
  mort.polygon <- spTransform(mort.polygon,master.proj)
  flight.polygon <- spTransform(flight.polygon,master.proj)
  
  ## Make sure at least one of the damage types was mortality (because defoliation is an option too), and remove those that aren't
  mortality <- (mort.polygon$DMG_TYPE1 == 2) | (mort.polygon$DMG_TYPE2 == 2) | (mort.polygon$DMG_TYPE3 == 2)
  mort.polygon <- mort.polygon[mortality,]
  
  ## Remove the one outlier (the only point in all years that is >1200 TPA). TPA is trees per acre, a field reported in the aerial surveys
  mort.polygon <- mort.polygon[mort.polygon$TPA1 < 1200,]
  
  #Add column to flight shapefile that indicates that it was flown (every polygon will have the value 1 for this attribute; useful for rasterization later)
  flight.polygon$FLOWN1 <- 1
  
  ## Eliminate polygons where mortality was attributed to non drought-related agents/causes
  fire <- 30000:30005
  animals <- 41000:42900
  physical <- c(50004:50006,50011:50020)
  human <- c(70000:70011,70013:71001)
  non.drought <- c(fire,animals,physical,human)
  mort.polygon <- mort.polygon[!((mort.polygon$DCA1 %in% non.drought) | (mort.polygon$DCA2 %in% non.drought) | (mort.polygon$DCA3 %in% non.drought)),]
  
  # Sum the three TPA columns wherever they are positive for total dead trees per acre
  tpa.cols <- cbind(mort.polygon$TPA1,mort.polygon$TPA2,mort.polygon$TPA3)
  tpa.cols[tpa.cols <= 0] <- NA # if any values are 0 or less, set to NA so they are disregarded
  mort.polygon$TPA.tot <- rowSums(tpa.cols,na.rm=TRUE)

  ## Eliminate "background mortality" using filters based on personal communication with Zachary Heath and Jeffrey Moore  
  # Eliminate polygons where total TPA == 1 (when TPA is less than 1 it is usually from large polygons with a small fixed number of trees, which are likely not "background" mortality)
  mort.polygon <- mort.polygon[mort.polygon$TPA.tot != 1,]
  # Eliminate polygons where the number of trees is 3 or less
  mort.polygon <- mort.polygon[which(as.numeric(mort.polygon$NO_TREES1) > 3) , ]
  
  ## When multiple host species present, divide TPA by the number of hosts
  # Tally number of hosts listed
  host.cols <- cbind((mort.polygon$HOST1), (mort.polygon$HOST2), (mort.polygon$HOST3))
  host.cols[host.cols < 1] <- NA
  n.host <- rowSums(host.cols > 0,na.rm=TRUE)
  
  # Divide TPA by number of hosts
  mort.polygon$TPA.split <- mort.polygon$TPA.tot/n.host
  mort.polygon$TPA.split[mort.polygon$TPA.split == Inf] <- 0 # because sometimes there are zero hosts with a host id > 0
  # Now TPA.split holds the average number of dead trees PER HOST
  
  
  
  #### Write flight path polygon to shapefile for rasterization
  # Write it
  writeOGR(flight.polygon, "features/ADS-intermediate-products","flightpolygon",driver="ESRI Shapefile",overwrite=TRUE)

  # Rasterize it to the master.fine grid
  flightraster <- gdal_rasterize("features/ADS-intermediate-products/flightpolygon.shp","features/ADS-intermediate-products/flightrasterfocal_test.tif",
                                 a="FLOWN1", tr=master.fine.res, te=master.fine.extent,
                                 l="flightpolygon",a_nodata=NA,verbose=TRUE,output_Raster=TRUE)
  
  # Aggregate flight raster to the scale for statistical analysis
  # Note that this operation sets all coarse cells that are partially outside the flight path to NA
  flightraster.coarse <- aggregate(flightraster,agg.factor,fun=mean,na.rm=FALSE)
  
  rm(flight.polygon)
  

  # For each species
  for(j in 1:(length(species.set)+1)){

    #in the last iteration through this loop, rasterize ALL mortality
    if(j == length(species.set)+1) {
      all.sp = TRUE
    } else {
      all.sp = FALSE
    }
    
    
    # Define tree species set 
    if(all.sp) {
      focal.sp <- -998
    } else {
      focal.sp <- species.set[[j]] 
    }
    
    # if(focal.sp == 202) {
    #   next()
    # }
    
    ## For each mortality polygon, convert polygon TPA to TPA of focal species only (i.e., set TPA 0 if polygon does not contain any of focal species)
    # See if focal species is present
    if(all.sp) {
      n.host.match.focal <- n.host
    } else {
      n.host.match.focal <- (mort.polygon$HOST1 %in% focal.sp) + (mort.polygon$HOST2 %in% focal.sp) + (mort.polygon$HOST3 %in% focal.sp)
    }
    
    if(sum(n.host.match.focal) == 0) {
      cat("\rNo mortality polygons found for species",focal.sp,", year",year)
      next()
    }
    
    # Multiply the average per-host TPA by the number of hosts that are focal species
    mort.polygon$tpa.focal <- mort.polygon$TPA.split * n.host.match.focal
    mort.polygon$tpa.focal <- ifelse(is.na(mort.polygon$tpa.focal),0,mort.polygon$tpa.focal) # if number of hosts was set to NA (meaning no matching hosts), set it to 0 so it can be added
    
    # Eliminate polygons with none of the current focal species
    if(sum(mort.polygon$tpa.focal) > 0) {
      mort.polygon.focal <- mort.polygon[mort.polygon$tpa.focal > 0,]
    } else {
      mort.polygon.focal <- mort.polygon[0,]
    }
    
    # Sort the polygons so the polygons with larger mortality density come later
    mort.polygon.focal <- mort.polygon.focal[order(mort.polygon.focal$TPA1),]
    
    ## Rasterize mortality polygon; define unobserved cells, i.e. not in flight path, as NA
    # Write mortality polygon of focal species to shapefile for next step (rasterization)
    writeOGR(mort.polygon.focal, "features/ADS-intermediate-products","mortpolygonfocal", driver="ESRI Shapefile",overwrite=TRUE)
    
    # Rasterize
    mortraster <- gdal_rasterize("features/ADS-intermediate-products/mortpolygonfocal.shp","features/ADS-intermediate-products/mortrasterfocal.tif",
                                 a="tpa_focal",tr=master.fine.res, te=master.fine.extent,
                                 l="mortpolygonfocal",verbose=TRUE,output_Raster=TRUE)
    
    # for raster cells that did not have any overlapping polygons, set mortality value to 0
    #a <- is.na(mortraster)
    mortraster[is.na(mortraster)] <- 0
    
    # Aggregate fine-scale raster to the scale for statistical analysis
    mortraster.coarse <- aggregate(mortraster,agg.factor,fun=mean,na.rm=FALSE)
    
    # Set aggregated mortality cells that are at all outside flight path to NA (because cells in flightraster.coarse that are outside the flight path are NA)
    mort.flight <- mortraster.coarse * flightraster.coarse
    
    # Divide mortality density by proportion of cell that is land to obtain the on-land mortality density
    mort.flight <- mort.flight / propland


    #### Prepare for statistical analysis
    
    # Stack rasters, name layers, and write to external file
    raster_stack <- stack(mort.flight)
    if(all.sp) {
      spgrp.text <- "ALL"
    } else {
      spgrp.text <- sprintf("%03d",focal.sp)[1] # add leading zeros
    }
    
    layer.names <- paste("Y",year,".sp",spgrp.text,".",c("mort.tpa"),sep="")
    names(raster_stack) <- layer.names
    
    writeRaster(raster_stack, file=paste0("features/ADS-rasterized/Y",year,"_","sp",spgrp.text,".tif",sep=""),overwrite=TRUE)
    cat("\rFinished Year",year,", species group",j)
  }
}