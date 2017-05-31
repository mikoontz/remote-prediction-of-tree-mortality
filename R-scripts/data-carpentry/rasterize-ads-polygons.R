#-----------------------------------------------------#
#### 1. Load packages and define utility functions ####
#-----------------------------------------------------#

# Specify packages
p <- c("sp","rgdal","raster","ndcf4","gdalUtils","spatial.tools","plyr","reshape2","lme4","ggplot2","rgeos","data.table")

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

# Load project boundary and raster template
project.area <- shapefile("features/so-sierra-subset-mask/so-sierra-subset-mask.shp")
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
  # Survey shapefiles should be in ESRI shapefile format, in the folder "Survey shapefiles" which is a subfolder of the working directory, with naming "ADS_2009.shp" etc.
  # There should also be flown area shapefiles (same format) in the same subfolder, one file for each year, with naming "Flown_area_2009.shp" etc.
  # Annual survey spatial data accessible here: http://www.fs.usda.gov/detail/r5/forest-grasslandhealth/?cid=fsbdev3_046696 (needs to be converted from geodatabases into shapefiles before opening in this script)
  mort.polygon <- readOGR("features/ADS-shapefiles",mort.file)
  flight.polygon <- readOGR("features/ADS-shapefiles",flight.file)
  
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
  mort.polygon <- mort.polygon[mort.polygon$NO_TREES1 > 3 , ]
  
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
  for(j in 1:length(species.set)){

    # Define tree species set 
    focal.sp <- species.set[[j]] 
    
    # if(focal.sp == 202) {
    #   next()
    # }
    
    ## For each mortality polygon, convert polygon TPA to TPA of focal species only (i.e., set TPA 0 if polygon does not contain any of focal species)
    # See if focal species is present
    n.host.match.focal <- (mort.polygon$HOST1 %in% focal.sp) + (mort.polygon$HOST2 %in% focal.sp) + (mort.polygon$HOST3 %in% focal.sp)
    
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
    spgrp.text <- sprintf("%03d",focal.sp)[1] # add leading zeros
    
    layer.names <- paste("Y",year,".sp",spgrp.text,".",c("mort.tpa"),sep="")
    names(raster_stack) <- layer.names
    
    writeRaster(raster_stack, file=paste0("features/ADS-rasterized/Y",year,"_","sp",spgrp.text,".tif",sep=""),overwrite=TRUE)
    cat("\rFinished Year",year,", species group",j)
  }
}



#####!!!! Code below this line not yet prepared for remote detection study !!!!#####
## It is just copy-pasted from EcoLett script


#-----------------------------------------------------#
####4. Forest structure data ##########################
#-----------------------------------------------------#

## Open layers with trees per hectare and basal area to use as predictor variables and to know where trees were present
# Useful because it's not interesting to report no mortality in a place there there were no trees

# The two layers used in this section are the "Live trees per hectare" and "Live basal area" layers from the GNN dataset: http://lemma.forestry.oregonstate.edu/data/structure-maps
# These maps need to be pre-processed into one layer for each of the two variables (with file names as below), and aggregated to a coarser scale so they are roughly the same resolution as the master.coarse raster used in this script
# This section syncs them so they are on exactly the same grid (origin and resolution) as the coarse mortality data grid.
# It then writes them to a different sub-folder ("RasterOutput") for later incorporation into the final raster stack used to create the data frame

s <- raster("DistrRastAgg/live_tph_all_agg.grd")
s <- spatial_sync_raster(s,master.coarse,method="bilinear")
names(s) <- "live.tph"
writeRaster(s,paste("RasterOutput/sp_all_distr_tph.grd",sep=""),overwrite=TRUE)

s <- raster("DistrRastAgg/live_bah_all_agg.grd")
s <- spatial_sync_raster(s,master.coarse,method="bilinear")
names(s) <- "live.bah"
writeRaster(s,paste("RasterOutput/sp_all_distr_bah.grd",sep=""),overwrite=TRUE)


#-----------------------------------------------------#
####5. Distribution mask ##############################
#-----------------------------------------------------#

# This is to exclude cells that are outside the forested area in California 
s <- raster("DistrRastAgg/live_tph_all_agg.grd") # Raster from the same source as in the above section
 
# Set it to 1s everywhere it is not NA (that is, everywhere there are any trees), and 0s everywhere it is NA (that is, everywhere there are no trees)
s[!is.na(s[])] <- 1
s[is.na(s[])] <- 0 

# Sync it to exactly the same grid (origin and resolution) as the coarse mortality data grid
s <- spatial_sync_raster(s,master.coarse,method="bilinear")

# Only keep areas where it is exactly 1 (all synced cells were completely within original cells); everything else -> 0. That is, do not use any coarse raster cells where any part of the coarse cell was outside an area that had any trees.
s[is.na(s[])] <- 0 
s[s[]<0.999] <- 0
s[s[]>=0.999] <- 1

names(s) <- "ca.mask"
writeRaster(s,paste("RasterOutput/ca_mask.grd",sep=""),overwrite=TRUE)

#-----------------------------------------------------#
####6. Climate data ###################################
#-----------------------------------------------------#

### Open and prep climate layers ###

# All of the climate layers used here are pre-processed to some extent. The layers used here are all at ~4 km resolution (must be aggregated in advance if they are supplied in a lower resolution.)
# For each climate variable, the rasters are supplied in the format of a "RasterStack" (in .grd format), with one layer for each year, from 1981 to 2015, named, "climate.variable.year", where "climate.variable" matches the name of the entire raster stack, and "year" is a value between 1981 and 2015.
# For example, "ppt.ann.1981","ppt.ann.1982", etc., and "def.ann.1981", "def.ann.1982", etc.
# The CWD and AET values were calculated following Willmott et al. (1985) Climatology of the terrestrial seasonal water cycle. Journal of Climatology. 5: 589-606, using monthly values of temperature (from TopoWx), and precipitation (from PRISM) for each year from 1980 to 2015, with output water balance values summarized annually from 1981 to 2015
# The variables below refer to the following climate variables:
#  tmean.ann: annual average of monthly minimum and monthly maximum temperatures in degrees Celsius from TopoWx dataset: http://www.ntsg.umt.edu/project/TopoWx
#  ppt.ann: total annual precipitation in millimeters from PRISM dataset: http://prism.oregonstate.edu
#  aet.ann: total annual actual evapotranspiration in millimeters modeled according to Willmott method (cited above)
#  def.ann: total annual climatc water deficit in millimeters modeled according to Willmott methods (cited above)
#  snow.ann: number of months per year that snow cover exists. Calculated according to this method: http://elib.dlr.de/95803/ (the dataset was provided to us by the authors of this method upon request) -- not using for the present analysis but included in script for potential future extension of analysis


tmean.ann <- brick("Climate layers/TopoWx_tmean_wateryear/TopoWx_tmean_wateryear.grd")
ppt.ann <- brick("Climate layers/PRISM_ppt_wateryear/PRISM_ppt_wateryear.grd")
aet.ann <- brick("Climate layers/ThornWil_AET_wateryear/ThornWil_AET_wateryear.grd")
def.ann <- brick("Climate layers/ThornWil_Def_wateryear/ThornWil_Def_wateryear.grd")

# These layers allow for including additional climate layers (e.g. snow extent, Dobrowski et al. 2013 water balance and BCM water balance) as predictors of mortality but were not ultimately used in our analysis. Also presented as raster stacks of annual values from 1981 to 2015 as above.
snow.ann <- brick("Climate layers/snowpack/CA_snowpack_4km.grd")
dob.pub.aet <- raster("Climate layers/published wb layers/dob_pub_aet.grd")
dob.pub.def <- raster("Climate layers/published wb layers/dob_pub_def.grd")
bcm.aet <- raster("Climate layers/published wb layers/bcm_aet.grd")
bcm.def <- raster("Climate layers/published wb layers/bcm_def.grd")

# Set projections appropriately so rasters can be synced and stacked
projection(tmean.ann) <- projection(ppt.ann)
ppt.temp.rast <- stack(tmean.ann,ppt.ann)
wb.rast <- stack(aet.ann,def.ann)

# Interpolate these ~4km rasters down to the 3.5 km grid used for statistical analysis (master.coarse), and stack them all together into one large raster stack object
ppt.temp.rast <- spatial_sync_raster(ppt.temp.rast,master.coarse,method="bilinear")
wb.rast <- spatial_sync_raster(wb.rast,master.coarse,method="bilinear")
climate.rast <- stack(ppt.temp.rast,wb.rast)

snow.sync <- spatial_sync_raster(snow.ann,climate.rast,method="bilinear")
dob.pub <- stack(dob.pub.aet,dob.pub.def)
bcm <- stack(bcm.aet,bcm.def)
dob.pub.sync <- spatial_sync_raster(dob.pub,climate.rast,method="bilinear")
bcm.sync <- spatial_sync_raster(bcm,climate.rast,method="bilinear")

climate.rast <- stack(climate.rast,snow.sync,dob.pub.sync,bcm.sync)

# Write the large stack of annual values from many climate variables to file to load later
writeRaster(climate.rast,"RasterOutput/climate.grd",overwrite=TRUE)
writeRaster(climate.rast,"RasterOutput/climate.tiff",overwrite=TRUE)


#-----------------------------------------------------#
####7. Jepson Ecoregions ##############################
#-----------------------------------------------------#

# Load Jepson ecoregions for potentially breaking analyisis up by region (not used for published analysis, but provided here for facilitating future extension)
jepson <- readOGR(dsn="jepson ecoregions/coarse good/jepson.shp",layer="jepson") # Jepson ecoregions shapefile, provided by email from the Jepson Herbarium in response to an email request.

jeps.rast <- gdal_rasterize("../GIS/jepson ecoregions/coarse good/jepson.shp","../GIS/jepson ecoregions/coarse good/jeps_rast_test.tif",
                            a="REGION", tr=master.coarse.res, te=master.coarse.extent,
                            l="jepson",verbose=TRUE,output_Raster=TRUE)

names(jeps.rast) <- "jepson"

writeRaster(jeps.rast,"RasterOutput/jepson.grd",overwrite=TRUE)



#-----------------------------------------------------#
####8. Convert all raster layers into a data frame ####
#-----------------------------------------------------#


# List all .grd files from RasterOutput folder
RasterOutput.files <- list.files(path="RasterOutput",pattern=c(".grd"),full.names=TRUE)

# Pull them all together into a stack
raster.output.uncropped <- stack(RasterOutput.files)

## Drop all cells that fall outside CA (because live tree BA does not go outside CA, and climate data in ocean not dependable)
# Define California Albers projection
albers.proj <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
# Open california shapefile
calif <- readOGR("California shapefile","statep010")  # Open a polygon shapefile with a single feature, representing the state of california, in a geographic (lat/long) projection. For example you could use the following dataset (but before loading the file in this script, remove all states except California): https://catalog.data.gov/dataset/usgs-small-scale-dataset-state-boundaries-of-the-united-states-200506-shapefile
calif <- spTransform(calif,albers.proj) #project to albers (meter units)

# Mask all the raster output (set cells that fall outside the state to NA)
raster.output <- mask(raster.output.uncropped,calif)

# Write the resulting file
writeRaster(raster.output,"raster_output.tif",overwrite=TRUE)
writeRaster(raster.output.uncropped,"raster_output_uncropped.tif",overwrite=TRUE)

# Convert raster stack to data.frame
d <- as.data.frame(raster.output) # NAs in the mortality columns mean cells were outside flightpath

# Add in the spatial position of each cell
cells <- 1:ncell(master.coarse)
pos <- rowColFromCell(master.coarse,cells)
d$pos.x <- pos[,1]
d$pos.y <- pos[,2]

# Get California Albers coordinates for each cell too, in m
pos.alb <- xyFromCell(master.coarse,cells)
d$alb.x <- round(pos.alb[,1])
d$alb.y <- round(pos.alb[,2])

# Remove cells that are outside the boundaries of the climate layer (can do this based on any year since the layers extents are the same each year)
d <- d[!is.na(d$tmean.ann.2010),]

# Set mortality to NA when distr = 0 or NA, unless mort > 0 in any year; this is so we don't consider (lack of) mortality in areas where there were no trees in the first place
# Also remove cells that are outside the boundaries of the distribution layer

dist.col <- list() # List to hold data frame indicating where trees were present. Currently only will have one element, but could have multiple in a future extension (e.g., to look at individual species)
mort.col <- list() # List to hold data frame indicating where mortality was present. Currently only will have one element, but could have multiple in a future extension (e.g., to look at individual species)

i <- 1
search <- "live.tph" # Find the live trees per hectare column
dist.col[[i]] <- grep(search,names(d))  

spgrp.text <- sprintf("%02d",i) # add leading zero
search <- paste("spgrp",spgrp.text,".mort",sep="") # Find the mortality amount columns for current selected species group
mort.col[[i]] <- grep(search,names(d))

total.mort <- rowSums(d[,mort.col[[i]]],na.rm=TRUE) # Look at the mortality across all species in the group (in current implementation there isjust be one column representing mortality regardless of species)

## If there is either mortality or live trees there (or both), consider it a cell that should be incorporated in the analysis.
dist.present <- ifelse( (total.mort > 0) | (d[,dist.col[[i]] ] > 0) ,1,NA)
d[, mort.col[[i]] ] <- d[, mort.col[[i]] ] * dist.present

## Melt data frame so that all mortality variables are in a single column
mort.cols <- grep("mort.tpa",names(d))
d.melt <- melt(d,measure.vars=mort.cols,variable.name="ID",value.name="mort.tpa")
d.melt$mort.bin <- fn.bin(d.melt$mort.tpa) # additionally add a 0/1 column indicating whether any mortality is present

# Remove rows where mortality is NA (outside of flight path or species distribution)
d.melt <- d.melt[!is.na(d.melt$mort.tpa),]

# Add other necessary columns
d.melt$year <- substr(d.melt$ID,2,5) #add column for year
d.melt$spgroup <- substr(d.melt$ID,7,13) #add column for species group


#-----------------------------------------------------#
####9. Derive climate variables for data frame ########
#-----------------------------------------------------#

# E.g., calculae normals and annual anomalies from annual time-series

# Calculate climate normal and standard deviation
normal.years <- 1981:2015
tmean.names <- paste("tmean.ann.",normal.years,sep="")
ppt.names <- paste("ppt.ann.",normal.years,sep="")
aet.names <- paste("aet.ann.",normal.years,sep="")
def.names <- paste("def.ann.",normal.years,sep="")
snow.names <- paste("snow.late.",2001:2015,sep="")

tmean.cols <- grep("tmean.ann.",names(d.melt))
ppt.cols <- grep("ppt.ann",names(d.melt))
aet.cols <- grep("aet.ann",names(d.melt))
def.cols <- grep("def.ann",names(d.melt))
snow.cols <- grep("snow.late",names(d.melt))

# Climatic normals
d.melt$Tnorm <- rowMeans(d.melt[,tmean.names],na.rm=TRUE) # Average temperature across all years
d.melt$Pnorm <- rowMeans(d.melt[,ppt.names],na.rm=TRUE) # Average precip across all years
d.melt$AETnorm <- rowMeans(d.melt[,aet.names],na.rm=TRUE) # Average AET across all years
d.melt$Defnorm <- rowMeans(d.melt[,def.names],na.rm=TRUE) # Average CWD across all years
d.melt$Snownorm <- rowMeans(d.melt[,snow.names],na.rm=TRUE) # Average snow across all years

# Standard deviations of each year's value from the long-term climatic mean
d.melt$Tsd <- apply(d.melt[,tmean.names],1,sd)
d.melt$Psd <- apply(d.melt[,ppt.names],1,sd)
d.melt$AETsd <- apply(d.melt[,aet.names],1,sd)
d.melt$Defsd <- apply(d.melt[,def.names],1,sd)
d.melt$Snowsd <- apply(d.melt[,snow.names],1,sd)


#For each year, calculate that year's weather values and the values X years prior, where x is 1 through 5
v.year <- 2005:2015 #
for(year in v.year) {
  
  tmean.colname <- paste("tmean.ann.",year,sep="")
  ppt.colname <- paste("ppt.ann.",year,sep="")
  aet.colname <- paste("aet.ann.",year,sep="")
  def.colname <- paste("def.ann.",year,sep="")
  snow.colname <- paste("snow.late.",year,sep="")
  
  tmean.col <- grep(tmean.colname,names(d.melt))
  ppt.col <- grep(ppt.colname,names(d.melt))
  aet.col <- grep(aet.colname,names(d.melt))
  def.col <- grep(def.colname,names(d.melt))
  snow.col <- grep(snow.colname,names(d.melt))
  
  # Temperature; T0 means current year, T1 means one year prior, etc.
  d.melt[d.melt$year==year,"T0"] <- d.melt[d.melt$year==year,tmean.col]
  d.melt[d.melt$year==year,"T1"] <- d.melt[d.melt$year==year,tmean.col-1]
  d.melt[d.melt$year==year,"T2"] <- d.melt[d.melt$year==year,tmean.col-2]
  d.melt[d.melt$year==year,"T3"] <- d.melt[d.melt$year==year,tmean.col-3]
  d.melt[d.melt$year==year,"T4"] <- d.melt[d.melt$year==year,tmean.col-4]
  d.melt[d.melt$year==year,"T5"] <- d.melt[d.melt$year==year,tmean.col-5]
  
  # Precipitation
  d.melt[d.melt$year==year,"P0"] <- d.melt[d.melt$year==year,ppt.col]
  d.melt[d.melt$year==year,"P1"] <- d.melt[d.melt$year==year,ppt.col-1]
  d.melt[d.melt$year==year,"P2"] <- d.melt[d.melt$year==year,ppt.col-2]
  d.melt[d.melt$year==year,"P3"] <- d.melt[d.melt$year==year,ppt.col-3]
  d.melt[d.melt$year==year,"P4"] <- d.melt[d.melt$year==year,ppt.col-4]
  d.melt[d.melt$year==year,"P5"] <- d.melt[d.melt$year==year,ppt.col-5]
  
  # AET
  d.melt[d.melt$year==year,"AET0"] <- d.melt[d.melt$year==year,aet.col]
  d.melt[d.melt$year==year,"AET1"] <- d.melt[d.melt$year==year,aet.col-1]
  d.melt[d.melt$year==year,"AET2"] <- d.melt[d.melt$year==year,aet.col-2]
  d.melt[d.melt$year==year,"AET3"] <- d.melt[d.melt$year==year,aet.col-3]
  d.melt[d.melt$year==year,"AET4"] <- d.melt[d.melt$year==year,aet.col-4]
  d.melt[d.melt$year==year,"AET5"] <- d.melt[d.melt$year==year,aet.col-5]
  
  # CWD
  d.melt[d.melt$year==year,"Def0"] <- d.melt[d.melt$year==year,def.col]
  d.melt[d.melt$year==year,"Def1"] <- d.melt[d.melt$year==year,def.col-1]
  d.melt[d.melt$year==year,"Def2"] <- d.melt[d.melt$year==year,def.col-2]
  d.melt[d.melt$year==year,"Def3"] <- d.melt[d.melt$year==year,def.col-3]
  d.melt[d.melt$year==year,"Def4"] <- d.melt[d.melt$year==year,def.col-4]
  d.melt[d.melt$year==year,"Def5"] <- d.melt[d.melt$year==year,def.col-5]
  
  # Snow cover (not used in our analysis)
  d.melt[d.melt$year==year,"Snow0"] <- d.melt[d.melt$year==year,snow.col]
  d.melt[d.melt$year==year,"Snow1"] <- d.melt[d.melt$year==year,snow.col-1]
  d.melt[d.melt$year==year,"Snow2"] <- d.melt[d.melt$year==year,snow.col-2]
  d.melt[d.melt$year==year,"Snow3"] <- d.melt[d.melt$year==year,snow.col-3]
  d.melt[d.melt$year==year,"Snow4"] <- d.melt[d.melt$year==year,snow.col-4]
  d.melt[d.melt$year==year,"Snow5"] <- d.melt[d.melt$year==year,snow.col-5]
  
}



#Calculate Temp z.score (number of standard deviations a given year's climate value is from the mean)
d.melt$Tz0 <- (d.melt$T0-d.melt$Tnorm)/d.melt$Tsd
d.melt$Tz1 <- (d.melt$T1-d.melt$Tnorm)/d.melt$Tsd
d.melt$Tz2 <- (d.melt$T2-d.melt$Tnorm)/d.melt$Tsd
d.melt$Tz3 <- (d.melt$T3-d.melt$Tnorm)/d.melt$Tsd
d.melt$Tz4 <- (d.melt$T4-d.melt$Tnorm)/d.melt$Tsd
d.melt$Tz5 <- (d.melt$T5-d.melt$Tnorm)/d.melt$Tsd

#Calculate Precip z.score
d.melt$Pz0 <- (d.melt$P0-d.melt$Pnorm)/d.melt$Psd
d.melt$Pz1 <- (d.melt$P1-d.melt$Pnorm)/d.melt$Psd
d.melt$Pz2 <- (d.melt$P2-d.melt$Pnorm)/d.melt$Psd
d.melt$Pz3 <- (d.melt$P3-d.melt$Pnorm)/d.melt$Psd
d.melt$Pz4 <- (d.melt$P4-d.melt$Pnorm)/d.melt$Psd
d.melt$Pz5 <- (d.melt$P5-d.melt$Pnorm)/d.melt$Psd

#Calculate AET z.score
d.melt$AETz0 <- (d.melt$AET0-d.melt$AETnorm)/d.melt$AETsd
d.melt$AETz1 <- (d.melt$AET1-d.melt$AETnorm)/d.melt$AETsd
d.melt$AETz2 <- (d.melt$AET2-d.melt$AETnorm)/d.melt$AETsd
d.melt$AETz3 <- (d.melt$AET3-d.melt$AETnorm)/d.melt$AETsd
d.melt$AETz4 <- (d.melt$AET4-d.melt$AETnorm)/d.melt$AETsd
d.melt$AETz5 <- (d.melt$AET5-d.melt$AETnorm)/d.melt$AETsd

#Calculate Deficit z.score
d.melt$Defz0 <- (d.melt$Def0-d.melt$Defnorm)/d.melt$Defsd
d.melt$Defz1 <- (d.melt$Def1-d.melt$Defnorm)/d.melt$Defsd
d.melt$Defz2 <- (d.melt$Def2-d.melt$Defnorm)/d.melt$Defsd
d.melt$Defz3 <- (d.melt$Def3-d.melt$Defnorm)/d.melt$Defsd
d.melt$Defz4 <- (d.melt$Def4-d.melt$Defnorm)/d.melt$Defsd
d.melt$Defz5 <- (d.melt$Def5-d.melt$Defnorm)/d.melt$Defsd

#Calculate Snow z.score
d.melt$Snowz0 <- (d.melt$Snow0-d.melt$Snownorm)/d.melt$Snowsd
d.melt$Snowz1 <- (d.melt$Snow1-d.melt$Snownorm)/d.melt$Snowsd
d.melt$Snowz2 <- (d.melt$Snow2-d.melt$Snownorm)/d.melt$Snowsd
d.melt$Snowz3 <- (d.melt$Snow3-d.melt$Snownorm)/d.melt$Snowsd
d.melt$Snowz4 <- (d.melt$Snow4-d.melt$Snownorm)/d.melt$Snowsd
d.melt$Snowz5 <- (d.melt$Snow5-d.melt$Snownorm)/d.melt$Snowsd

# Remove raw climate data from dataframe
d.melt <- d.melt[,-c(tmean.cols,ppt.cols,aet.cols,def.cols,snow.cols)]

# Only look at cells that were flown and were within forested area
d.melt <- d.melt[d.melt$ca.mask == 1,]

# Calculate dead trees per hectare based on dead trees per acre
d.melt$mort.tph <- d.melt$mort.tpa * (1/0.405)

# Calculate the proportion of trees that died
d.melt$mort.prop <- d.melt$mort.tph/d.melt$live.tph

# Screen out "non-forest" plots where live.bah < 2. Also  screen out ~15 outliers with extremely high modeled live basal area.
d.melt=d.melt[between(d.melt$live.bah,1.999,60.1),]

# Only work with the (relevant) numeric cols
d.reduced=d.melt[,c("live.bah","Defnorm","alb.x","alb.y")]
d.unique <- unique(d.reduced)


# Derive deficit quantile (FDCI)
# This determines what quantile the CWD value for a given pixel falls into, within the full range of possible CWD values for pixels that have similar basal area values to the pixel in question (within 2.5 m2/ha)
# The parameter of this function is a data frame containing "live.bah" and "Defnorm" variables
def_quant_calc=function(x){  
  ba=x[which(names(x)=="live.bah")]
  def=x[which(names(x)=="Defnorm")]
  ba.range=c((ba-2.5),(ba+2.5))
  def.range=d.unique[between(d.unique$live.bah,ba.range[1],ba.range[2]),"Defnorm"]
  return(ecdf(def.range)(def))
}
d.melt$Defquant=apply(d.reduced,1,def_quant_calc) #Calculate the deficit quantiles. This takes about 5 minutes to run.


## Write table to external file to use for analysis
# In this table, there is one row per grid cell, per year surveyed. Not all cells were surveyed in all years.
# For a description of the values represented by each variable, see the data description accompanying the data file in the the same Dryad repository as this script.
write.table(d.melt,file="Young_et_al_Data.csv",sep=",",col.names=TRUE,row.names=FALSE)