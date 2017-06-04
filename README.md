# Overview
Using time series analysis of satellite imagery to predict Sierra Nevada, California tree mortality due to the 2012-2016 drought.

# Data Sources
## MODIS Satellite Imagery
[MODIS/006/MOD13Q1 -- MODIS Satellite Imagery -- Vegetation Indices 16-Day Product v006](https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1)

Accessed via Google Earth Engine. Template generated and exported in order to match other data products to it (see create-raster-template.js script). Final masking of the Sierra Nevada Ecoregion done using `R` (see create-raster-template.R script). Extraction of time series done with `python` (pending).

## Sierra Nevada Ecoregion according to Jepson "Geographic Subdivisions of California"

Jepson Flora Project (eds.) 2016. Jepson eFlora, (http://ucjeps.berkeley.edu/eflora/)[http://ucjeps.berkeley.edu/eflora/] [accessed on Mar 07, 2016]

For more information, and to request the GIS file related to the subdivisions, visit (http://ucjeps.berkeley.edu/eflora/geography.html)[http://ucjeps.berkeley.edu/eflora/geography.html].

We derived our map of the Sierra Nevada by subsetting the Geographic Subdivisions of California mpa to include the following subdivisions:
+ Central High Sierra Nevada District
+ Central Sierra Nevada Foothills District
+ Northern Sierra Nevada Foothills District
+ Southern Sierra Nevada Foothills District
+ Northern High Sierra Nevada District 
+ Southern High Sierra Nevada District 
+ Tehachapi Mountain Area Subregion

The derived product is found here: (https://www.google.com/fusiontables/DataSource?docid=1vdDUTu09Rkw5qKR_DSfmFX-b_7kqy4E-pjxg9Sq6)[https://www.google.com/fusiontables/DataSource?docid=1vdDUTu09Rkw5qKR_DSfmFX-b_7kqy4E-pjxg9Sq6]

## CalVeg Vegetation Classification
(https://www.fs.usda.gov/detail/r5/landmanagement/resourcemanagement/?cid=stelprdb5347192)

[Wildlife Habitat Relationship Lifeform description](https://www.fs.fed.us/r5/rsl/projects/classification/cv-cwhr-xwalk.html)

[Wildlife Habitat Relationship Type description](http://frap.fire.ca.gov/projects/frap_veg/classification)

We downloaded the geodatabase files for Zone 3, North Sierra and Zone 4, South Sierra, imported them into `R` using the `sf` package, subsetted them by Wildlife Habitat Relationship Types that included conifers, and rasterized them to the MODIS raster template from Google Earth Engine (see the rasterize-calveg-forest-polygons.R script).

# Data Workflow

1. Jepson Geographic Subdivisions of California (jepcodes-v7.kmz/EPSG:4326)

1. Open in Google Earth, save as .kml (jepcodes-v7.kml/EPSG:4326)

1. Read into R (create-raster-template.R), scrape HTML to get attribute table, subset to desired Sierra Nevada regions, export as .kml (same projection) and .shp (transform to California Albers (SierraEcoregions_Jepson.kml/EPSG:4326; SierraEcoregions_Jepson.shp/EPSG:3310)

1. Upload SierraEcoregions_Jepson.kml/EPSG:4326 to Fusion Table (ft:1vdDUTu09Rkw5qKR_DSfmFX-b_7kqy4E-pjxg9Sq6)

1. Read Fusion Table into Earth Engine, get a single MODIS EVI image at 250m spatial resolution, clip to the SierraEcoregion_Jepson.kml/EPSG:4326 polygon, export to Google Drive as sierra-nevada-250m-evi-template.tif/EPSG:3310

1. Download from Google Drive to local repository (remote-prediction-of-tree-mortality/features/sierra-nevada-250m-evi-template.tif/EPSG:3310)



