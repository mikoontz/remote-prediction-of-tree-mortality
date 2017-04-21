# Overview
Using time series analysis of satellite imagery to predict Sierra Nevada, California tree mortality due to the 2012-2016 drought.

# Data Sources
[MODIS/006/MOD13Q1 -- MODIS Satellite Imagery -- Vegetation Indices 16-Day Product v006](https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1)

Accessed via Google Earth Engine. Template generated and exported in order to match other data products to it (see create-raster-template.js script). Final masking of the Sierra Nevada Ecoregion done using `R` (see create-raster-template.R script). Extraction of time series done with `python` (pending).

Sierra Nevada TNC Ecoregion

Need more information on the original source. 

We converted this vector shapefile from ESRI Shapefile to a KML using QGIS so that we could upload it as a Google Fusion Table for use in Google Earth Engine. The Fusion Table is publically available [here](https://www.google.com/fusiontables/DataSource?docid=1VX1Uny1uAIYdgfuvsjLJVw1wMkhmSZq7fCQ0gowV).

[CalVeg Vegetation Classification](https://www.fs.usda.gov/detail/r5/landmanagement/resourcemanagement/?cid=stelprdb5347192)
[Wildlife Habitat Relationship Lifeform description](https://www.fs.fed.us/r5/rsl/projects/classification/cv-cwhr-xwalk.html)
[Wildlife Habitat Relationship Type description](http://frap.fire.ca.gov/projects/frap_veg/classification)

We downloaded the geodatabase files for Zone 3, North Sierra and Zone 4, South Sierra, imported them into `R` using the `sf` package, subsetted them by Wildlife Habitat Relationship Types that included conifers, and rasterized them to the MODIS raster template from Google Earth Engine (see the rasterize-calveg-forest-polygons.R script).
