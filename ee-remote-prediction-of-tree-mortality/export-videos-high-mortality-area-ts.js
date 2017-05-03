/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var sn = ee.FeatureCollection("ft:1VX1Uny1uAIYdgfuvsjLJVw1wMkhmSZq7fCQ0gowV"),
    modis = ee.ImageCollection("MODIS/006/MOD13Q1"),
    conifer_forest_v1 = ee.Image("users/mkoontz/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-type-no-mask"),
    geo = ee.FeatureCollection("ft:13cCOGKaj6ae1GScp4k6QsZrUo4IuhFNurz9EPL1N"),
    l8 = ee.ImageCollection("LANDSAT/LC8_SR"),
    l7 = ee.ImageCollection("LANDSAT/LE7_SR"),
    conifer_forest = ee.Image("users/mkoontz/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-type-no-mask-full-cell");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
Map.addLayer(geo); // An area of high mortality according to Young et al 2017

// For the MODIS 250m Vegetation Indices product
// This function will mask each image in the image collection one at a time.
// For each image, each pixel that is non-forest and/or of bad quality will be
// masked.
// Args: img, an image
// Returns: an image with some "masking" 

var quality_forest_mask = function(img) {
  var quality_img = img
                      .select("SummaryQA")
                      .rename("quality_pixel");
                      
  // The first video also masked out pixels of "marginal data" (when
  // the SummaryQA band was equal to 1
  // This mask includes "good data" and the "marginal data"
  // Comment out the line .or(quality_img.eq(1)) to get a more strict
  // mask, and then export the modis video called
  // "sn-high-mortality-area-evi-ts-modis-forest-quality-mask"
  var quality_pixels = quality_img
                        .eq(0)
                        .or(quality_img.eq(1));

  var clean_img = img
                    .select([
                    "sur_refl_b01", 
                    "sur_refl_b02", 
                    "sur_refl_b03"])
                    .updateMask(conifer_forest);
  
  clean_img = clean_img
                 .divide(10000)
                 .multiply(255)
                 .toUint8()
                 .clip(geo);

  return clean_img;

};


// The masked modis vegetation indices dataset is created by mapping the masking
// function to each image.
var m_modis = modis
                .filterDate("2012-01-01", "2017-04-01")
                .map(quality_forest_mask);

//   Export.video.toDrive({
//       collection: m_modis,
//       description: 'sn-high-mortality-area-evi-ts-modis-lower-quality-threshold-mask',
//       folder: 'ee/sierra-nevada-forest-quality-mask-time-series-videos',
//       framesPerSecond: 10,
//       scale: 250,
//       region: geo,
//       crs: 'EPSG:3310'
//   });

// var landsat_7_plus = function(img) {
//   return img
//             .select("B3", "B2", "B1")
//             .divide(1000)
//             .multiply(255)
//             .toUint8()
//             .clip(geo);
// };

// var l7p = l7
//           .filterDate("2012-04-01", "2017-04-15")
//           .filterBounds(geo)
//           .map(landsat_plus);                  

//   Export.video.toDrive({
//       collection: l7p,
//       description: 'sn-high-mortality-area-rgb-ts-landsat-7-no-mask',
//       folder: 'ee/sierra-nevada-forest-quality-mask-time-series-videos',
//       framesPerSecond: 10,
//       scale: 30,
//       region: geo,
//       crs: 'EPSG:3310'
//   });


// var landsat_8_plus = function(img) {
//   return img
//             .select("B4", "B3", "B2")
//             .divide(1000)
//             .multiply(255)
//             .toUint8()
//             .clip(geo);
// };

// var l8p = l8
//           .filterDate("2012-04-01", "2017-04-15")
//           .filterBounds(geo)
//           .map(landsat_plus);                  

//   Export.video.toDrive({
//       collection: l8p,
//       description: 'sn-high-mortality-area-rgb-ts-landsat-8-no-mask',
//       folder: 'ee/sierra-nevada-forest-quality-mask-time-series-videos',
//       framesPerSecond: 10,
//       scale: 30,
//       region: geo,
//       crs: 'EPSG:3310'
//   });