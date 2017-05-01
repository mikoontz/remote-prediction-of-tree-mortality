/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var sn = ee.FeatureCollection("ft:1VX1Uny1uAIYdgfuvsjLJVw1wMkhmSZq7fCQ0gowV"),
    modis = ee.ImageCollection("MODIS/006/MOD13Q1"),
    conifer_forest = ee.Image("users/mkoontz/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-type-no-mask");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// This function will mask each image in the image collection one at a time.
// For each image, each pixel that is non-forest and/or of bad quality will be
// masked.
// Args: img, an image
// Returns: an image with some "masking" 

var quality_forest_mask = function(img) {
  var quality_pixels = img
                  .select(["SummaryQA"])
                  .rename(["quality_pixel"])
                  .eq(0);

  var clean_img = img
                    .select([
                    "sur_refl_b01", 
                    "sur_refl_b02", 
                    "sur_refl_b03"])
                    .updateMask(quality_pixels)
                    .updateMask(conifer_forest);
  
  clean_img = clean_img
                 .divide(10000)
                 .multiply(255)
                 .toUint8()
                 .clip(sn);

  return clean_img;

  
};

// The masked modis vegetation indices dataset is created by mapping the masking
// function to each image.
var m_modis = modis.map(quality_forest_mask);

  Export.video.toDrive({
      collection: m_modis,
      description: 'sn-whole-evi-ts-modis-forest-quality-mask',
      folder: 'ee/sierra-nevada-forested-quality-pixel-time-series-videos',
      framesPerSecond: 10,
      scale: 250,
      region: sn,
      crs: 'EPSG:3310'
  });