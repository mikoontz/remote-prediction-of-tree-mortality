/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var sn = ee.FeatureCollection("ft:1VX1Uny1uAIYdgfuvsjLJVw1wMkhmSZq7fCQ0gowV"),
    modis = ee.ImageCollection("MODIS/006/MOD13Q1"),
    conifer_forest = ee.Image("users/mkoontz/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-type-no-mask");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// This function will mask each image in the image collection one at a time.
// For each image, each pixel that is non-forest and/or of bad quality will be
// masked.
// Args: img, an image
// Returns: an image with some masking

var quality_forest_mask = function(img) {
  var the_mask = img
                    .select(["SummaryQA"]) // Look at the SummaryQA band
                    .eq(0)                 // Return 1 if value is 0 (good quality)
                    .and(conifer_forest);  // ... AND if pixel is forested. 
  
  return img.updateMask(the_mask);         // Mask pixels with a value of 0
};

// The masked modis vegetation indices dataset is created by mapping the masking
// function to each image.
var m_modis = modis.map(quality_forest_mask);

Map.addLayer(ee.Image(m_modis.first()));
Map.addLayer(ee.Image(modis.first()).clip(sn));
Map.centerObject(sn);