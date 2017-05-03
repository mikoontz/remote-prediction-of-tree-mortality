/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var sn = ee.FeatureCollection("ft:1VX1Uny1uAIYdgfuvsjLJVw1wMkhmSZq7fCQ0gowV"),
    modis = ee.ImageCollection("MODIS/006/MOD13Q1"),
    conifer_forest_v1 = ee.Image("users/mkoontz/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-type-no-mask"),
    conifer_forest = ee.Image("users/mkoontz/sierra-nevada-250m-calveg-conifer-forested-pixels-by-whr-type-no-mask-full-cell");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// This function will mask each image in the image collection one at a time.
// For each image, each pixel that is non-forest and/or of bad quality will be
// masked.
// Args: img, an image
// Returns: an image with some "masking" which is really just filling all bands in with 1

var quality_forest_mask = function(img) {
  var date = ee.Date(img.get('system:time_start')).format("YYYYMMdd");
  date = ee.Number.parse(date);
  
  date = ee.Image(ee.Number(img.get('system:time_start')))
              .rename('date')
              .divide(1000);
  
  var quality_pixels = img
                  .select(["SummaryQA"])
                  .rename(["quality_pixel"])
                  .eq(0);

  var clean_img = img
                    .select(["NDVI", 
                    "EVI", 
                    "sur_refl_b01", 
                    "sur_refl_b02", 
                    "sur_refl_b03", 
                    "sur_refl_b07",
                    "ViewZenith",
                    "SolarZenith",
                    "RelativeAzimuth",
                    "DayOfYear",
                    "SummaryQA",
                    "DetailedQA"])
                    .addBands(date)
                    .updateMask(quality_pixels)
                    .updateMask(conifer_forest)
                    .toInt32();

  // In order to export a "masked" image, fill the masked pixels with a silly number
  // All "masked" pixels will have a value of 1 (or -9999 for the commented out code)
  // I switched it to 1 in an effort to save space?

  // First, generate a dummy image with the same number of bands as the MODIS image
  // var filled_img = ee.Image([ -9999, -9999, -9999, 
  //                           -9999, -9999, -9999, 
  //                           -9999, -9999, -9999, 
  //                           -9999, -9999, -9999])
  var filled_img = ee.Image([ 1, 1, 1, 
                              1, 1, 1, 
                              1, 1, 1, 
                              1, 1, 1, 1])
                    .rename(["NDVI", 
                    "EVI", 
                    "sur_refl_b01", 
                    "sur_refl_b02", 
                    "sur_refl_b03", 
                    "sur_refl_b07",
                    "ViewZenith",
                    "SolarZenith",
                    "RelativeAzimuth",
                    "DayOfYear",
                    "SummaryQA",
                    "DetailedQA",
                    "date"]);

  // Where there is a mask in the clean image, replace the value of the clean image with the 
  // value from the dummy image (with all 1 (or -9999 from the commented out code))
  var export_img = filled_img
                  .where(clean_img.mask(), clean_img)
                  .clip(sn.geometry().bounds());
  
  export_img = ee.Image(export_img.copyProperties({
    source: clean_img,
    properties: ["system:time_start"]
    }));
  
  return export_img;
};

// The masked modis vegetation indices dataset is created by mapping the masking
// function to each image.
var m_modis = modis.map(quality_forest_mask);

// Visualize
var viridis = ["440154", "482878", "3E4A89", "31688E", "26828E", "1F9E89", "35B779", "6DCD59", "B4DE2C", "FDE725"];

var img = ee.Image(m_modis.first());

Map.addLayer(img.select(["EVI"]), {min: 0, max: 3000, palette: viridis}, "Export img");
Map.addLayer(conifer_forest, {}, "Conservative forest mask");
Map.addLayer(conifer_forest_v1, {}, "Anti-conservative forest mask");

Map.centerObject(sn);

// There are 391 images
// print(m_modis.size().getInfo());

// var num_imgs = 391;
var first_img = 0;
var num_imgs = 391; // in the next round, first_img should equal this value; It will be 391 in the last round.

// Create a list with num_imgs empty list elements
var metadata = ee.List.sequence({
  start: 0,
  end: 390,
  step: 1
});

// Loop through all images, add the date to the appropriate list element, and export image
for (var i = first_img; i < num_imgs; i++) {
  var img = ee.Image(m_modis.toList(1, i).get(0));
  var date = ee.Date(img.get('system:time_start')).format("YYYYMMdd");
  var dict = ee.Dictionary({
    "image_num": i,
    "date": date
  });
  var ftr = ee.Feature(null, dict);
  metadata = metadata.set(i, ftr);
  
  Export.image.toDrive({
      image: img,
      description: 'sn-whole-ts-modis-forest-quality-mask-' + i,
      folder: 'ee/sierra-nevada-forest-quality-mask-modis-time-series',
      scale: 250,
      region: sn,
      crs: 'EPSG:3310'
  });
}

metadata = ee.FeatureCollection(metadata);

Export.table.toDrive({
  collection: metadata,
  description: "metadata",
  folder: "ee/sierra-nevada-forest-quality-mask-modis-time-series",
  fileNamePrefix: "metadata",
  fileFormat: "CSV"
});
