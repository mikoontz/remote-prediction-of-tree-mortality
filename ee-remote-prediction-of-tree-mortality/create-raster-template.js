/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var modis = ee.ImageCollection("MODIS/006/MOD13Q1"),
    sn = ee.FeatureCollection("ft:1vdDUTu09Rkw5qKR_DSfmFX-b_7kqy4E-pjxg9Sq6");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var evi = modis.filterDate("2017-01-01", "2017-04-01").
                select(["EVI"]).
                mosaic().
                clip(sn);

var template_img = evi.
                    mask();
                    
Map.addLayer(template_img);

Map.addLayer(evi, {min: -2000, max: 10000, palette: ["ff0000", "00ff00"]});
Map.centerObject(sn);

// Input is EPSG:4326 by default for .kml files upon import,
// so that's what we start with, but we want an exported
// product that is best for our California focus
// Export as California Albers projection

Export.image.toDrive({image: template_img,
                      description: "sierra-nevada-250m-evi-template",
                      fileNamePrefix: "sierra-nevada-250m-evi-template",
                      scale: 250,
                      region: sn.geometry().getInfo()['coordinates'],
                      crs: 'EPSG:3310',
                      folder: "ee"
});