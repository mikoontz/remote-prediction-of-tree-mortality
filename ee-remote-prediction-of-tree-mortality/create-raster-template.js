/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var sn = ee.FeatureCollection("ft:1VX1Uny1uAIYdgfuvsjLJVw1wMkhmSZq7fCQ0gowV"),
    modis = ee.ImageCollection("MODIS/006/MOD13Q1");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var evi = modis.filterDate("2017-01-01", "2017-04-01").
                select(["EVI"]).
                mosaic().
                clip(sn);
                
Map.addLayer(evi, {min: -2000, max: 10000, palette: ["ff0000", "00ff00"]});
Map.centerObject(sn);

print(evi.projection()); // Confirm projection

Export.image.toDrive({image: evi,
                      description: "sierra-nevada-250m-evi-template",
                      fileNamePrefix: "sierra-nevada-250m-evi-template",
                      scale: 250,
                      region: sn,
                      folder: "ee",
                      skipEmptyTiles: true
});
