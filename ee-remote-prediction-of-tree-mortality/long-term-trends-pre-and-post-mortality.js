/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var jepcodes = ee.FeatureCollection("ft:1AE_cZiXvH2L7cXPsgvKse4y5-1cNDdt3kyCytHdV"),
    modis = ee.ImageCollection("MODIS/MOD13Q1"),
    sn = ee.FeatureCollection("ft:1VX1Uny1uAIYdgfuvsjLJVw1wMkhmSZq7fCQ0gowV"),
    nlcd = ee.Image("USGS/NLCD/NLCD2011");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var extent = modis.filterDate("2017-01-01", "2017-04-01").
                mosaic().
                clip(sn);

var quality_mask = function(img) {
  var the_mask = img.select(["SummaryQA"]);
  return img.mask(the_mask);
};

var landcover = nlcd.
                select(["landcover"]);

var the_forest_mask = landcover.
                        gte(41).
                        add(landcover.lte(43)).
                        lt(2);

var forest_mask = function(img) {
  return img.mask(the_forest_mask);
};

var m_modis = modis.map(quality_mask).
                    map(forest_mask);

var pre_mortality = m_modis.filterDate("2000-01-01", "2012-05-31").
                    select(["EVI"]).
                    formaTrend();
                    
var post_mortality = m_modis.filterDate("2012-05-31", "2017-04-20").
                    select(["EVI"]).
                    formaTrend();

var viridis = ["440154", "482878", "3E4A89", "31688E", "26828E", "1F9E89", "35B779", "6DCD59", "B4DE2C", "FDE725"];

Map.addLayer(pre_mortality.select(["long-trend"]).clip(sn), {min: -5, max: 5, palette: viridis}, "Pre-mortality trend");
Map.addLayer(pre_mortality.select(["long-tstat"]).clip(sn), {}, "Pre T Stat");
Map.addLayer(post_mortality.select(["long-trend"]).clip(sn), {min: -5, max: 5, palette: viridis}, "Post-mortality trend");
Map.addLayer(post_mortality.select(["long-tstat"]).clip(sn), {}, "Post T Stat");

Map.centerObject(sn);

// Export.image.toDrive({image: trend.select(["long-trend"]).clip(sn),
//                       description: "long-term-trend-export",
//                       folder: "ee",
//                       fileNamePrefix: "long-term-EVI-trend-sierra-nevada",
//                       region: sn,
//                       scale: 250});
                      
