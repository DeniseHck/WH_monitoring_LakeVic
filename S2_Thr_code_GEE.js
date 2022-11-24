/*
BEFORE STARTING:
//--------------------------------------------------------------------------------
Import assets:
- a Point (34.53, -0.27)
- imageVisParam2 (image visualisation parameters)
- Winam Gulf boundary
*/

//SENTINEL CLOUD MASKING
//--------------------------------------------------------------------------------
function maskS2clouds(image) {
    var qa = image.select('QA60');
    //Bits 10 and 11 are clouds and cirrus, respectively
    var cloudBitMask = 1 << 10;
    var cirrusBitMask = 1 << 11;
    //Both flags should be set to zero, indicating clear conditions
    var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(qa.bitwiseAnd(cirrusBitMask).eq(0));
    //Return masked and scaled data, without the QA bands
    return image.updateMask(mask).divide(10000)/*.select('B.*')*/.copyProperties(image, ["system:time_start"]);
}


//FILTERING IMAGE COLLECTION and adding CLOUD MASK
//---------------------------------------------------------------------------------
var collection = ee.ImageCollection("COPERNICUS/S2_SR")
.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 5))
.filterBounds(table)
.filterBounds(point)
.filterDate('2018-01-01', '2021-12-31')
.filter(ee.Filter.neq('system:index','20200216T075011_20200216T080428_T36MXE')) //Eliminate 16 feb 2021 due to clouds
.filter(ee.Filter.neq('system:index','20200725T074621_20200725T080832_T36MXE'))
.sort("system:time_start")
.map(maskS2clouds);


// a test
var addTime = function(image) {
  return image.addBands(image.getNumber("system:time_start"));
};

//count & print the number of images in this collection.
//----------------------------------------------------------------------------------
var count = collection.size();
print('Count: ', count);


// Get the date range of images in the collection.
//----------------------------------------------------------------------------------
var range = collection.reduceColumns(ee.Reducer.minMax(), ["system:time_start"]);
print('Date range: ', ee.Date(range.get('min')), ee.Date(range.get('max')));


//IMAGE STYLE
//----------------------------------------------------------------------------------
var trueColor = {
  bands: ['B4', 'B3', 'B2'],
  min: 0,
  max: 1200,
};


var colors = ["000000","ffff00", "00cc00"];


//IMAGE PROCESSING
//----------------------------------------------------------------------------------------------

var addClass = function(image) {
  // define bands of interest
  var nir = image.select('B8');
  var green = image.select('B3');
  var red = image.select('B4');
  var swir = image.select('B11');
  var nirA = image.select('B8A');//20m
  var red5 = image.select('B5'); //20m
  var red6 = image.select('B6'); //20m
  var red7 = image.select('B7'); //20m

  var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
  var reNDVI = nirA.subtract(red6).divide(nirA.add(red6)).rename("RENDVI");
  var afai = nir.subtract(red).add(swir.subtract(red)).multiply(0.5).rename("AFAI");

  //NDVI mask (<0 = veg)
  var threshold = ndvi.gt(0);
  var mask = threshold.updateMask(threshold);
  var threshold2 = reNDVI.gt(0.04);
  var mask2 =  threshold2.updateMask(threshold2);



  //reNDVI reclass
  var reclassRENDVI = ee.Image(2)
    .where(reNDVI.gt(-2).and(reNDVI.lte(-0.3)), 0) //gt: greater than, lte:less than or equal to
    .where(reNDVI.gt(-0.3).and(reNDVI.lte(0.04)), 2)
    .where(reNDVI.gt(0.04), 1)
    .updateMask(mask).clip(table)
    .updateMask(mask2)
    .rename("RECLASSIFIED");

  var threshold3 = reNDVI.lte(0.04);
  var mask3 =  threshold3.updateMask(threshold3);
  var threshold4 = afai.gt(0);
  var mask4=threshold4.updateMask(threshold4);

  var reclassRENDVI_cyano = ee.Image(2)
    .where(reNDVI.gt(-2).and(reNDVI.lte(-0.3)), 0) //gt: greater than, lte:less than or equal to
    .where(reNDVI.gt(-0.3).and(reNDVI.lte(0.04)), 1)
    .where(reNDVI.gt(0.04), 1)
    .updateMask(mask).clip(table)
    .updateMask(mask3)
    .rename("RECLASSIFIED_cyano");


  return image.addBands(reclassRENDVI).addBands(reclassRENDVI_cyano);
};


//MEASURE WH IN KM2
//----------------------------------------------------
var addMeasure = function(img) {
  var WH_cover = img.select("RECLASSIFIED").eq(1);
  var cyano = img.select("RECLASSIFIED").eq(2);

  var areaImage_WH = WH_cover.multiply(ee.Image.pixelArea()).rename("areaImage_WH").divide(1e6);
  var areaImage_CY = cyano.multiply(ee.Image.pixelArea()).rename("areaImage_cyano").divide(1e6);

  //adding area of WH as band
  var test = img.addBands(areaImage_WH);

  var statsWH = areaImage_WH.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: table.geometry(),
    scale: 20,
    maxPixels: 1e10
  });


  var statscyano = areaImage_CY.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: table.geometry(),
    scale: 20,
    maxPixels: 1e10
  });

  return img.set(statsWH).set(statscyano);
}


//Mapped reclass to whole collection
var mapped = collection.map(addClass).map(addMeasure);
print(mapped);


var image1 = mapped.first();
var nat1 = image1;
var classified1 = image1.select("RECLASSIFIED");
var classified2 = image1.select("RECLASSIFIED_cyano");

//Timeseries
//---------------------------------------------------------------------

//use seriesByRegion instead - it works

var WHtimeseries1 = ui.Chart.image.seriesByRegion(
  mapped,               //the collection
  table.geometry(),     //the aoi
  ee.Reducer.sum(),     //reducer type
  "RECLASSIFIED",       //band selected
  20)                  //scale


var plotWH = WHtimeseries1
.setSeriesNames(["Water Hyacinth"])
.setChartType("LineChart")
//.setSeriesNames(["WH", "cyanobacteria"])
.setOptions({
  interpolateNulls: true,
  lineWidth: 1,
  colors: ['violet'],
  pointSize: 3,
  title: "Water Hyacinth evolution 2018-2022 from Sentinel-2 data",
  hAxis: {title: "Date"},
  vAxis: {title: "Water Hyacinth pixels"}
});

print(plotWH);


var WHtimeseries2 = ui.Chart.image.seriesByRegion(
  mapped,               //the collection
  table.geometry(),     //the aoi
  ee.Reducer.sum(),     //reducer type
  "RECLASSIFIED_cyano",       //band selected
  20)                  //scale

var plotCyano = WHtimeseries2
.setSeriesNames(["cyanobacteria"])
.setChartType("LineChart")
//.setSeriesNames(["WH", "cyanobacteria"])
.setOptions({
  interpolateNulls: true,
  lineWidth: 1,
  colors: ['greenyellow'],
  pointSize: 3,
  title: "Algal bloom evolution 2018-2022 from Sentinel-2 data",
  hAxis: {title: "Date"},
  vAxis: {title: "Cyanobacteria pixels"}
});

print(plotCyano);


//Export.table.toDrive({collection: mapped});


//make a video
//----------------------------------------------------------

var plotted = mapped.select("RECLASSIFIED");
var projection = image1.select("B3").projection().getInfo();

/*
// Define arguments for animation function parameters.
var videoArgs = {
  dimensions: 768,
  region: table,
  framesPerSecond: 7,
  crs: projection.crs,
  crsTransform: projection.transform("EPSG:4326"),
  min: -40.0,
  max: 35.0,
};

print(ui.Thumbnail(plotted, videoArgs));
*/

//print(ui.Thumbnail(plotted, videoArgs));

// Alternatively, print a URL that will produce the animation when accessed.
//print(plotted.getVideoThumbURL(videoArgs));


//DATAVIZ
//---------------------------------------------

Map.addLayer(nat1, imageVisParam2, 'Natural color');
Map.addLayer(classified1, {palette:colors}, 'Water hyacinth');
Map.addLayer(classified2, {palette:"greenyellow"}, 'cyanobacteria')


//MAKE CHART INTERACTIVE
//----------------------------------------------------


var label = ui.Label('Click a point on the chart to show the image for that date.');
Map.add(label);


//Create callbakc function that adds image to the map coresponding with clicked data point on chart
plotWH.onClick(function(xValue, yValue, seriesName) {
    if (!xValue) return;  // Selection was cleared.

    // Show the image for the clicked date.
    var equalDate = ee.Filter.equals('system:time_start', xValue);
    //Find image coresponding with clicked data and clip water classification to roi
    var classification = ee.Image(mapped.filter(equalDate).first()).select("RECLASSIFIED");
    var classification_cyano = ee.Image(mapped.filter(equalDate).first()).select("RECLASSIFIED_cyano");
    var S2image = ee.Image(mapped.filter(equalDate).first());

    //Make map layer based on SAR image, reset the map layers, and add this new layer
    var Layer = ui.Map.Layer(S2image, imageVisParam2, 'TrueColor');

    Map.layers().reset([Layer]);
    var visParamsWH = {
      min: 0,
      max: 1,
      palette: ['violet']
    }

     var visParamsCyano = {
      min: 0,
      max: 1,
      palette: ['greenyellow']
    }
    //Add water classification on top of SAR image
    Map.addLayer(classification,visParamsWH,'Water Hyacinth')
    Map.addLayer(classification_cyano, visParamsCyano, "cyanobacteria")



    // Show a label with the date on the map.
    label.setValue((new Date(xValue)).toUTCString());
  });
/*
plotCyano.onClick(function(xValue, yValue, seriesName) {
    if (!xValue) return;  // Selection was cleared.

    // Show the image for the clicked date.
    var equalDate = ee.Filter.equals('system:time_start', xValue);
    //Find image coresponding with clicked data and clip water classification to roi
    var classification_cyano = ee.Image(mapped.filter(equalDate).first()).select("RECLASSIFIED_cyano");
    var classification = ee.Image(mapped.filter(equalDate).first()).select("RECLASSIFIED");
    var S2image = ee.Image(mapped.filter(equalDate).first());

    //Make map layer based on SAR image, reset the map layers, and add this new layer
    var Layer = ui.Map.Layer(S2image, imageVisParam2, 'TrueColor');

    Map.layers().reset([Layer]);
    var visParamsWH = {
      min: 0,
      max: 1,
      palette: ['violet']
    }

     var visParamsCyano = {
      min: 0,
      max: 1,
      palette: ['orange']
    }
    //Add water classification on top of SAR image
    Map.addLayer(classification_cyano,visParamsCyano,'Cyanobacteria')
    Map.addLayer(classification,visParamsWH,'Water hyacinth')


    // Show a label with the date on the map.
    label.setValue((new Date(xValue)).toUTCString());
  });
*/


/*
Map.addLayer(nat2, trueColor, 'NAT1');
Map.addLayer(classified2, {palette:"blue"});

Map.addLayer(nat3, trueColor, 'NAT1');
Map.addLayer(classified3, {palette:"red"});
*/



/*
//IMAGE PROCESSING
//---------------------------------------------------------------------------------


//COLOR PALETTE
//---------------------------------------------------------------------

*/
var colors = ["ee82ee", "adff2f"];   //undefined errod
var names = ["water hyacinth", "cyanobacteria"];


//ADD LEGEND
//---------------------------------------------------------------------

// set position of panel
var legend = ui.Panel({
  style: {
    position: 'bottom-right',
    padding: '8px 15px'
  }
});

// Create legend title
var legendTitle = ui.Label({
  value: 'Key',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});

// Add the title to the panel
legend.add(legendTitle);

// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {

      // Create the label that is actually the colored box.
      var colorBox = ui.Label({
        style: {
          backgroundColor: '#'+ color,
          // Use padding to give the box height and width.
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });

      // Create the label filled with the description text.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });

      // return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};


// Add color and and names
for (var i = 0; i< 2;i++) {
  legend.add(makeRow(colors[i], names[i]));
  }


//DISPLAY RESULTS
//-----------------------------------------------------------------------------------

// add legend to map (alternatively you can also print the legend to the console)
Map.add(legend);

/*
//EXPORT - Export a cloud-optimized GeoTIFF.
//-------------------------------------------------------------------------------------



Export.image.toDrive({
  image: nat1.select("B4", "B3", "B2"),
  description: '2022_03_07_S2_TrueColor',
  crs: "EPSG:32636",
  //crs: projection.crs,
  //crsTransform: projection.transform,
  scale: 20,
  region: table
});

*/
Export.image.toDrive({
  image: classified1,
  description: '2021_12_22_ThrClass',
  crs: "EPSG:32636",
  //crs: projection.crs,
  //crsTransform: projection.transform,
  scale: 20,
  region: table
});
