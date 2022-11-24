/*
BEFORE STARTING:
---------------------------------------------------------------------------------
Import assets:
- point (43.56, -0.25)
- Winam Gulf boundary
- sample points (WH, algae, water) x RF
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
    return image.updateMask(mask).divide(10000).copyProperties(image, ["system:time_start"]);
}


//FILTERING IMAGE COLLECTION and adding CLOUD MASK
//---------------------------------------------------------------------------------
var collection = ee.ImageCollection("COPERNICUS/S2_SR")
.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 5))
.filterBounds(table)
.filterBounds(point)
.filterDate('2018-01-01', '2022-04-01')
.filter(ee.Filter.neq('system:index','20200216T075011_20200216T080428_T36MXE')) //Eliminate 16 feb 2021 due to clouds
.filter(ee.Filter.neq('system:index','20200725T074621_20200725T080832_T36MXE'))
.sort("system:time_start")
.map(maskS2clouds);


var image = collection.first();

// a test
var addTime = function(image) {
  return image.addBands(image.getNumber("system:time_start"));
};


//CLASSIFICATION (superivised)
//from https://www.youtube.com/watch?v=moFOIpm3JcI
//------------------------------------

var bands = ['B2', 'B3', 'B4', 'B5', 'B7', 'B8A', 'B11', 'B12']; //After Ongore et al. 2018

var points = table2.merge(table3).merge(table5).merge(table6).merge(sparsearse);



//one for training one for validating
var sample = points.randomColumn();
var trainingsample = sample.filter("random <= 0.8"); //change number to get more or less
var validationsample = sample.filter("random > 0.8");

print(trainingsample, "training sample");
print(validationsample, "validation sample");

/*
to export sampling points to assets for later cloud use

Export.table.toAsset({
  collection:sparsearse,
  description: "Sample points wh x RF",
  assetId: "WH"

});

// Export the FeatureCollection to a KML file.

Export.table.toDrive({
  collection: points,
  description:'RF_samples',
  fileFormat: 'SHP'
});
*/

// Sample rasters to create training data
var training = image.sampleRegions({
  collection: trainingsample,
  properties: ['Class'],
  scale: 20
});

print(training, "training data band values");


var validation = image.sampleRegions({
  collection: validationsample,
  properties: ['Class'],
  scale:20
});


//RF classifier model building
var RFclassifier = ee.Classifier.smileRandomForest(500).train(training, "Class");


var classfunc = function(image) {
  var nir = image.select('B8');
  var red = image.select('B4');

  var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');

  var threshold = ndvi.gt(0.2);
  var mask = threshold.updateMask(threshold);

  var classified = image.classify(RFclassifier).updateMask(mask);
  var nl2012 = image.addBands(classified);
  var l2022 = nl2012.select("classification")

  return image.addBands(classified);
}

var mapped = collection.map(classfunc);


//if low accuracy, reiterate the process changing trees

//CHECK ACCURACY
var trainAccuracy = RFclassifier.confusionMatrix();
//print("training error matrix", trainAccuracy);
//print("Training overall accuracy", trainAccuracy.accuracy());

//Get confusion matrix and overall accuracy for validation sample
validation = validation.classify(RFclassifier);

var validationAccuracy = validation.errorMatrix("Class", "classification");
//print("validation error matrix", validationAccuracy);
//print("Validation Accuracy", validationAccuracy.accuracy());


//Variable importance
//--------------------------

var explain = RFclassifier.explain();
print(explain);

var variable_importance = ee.Feature(null, ee.Dictionary(explain).get("importance"));
var chartTitle = "Random forest: bands of major importance";
var chart = ui.Chart.feature.byProperty(variable_importance)
  .setChartType("BarChart")
  .setOptions({
    title:chartTitle
  });

print(chart)

//IMAGE STYLE
//----------------------------------------------------------------------------------
var trueColor = {
  bands: ['B4', 'B3', 'B2'],
  min: 0.03,
  max: 0.112,
};

var ndviParams = {min: -1, max: 1, palette: ['black','lightgreen']};
var styling = {color:'magenta', fillColor: '00000000'};

var palette = [
  'violet',
  'greenyellow',
  'black'
];

Map.addLayer(image.clip(table), trueColor, 'RGB_least_cloud');

var classified = mapped.first().select("classification");


Map.addLayer(classified.clip(table), {palette:palette, min:1, max:3}, "CLASSIFIED")


var addWHonly = function(img) {
  var WH_cover = img.select("classification").eq(1).rename("WH");
  var cyano = img.select("classification").eq(2).rename("Algae");

  return img.addBands(WH_cover).addBands(cyano);
}

var mapled = mapped.map(addWHonly)
print(mapled);

// Define the chart and print it to the console.
var chartWH =
    ui.Chart.image
        .seriesByRegion({
          imageCollection: mapled,
          band: 'WH',
          regions: table, //table works
          reducer: ee.Reducer.sum(),
          scale: 20,
          seriesProperty: 'Class',
          xProperty: 'system:time_start'
        })
        .setOptions({
          title: 'Water Hyacinth in Winam Gulf 2018-2022',
          hAxis: {title: 'Date', titleTextStyle: {italic: false, bold: true}},
          vAxis: {
            title: 'Pixels',
            titleTextStyle: {italic: false, bold: true}
          },
          lineWidth: 5,
          colors: ['violet'],
        });
print(chartWH);

// Define the chart and print it to the console.
var chartCyano =
    ui.Chart.image
        .seriesByRegion({
          imageCollection: mapled,
          band: 'Algae',
          regions: table, //table works
          reducer: ee.Reducer.sum(),
          scale: 20,
          seriesProperty: 'Class',
          xProperty: 'system:time_start'
        })
        .setOptions({
          title: 'Water Hyacinth in Winam Gulf 2018-2022',
          hAxis: {title: 'Date', titleTextStyle: {italic: false, bold: true}},
          vAxis: {
            title: 'Pixels',
            titleTextStyle: {italic: false, bold: true}
          },
          lineWidth: 5,
          colors: ['greenyellow'],
        });
print(chartCyano);

var label = ui.Label('Click a point on the chart to show the image for that date.');
Map.add(label);


chartWH.onClick(function(xValue, yValue, seriesName) {
    if (!xValue) return;  // Selection was cleared.

    // Show the image for the clicked date.
    var equalDate = ee.Filter.equals('system:time_start', xValue);
    //Find image coresponding with clicked data and clip water classification to roi
    var classification = ee.Image(mapled.filter(equalDate).first()).select("classification");
    //var classification_cyano = ee.Image(mapped.filter(equalDate).first()).select("RECLASSIFIED_cyano");
    var S2image = ee.Image(mapled.filter(equalDate).first());

    //Make map layer based on SAR image, reset the map layers, and add this new layer
    var Layer = ui.Map.Layer(S2image, trueColor, 'TrueColor');

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
    Map.addLayer(classification.clip(table),{palette:palette, min:1, max:3},'Water Hyacinth')
    //Map.addLayer(classification_cyano, visParamsCyano, "cyanobacteria")



    // Show a label with the date on the map.
    label.setValue((new Date(xValue)).toUTCString());
  });

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

Map.add(legend);
Map.setOptions("SATELLITE");

/*
Map.addLayer(table2)
Map.addLayer(table3)
Map.addLayer(table4)
*/


//EXPORT - Export a cloud-optimized GeoTIFF.
//-------------------------------------------------------------------------------------


/*
Export.image.toDrive({
  image: nat1.select("B4", "B3", "B2"),
  description: '2022_03_07_S2_TrueColor',
  crs: "EPSG:32636",
  //crs: projection.crs,
  //crsTransform: projection.transform,
  scale: 20,
  region: table
});

var img1 = mapped.filterDate('2021-12-22', '2022-04-01').first().clip(table);
Map.addLayer(img1.select("classification"), {palette:palette, min:1, max:3})

Export.image.toDrive({
  image: img1.select("B4", "B3", "B2"),
  description: '2021_12_22_S2_TrueColor',
  crs: "EPSG:32636",
  //crs: projection.crs,
  //crsTransform: projection.transform,
  scale: 20,
  region: table
});


Export.image.toDrive({
  image: img1.select("classification"),
  description: '2021_12_22_S2_Classified_RF',
  crs: "EPSG:32636",
  //crs: projection.crs,
  //crsTransform: projection.transform,
  scale: 20,
  region: table
});

*/
