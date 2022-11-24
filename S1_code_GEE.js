/*
BEFORE BEGINNING:
-----------------------------------------------------------------------------------------------
- import Winam Gulf Boundary assets ("table")

*/


// Load the Sentinel-1 ImageCollection, filter to Jun-Sep 2020 observations.
var sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD')
                    .filterDate('2018-05-15', '2018-05-31')
                    .filterBounds(table)
                    .filter(ee.Filter.eq('instrumentMode', 'IW'))
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'));

//var projection = sentinel1.projection().getInfo();

var vvVhIw = sentinel1
  // Filter to get images with VV and VH dual polarization.
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))

  // Filter to get images collected in interferometric wide swath mode.
  .filter(ee.Filter.eq('instrumentMode', 'IW'));

//print(sentinel1);

// Separate ascending and descending orbit images into distinct collections.
var vvVhIwAsc = vvVhIw.filter(
  ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));

var vvVhIwDesc = vvVhIw.filter(
  ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));

// Calculate temporal means for various observations to use for visualization.
// Mean VV for combined ascending and descending image collections.
var vvIwAscDesc = vvVhIwAsc.merge(vvVhIwDesc).select('VV').mean();

var clipped = vvIwAscDesc/*.clip(table)*/;


// Lee filter
//---------------------------------------------------------------------------//
/** Lee Filter applied to one image. It is implemented as described in
 J. S. Lee, “Digital image enhancement and noise filtering by use of local statistics,”
 IEEE Pattern Anal. Machine Intell., vol. PAMI-2, pp. 165–168, Mar. 1980.*/
 /*
var leefilter = function(image,KERNEL_SIZE) {
        var bandNames = image.bandNames().remove('angle');
        //S1-GRD images are multilooked 5 times in range
        var enl = 5
        // Compute the speckle standard deviation
        var eta = 1.0/Math.sqrt(enl);
        eta = ee.Image.constant(eta);

        // MMSE estimator
        // Neighbourhood mean and variance
        var oneImg = ee.Image.constant(1);

        var reducers = ee.Reducer.mean().combine({
                      reducer2: ee.Reducer.variance(),
                      sharedInputs: true
                      });
        var stats = image.select(bandNames).reduceNeighborhood({reducer: reducers,kernel: ee.Kernel.square(KERNEL_SIZE/2,'pixels'), optimization: 'window'})
        var meanBand = bandNames.map(function(bandName){return ee.String(bandName).cat('_mean')});
        var varBand = bandNames.map(function(bandName){return ee.String(bandName).cat('_variance')});

        var z_bar = stats.select(meanBand);
        var varz = stats.select(varBand);

        // Estimate weight
        var varx = (varz.subtract(z_bar.pow(2).multiply(eta.pow(2)))).divide(oneImg.add(eta.pow(2)));
        var b = varx.divide(varz);

        //if b is negative set it to zero
        var new_b = b.where(b.lt(0), 0)
        var output = oneImg.subtract(new_b).multiply(z_bar.abs()).add(new_b.multiply(image.select(bandNames)));
        output = output.rename(bandNames);
        return image.addBands(output, null, true);
  }
*/

// Refined Lee filter
//---------------------------------------------------------------------------//
/** This filter is modified from the implementation by Guido Lemoine
 * Source: Lemoine et al.; https://code.earthengine.google.com/5d1ed0a0f0417f098fdfd2fa137c3d0c */

 var refinedLee = function(image) {

    var bandNames = image.bandNames().remove('angle');

    var result = ee.ImageCollection(bandNames.map(function(b){
    var img = image.select([b]);

    // img must be linear, i.e. not in dB!
    // Set up 3x3 kernels
    var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
    var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);

    var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
    var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);

    // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
    var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);

    var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);

    // Calculate mean and variance for the sampled windows and store as 9 bands
    var sample_mean = mean3.neighborhoodToBands(sample_kernel);
    var sample_var = variance3.neighborhoodToBands(sample_kernel);

    // Determine the 4 gradients for the sampled windows
    var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
    gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
    gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
    gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());

    // And find the maximum gradient amongst gradient bands
    var max_gradient = gradients.reduce(ee.Reducer.max());

    // Create a mask for band pixels that are the maximum gradient
    var gradmask = gradients.eq(max_gradient);

    // duplicate gradmask bands: each gradient represents 2 directions
    gradmask = gradmask.addBands(gradmask);

    // Determine the 8 directions
    var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
    directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
    directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
    directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
    // The next 4 are the not() of the previous 4
    directions = directions.addBands(directions.select(0).not().multiply(5));
    directions = directions.addBands(directions.select(1).not().multiply(6));
    directions = directions.addBands(directions.select(2).not().multiply(7));
    directions = directions.addBands(directions.select(3).not().multiply(8));

    // Mask all values that are not 1-8
    directions = directions.updateMask(gradmask);

    // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
    directions = directions.reduce(ee.Reducer.sum());

    //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
    //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

    var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));

    // Calculate localNoiseVariance
    var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

    // Set up the 7*7 kernels for directional statistics
    var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));

    var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0],
      [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

    var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
    var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);

    // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
    var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
    var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));

    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));

    // and add the bands for rotated kernels
    for (var i=1; i<4; i++) {
      dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
      dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
      dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
      dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    }

    // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
    dir_mean = dir_mean.reduce(ee.Reducer.sum());
    dir_var = dir_var.reduce(ee.Reducer.sum());

    // A finally generate the filtered value
    var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));

    var b = varX.divide(dir_var);

    return dir_mean.add(b.multiply(img.subtract(dir_mean)))
      .arrayProject([0])
      // Get a multi-band image bands.
      .arrayFlatten([['sum']])
      .float();
  })).toBands().rename(bandNames).copyProperties(image);
  return image.addBands(result, null, true)
  }



//LEGEND
//-------------------------------------------------------------------------------

//COLOR PALETTE
//---------------------------------------------------------------------
var colors = ["00cc00","000000",];
//var colors = ["309e48", "ffad16","fff80c"]
var names = ["floating vegetation", "open water"];


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
  value: '',
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

var filtered = refinedLee(clipped, 3);


//IMAGE THRESHOLDING
//---------------------------------------------------------------
//var wh_unfiltered = clipped.lt(10).rename('Water')
//var wh = filtered.gt(10).rename("WH")  //Identify all pixels below threshold and set them equal to 1. All other pixels set to 0
var wh =  ee.Image(1)
          .where(filtered.gt(-100).and(filtered.lte(-16)), 2)
          .where(filtered.gt(-16).and(filtered.lte(100)), 1);


//try to convert to array, sort, get cummulative sum
//convert back to image, divide by max value, and get values great than .9


Map.setCenter(34.5, -0.25, 10);

var leecolors = ["black", "white"];

// Compute the histogram of the  band.  The mean and variance are only FYI.
var histogram = filtered.reduceRegion({
  reducer: ee.Reducer.histogram(255, 2)
      .combine('mean', null, true)
      .combine('variance', null, true),
  geometry: table,
  scale: 20,
  bestEffort: true
});




//print(Chart.image.histogram(filtered, table, 20));





/*
var glcm = wh.glcmTexture({size: 4});
var contrast = glcm.select('N_contrast');
Map.addLayer(contrast,
             {min: 0, max: 1500, palette: ['0000CC', 'CC0000']},
             'contrast');
   */


//VISUALISE DATA
//-----------------------------------------------------------------------

Map.addLayer(clipped.clip(table), {min:-20, max:-5, palette: ["black", "white"]});
//Map.addLayer(filtered, {min:-20, max:-6, palette: ["black", "white"]});


//Map.addLayer(wh_unfiltered, {palette:colors});
Map.addLayer(filtered.clip(table), {min:-20, max:-5, palette:leecolors},"refinedLee" );

Map.addLayer(wh.clip(table), {palette: colors},"Thresholded" );
Map.add(legend);


//MEASURE & PLOT
//--------------------------------------------------------------------------

var WH_cover = wh.eq(1).clip(table);

var areaImage = WH_cover.multiply(ee.Image.pixelArea())

print(WH_cover)

var area = areaImage.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: table.geometry(),
  scale: 20,
  maxPixels: 1e10
  });

var WHSqKm = ee.Number(
  area.get('constant')).divide(1e6)
print(WHSqKm);

Map.setOptions("SATELLITE");
// Export the image, specifying the CRS, transform, and region.
var projection = filtered.projection().getInfo();

Export.image.toDrive({
  image: filtered,
  description: '2018_yearly_avg',
  crs: "EPSG:32636",
  //crs: projection.crs,
  //crsTransform: projection.transform,
  scale: 10,
  region: table
});

/*
Export.image.toDrive({
  image: filtered,
  description: '2021_mean_S1',
  crs: projection.crs,
  crsTransform: projection.transform,
  region: table
});

*/
