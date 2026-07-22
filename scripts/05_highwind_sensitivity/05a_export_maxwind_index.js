/**
 * Forest Sensitivity Analysis Pipeline — Script 5a
 * High Windspeed Index Export (Annual Max Hourly Windspeed, Hours > Threshold, Mean > Threshold)
 *
 * Computes three quantities per pixel per year and exports as a single
 * multiband asset — three bands per year:
 *
 * WSmax_{year}     = maximum hourly windspeed within the year (ERA5-Land)
 * WShoursGT_{year} = total hours where windspeed > WIND_THRESHOLD
 * WSmeanGT_{year}  = mean windspeed during the hours it exceeded WIND_THRESHOLD
 *
 * Windspeed computed from ERA5-Land hourly u/v 10m wind components:
 * windspeed = sqrt(u_component_of_wind_10m^2 + v_component_of_wind_10m^2)
 *
 * Output asset bands (e.g., for 2004-2022 = 19 years * 3 = 57 bands):
 * WSmax_2004, WShoursGT_2004, WSmeanGT_2004, ... 
 *
 * Requires: nothing — only public datasets (ERA5-Land Hourly)
 */

//===========================================================================
//                          1. CONFIGURATION
//===========================================================================

var STATE_NAME      = 'Andhra Pradesh';
var START_YEAR      = 2004;
var END_YEAR        = 2024;

// Set your wind speed threshold here (in m/s)
var WIND_THRESHOLD  = 10.0; 

var OUTPUT_ASSET_ID = 'projects/sura-496709/assets/AP_Wind_Index';
var OUTPUT_DESC     = 'AP_Wind_Index';

//===========================================================================
//                          2. AOI
//===========================================================================

var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
            .geometry();

Map.centerObject(aoi, 7);

//===========================================================================
//                    3. HOURLY WINDSPEED FROM U/V COMPONENTS
//===========================================================================

var era5Hourly = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')
                    .filterBounds(aoi)
                    .filterDate('2000-01-01', ee.Date.fromYMD(END_YEAR, 12, 31))
                    .select(['u_component_of_wind_10m', 'v_component_of_wind_10m']);

var proj = era5Hourly.first().projection();

var toWindSpeed = function(img) {
  var ws = img.select('u_component_of_wind_10m').pow(2)
              .add(img.select('v_component_of_wind_10m').pow(2))
              .sqrt()
              .rename('windspeed');
  return ws.copyProperties(img, ['system:time_start']);
};

var windSpeedCol = era5Hourly.map(toWindSpeed);

//===========================================================================
//                    4. ANNUAL METRICS (Max, Hours > Thresh, Mean > Thresh)
//===========================================================================

var years = ee.List.sequence(START_YEAR, END_YEAR);

var annualWSMetrics = ee.ImageCollection(years.map(function(y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  // Note: filterDate is exclusive on the end date. Using y+1 ensures Dec 31 is included.
  var end   = ee.Date.fromYMD(ee.Number(y).add(1), 1, 1);

  var yearCol = windSpeedCol.filterDate(start, end);

  // 1. Max Wind Speed
  var wsMax = yearCol.max()
                     .rename('WSmax');

  // 2. Number of hours wind speed > threshold
  var wsHoursGT = yearCol.map(function(img) {
                           return img.gt(WIND_THRESHOLD).rename('WShoursGT');
                         })
                         .sum();

  // 3. Mean wind speed when > threshold
  var wsMeanGT = yearCol.map(function(img) {
                          // Mask out pixels below the threshold before calculating mean
                          return img.updateMask(img.gt(WIND_THRESHOLD));
                        })
                        .mean()
                        .rename('WSmeanGT');

  // Combine all three into a single image for the year
  return ee.Image([wsMax, wsHoursGT, wsMeanGT])
           .setDefaultProjection(proj)
           .set('year', y);
}));

//===========================================================================
//                    5. STACK INTO SINGLE MULTIBAND IMAGE
//===========================================================================

// Build one image with bands for each year
var outputImage = ee.Image([]);

years.evaluate(function(yearList) {
  yearList.forEach(function(y) {
    var yearImg = annualWSMetrics.filter(ee.Filter.eq('year', y)).first();
    
    // Rename bands to include the year
    var wsMaxBand   = yearImg.select('WSmax').rename('WSmax_' + y);
    var wsHoursBand = yearImg.select('WShoursGT').rename('WShoursGT_' + y);
    var wsMeanBand  = yearImg.select('WSmeanGT').rename('WSmeanGT_' + y);

    outputImage = outputImage.addBands([wsMaxBand, wsHoursBand, wsMeanBand]);
  });

  // Quick preview (visualizing the max wind speed band for 2019)
  Map.addLayer(
    annualWSMetrics.filter(ee.Filter.eq('year', 2019)).first().select('WSmax'),
    {min: 0, max: 25, palette: ['white', 'orange', 'darkred']},
    'WSmax 2019 (preview)'
  );

  // Export
  Export.image.toAsset({
    image       : outputImage.clip(aoi),
    description : OUTPUT_DESC,
    assetId     : OUTPUT_ASSET_ID,
    region      : aoi,
    scale       : 11132,   // ERA5-Land native resolution ~11.1km
    crs         : 'EPSG:4326',
    maxPixels   : 1e13
  });

  print('✅ Export task ready — go to Tasks tab and click Run.');
  print('Output bands: WSmax_2004, WShoursGT_2004, WSmeanGT_2004 ... WSmax_2022, WShoursGT_2022, WSmeanGT_2022');
});
