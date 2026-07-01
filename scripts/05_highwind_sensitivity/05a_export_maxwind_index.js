/**
 * Forest Sensitivity Analysis Pipeline — Script 5a
 * High Windspeed Index Export (Annual Max Hourly Windspeed)
 *
 * Computes one quantity per pixel per year and exports as a single
 * multiband asset — one band per year:
 *
 * WSmax_{year} = maximum hourly windspeed within the year (ERA5-Land)
 *
 * Windspeed computed from ERA5-Land hourly u/v 10m wind components:
 * windspeed = sqrt(u_component_of_wind_10m^2 + v_component_of_wind_10m^2)
 *
 * Output asset bands:
 * WSmax_2004, WSmax_2005, ..., WSmax_2022      (19 bands)
 *
 * Requires: nothing — only public datasets (ERA5-Land Hourly)
 */

//===========================================================================
//                          1. CONFIGURATION
//===========================================================================

var STATE_NAME      = 'Andhra Pradesh';
var START_YEAR      = 2004;
var END_YEAR        = 2022;

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
                    .filterDate('2000-01-01', '2023-12-31')
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
//                    4. ANNUAL MAX WINDSPEED (WSmax)
//===========================================================================

var years = ee.List.sequence(START_YEAR, END_YEAR);

var annualWSmax = ee.ImageCollection(years.map(function(y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var end   = ee.Date.fromYMD(y, 12, 31);

  var wsMax = windSpeedCol.filterDate(start, end)
                          .max()
                          .setDefaultProjection(proj)
                          .rename('WSmax')
                          .set('year', y);
  return wsMax;
}));

//===========================================================================
//                    5. STACK INTO SINGLE MULTIBAND IMAGE
//===========================================================================

// Build one image with 19 named bands: WSmax_2004 ... WSmax_2022

var outputImage = ee.Image([]);

years.evaluate(function(yearList) {
  yearList.forEach(function(y) {
    var wsBand = annualWSmax.filter(ee.Filter.eq('year', y)).first()
                            .rename('WSmax_' + y);
    outputImage = outputImage.addBands(wsBand);
  });

  // Quick preview
  Map.addLayer(
    annualWSmax.filter(ee.Filter.eq('year', 2019)).first(),
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
  print('Output bands: WSmax_2004...WSmax_2022');
});