/**
 * Forest Sensitivity Analysis Pipeline — Script 3a
 * Heavy Rainfall Index Export
 *
 * Computes two quantities per pixel per year and exports as a single
 * multiband asset — one band per year for each quantity:
 *
 *   Hm_{year}     = annual sum of precipitation on heavy days
 *   zScore_{year} = z-score of Hm relative to the full period mean/stddev
 *
 * Heavy day definition: daily precipitation > long-term 95th percentile (CHIRPS)
 * Z-score computed across all years in the period.
 *
 * Output asset bands:
 *   Hm_2004, Hm_2005, ..., Hm_2023      (19 bands)
 *   zScore_2004, ..., zScore_2023        (19 bands)
 *   Total: 38 bands
 *
 * This asset is the direct input to Script 3b, analogous to how
 * SPEI assets are the input to Script 2 (drought).
 *
 * Requires: nothing — only public datasets (CHIRPS)
 */

//===========================================================================
//                          1. CONFIGURATION
//===========================================================================

var STATE_NAME      = 'Madhya Pradesh';
var START_YEAR      = 2004;
var END_YEAR        = 2022;

var OUTPUT_ASSET_ID = 'projects/cs5-pushkinmangla/assets/MP_Rain_Index';
var OUTPUT_DESC     = 'MP_Rain_Index';

//===========================================================================
//                          2. AOI
//===========================================================================

var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
            .geometry();

Map.centerObject(aoi, 7);

//===========================================================================
//                    3. HEAVY RAINFALL INDEX (Hm)
//===========================================================================

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
               .filterBounds(aoi)
               .filterDate('2000-01-01', '2023-12-31')
               .select('precipitation');

var proj = chirps.first().projection();

// Long-term 95th percentile — defines what counts as a heavy day
var p95 = chirps.reduce(ee.Reducer.percentile([95]))
                .setDefaultProjection(proj)
                .rename('p95');

// Annual heavy rain sum per year
var years = ee.List.sequence(START_YEAR, END_YEAR);

var annualHm = ee.ImageCollection(years.map(function(y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var end   = ee.Date.fromYMD(y, 12, 31);

  var heavySum = chirps.filterDate(start, end)
                       .map(function(img) {
                         return img.multiply(img.gt(p95));
                       })
                       .sum()
                       .setDefaultProjection(proj)
                       .rename('Hm')
                       .set('year', y);
  return heavySum;
}));

//===========================================================================
//                    4. Z-SCORE ACROSS ALL YEARS
//===========================================================================

var hmMean   = annualHm.mean().rename('Hm_mean');
var hmStdDev = annualHm.reduce(ee.Reducer.stdDev()).rename('Hm_stdDev');

var annualZScore = ee.ImageCollection(years.map(function(y) {
  var year = ee.Number(y);
  var hm   = annualHm.filter(ee.Filter.eq('year', year)).first();
  var z    = hm.subtract(hmMean).divide(hmStdDev).rename('zScore');
  return z.set('year', year);
}));

//===========================================================================
//                    5. STACK INTO SINGLE MULTIBAND IMAGE
//===========================================================================

// Build one image with 38 named bands:
// Hm_2004 ... Hm_2022, zScore_2004 ... zScore_2022

var outputImage = ee.Image([]);

years.evaluate(function(yearList) {
  yearList.forEach(function(y) {
    var hmBand = annualHm.filter(ee.Filter.eq('year', y)).first()
                         .rename('Hm_' + y);
    var zBand  = annualZScore.filter(ee.Filter.eq('year', y)).first()
                             .rename('zScore_' + y);
    outputImage = outputImage.addBands(hmBand).addBands(zBand);
  });

  // Quick preview
  Map.addLayer(
    annualZScore.filter(ee.Filter.eq('year', 2019)).first(),
    {min: -2, max: 2, palette: ['white','blue','darkblue']},
    'Z-score 2019 (preview)'
  );

  // Export
  Export.image.toAsset({
    image       : outputImage.clip(aoi),
    description : OUTPUT_DESC,
    assetId     : OUTPUT_ASSET_ID,
    region      : aoi,
    scale       : 5566,   // CHIRPS native resolution ~5.5km
    crs         : 'EPSG:4326',
    maxPixels   : 1e13
  });

  print('✅ Export task ready — go to Tasks tab and click Run.');
  print('Output bands: Hm_2004...Hm_2022, zScore_2004...zScore_2022');
});
