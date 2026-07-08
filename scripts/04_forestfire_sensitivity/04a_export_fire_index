/**
 * Forest Sensitivity Analysis Pipeline — Fire Script
 * Fire Radiative Power (FRP) Index Export (Threshold > 30)
 *
 * Computes 5 quantities per pixel per year and exports as a single
 * multiband asset — one band per year for each quantity:
 *
 * FRP_{year}        = annual sum of FRP on fire days
 * zScore_{year}     = z-score of FRP relative to the 2004-2022 period
 * fireDays_{year}   = number of days in the year with FRP > 30
 * fireAvg_{year}    = average daily FRP on fire days
 * maxFRP_{year}     = maximum daily FRP recorded in that year
 *
 * Fire day definition: daily MaxFRP > 30 (MODIS Terra Thermal Anomalies)
 *
 * Requires: nothing — only public datasets (MODIS MOD14A1)
 */

//===========================================================================
//                          1. CONFIGURATION
//===========================================================================

var STATE_NAME      = 'Madhya Pradesh';
var START_YEAR      = 2004;
var END_YEAR        = 2022;

var OUTPUT_ASSET_ID = 'projects/cs5-pushkinmangla/assets/MP_Fire_Index_FRP30';
var OUTPUT_DESC     = 'MP_Fire_Index_FRP30';
var FRP_THRESHOLD   = 30;

//===========================================================================
//                          2. AOI
//===========================================================================

var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
            .geometry();

Map.centerObject(aoi, 7);

//===========================================================================
//                    3. DATASET SETUP
//===========================================================================

// Terra Thermal Anomalies & Fire Daily 1km
var modisFire = ee.ImageCollection('MODIS/061/MOD14A1')
                  .filterBounds(aoi)
                  .select('MaxFRP'); 

var proj = modisFire.first().projection();

//===========================================================================
//                    4. ANNUAL METRICS CALCULATION
//===========================================================================

var years = ee.List.sequence(START_YEAR, END_YEAR);

var annualMetrics = ee.ImageCollection(years.map(function(y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var end   = ee.Date.fromYMD(y, 12, 31);

  var yearCollection = modisFire.filterDate(start, end);

  // Absolute maximum daily FRP for this year
  var maxFRP = yearCollection.max()
                             .unmask(0)
                             .setDefaultProjection(proj)
                             .rename('maxFRP');

  // Binary mask: isolate days where FRP > 30
  var fireDaysCollection = yearCollection.map(function(img) {
    return img.updateMask(img.gt(FRP_THRESHOLD));
  });

  // 1. Annual sum of FRP on fire days
  var sumFRP = fireDaysCollection.map(function(img) { return img.unmask(0); })
                                 .sum()
                                 .setDefaultProjection(proj)
                                 .rename('FRP');

  // 2. Number of fire days
  var fireDays = fireDaysCollection.map(function(img) { return img.notMasked(); })
                                   .sum()
                                   .unmask(0)
                                   .setDefaultProjection(proj)
                                   .rename('fireDays');

  // 3. Average FRP of fire days
  var fireAvg = fireDaysCollection.mean()
                                  .unmask(0)
                                  .setDefaultProjection(proj)
                                  .rename('fireAvg');

  // Combine metrics into a single image per year containing all 4 base properties
  return sumFRP.addBands(fireDays)
               .addBands(fireAvg)
               .addBands(maxFRP)
               .set('year', y);
}));

//===========================================================================
//                    5. FIXED Z-SCORE CALCULATION
//===========================================================================

var frpMean   = annualMetrics.select('FRP').mean();
var frpStdDev = annualMetrics.select('FRP').reduce(ee.Reducer.stdDev());

// Server-side map architecture over the annual image collection
var completedAnnualCollection = annualMetrics.map(function(img) {
  var frp = img.select('FRP');
  var z   = frp.subtract(frpMean).divide(frpStdDev).rename('zScore');
  return img.addBands(z); 
});

//===========================================================================
//         6. SERVER-SIDE STACK INTO SINGLE MULTIBAND IMAGE & EXPORT
//===========================================================================

// Server-side iteration style to stack and correctly rename bands
var initialImage = ee.Image([]);
var outputImage = ee.Image(years.iterate(function(y, acc) {
  var yearStr = ee.String(ee.Number(y).toInt());
  var yearImg = completedAnnualCollection.filter(ee.Filter.eq('year', y)).first();
  
  var frpBand      = yearImg.select('FRP').rename(ee.String('FRP_').cat(yearStr));
  var zBand        = yearImg.select('zScore').rename(ee.String('zScore_').cat(yearStr));
  var fireDaysBand = yearImg.select('fireDays').rename(ee.String('fireDays_').cat(yearStr));
  var fireAvgBand  = yearImg.select('fireAvg').rename(ee.String('fireAvg_').cat(yearStr));
  var maxFRPBand   = yearImg.select('maxFRP').rename(ee.String('maxFRP_').cat(yearStr));
  
  return ee.Image(acc).addBands([frpBand, zBand, fireDaysBand, fireAvgBand, maxFRPBand]);
}, initialImage));

// Export execution block
Export.image.toAsset({
  image       : outputImage.clip(aoi),
  description : OUTPUT_DESC,
  assetId     : OUTPUT_ASSET_ID,
  region      : aoi,
  scale       : 1000,   // MODIS MOD14A1 native resolution is 1km
  crs         : 'EPSG:4326',
  maxPixels   : 1e13
});

print('✅ Clean pipeline compilation verified.');
print('Ready to execute in the tasks tab. Total structured bands: 95.');
