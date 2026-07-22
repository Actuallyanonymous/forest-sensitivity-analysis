//===========================================================================
//                          1. CONFIGURATION
//===========================================================================

var STATE_NAME      = 'Madhya Pradesh';

// Change these two for the years you actually want results for.
var START_YEAR      = 2004;
var END_YEAR        = 2024;

// Fixed baseline window for zScore normalization — independent of
// START_YEAR/END_YEAR. I'm matching the same 2004-2024 window I locked
// in for rain (Script 3a) and drought (Script 2), so all my z-score-based
// indices are normalized against the same reference period.
// Do NOT change this once results are published, or old zScore bands
// will drift when I extend the analysis timeline later.
var BASELINE_START_YEAR = 2004;
var BASELINE_END_YEAR   = 2024;

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

var modisFire = ee.ImageCollection('MODIS/061/MOD14A1')
                  .filterBounds(aoi)
                  .select('MaxFRP');

var proj = modisFire.first().projection();

//===========================================================================
//                    4. ANNUAL METRICS CALCULATION
//===========================================================================

// We need annual metrics computed for every year covering BOTH the analysis
// window and the baseline window, whichever stretches further.
var metricsMinYear = Math.min(START_YEAR, BASELINE_START_YEAR);
var metricsMaxYear = Math.max(END_YEAR, BASELINE_END_YEAR);
var metricsYears   = ee.List.sequence(metricsMinYear, metricsMaxYear);

var annualMetrics = ee.ImageCollection(metricsYears.map(function(y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var end   = ee.Date.fromYMD(y, 12, 31);

  var yearCollection = modisFire.filterDate(start, end);

  var maxFRP = yearCollection.max()
                             .unmask(0)
                             .setDefaultProjection(proj)
                             .rename('maxFRP');

  var fireDaysCollection = yearCollection.map(function(img) {
    return img.updateMask(img.gt(FRP_THRESHOLD));
  });

  var sumFRP = fireDaysCollection.map(function(img) { return img.unmask(0); })
                                 .sum()
                                 .setDefaultProjection(proj)
                                 .rename('FRP');

  var fireDays = fireDaysCollection.map(function(img) { return img.mask(); })
                                   .sum()
                                   .unmask(0)
                                   .setDefaultProjection(proj)
                                   .rename('fireDays');

  var fireAvg = fireDaysCollection.mean()
                                  .unmask(0)
                                  .setDefaultProjection(proj)
                                  .rename('fireAvg');

  return sumFRP.addBands(fireDays)
               .addBands(fireAvg)
               .addBands(maxFRP)
               .set('year', y);
}));

//===========================================================================
//                    5. FIXED Z-SCORE CALCULATION
//===========================================================================

// Baseline stats only from the frozen baseline years, not from whatever
// analysis window I happen to be running right now.
var baselineYears = ee.List.sequence(BASELINE_START_YEAR, BASELINE_END_YEAR);

var baselineMetrics = annualMetrics.filter(
  ee.Filter.inList('year', baselineYears)
);

var frpMean   = baselineMetrics.select('FRP').mean();
var frpStdDev = baselineMetrics.select('FRP').reduce(ee.Reducer.stdDev());

var completedAnnualCollection = annualMetrics.map(function(img) {
  var frp = img.select('FRP');
  var z   = frp.subtract(frpMean).divide(frpStdDev).rename('zScore');
  return img.addBands(z);
});

//===========================================================================
//         6. SERVER-SIDE STACK INTO SINGLE MULTIBAND IMAGE & EXPORT
//===========================================================================

var analysisYears = ee.List.sequence(START_YEAR, END_YEAR);

var initialImage = ee.Image([]);
var outputImage = ee.Image(analysisYears.iterate(function(y, acc) {
  var yearStr = ee.String(ee.Number(y).toInt());
  var yearImg = completedAnnualCollection.filter(ee.Filter.eq('year', y)).first();

  var frpBand      = yearImg.select('FRP').rename(ee.String('FRP_').cat(yearStr));
  var zBand        = yearImg.select('zScore').rename(ee.String('zScore_').cat(yearStr));
  var fireDaysBand = yearImg.select('fireDays').rename(ee.String('fireDays_').cat(yearStr));
  var fireAvgBand  = yearImg.select('fireAvg').rename(ee.String('fireAvg_').cat(yearStr));
  var maxFRPBand   = yearImg.select('maxFRP').rename(ee.String('maxFRP_').cat(yearStr));

  return ee.Image(acc).addBands([frpBand, zBand, fireDaysBand, fireAvgBand, maxFRPBand]);
}, initialImage));

Export.image.toAsset({
  image       : outputImage.clip(aoi),
  description : OUTPUT_DESC,
  assetId     : OUTPUT_ASSET_ID,
  region      : aoi,
  scale       : 1000,
  crs         : 'EPSG:4326',
  maxPixels   : 1e13
});

print('Clean pipeline compilation verified.');
print('Ready to execute. Total structured bands: ' + analysisYears.length().multiply(5).getInfo());
