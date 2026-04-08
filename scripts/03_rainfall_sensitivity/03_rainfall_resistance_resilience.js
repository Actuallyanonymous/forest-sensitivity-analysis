/**
 * Forest Sensitivity Analysis Pipeline — Script 3
 * Heavy Rainfall Resistance & Resilience (Signed/Directional)
 *
 * For each forest pixel, computes mean signed resistance and resilience
 * across all anomalous rainfall years (z-score > threshold).
 *
 * Signed metrics:
 *   Resistance = Yn_bar / |Ye - Yn_bar| × sign(Ye - Yn_bar)
 *   Resilience = |Ye - Yn_bar| / |Y(e+1) - Yn_bar| × sign(Y(e+1) - Yn_bar)
 *
 * Positive resistance = forest maintained/grew during heavy rain
 * Negative resistance = forest declined during heavy rain
 * Positive resilience = forest recovered after heavy rain
 * Negative resilience = forest continued declining after heavy rain
 *
 * Rainfall anomaly detection:
 *   Heavy day = daily precipitation > long-term 95th percentile (CHIRPS)
 *   Anomalous year = z-score of annual heavy rain sum > Z_THRESHOLD
 *
 * Requires:
 *   - Forest mask asset from Script 1
 */

// CONFIGURATION :=


var TREE_COVER_ASSET = 'projects/cs5-pushkinmangla/assets/MP_Hybrid_Tree_Period_2003_2022';

var OUTPUT_ASSET_ID  = 'projects/cs5-pushkinmangla/assets/MP_Rain_Metrics';
var OUTPUT_DESC      = 'MP_Rain_Metrics';

var STATE_NAME       = 'Madhya Pradesh';
var START_YEAR       = 2004;
var END_YEAR         = 2022;
var Z_THRESHOLD      = 1.0;   // z-score above this = anomalous rainfall year


//  AOI :=


var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
            .geometry();

Map.centerObject(aoi, 7);


//  LOAD ASSETS :=

var treeMeta      = ee.Image(TREE_COVER_ASSET);
var startYearTree = treeMeta.select('start_year');
var endYearTree   = treeMeta.select('end_year');


// HEAVY RAINFALL ANOMALY DETECTION :=


var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
               .filterBounds(aoi)
               .filterDate('2000-01-01', '2023-12-31')
               .select('precipitation');

var proj = chirps.first().projection();

// Long-term 95th percentile threshold for heavy day definition
var p95 = chirps.reduce(ee.Reducer.percentile([95]))
                .setDefaultProjection(proj)
                .rename('p95');

// Annual sum of heavy-day precipitation
var years = ee.List.sequence(START_YEAR, END_YEAR + 1);

var annualHm = ee.ImageCollection(years.map(function(y) {
  var start    = ee.Date.fromYMD(y, 1, 1);
  var end      = ee.Date.fromYMD(y, 12, 31);
  var heavySum = chirps.filterDate(start, end)
                       .map(function(img) {
                         return img.multiply(img.gt(p95));
                       }).sum()
                        .setDefaultProjection(proj)
                        .rename('Hm')
                        .set('year', y);
  return heavySum;
}));

// Z-score of annual heavy rain sum across the period
var hmStats = annualHm.reduce(
  ee.Reducer.mean().combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true})
);

var anomalyCol = annualHm.map(function(img) {
  var zScore = img.subtract(hmStats.select('Hm_mean'))
                  .divide(hmStats.select('Hm_stdDev'))
                  .rename('zScore');
  return img.addBands(zScore).set('year', img.get('year'));
});


// LANDSAT NDVI :=


var maskClouds = function(image) {
  var qa   = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0)
               .and(qa.bitwiseAnd(1 << 4).eq(0));
  return image.updateMask(mask);
};

var getAnnualNDVI = function(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end   = ee.Date.fromYMD(year, 12, 31);

  var l89 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
              .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
              .filterDate(start, end).filterBounds(aoi)
              .map(maskClouds)
              .map(function(img) {
                return img.normalizedDifference(['SR_B5','SR_B4']).rename('ndvi');
              });

  var l57 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
              .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
              .filterDate(start, end).filterBounds(aoi)
              .map(maskClouds)
              .map(function(img) {
                return img.normalizedDifference(['SR_B4','SR_B3']).rename('ndvi');
              });

  return l89.merge(l57).median().set('year', year).rename('ndvi');
};

var ndviCol = ee.ImageCollection(years.map(getAnnualNDVI));


// BASELINE NDVI (Yn_bar) :=
// Mean NDVI across non-anomalous years only


var analysisYears = ee.List.sequence(START_YEAR, END_YEAR);

var Yn_bar = ee.ImageCollection(analysisYears.map(function(y) {
  var year    = ee.Number(y);
  var ndvi    = ee.Image(ndviCol.filter(ee.Filter.eq('year', year)).first());
  var anomaly = ee.Image(anomalyCol.filter(ee.Filter.eq('year', year)).first())
                  .reproject({crs: ndvi.projection(), scale: 30});

  var isNormal  = anomaly.select('zScore').abs().lt(Z_THRESHOLD);
  var isForest  = startYearTree.lte(year).and(endYearTree.gte(year));

  return ndvi.updateMask(isNormal.and(isForest)).set('year', year);
})).mean().rename('ndvi_baseline');

//===========================================================================
//                    7. SIGNED RESISTANCE & RESILIENCE
//===========================================================================

var metricsCol = ee.ImageCollection(analysisYears.map(function(y) {
  var year    = ee.Number(y);

  var ndviYe  = ee.Image(ndviCol.filter(ee.Filter.eq('year', year)).first());
  var anomaly = ee.Image(anomalyCol.filter(ee.Filter.eq('year', year)).first())
                  .reproject({crs: ndviYe.projection(), scale: 30});

  // Only compute on forest pixels during anomalous rainfall years
  var isAnomalous = anomaly.select('zScore').gt(Z_THRESHOLD);
  var isForest    = startYearTree.lte(year).and(endYearTree.gte(year));
  var mask        = isAnomalous.and(isForest);

  var diffRaw    = ndviYe.subtract(Yn_bar);
  var diffAbs    = diffRaw.abs();

  // Signed resistance — positive means forest held up or grew
  var resistance = Yn_bar.divide(diffAbs)
                         .multiply(diffRaw.signum())
                         .rename('resistance');

  // Signed resilience — positive means recovery, negative means continued decline
  var ndviNext   = ee.Image(ndviCol.filter(ee.Filter.eq('year', year.add(1))).first());
  var diffNext   = ndviNext.subtract(Yn_bar);
  var resilience = diffAbs.divide(diffNext.abs())
                          .multiply(diffNext.signum())
                          .rename('resilience');

  return ee.Image.cat([resistance, resilience])
    .updateMask(mask)
    .set('year', year);
}));


// AGGREGATE & EXPORT :=


var meanResist = metricsCol.select('resistance').mean().clip(aoi);
var meanResil  = metricsCol.select('resilience').mean().clip(aoi);

var finalOutput = meanResist.rename('resistance')
                            .addBands(meanResil.rename('resilience'));

// Preview
var visParams = {
  min: -3, max: 3,
  palette: ['8b0000','ff0000','ffffff','00ff00','006400']
};
Map.addLayer(meanResist, visParams, 'Rainfall resistance');
Map.addLayer(meanResil,  visParams, 'Rainfall resilience');

Export.image.toAsset({
  image       : finalOutput,
  description : OUTPUT_DESC,
  assetId     : OUTPUT_ASSET_ID,
  region      : aoi,
  scale       : 30,
  maxPixels   : 1e13
});
