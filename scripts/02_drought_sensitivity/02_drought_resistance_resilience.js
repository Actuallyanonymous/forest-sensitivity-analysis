/**
 * Forest Sensitivity Analysis Pipeline — Script 2
 * Drought Resistance & Resilience
 *
 * For each forest pixel, computes mean resistance and resilience
 * across all drought years (SPEI-12 < threshold).
 *
 * Resistance  = Yn_bar / |Ye - Yn_bar|
 * Resilience  = |Ye - Yn_bar| / |Y(e+1) - Yn_bar|
 *
 * Where:
 *   Yn_bar = mean NDVI across non-drought years (baseline)
 *   Ye     = NDVI during drought year
 *   Y(e+1) = NDVI the year after drought
 *
 * Requires:
 *   - Forest mask asset from Script 1
 *   - SPEI-12 assets from spei-drought-analysis-pipeline
 */

// CONFIGURATION :=

var TREE_COVER_ASSET  = 'projects/cs5-pushkinmangla/assets/MP_Hybrid_Tree_Period_2003_2022';
var OUTPUT_ASSET_ID   = 'projects/cs5-pushkinmangla/assets/MP_Drought_Metrics';
var OUTPUT_DESC       = 'MP_Drought_Metrics';
var STATE_NAME        = 'Madhya Pradesh';
var START_YEAR        = 2004;
var END_YEAR          = 2022;
var DROUGHT_THRESHOLD = -1.0;   // SPEI-12 below this = drought year

// AOI :=

var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
            .geometry();

Map.centerObject(aoi, 7);

// Loading the assets :=

var treeMeta  = ee.Image(TREE_COVER_ASSET);
var startYear = treeMeta.select('start_year');
var endYear   = treeMeta.select('end_year');

// Load SPEI-12 collection from single multiband asset
var SPEI12_ASSET = 'projects/cs5-pushkinmangla/assets/SPEI12_Madhya_Pradesh';
var spei12_raw   = ee.Image(SPEI12_ASSET);
var spei12_bandnames = [];
for (var yn = 2004; yn <= 2023; yn++) {
  spei12_bandnames.push('y' + yn);
}
var spei12_named = spei12_raw.rename(spei12_bandnames);

// Build per-year SPEI collection
var speiImages = [];
for (var y = START_YEAR; y <= END_YEAR; y++) {
  speiImages.push(
    spei12_named.select('y' + y)
      .rename('spei')
      .set('year', y)
  );
}
var speiCol = ee.ImageCollection(speiImages);

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

// Load NDVI for START_YEAR to END_YEAR+1 (need next year for resilience)
var ndviYears = ee.List.sequence(START_YEAR, END_YEAR + 1);
var ndviCol   = ee.ImageCollection(ndviYears.map(getAnnualNDVI));

// BASELINE NDVI (Yn_bar) :=
// Mean NDVI across non-drought years only
// Uses simple masked ImageCollection mean — we are trying to avoid GEE array scalign issues here

var analysisYears = ee.List.sequence(START_YEAR, END_YEAR);

var ndviNonDrought = ee.ImageCollection(analysisYears.map(function(y) {
  var year = ee.Number(y);
  var ndvi = ndviCol.filter(ee.Filter.eq('year', year)).first();
  var spei = speiCol.filter(ee.Filter.eq('year', year)).first()
               .resample('bilinear')
               .reproject({crs: ndvi.projection(), scale: 30});
  var isNonDrought = spei.gte(DROUGHT_THRESHOLD);
  return ndvi.updateMask(isNonDrought).set('year', year);
}));

var Yn_bar = ndviNonDrought.mean().rename('ndvi_baseline');

// RESISTANCE & RESILIENCE :=

var metricsCol = ee.ImageCollection(analysisYears.map(function(y) {
  var year = ee.Number(y);

  var ndviYe = ndviCol.filter(ee.Filter.eq('year', year)).first();
  var speiYe = speiCol.filter(ee.Filter.eq('year', year)).first()
                 .resample('bilinear')
                 .reproject({crs: ndviYe.projection(), scale: 30});

  // Only compute on forest pixels during drought years
  var isForest  = startYear.lte(year).and(endYear.gte(year));
  var isDrought = speiYe.lt(DROUGHT_THRESHOLD);
  var mask      = isForest.and(isDrought);

  var diff       = ndviYe.subtract(Yn_bar).abs().max(1e-6);
  var resistance = Yn_bar.divide(diff).rename('resistance');

  var ndviNext = ndviCol.filter(ee.Filter.eq('year', year.add(1))).first();
  var diffNext = ndviNext.subtract(Yn_bar).abs().max(1e-6);
  var resilience = diff.divide(diffNext).rename('resilience');

  return ee.Image.cat([resistance, resilience])
    .updateMask(mask)
    .set('year', year);
}));

// AGGREGATE & EXPORT :=

var meanResistance = metricsCol.select('resistance').mean().clip(aoi);
var meanResilience = metricsCol.select('resilience').mean().clip(aoi);

var finalOutput = meanResistance.rename('resistance')
                                .addBands(meanResilience.rename('resilience'));

// Preview
var visParams = {min: 0, max: 10, palette: ['red','yellow','green','blue']};
Map.addLayer(meanResistance, visParams, 'Drought resistance');
Map.addLayer(meanResilience, visParams, 'Drought resilience');

Export.image.toAsset({
  image       : finalOutput,
  description : OUTPUT_DESC,
  assetId     : OUTPUT_ASSET_ID,
  region      : aoi,
  scale       : 30,
  maxPixels   : 1e13
});
