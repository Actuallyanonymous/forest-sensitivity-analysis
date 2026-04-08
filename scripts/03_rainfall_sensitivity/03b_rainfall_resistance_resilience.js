/**
 * Forest Sensitivity Analysis Pipeline — Script 3b
 * Heavy Rainfall Resistance & Resilience
 *
 * Loads the rainfall index asset from Script 3a, applies z-score
 * threshold to identify anomalous years, then computes signed
 * resistance and resilience for forest pixels.
 *
 * Signed metrics:
 *   Resistance = Yn_bar / |Ye - Yn_bar| × sign(Ye - Yn_bar)
 *   Resilience = |Ye - Yn_bar| / |Y(e+1) - Yn_bar| × sign(Y(e+1) - Yn_bar)
 *
 * Positive = forest maintained/recovered
 * Negative = forest declined/continued declining
 *
 * Requires:
 *   - Forest mask asset (Script 1)
 *   - Rainfall index asset (Script 3a)
 */

// COnfiguration :=

var TREE_COVER_ASSET  = 'projects/cs5-pushkinmangla/assets/MP_Hybrid_Tree_Period_2003_2022';
var RAIN_INDEX_ASSET  = 'projects/cs5-pushkinmangla/assets/MP_Rain_Index';

var OUTPUT_ASSET_ID   = 'projects/cs5-pushkinmangla/assets/MP_Rain_Metrics';
var OUTPUT_DESC       = 'MP_Rain_Metrics';

var STATE_NAME        = 'Madhya Pradesh';
var START_YEAR        = 2004;
var END_YEAR          = 2022;
var Z_THRESHOLD       = 1.0;   // z-score above this = anomalous rainfall year

//Choosing the AOI :=

var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
            .geometry();

Map.centerObject(aoi, 7);

//Loading the assets :=

var treeMeta      = ee.Image(TREE_COVER_ASSET);
var startYearTree = treeMeta.select('start_year');
var endYearTree   = treeMeta.select('end_year');

var rainIndex     = ee.Image(RAIN_INDEX_ASSET);

// Reconstruct per-year collections from the multiband asset
var analysisYears = ee.List.sequence(START_YEAR, END_YEAR);

var hmCol = ee.ImageCollection(
  ee.List.sequence(START_YEAR, END_YEAR).map(function(y) {
    var year = ee.Number(y).format('%.0f');
    return rainIndex.select(ee.String('Hm_').cat(year))
                    .rename('Hm')
                    .set('year', ee.Number(y));
  })
);

var zScoreCol = ee.ImageCollection(
  ee.List.sequence(START_YEAR, END_YEAR).map(function(y) {
    var year = ee.Number(y).format('%.0f');
    return rainIndex.select(ee.String('zScore_').cat(year))
                    .rename('zScore')
                    .set('year', ee.Number(y));
  })
);

//tkaing the ndvi from landsat :=

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

// Need END_YEAR+1 for resilience (next year NDVI)
var ndviCol = ee.ImageCollection(
  ee.List.sequence(START_YEAR, END_YEAR + 1).map(getAnnualNDVI)
);


// BASELINE NDVI (Yn_bar)
// Mean NDVI across non-anomalous years only


var Yn_bar = ee.ImageCollection(analysisYears.map(function(y) {
  var year   = ee.Number(y);
  var ndvi   = ee.Image(ndviCol.filter(ee.Filter.eq('year', year)).first());
  var zScore = ee.Image(zScoreCol.filter(ee.Filter.eq('year', year)).first())
                 .reproject({crs: ndvi.projection(), scale: 30});

  var isNormal = zScore.select('zScore').abs().lt(Z_THRESHOLD);
  var isForest = startYearTree.lte(year).and(endYearTree.gte(year));

  return ndvi.updateMask(isNormal.and(isForest)).set('year', year);
})).mean().rename('ndvi_baseline');


// SIGNED RESISTANCE & RESILIENCE


var metricsCol = ee.ImageCollection(analysisYears.map(function(y) {
  var year   = ee.Number(y);

  var ndviYe = ee.Image(ndviCol.filter(ee.Filter.eq('year', year)).first());
  var zScore = ee.Image(zScoreCol.filter(ee.Filter.eq('year', year)).first())
                 .reproject({crs: ndviYe.projection(), scale: 30});

  // Only compute on forest pixels during anomalous rainfall years
  var isAnomalous = zScore.select('zScore').gt(Z_THRESHOLD);
  var isForest    = startYearTree.lte(year).and(endYearTree.gte(year));
  var mask        = isAnomalous.and(isForest);

  var diffRaw    = ndviYe.subtract(Yn_bar);
  var diffAbs    = diffRaw.abs();

  var resistance = Yn_bar.divide(diffAbs)
                         .multiply(diffRaw.signum())
                         .rename('resistance');

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
