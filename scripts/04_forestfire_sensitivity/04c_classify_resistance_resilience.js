/**
 * Fire — Script 4c
 * Keep per-year values from -20 to 20, then mean what remains.
 * Values outside that range are masked and do not enter the mean.
 * Exact ±20 is kept. Output bands are resistance and resilience only.
 * Resistance: all fire years. Resilience: eligible_YYYY only (Ye < 0.95*Yn).
 */

var YEARLY_ASSET = 'projects/cs5-pushkinmangla/assets/MP_Fire_Metrics_Harmonized_kNDVI';
var OUTPUT_ASSET_ID = 'projects/cs5-pushkinmangla/assets/MP_Fire_Metrics_Filtered';
var OUTPUT_DESC = 'MP_Fire_Metrics_Filtered';
var STATE_NAME = 'Madhya Pradesh';
var START_YEAR = 2004;
var END_YEAR = 2024;

var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
  .geometry();
var img = ee.Image(YEARLY_ASSET);

var keepInRange = function(metric) {
  return metric.updateMask(metric.gte(-20).and(metric.lte(20)));
};

var resistList = [];
var resilList = [];
for (var y = START_YEAR; y <= END_YEAR; y++) {
  resistList.push(keepInRange(img.select('resistance_' + y).rename('resistance')));
  resilList.push(keepInRange(
    img.select('resilience_' + y).rename('resilience')
      .updateMask(img.select('eligible_' + y))
  ));
}

var meanResistance = ee.ImageCollection(resistList).mean().rename('resistance');
var meanResilience = ee.ImageCollection(resilList).mean().rename('resilience');

var finalOutput = meanResistance.addBands(meanResilience).clip(aoi);

Export.image.toAsset({
  image: finalOutput,
  description: OUTPUT_DESC,
  assetId: OUTPUT_ASSET_ID,
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});
