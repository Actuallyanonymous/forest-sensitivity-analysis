/**
 * Wind — Script 5c
 * Mean the per-year bands from Script 5b, then class each pixel.
 * Resistance: all high-wind years. Resilience: eligible_YYYY only (Ye < 0.95*Yn).
 * Class codes: 1 = low, 2 = medium, 3 = high. Cuts at ±20.
 */

var YEARLY_ASSET = 'projects/sura-496709/assets/AP_Wind_Metrics_Harmonized_kNDVI';
var OUTPUT_ASSET_ID = 'projects/sura-496709/assets/AP_Wind_Metrics_Classified';
var OUTPUT_DESC = 'AP_Wind_Metrics_Classified';
var STATE_NAME = 'Andhra Pradesh';
var START_YEAR = 2004;
var END_YEAR = 2024;

var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
  .geometry();
var img = ee.Image(YEARLY_ASSET);

var resistList = [];
var resilList = [];
for (var y = START_YEAR; y <= END_YEAR; y++) {
  resistList.push(img.select('resistance_' + y));
  resilList.push(img.select('resilience_' + y).updateMask(img.select('eligible_' + y)));
}

var meanResistance = ee.ImageCollection(resistList).mean().rename('resistance');
var meanResilience = ee.ImageCollection(resilList).mean().rename('resilience');

var classOf = function(metric) {
  var low = metric.lt(-20);
  var high = metric.gt(20);
  var medium = metric.gte(-20).and(metric.lte(20));
  return low.multiply(1).add(medium.multiply(2)).add(high.multiply(3))
    .toInt16().updateMask(metric.mask());
};

var finalOutput = meanResistance.addBands(meanResilience)
  .addBands(classOf(meanResistance).rename('resistance_class'))
  .addBands(classOf(meanResilience).rename('resilience_class'))
  .clip(aoi);

Export.image.toAsset({
  image: finalOutput,
  description: OUTPUT_DESC,
  assetId: OUTPUT_ASSET_ID,
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});
