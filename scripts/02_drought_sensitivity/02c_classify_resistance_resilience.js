/**
 * Drought — Script 2c
 * Mean the yearly class numbers from Script 2, then round to a class.
 *
 * Resistance uses every scored drought year. Classes 1, 2, 3, 4.
 * Resilience leaves out class 0. If every year is class 0, the pixel is 0.
 * A mean ending in .5 rounds up.
 * Script 2 does not score a year unless the next year's kNDVI exists.
 */

var YEARLY_ASSET = 'projects/cs5-pushkinmangla/assets/MP_Drought_Metrics_SPEI3_JAS';
var OUTPUT_ASSET_ID = 'projects/cs5-pushkinmangla/assets/MP_Drought_Metrics_ClassMean';
var OUTPUT_DESC = 'MP_Drought_Metrics_ClassMean';
var STATE_NAME = 'Madhya Pradesh';
var START_YEAR = 2004;
var END_YEAR = 2022;

var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
  .geometry();

var img = ee.Image(YEARLY_ASSET);

var roundHalfUp = function(image) {
  return image.add(0.5).floor().toInt16();
};

var resistList = [];
var resilList = [];
for (var y = START_YEAR; y <= END_YEAR; y++) {
  resistList.push(img.select('resistance_class_' + y).rename('class'));
  var resilienceYear = img.select('resilience_class_' + y).rename('class');
  resilList.push(resilienceYear.updateMask(resilienceYear.gt(0)));
}

var resistanceClass = roundHalfUp(ee.ImageCollection(resistList).mean())
  .rename('resistance_class');
var resilienceClass = roundHalfUp(ee.ImageCollection(resilList).mean())
  .unmask(0)
  .updateMask(resistanceClass.mask())
  .rename('resilience_class');

var finalOutput = resistanceClass.addBands(resilienceClass).clip(aoi);

Map.centerObject(aoi, 7);
Map.addLayer(finalOutput.select('resistance_class'),
  {min: 1, max: 4, palette: ['8B0000', 'A6DBA0', '1B7837', '2166AC']},
  'Resistance class 1 deep drop, 2 within 5% below, 3 within 5% above, 4 at least 5% above');
Map.addLayer(finalOutput.select('resilience_class'),
  {min: 0, max: 5, palette: ['BDBDBD', '8B0000', 'FFFFFF', 'A6DBA0', '1B7837', '2166AC']},
  'Resilience class 0 not computed, 1 no recovery, 5 at least 5% above',
  false);

Export.image.toAsset({
  image: finalOutput,
  description: OUTPUT_DESC,
  assetId: OUTPUT_ASSET_ID,
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});
