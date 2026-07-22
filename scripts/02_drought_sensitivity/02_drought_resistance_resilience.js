/**
 * Forest Sensitivity Analysis Pipeline — Script 2
 * Drought Resistance & Resilience (Harmonized kNDVI + Signed Formulas)
 *
 * For each forest pixel, computes mean resistance and resilience
 * across all drought years (SPEI-12 < threshold).
 *
 * Harmonization: Transforms Landsat 8/9 (OLI) to Landsat 5/7 (ETM+)
 * equivalent before computing kNDVI to eliminate sensor-shift bias.
 *
 * Signed resistance (both +ve and -ve events):
 * Resistance = Yn_bar / |Ye - Yn_bar| × sign(Ye - Yn_bar)
 *
 * Resilience computed ONLY when Ye < Yn_bar (negative effect years):
 * Resilience = |Ye - Yn_bar| / |Ye+1 - Yn_bar| × sign(Ye+1 - Yn_bar)
 *
 * Requires:
 * - Forest mask asset from Script 1
 * - SPEI-12 assets from spei-drought-analysis-pipeline
 */

// CONFIGURATION :=

var TREE_COVER_ASSET  = 'projects/cs5-pushkinmangla/assets/MP_Hybrid_Tree_Period_2003_2022';
var OUTPUT_ASSET_ID   = 'projects/cs5-pushkinmangla/assets/MP_Drought_Metrics_Harmonized_kNDVI';
var OUTPUT_DESC       = 'MP_Drought_Metrics_Harmonized_kNDVI';
var STATE_NAME        = 'Madhya Pradesh';
var START_YEAR        = 2004;
var END_YEAR          = 2022;
var DROUGHT_THRESHOLD = -1.0;   // SPEI-12 below this = drought year

// Fixed baseline window. This is independent of analysis START_YEAR/END_YEAR.
// Please don't EVER change this once results are published, or old outputs will change when the pipeline timeline is extended.
// This helps in fixing the Yn_bar (the average of the non-drought year NDVI) to a constant value. 
var BASELINE_START_YEAR = 2004;  // SPEI has no data before 2004
var BASELINE_END_YEAR   = 2024;

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

// Building the per-year SPEI collection here. 
var speiMinYear = Math.min(START_YEAR, BASELINE_START_YEAR);
var speiMaxYear = Math.max(END_YEAR, BASELINE_END_YEAR);

var speiImages = [];
for (var y = speiMinYear; y <= speiMaxYear; y++) {
  speiImages.push(
    spei12_named.select('y' + y)
      .rename('spei')
      .set('year', y)
  );
}
var speiCol = ee.ImageCollection(speiImages);

// LANDSAT HARMONIZATION & kNDVI :=

// Chastain et al. coefficients (OLI to ETM+)
var chastainBandNames = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
var oliETMSlopes      = ee.Image.constant([1.03501, 1.00921, 1.01991, 1.14061, 1.04351, 1.05271]);
var oliETMIntercepts  = ee.Image.constant([-0.0055, -0.0008, -0.0021, -0.0163, -0.0045, 0.00261]);

// Pre-process Landsat 5/7 (Baseline)
var prepL57 = function(image) {
  var qa   = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  
  // Apply mask, select optical bands, and apply Collection 2 scale factors
  var scaled = image.updateMask(mask)
                    .select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'])
                    .multiply(0.0000275).add(-0.2);
                    
  return scaled.rename(chastainBandNames).copyProperties(image, ["system:time_start"]);
};

// Pre-process Landsat 8/9 and Harmonize to ETM+
var prepL89 = function(image) {
  var qa   = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  
  // Apply mask, select optical bands, and apply Collection 2 scale factors
  var scaled = image.updateMask(mask)
                    .select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'])
                    .multiply(0.0000275).add(-0.2)
                    .rename(chastainBandNames);
                    
  // Apply Chastain regression model (OLI -> ETM+)
  var harmonized = scaled.multiply(oliETMSlopes).add(oliETMIntercepts);
  
  return harmonized.copyProperties(image, ["system:time_start"]);
};

// Calculate annual median kNDVI
var getAnnualKNDVI = function(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end   = ee.Date.fromYMD(year, 12, 31);

  var l89 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
              .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
              .filterDate(start, end).filterBounds(aoi)
              .map(prepL89)
              .map(function(img) {
                return img.normalizedDifference(['NIR','RED']).pow(2).tanh().rename('kndvi');
              });

  var l57 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
              .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
              .filterDate(start, end).filterBounds(aoi)
              .map(prepL57)
              .map(function(img) {
                return img.normalizedDifference(['NIR','RED']).pow(2).tanh().rename('kndvi');
              });

  return l89.merge(l57).median().set('year', year).rename('kndvi');
};

// Load kNDVI for START_YEAR to END_YEAR+1 (need next year for resilience)
//Here this is also changed based on the baseline years. 
var kndviMinYear = Math.min(START_YEAR, BASELINE_START_YEAR);
var kndviMaxYear = Math.max(END_YEAR + 1, BASELINE_END_YEAR);

var kndviYears = ee.List.sequence(kndviMinYear, kndviMaxYear);
var kndviCol   = ee.ImageCollection(kndviYears.map(getAnnualKNDVI));

// BASELINE kNDVI (Yn_bar) :=
// Mean kNDVI across non-drought years only

//the analaysis years will still remain the same , so if some internal year widtch is given like 2010-2018 for example
//then too the baseline yn_bar would be same of the bigger normalization. But the analysis will be resulting only of the analyssi years.
var analysisYears = ee.List.sequence(START_YEAR, END_YEAR);

var baselineYears = ee.List.sequence(BASELINE_START_YEAR, BASELINE_END_YEAR);

var kndviNonDrought = ee.ImageCollection(baselineYears.map(function(y) {
  var year  = ee.Number(y);
  var kndvi = kndviCol.filter(ee.Filter.eq('year', year)).first();
  var spei  = speiCol.filter(ee.Filter.eq('year', year)).first()
                 .resample('bilinear')
                 .reproject({crs: kndvi.projection(), scale: 30});
  var isNonDrought = spei.gte(DROUGHT_THRESHOLD);
  return kndvi.updateMask(isNonDrought).set('year', year);
}));

var Yn_bar = kndviNonDrought.mean().rename('kndvi_baseline');

// SIGNED RESISTANCE & RESILIENCE :=

var metricsCol = ee.ImageCollection(analysisYears.map(function(y) {
  var year = ee.Number(y);

  var kndviYe = kndviCol.filter(ee.Filter.eq('year', year)).first();
  var speiYe  = speiCol.filter(ee.Filter.eq('year', year)).first()
                  .resample('bilinear')
                  .reproject({crs: kndviYe.projection(), scale: 30});

  // Only compute on forest pixels during drought years
  var isForest  = startYear.lte(year).and(endYear.gte(year));
  var isDrought = speiYe.lt(DROUGHT_THRESHOLD);
  var eventMask = isForest.and(isDrought);

  var diffRaw = kndviYe.subtract(Yn_bar);
  var diffAbs = diffRaw.abs().max(1e-6);

  // Resistance: signed, computed for ALL drought years
  var resistance = Yn_bar.divide(diffAbs)
                         .multiply(diffRaw.signum())
                         .rename('resistance')
                         .updateMask(eventMask);

  // Resilience: ONLY computed when Ye < Yn_bar (negative effect years)
  var isNegativeEffect = kndviYe.lt(Yn_bar);
  var resilMask        = eventMask.and(isNegativeEffect);

  var kndviNext   = kndviCol.filter(ee.Filter.eq('year', year.add(1))).first();
  var diffNext    = kndviNext.subtract(Yn_bar);
  var diffNextAbs = diffNext.abs().max(1e-6);

  var resilience = diffAbs.divide(diffNextAbs)
                          .multiply(diffNext.signum())
                          .rename('resilience')
                          .updateMask(resilMask);

  return ee.Image.cat([resistance, resilience]).set('year', year);
}));

// AGGREGATE & EXPORT :=

var meanResistance = metricsCol.select('resistance').mean().clip(aoi);
var meanResilience = metricsCol.select('resilience').mean().clip(aoi);

var finalOutput = meanResistance.rename('resistance')
                                .addBands(meanResilience.rename('resilience'));

// Preview
var visResist = {min: -3, max: 3,
  palette: ['8b0000','ff0000','ffffff','00ff00','006400']};
var visResil  = {min: -3, max: 3,
  palette: ['006400','ffffff','8b0000','ffffff','004d00']};
Map.addLayer(meanResistance, visResist, 'Drought resistance (Harmonized)');
Map.addLayer(meanResilience, visResil,  'Drought resilience (Harmonized)');

Export.image.toAsset({
  image       : finalOutput,
  description : OUTPUT_DESC,
  assetId     : OUTPUT_ASSET_ID,
  region      : aoi,
  scale       : 30,
  maxPixels   : 1e13
});
