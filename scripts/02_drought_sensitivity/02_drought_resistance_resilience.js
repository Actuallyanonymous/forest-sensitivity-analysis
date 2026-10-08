/**
 * Forest Sensitivity Analysis Pipeline — Script 2
 * Drought Resistance & Resilience (Harmonized kNDVI + Signed Formulas)
 *
 * For each forest pixel and drought year, stores signed resistance and
 * resilience. Does not average across years. Script 02c aggregates and
 * classifies those yearly bands.
 *
 * Harmonization: Transforms Landsat 8/9 (OLI) to Landsat 5/7 (ETM+)
 * equivalent before computing kNDVI to eliminate sensor-shift bias.
 *
 * Signed resistance (both +ve and -ve events):
 * Resistance = Yn_bar / |Ye - Yn_bar| × sign(Ye - Yn_bar)
 *
 * Resilience = |Ye - Yn_bar| / |Ye+1 - Yn_bar| × sign(Ye+1 - Yn_bar)
 * Stored for every drought year. The Ye < 0.95*Yn filter is applied in 02c.
 * Band eligible_YYYY = 1 when the pixel is in drought and Ye < 0.95*Yn.
 *
 * Drought year: SPEI-3 ending September (July–August–September) < -1.
 *
 * Requires:
 * - Forest mask asset from Script 1
 * - SPEI-3 asset from the SPEI pipeline (band yYYYY_m07_09)
 */

// CONFIGURATION :=

var TREE_COVER_ASSET  = 'projects/cs5-pushkinmangla/assets/MP_Hybrid_Tree_Period_2003_2022';
var SPEI3_ASSET       = 'projects/cs5-pushkinmangla/assets/SPEI/SPEI3_Madhya_Pradesh';
var OUTPUT_ASSET_ID   = 'projects/cs5-pushkinmangla/assets/MP_Drought_Metrics_SPEI3_JAS';
var OUTPUT_DESC       = 'MP_Drought_Metrics_SPEI3_JAS';
var STATE_NAME        = 'Madhya Pradesh';
var START_YEAR        = 2004;
var END_YEAR          = 2022;
var DROUGHT_THRESHOLD = -1.0;   // SPEI-3 JAS below this = drought year

// Fixed baseline window. This is independent of analysis START_YEAR/END_YEAR.
// Please don't EVER change this once results are published, or old outputs will change when the pipeline timeline is extended.
// This helps in fixing the Yn_bar (the average of the non-drought year NDVI) to a constant value.
// SPEI-3 in the current file runs through 2023 (y2004_m07_09 ... y2023_m07_09). Do not set the end to 2024 until that band exists.
var BASELINE_START_YEAR = 2004;  // SPEI has no data before 2004
var BASELINE_END_YEAR   = 2023;

// AOI :=

var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
            .geometry();

Map.centerObject(aoi, 7);

// Loading the assets :=

var treeMeta  = ee.Image(TREE_COVER_ASSET);
var startYear = treeMeta.select('start_year');
var endYear   = treeMeta.select('end_year');

var spei3_raw = ee.Image(SPEI3_ASSET);

var speiMinYear = Math.min(START_YEAR, BASELINE_START_YEAR);
var speiMaxYear = Math.max(END_YEAR, BASELINE_END_YEAR);

var speiImages = [];
for (var y = speiMinYear; y <= speiMaxYear; y++) {
  // SPEI-3 ending September = July–August–September
  speiImages.push(
    spei3_raw.select('y' + y + '_m07_09')
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

  // Resilience for every drought year. 02c keeps only Ye < 0.95*Yn.
  var kndviNext   = kndviCol.filter(ee.Filter.eq('year', year.add(1))).first();
  var diffNext    = kndviNext.subtract(Yn_bar);
  var diffNextAbs = diffNext.abs().max(1e-6);

  var resilience = diffAbs.divide(diffNextAbs)
                          .multiply(diffNext.signum())
                          .rename('resilience')
                          .updateMask(eventMask);

  var eligible = eventMask.and(kndviYe.lt(Yn_bar.multiply(0.95)))
                          .rename('eligible');

  return ee.Image.cat([resistance, resilience, eligible]).set('year', year);
}));

// One band per year. Names: resistance_YYYY, resilience_YYYY, eligible_YYYY.

var yearlyBands = [];
for (var yy = START_YEAR; yy <= END_YEAR; yy++) {
  var yearImg = ee.Image(metricsCol.filter(ee.Filter.eq('year', yy)).first());
  yearlyBands.push(yearImg.select('resistance').rename('resistance_' + yy));
  yearlyBands.push(yearImg.select('resilience').rename('resilience_' + yy));
  yearlyBands.push(yearImg.select('eligible').rename('eligible_' + yy));
}
var finalOutput = ee.Image.cat(yearlyBands).clip(aoi);

Map.addLayer(finalOutput.select('resistance_' + END_YEAR),
  {min: -3, max: 3, palette: ['8b0000','ff0000','ffffff','00ff00','006400']},
  'Drought resistance ' + END_YEAR);

Export.image.toAsset({
  image       : finalOutput,
  description : OUTPUT_DESC,
  assetId     : OUTPUT_ASSET_ID,
  region      : aoi,
  scale       : 30,
  maxPixels   : 1e13
});
