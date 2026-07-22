/**
 * Forest Sensitivity Analysis Pipeline — Script 5b
 * High Windspeed Resistance & Resilience (Harmonized kNDVI + Signed Formulas)
 *
 * For each forest pixel, computes mean resistance and resilience
 * across all high-windspeed years (WSmax > threshold).
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
 * - Windspeed index asset from Script 5a
 */

// CONFIGURATION :=

var TREE_COVER_ASSET  = 'projects/valiant-complex-468922-j7/assets/AP_Hybrid_Tree_Period_2003_2022';
var WIND_INDEX_ASSET  = 'projects/sura-496709/assets/AP_Wind_Index';

var OUTPUT_ASSET_ID   = 'projects/sura-496709/assets/AP_Wind_Metrics_Harmonized_kNDVI';
var OUTPUT_DESC       = 'AP_Wind_Metrics_Harmonized_kNDVI';

var STATE_NAME        = 'Andhra Pradesh';

// These are the years I actually want windspeed-resistance results FOR.
// If I want to extend my analysis later, I just change END_YEAR — nothing
// else here needs touching.
var START_YEAR        = 2004;
var END_YEAR          = 2024;
var WIND_THRESHOLD    = 15;

// Same fix as everywhere else in the pipeline (drought, rain, fire). Yn_bar
// is my baseline — "what kNDVI normally looks like in a non-high-wind
// year." If I compute it over analysisYears directly, then every time I
// extend END_YEAR in the future, Yn_bar shifts a little, and that quietly
// rewrites all my past resistance/resilience results too. So I'm freezing
// this window separately.
//
// IMPORTANT: I'm capping this at 2022, not 2024 like my other baselines
// (drought/rain/fire), because Script 5a — the windspeed index this script
// depends on — has only been exported through 2022 so far. If I ever
// extend Script 5a's END_YEAR to 2024 and re-export AP_Wind_Index, I can
// bump this to 2024 to match, matching how I did it for rain (3a -> 3b).
// Until then, trying to read WSmax_2023 or WSmax_2024 here would fail,
// since those bands don't exist yet in AP_Wind_Index.
var BASELINE_START_YEAR = 2004;
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

// Load windspeed index from single multiband asset (Script 5a output)
var windIndex_raw = ee.Image(WIND_INDEX_ASSET);

// I need WSmax bands covering BOTH my analysis window and my baseline
// window — whichever stretches further. Right now they're the same range
// (2004-2022), so this doesn't change anything today, but it protects me
// once the two windows stop lining up in the future.
var wsMinYear = Math.min(START_YEAR, BASELINE_START_YEAR);
var wsMaxYear = Math.max(END_YEAR, BASELINE_END_YEAR);

var wsImages = [];
for (var y = wsMinYear; y <= wsMaxYear; y++) {
  wsImages.push(
    windIndex_raw.select('WSmax_' + y)
      .rename('windspeed')
      .set('year', y)
  );
}
var wsCol = ee.ImageCollection(wsImages);

// LANDSAT HARMONIZATION & kNDVI :=

// Chastain et al. coefficients (OLI to ETM+)
var chastainBandNames = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
var oliETMSlopes      = ee.Image.constant([1.03501, 1.00921, 1.01991, 1.14061, 1.04351, 1.05271]);
var oliETMIntercepts  = ee.Image.constant([-0.0055, -0.0008, -0.0021, -0.0163, -0.0045, 0.00261]);

// Pre-process Landsat 5/7 (Baseline)
var prepL57 = function(image) {
  var qa   = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));

  var scaled = image.updateMask(mask)
                    .select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'])
                    .multiply(0.0000275).add(-0.2);

  return scaled.rename(chastainBandNames).copyProperties(image, ["system:time_start"]);
};

// Pre-process Landsat 8/9 and Harmonize to ETM+
var prepL89 = function(image) {
  var qa   = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));

  var scaled = image.updateMask(mask)
                    .select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'])
                    .multiply(0.0000275).add(-0.2)
                    .rename(chastainBandNames);

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

// I need kNDVI images covering whichever is bigger: my baseline window, or
// my analysis window PLUS ONE year (since resilience for my last analysis
// year needs next-year kNDVI to compare against).
var kndviMinYear = Math.min(START_YEAR, BASELINE_START_YEAR);
var kndviMaxYear = Math.max(END_YEAR + 1, BASELINE_END_YEAR);

var kndviYears = ee.List.sequence(kndviMinYear, kndviMaxYear);
var kndviCol   = ee.ImageCollection(kndviYears.map(getAnnualKNDVI));

// BASELINE kNDVI (Yn_bar) :=
// Mean kNDVI across non-high-wind years only

// This is my actual analysis window — the years I want results FOR.
// Untouched, exactly as before.
var analysisYears = ee.List.sequence(START_YEAR, END_YEAR);

// This is my frozen baseline window — the years I use to WORK OUT what
// Yn_bar (my "normal" kNDVI reference) is, kept separate from
// analysisYears on purpose.
var baselineYears = ee.List.sequence(BASELINE_START_YEAR, BASELINE_END_YEAR);

var kndviNonEvent = ee.ImageCollection(baselineYears.map(function(y) {
  var year  = ee.Number(y);
  var kndvi = kndviCol.filter(ee.Filter.eq('year', year)).first();
  var ws    = wsCol.filter(ee.Filter.eq('year', year)).first()
                 .resample('bilinear')
                 .reproject({crs: kndvi.projection(), scale: 30});
  var isNonEvent = ws.lte(WIND_THRESHOLD);
  return kndvi.updateMask(isNonEvent).set('year', year);
}));

var Yn_bar = kndviNonEvent.mean().rename('kndvi_baseline');

// SIGNED RESISTANCE & RESILIENCE :=

// Unchanged — still only runs over analysisYears, just uses the frozen
// Yn_bar from above instead of one that would silently drift.
var metricsCol = ee.ImageCollection(analysisYears.map(function(y) {
  var year = ee.Number(y);

  var kndviYe = kndviCol.filter(ee.Filter.eq('year', year)).first();
  var wsYe    = wsCol.filter(ee.Filter.eq('year', year)).first()
                  .resample('bilinear')
                  .reproject({crs: kndviYe.projection(), scale: 30});

  // Only compute on forest pixels during high-windspeed years
  // Flag ANY year where WSmax crosses the threshold, however briefly
  var isForest    = startYear.lte(year).and(endYear.gte(year));
  var isHighWind  = wsYe.gt(WIND_THRESHOLD);
  var eventMask   = isForest.and(isHighWind);

  var diffRaw = kndviYe.subtract(Yn_bar);
  var diffAbs = diffRaw.abs().max(1e-6);

  // Resistance: signed, computed for ALL high-wind years
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
Map.addLayer(meanResistance, visResist, 'Windspeed resistance (Harmonized)');
Map.addLayer(meanResilience, visResil,  'Windspeed resilience (Harmonized)');

Export.image.toAsset({
  image       : finalOutput,
  description : OUTPUT_DESC,
  assetId     : OUTPUT_ASSET_ID,
  region      : aoi,
  scale       : 30,
  maxPixels   : 1e13
});
