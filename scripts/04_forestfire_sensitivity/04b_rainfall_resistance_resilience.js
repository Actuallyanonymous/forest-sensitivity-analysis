/**
 * Forest Sensitivity Analysis Pipeline — Fire Resistance & Resilience
 * Fire Shock Resistance & Resilience (Harmonized kNDVI + Signed Formulas)
 *
 * Signed resistance (both +ve and -ve events):
 * Resistance = Yn_bar / |Ye - Yn_bar| × sign(Ye - Yn_bar)
 *
 * Resilience computed ONLY when Ye < Yn_bar (negative effect years):
 * Resilience = |Ye - Yn_bar| / |Ye+1 - Yn_bar| × sign(Ye+1 - Yn_bar)
 *
 * Harmonization: Transforms Landsat 8/9 (OLI) to Landsat 5/7 (ETM+)
 * equivalent before computing kNDVI to eliminate sensor-shift bias.
 *
 * Requires:
 * - Forest mask asset (Script 1)
 * - Fire index asset (FRP > 30)
 */

// Configuration :=

var TREE_COVER_ASSET  = 'projects/cs5-pushkinmangla/assets/MP_Hybrid_Tree_Period_2003_2022';
var FIRE_INDEX_ASSET  = 'projects/cs5-pushkinmangla/assets/MP_Fire_Index_FRP30';

var OUTPUT_ASSET_ID   = 'projects/cs5-pushkinmangla/assets/MP_Fire_Metrics_Harmonized_kNDVI';
var OUTPUT_DESC       = 'MP_Fire_Metrics_Harmonized_kNDVI';

var STATE_NAME        = 'Madhya Pradesh';

// These are the years we actually want fire-resistance/resilience RESULTS
// for. If ever want to extend my analysis further (say through 2028),
// I just change END_YEAR here , nothing else needs touching.
var START_YEAR        = 2004;
var END_YEAR          = 2024;
var Z_THRESHOLD        = 1.0;

// Same fix I applied in the rain script (3b) and the drought script
// (Script 2). Yn_bar is my baseline — "what kNDVI normally looks like on
// a healthy, non-anomalous year." Before, I was computing Yn_bar over
// analysisYears directly, which right now HAPPENS to equal 2004-2024
// anyway, so it looked fine — but that's only a coincidence of my current
// settings. The moment I push END_YEAR further out in the future, Yn_bar
// would silently shift too, and that would quietly rewrite every past
// year's resistance/resilience numbers, including ones I've already
// published. So I'm freezing the baseline window here, completely
// separate from my analysis window, so that can never happen again.
//
// I'm using 2004-2024 to match the same baseline window I already locked
// in for drought (Script 2), rain (Script 3b), and the fire z-score
// itself (fire index script) — keeping all my baselines consistent with
// each other.
var BASELINE_START_YEAR = 2004;
var BASELINE_END_YEAR   = 2024;

// AOI :=

var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
            .geometry();

Map.centerObject(aoi, 7);

// Loading the assets :=

var treeMeta      = ee.Image(TREE_COVER_ASSET);
var startYearTree = treeMeta.select('start_year');
var endYearTree   = treeMeta.select('end_year');

var fireIndex = ee.Image(FIRE_INDEX_ASSET);

// I need zScore bands covering BOTH my analysis window and my baseline
// window — whichever stretches further in either direction. Right now
// they're the same range (2004-2024), so this doesn't change anything
// today, but it protects me for later when I extend END_YEAR and the two
// windows stop lining up.
var zMinYear = Math.min(START_YEAR, BASELINE_START_YEAR);
var zMaxYear = Math.max(END_YEAR, BASELINE_END_YEAR);

var zScoreCol_list = [];
for (var y = zMinYear; y <= zMaxYear; y++) {
  zScoreCol_list.push(
    fireIndex.select('zScore_' + y).rename('zScore').set('year', y)
  );
}
var zScoreCol = ee.ImageCollection(zScoreCol_list);

// LANDSAT HARMONIZATION & kNDVI :=

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

// Same idea as zScoreCol above — I need kNDVI images covering whichever
// is bigger: my baseline window, or my analysis window PLUS ONE year
// (because resilience for my very last analysis year needs next-year
// kNDVI to compare against).
var kndviMinYear = Math.min(START_YEAR, BASELINE_START_YEAR);
var kndviMaxYear = Math.max(END_YEAR + 1, BASELINE_END_YEAR);

var kndviYears = ee.List.sequence(kndviMinYear, kndviMaxYear);
var kndviCol   = ee.ImageCollection(kndviYears.map(getAnnualKNDVI));

// BASELINE kNDVI (Yn_bar) :=
// Mean kNDVI across non-anomalous years only

// This is my actual analysis window — the years I want results FOR.
// Untouched, exactly as before.
var analysisYears = ee.List.sequence(START_YEAR, END_YEAR);

// This is my frozen baseline window — the years I use to WORK OUT what
// Yn_bar (my "normal" kNDVI reference) is. On purpose, kept separate from
// analysisYears so extending my results later never shifts this.
var baselineYears = ee.List.sequence(BASELINE_START_YEAR, BASELINE_END_YEAR);

var Yn_bar = ee.ImageCollection(baselineYears.map(function(y) {
  var year   = ee.Number(y);
  var kndvi  = ee.Image(kndviCol.filter(ee.Filter.eq('year', year)).first());
  var zScore = ee.Image(zScoreCol.filter(ee.Filter.eq('year', year)).first())
                 .resample('bilinear')
                 .reproject({crs: kndvi.projection(), scale: 30});

  var isNormal = zScore.select('zScore').abs().lt(Z_THRESHOLD);
  var isForest = startYearTree.lte(year).and(endYearTree.gte(year));

  return kndvi.updateMask(isNormal.and(isForest)).set('year', year);
})).mean().rename('kndvi_baseline');

// SIGNED RESISTANCE & RESILIENCE :=

// This part is untouched — it still only runs over analysisYears (my
// actual START_YEAR..END_YEAR), it just now uses the frozen Yn_bar from
// above instead of one that would silently move whenever I change
// START_YEAR/END_YEAR.
var metricsCol = ee.ImageCollection(analysisYears.map(function(y) {
  var year   = ee.Number(y);

  var kndviYe = ee.Image(kndviCol.filter(ee.Filter.eq('year', year)).first());
  var zScore  = ee.Image(zScoreCol.filter(ee.Filter.eq('year', year)).first())
                  .resample('bilinear')
                  .reproject({crs: kndviYe.projection(), scale: 30});

  // Only compute on forest pixels during anomalous fire years
  var isAnomalous = zScore.select('zScore').gt(Z_THRESHOLD);
  var isForest    = startYearTree.lte(year).and(endYearTree.gte(year));
  var eventMask   = isAnomalous.and(isForest);

  var diffRaw = kndviYe.subtract(Yn_bar);
  var diffAbs = diffRaw.abs().max(1e-6);

  // Resistance: signed, computed for ALL anomalous years (both +ve and -ve)
  var resistance = Yn_bar.divide(diffAbs)
                         .multiply(diffRaw.signum())
                         .rename('resistance')
                         .updateMask(eventMask);

  // Resilience: ONLY computed when Ye < Yn_bar (negative effect years)
  var isNegativeEffect = kndviYe.lt(Yn_bar);
  var resilMask        = eventMask.and(isNegativeEffect);

  var kndviNext   = ee.Image(kndviCol.filter(ee.Filter.eq('year', year.add(1))).first());
  var diffNext    = kndviNext.subtract(Yn_bar);
  var diffNextAbs = diffNext.abs().max(1e-6);

  var resilience = diffAbs.divide(diffNextAbs)
                          .multiply(diffNext.signum())
                          .rename('resilience')
                          .updateMask(resilMask);

  return ee.Image.cat([resistance, resilience]).set('year', year);
}));

// AGGREGATE & EXPORT :=

var meanResist = metricsCol.select('resistance').mean().clip(aoi);
var meanResil  = metricsCol.select('resilience').mean().clip(aoi);

var finalOutput = meanResist.rename('resistance')
                            .addBands(meanResil.rename('resilience'));

// Preview
var visResist = {min: -3, max: 3,
  palette: ['8b0000','ff0000','ffffff','00ff00','006400']};
var visResil  = {min: -3, max: 3,
  palette: ['006400','ffffff','8b0000','ffffff','004d00']};

Map.addLayer(meanResist, visResist, 'Fire resistance (Harmonized kNDVI)');
Map.addLayer(meanResil,  visResil,  'Fire resilience (Harmonized kNDVI)');

Export.image.toAsset({
  image       : finalOutput,
  description : OUTPUT_DESC,
  assetId     : OUTPUT_ASSET_ID,
  region      : aoi,
  scale       : 30,
  maxPixels   : 1e13
});
