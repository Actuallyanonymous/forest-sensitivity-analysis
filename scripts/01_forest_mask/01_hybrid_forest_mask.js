/**
 * Forest Sensitivity Analysis Pipeline — Script 1
 * Hybrid 30m Agricultural-Year Tree Cover Mask + Contiguous Forest Period
 *
 * AGRICULTURAL YEAR CONVENTION:
 * Year `y` = Jul 1 of `y` -> Jun 30 of `y+1`, matching Scripts 2-5.
 *
 * Produces per-pixel: length (ag-years), start_year, end_year (ag-year
 * labels) of the most recent unbroken forest period (2004–2024).
 *
 * Sources, and how each maps onto an ag-year `y`:
 *   1. GLC-FCS30D (2003–2022, one classification per calendar year).
 *      An ag-year straddles calendar years `y` and `y+1`, so the ag-year
 *      signal is the average of those two bands (>= 0.5 says forest).
 *      When only one of the two bands exists, that band is used alone.
 *   2. Dynamic World (2015–present) — windowed to Jul 1(y) -> Jun 30(y+1).
 *   3. IndiaSat LULC (2017–present), Core-stack — assets are named
 *      pan_india_lulc_v3_{y}_{y+1}, already ag-year aligned.
 *
 * Union logic: majority vote among active datasets per ag-year.
 *   GLC-FCS30D: active while y <= 2022 (classes 51–92)
 *   Dynamic World: active from y >= 2015 (class 1 = Trees)
 *   IndiaSat LULC: active from y >= 2017 (class 6 = Trees)
 *
 * Temporal correction: ±2 ag-year window as used in other places too by the team.
 **/


//CONFIGURATION :=


// --- Change these two lines to run for any state or district ---
var ADMIN_LEVEL  = 'state';               // 'state' or 'district'
var STATE_NAME   = 'Orissa';
var DISTRICT_NAME = '';                   // only used if ADMIN_LEVEL = 'district'

var START_YEAR      = 2004;
var END_YEAR        = 2024;
var LULC_START_YEAR = 2017;
var TEMPORAL_WINDOW = 2;

var OUTPUT_ASSET_ID = 'projects/sura-496709/assets/Orissa_Hybrid_Tree_Period_2004_2024';
var OUTPUT_DESC     = 'Orissa_Hybrid_Tree_Period_2004_2024';


// AOI :=


var aoi;
if (ADMIN_LEVEL === 'district') {
  aoi = ee.FeatureCollection('FAO/GAUL/2015/level2')
          .filter(ee.Filter.and(
            ee.Filter.eq('ADM1_NAME', STATE_NAME),
            ee.Filter.eq('ADM2_NAME', DISTRICT_NAME)
          )).geometry();
} else {
  aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
          .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
          .geometry();
}

Map.centerObject(aoi, 7);
Map.addLayer(aoi, {color: 'black', width: 2}, STATE_NAME + ' AOI', true, 0.5);


//DATASET PREPARATION :=


// --- GLC-FCS30D ---
var glcMosaic = ee.ImageCollection('projects/sat-io/open-datasets/GLC-FCS30D/annual')
  .mosaic();

// --- Dynamic World ---
var dwCol = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1').filterBounds(aoi);

// --- IndiaSat LULC (already ag-year aligned) ---
var indiaSatList = [];
for (var ly = LULC_START_YEAR; ly <= END_YEAR; ly++) {
  indiaSatList.push(
    ee.Image('projects/corestack-datasets/assets/datasets/LULC_v3_river_basin/pan_india_lulc_v3_' + ly + '_' + (ly + 1))
      .set('year', ly)
  );
}
var indiaSatCol = ee.ImageCollection(indiaSatList);

var getIndiaSatMask = function(year) {
  var img = indiaSatCol.filter(ee.Filter.eq('year', year)).first();
  return ee.Image(ee.Algorithms.If(
    img,
    ee.Image(img).select('predicted_label')
      .eq(6)
      .unmask(0),
    ee.Image(0)
  )).rename('tree');
};

var glcBandForCalendarYear = function(calYear) {
  var bandName = ee.String('b').cat(calYear.subtract(1999).format('%.0f'));
  return ee.Image(ee.Algorithms.If(
    calYear.lte(2022),
    glcMosaic.select(bandName).gte(51).and(
      glcMosaic.select(bandName).lte(92)),
    ee.Image(0)
  )).rename('tree');
};


// HYBRID MASK GENERATION :=


var years = ee.List.sequence(START_YEAR, END_YEAR);

var annualTreeCoverMasks = ee.ImageCollection(years.map(function(year) {
  year = ee.Number(year);

  // GLC — combine the two calendar-year bands this ag-year straddles.
  var glcY  = glcBandForCalendarYear(year);
  var glcY1 = glcBandForCalendarYear(year.add(1));
  var glcAvailCount = ee.Number(ee.Algorithms.If(year.lte(2022), 1, 0))
    .add(ee.Algorithms.If(year.add(1).lte(2022), 1, 0));
  var glcMask = ee.Image(ee.Algorithms.If(
    glcAvailCount.gt(0),
    glcY.add(glcY1).divide(glcAvailCount.max(1)).gte(0.5),
    ee.Image(0)
  )).rename('tree');
  var proj = glcMask.projection();

  // Dynamic World — windowed to the ag-year (Jul 1(year) -> Jun 30(year+1))
  var dwStart = ee.Date.fromYMD(year, 7, 1);
  var dwEnd   = ee.Date.fromYMD(year.add(1), 6, 30).advance(1, 'day');
  var dwYear  = dwCol.filterDate(dwStart, dwEnd).select('label');
  var dwMask  = ee.Image(ee.Algorithms.If(
    dwYear.size().gt(0),
    dwYear.mode().eq(1).unmask(0),
    ee.Image(0)
  )).rename('tree').reproject({crs: proj, scale: 30});

  // IndiaSat — already ag-year aligned
  var indiaSatMask = getIndiaSatMask(year).reproject({crs: proj, scale: 30});

  // Majority vote
  var glcActive      = ee.Algorithms.If(year.lte(2022), 1, 0);
  var forestSum      = glcMask.unmask(0).add(dwMask.unmask(0)).add(indiaSatMask.unmask(0));
  var dwActive       = ee.Algorithms.If(year.gte(2015), 1, 0);
  var indiasatActive = ee.Algorithms.If(year.gte(LULC_START_YEAR), 1, 0);
  var activeDatasets = ee.Number(glcActive).add(dwActive).add(indiasatActive);
  var requiredVotes  = ee.Algorithms.If(activeDatasets.eq(1), 1, 2);
  var hybridMask     = forestSum.gte(ee.Number(requiredVotes));

  return hybridMask.set('year', year).rename('tree');
}));


// TEMPORAL CORRECTION :=


var correctedTreeCoverMasks = ee.ImageCollection(years.map(function(year) {
  year = ee.Number(year);
  var originalMask = annualTreeCoverMasks.filter(ee.Filter.eq('year', year)).first();

  var windowMasks = annualTreeCoverMasks.filter(ee.Filter.and(
    ee.Filter.neq('year', year),
    ee.Filter.gte('year', year.subtract(TEMPORAL_WINDOW)),
    ee.Filter.lte('year', year.add(TEMPORAL_WINDOW))
  ));

  var corrected = originalMask.unmask(0).where(windowMasks.max().eq(1), 1);
  return corrected.set('year', year).rename('tree');
}));


// CONTIGUOUS FOREST PERIOD (LENGTH, START, END) :=


var calculateConsecutive = function(currentImage, previousState) {
  var prevCount  = ee.Image(ee.List(previousState).get(0));
  var prevStop   = ee.Image(ee.List(previousState).get(1));
  var currentMask = currentImage.select('tree');
  var stillCounting = prevStop.not();
  var newCount   = prevCount.add(currentMask.multiply(stillCounting));
  var newStop    = prevStop.or(currentMask.not());
  return ee.List([newCount, newStop]);
};

// Iterate backwards so we get the MOST RECENT contiguous period
var reversedCollection = correctedTreeCoverMasks.sort('year', false);
var initialState = ee.List([ee.Image(0).byte(), ee.Image(0).byte()]);
var finalState   = ee.List(reversedCollection.iterate(calculateConsecutive, initialState));

var recentLength    = ee.Image(finalState.get(0)).rename('length');
var forestEndYear   = ee.Image(END_YEAR).multiply(recentLength.gt(0)).rename('end_year');
var forestStartYear = forestEndYear.subtract(recentLength).add(1)
                        .multiply(recentLength.gt(0)).rename('start_year');

var finalOutput = recentLength
  .addBands(forestStartYear)
  .addBands(forestEndYear)
  .clip(aoi);


// VISUALISE & EXPORT :=


var RUN_LENGTH_VIS = {
  min: 0, max: END_YEAR - START_YEAR + 1,
  palette: ['#fef0d9','#fdcc8a','#fc8d59','#e34a33','#b30000','#4d0000']
};

Map.addLayer(finalOutput.select('length'),     RUN_LENGTH_VIS,       'Forest period length');
Map.addLayer(finalOutput.select('start_year'), {min: START_YEAR, max: END_YEAR,
  palette:['#ffffd4','#fed98e','#fe9929','#d95f0e','#993404']},      'Forest start year');
Map.addLayer(finalOutput.select('end_year'),   {min: START_YEAR, max: END_YEAR,
  palette:['#f7fbff','#c6dbef','#6baed6','#2171b5','#08306b']},      'Forest end year');

Export.image.toAsset({
  image       : finalOutput,
  description : OUTPUT_DESC,
  assetId     : OUTPUT_ASSET_ID,
  region      : aoi,
  scale       : 30,
  maxPixels   : 1e13
});
