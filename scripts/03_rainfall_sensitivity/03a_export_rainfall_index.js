/**
 * Forest Sensitivity Analysis Pipeline — Script 3a 
 * Heavy Rainfall Index Export (Wet-Days Only Threshold + Extended Metrics)
 *
 * Computes 5 quantities per pixel per year and exports as a single
 * multiband asset — one band per year for each quantity:
 *
 * Hm_{year}         = annual sum of precipitation on heavy days
 * zScore_{year}     = z-score of Hm relative to the 2004-2022 period
 * heavyDays_{year}  = number of days in the year with heavy rainfall
 * heavyAvg_{year}   = average daily precipitation on heavy rainfall days
 * maxDay_{year}     = maximum daily precipitation recorded in that year
 *
 * Heavy day definition: daily precipitation > long-term 95th percentile 
 * of WET DAYS ONLY (> 1mm) (CHIRPS)
 *
 * This asset is the direct input to Script 3b, analogous to how
 * SPEI assets are the input to Script 2 (drought).
 *
 * Requires: nothing — only public datasets (CHIRPS)
 */



// 1. CONFIGURATION

var STATE_NAME      = 'Madhya Pradesh';

//change these START_YEAR & END_YEAR for changing the analysis timeline . don't touch the baseline years, if not sure how this works.
var START_YEAR      = 2004;
var END_YEAR        = 2024;


// Fixed baseline window for zScore normalization — independent of START_YEAR/END_YEAR.
// Do NOT change this once results are published, or old zScore bands will drift
// when the pipeline timeline is extended.
var BASELINE_START_YEAR = 2004;
var BASELINE_END_YEAR   = 2024;

var OUTPUT_ASSET_ID = 'projects/cs5-pushkinmangla/assets/MP_Rain_Index_Wet95_Final';
var OUTPUT_DESC     = 'MP_Rain_Index_Wet95_Final';



//  2. AOI

var aoi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', STATE_NAME))
            .geometry();

Map.centerObject(aoi, 7);



// 3. BASELINE HEAVY RAINFALL THRESHOLD

//extend this too if the pipeline is extended further. 
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
               .filterBounds(aoi)
               .filterDate('2000-01-01', ee.Date.fromYMD(BASELINE_END_YEAR, 12, 31))
               .select('precipitation');
//now here I have written BASELINE_END_YEAR so this won't be needed to change even if baseline year changes and peipeline extends, they'll go hand in hand . Yayy
var proj = chirps.first().projection();

// Long-term 95th percentile of WET DAYS ONLY (> 1mm)
var p95 = chirps.map(function(img) {
                  return img.updateMask(img.gt(1)); 
                })
                .reduce(ee.Reducer.percentile([95]))
                .setDefaultProjection(proj)
                .rename('p95');




// 4. ANNUAL METRICS CALCULATION

var metricsMinYear = Math.min(START_YEAR, BASELINE_START_YEAR);
var metricsMaxYear = Math.max(END_YEAR, BASELINE_END_YEAR);
var metricsYears   = ee.List.sequence(metricsMinYear, metricsMaxYear);

var annualMetrics = ee.ImageCollection(metricsYears.map(function(y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var end   = ee.Date.fromYMD(y, 12, 31);

  var yearCollection = chirps.filterDate(start, end);

  // Absolute maximum daily rainfall for this year
  var maxDay = yearCollection.max()
                             .unmask(0)
                             .setDefaultProjection(proj)
                             .rename('maxDay');

  // Binary mask: isolated heavy rain days
  var heavyRainCollection = yearCollection.map(function(img) {
    return img.updateMask(img.gt(p95));
  });

  // 1. Annual sum of heavy rainfall
  var hm = heavyRainCollection.map(function(img) { return img.unmask(0); })
                              .sum()
                              .setDefaultProjection(proj)
                              .rename('Hm');

  // 2. Number of heavy days
  var heavyDays = heavyRainCollection.map(function(img) { return img.mask(); })
                                     .sum()
                                     .unmask(0)
                                     .setDefaultProjection(proj)
                                     .rename('heavyDays');

  // 3. Average intensity of heavy days
  var heavyAvg = heavyRainCollection.mean()
                                    .unmask(0)
                                    .setDefaultProjection(proj)
                                    .rename('heavyAvg');

  // Combine metrics into a single image per year containing all 4 base properties
  return hm.addBands(heavyDays)
           .addBands(heavyAvg)
           .addBands(maxDay)
           .set('year', y);
}));





// 5. FIXED Z-SCORE CALCULATION

var baselineYears = ee.List.sequence(BASELINE_START_YEAR, BASELINE_END_YEAR);

var baselineMetrics = annualMetrics.filter(
  ee.Filter.inList('year', baselineYears)
);

var hmMean   = baselineMetrics.select('Hm').mean();
var hmStdDev = baselineMetrics.select('Hm').reduce(ee.Reducer.stdDev());

// Changed to server-side map architecture instead of the annual image collection, done for GEE optimisation.. as it's a better practice. 
var completedAnnualCollection = annualMetrics.map(function(img) {
  var hm = img.select('Hm');
  var z  = hm.subtract(hmMean).divide(hmStdDev).rename('zScore');
  return img.addBands(z); 
});




// 6. SERVER-SIDE STACK INTO SINGLE MULTIBAND IMAGE & EXPORT

var analysisYears = ee.List.sequence(START_YEAR, END_YEAR); // same as original `years`

// Changed server-side iteration style to stack and correctly rename bands
var initialImage = ee.Image([]);
var outputImage = ee.Image(analysisYears.iterate(function(y, acc) {
  var yearStr = ee.String(ee.Number(y).toInt());
  var yearImg = completedAnnualCollection.filter(ee.Filter.eq('year', y)).first();
  
  var hmBand        = yearImg.select('Hm').rename(ee.String('Hm_').cat(yearStr));
  var zBand         = yearImg.select('zScore').rename(ee.String('zScore_').cat(yearStr));
  var heavyDaysBand = yearImg.select('heavyDays').rename(ee.String('heavyDays_').cat(yearStr));
  var heavyAvgBand  = yearImg.select('heavyAvg').rename(ee.String('heavyAvg_').cat(yearStr));
  var maxDayBand    = yearImg.select('maxDay').rename(ee.String('maxDay_').cat(yearStr));
  
  return ee.Image(acc).addBands([hmBand, zBand, heavyDaysBand, heavyAvgBand, maxDayBand]);
}, initialImage));

// Export execution block
Export.image.toAsset({
  image       : outputImage.clip(aoi),
  description : OUTPUT_DESC,
  assetId     : OUTPUT_ASSET_ID,
  region      : aoi,
  scale       : 5566,   
  crs         : 'EPSG:4326',
  maxPixels   : 1e13
});

print('Clean pipeline compilation verified.');
print('Ready to execute. Total structured bands: ' + analysisYears.length().multiply(5).getInfo());
