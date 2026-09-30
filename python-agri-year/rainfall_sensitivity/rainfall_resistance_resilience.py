import ee
from datetime import date

from utilities.constants import AEZ
from utilities.gee_utils import (
    ee_initialize,
    is_gee_asset_exists,
    export_raster_asset_to_gee,
)


def generate_rainfall_resilience(
    aez, start_year=2004, end_year=None, gee_account_id=None
):
    """
    * Forest Sensitivity Analysis Pipeline — Script 3b
    * Heavy Rainfall Resistance & Resilience (Harmonized kNDVI + Signed Formulas)
    *
    * AGRICULTURAL YEAR CONVENTION:
    * Year `y` = Jul 1 of `y` -> Jun 30 of `y+1`, matching Script 3a.
    * kNDVI compositing windows, zScore band lookups, and the baseline
    * (Yn_bar) window are all on this convention.
    *
    * Tree-cover / forest mask asset (Script 1) is also ag-year now, so
    * startYearTree/endYearTree comparisons below against `year` are on
    * the same convention with no fuzz.
    *
    * For each forest pixel and heavy-rain ag-year, stores signed resistance
    * and resilience. Does not average across years.
    * classify_resistance_resilience.py averages those yearly bands and
    * assigns low / medium / high classes.
    *
    * Annual kNDVI is the median of scenes in the ag-year. Yn_bar is the
    * mean of those annual values over non-anomalous forest ag-years in the
    * fixed baseline window.
    *
    * Signed resistance (both +ve and -ve events):
    * Resistance = Yn_bar / |Ye - Yn_bar| × sign(Ye - Yn_bar)
    *
    * Resilience = |Ye - Yn_bar| / |Ye+1 - Yn_bar| × sign(Ye+1 - Yn_bar)
    * Stored for every heavy-rain year. The Ye < 0.95*Yn_bar filter is
    * applied at classification. Band eligible_YYYY = 1 when the pixel is
    * in a heavy-rain year and Ye < 0.95*Yn_bar.
    *
    * Output bands: resistance_YYYY, resilience_YYYY, eligible_YYYY.
    *
    * Harmonization: Transforms Landsat 8/9 (OLI) to Landsat 5/7 (ETM+)
    * equivalent before computing kNDVI to eliminate sensor-shift bias.
    *
    * Requires:
    * - Forest mask asset (Script 1)
    * - Rainfall index asset (Script 3a) — must have been exported with an
    *   end_year >= this script's zMaxYear, or the zScore_{y} band lookup
    *   below will fail.
    """

    ee_initialize(gee_account_id)
    TREE_COVER_ASSET = f"projects/corestack-datasets-alpha/assets/datasets/SPEI/Hybrid_Tree_AEZ_{aez}_{str(2004)}_{str(end_year)}"
    RAIN_INDEX_ASSET = (
        f"projects/corestack-datasets-alpha/assets/datasets/SPEI/rain_index_AEZ_{aez}"
    )

    OUTPUT_DESC = f"Rain_Metrics_{aez}"  # f"Rain_Metrics_AEZ_{aez}"
    OUTPUT_ASSET_ID = (
        f"projects/corestack-datasets-alpha/assets/datasets/SPEI/{OUTPUT_DESC}"
    )

    # Same asset id as the previous 2-band mean export. Delete that asset
    # before re-running; this script now writes per-year bands.
    if is_gee_asset_exists(OUTPUT_ASSET_ID):
        return None

    if end_year is None:
        raise ValueError(
            "end_year must be specified explicitly — agricultural-year pipelines "
            "cannot infer a safe default."
        )

    Z_THRESHOLD = 1.0

    # Same frozen-baseline reasoning as Script 3a and Script 2 — Yn_bar (the
    # "normal" kNDVI reference) must not shift every time the analysis
    # window is extended, or already-published resistance/resilience
    # numbers would silently drift.
    #
    # BASELINE_START_YEAR=2004 / BASELINE_END_YEAR=2024 mean the same
    # ag-years here as in Script 3a: Jul 2004 -> Jun 2025 as the frozen
    # normalization period. Kept identical on purpose so 3a's zScore
    # baseline and 3b's kNDVI baseline stay aligned on the same years.
    BASELINE_START_YEAR = 2004
    BASELINE_END_YEAR = 2024

    aoi = ee.FeatureCollection(AEZ).filter(ee.Filter.eq("ae_regcode", aez)).geometry()

    # zScoreCol needs bands up through zMaxYear from the 3a asset, and
    # kndviCol needs data one year further out than that for resilience on
    # the final analysis year (see completeness guard below).
    zMinYear = min(start_year, BASELINE_START_YEAR)
    zMaxYear = max(end_year, BASELINE_END_YEAR)

    kndviMinYear = min(start_year, BASELINE_START_YEAR)
    kndviMaxYear = max(end_year + 1, BASELINE_END_YEAR)

    # kndviMaxYear is the furthest single ag-year fetched anywhere in this
    # script; that ag-year isn't complete until Jun 30 of kndviMaxYear+1.
    # Refuse to silently compute on a partial year.
    required_complete_by = date(kndviMaxYear + 1, 6, 30)
    if required_complete_by > date.today():
        raise ValueError(
            f"Agricultural year {kndviMaxYear} (Jul {kndviMaxYear} -> "
            f"Jun {kndviMaxYear + 1}) is not yet complete as of "
            f"{date.today().isoformat()}. Reduce end_year."
        )

    treeMeta = ee.Image(TREE_COVER_ASSET)
    startYearTree = treeMeta.select("start_year")
    endYearTree = treeMeta.select("end_year")

    rainIndex = ee.Image(RAIN_INDEX_ASSET)

    zScoreCol_list = []
    for y in range(zMinYear, zMaxYear + 1):
        zScoreCol_list.append(
            rainIndex.select("zScore_" + str(y)).rename("zScore").set("year", y)
        )

    zScoreCol = ee.ImageCollection(zScoreCol_list)

    # LANDSAT HARMONIZATION & kNDVI :=

    chastainBandNames = ["BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2"]
    oliETMSlopes = ee.Image.constant(
        [1.03501, 1.00921, 1.01991, 1.14061, 1.04351, 1.05271]
    )
    oliETMIntercepts = ee.Image.constant(
        [-0.0055, -0.0008, -0.0021, -0.0163, -0.0045, 0.00261]
    )

    # Pre-process Landsat 5/7 (Baseline)
    def prepL57(image):
        qa = image.select("QA_PIXEL")
        mask = qa.bitwiseAnd(1 << 3).eq(0).And(qa.bitwiseAnd(1 << 4).eq(0))

        scaled = (
            image.updateMask(mask)
            .select(["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"])
            .multiply(0.0000275)
            .add(-0.2)
        )

        return scaled.rename(chastainBandNames).copyProperties(
            image, ["system:time_start"]
        )

    # Pre-process Landsat 8/9 and Harmonize to ETM+
    def prepL89(image):
        qa = image.select("QA_PIXEL")
        mask = qa.bitwiseAnd(1 << 3).eq(0).And(qa.bitwiseAnd(1 << 4).eq(0))

        scaled = (
            image.updateMask(mask)
            .select(["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"])
            .multiply(0.0000275)
            .add(-0.2)
            .rename(chastainBandNames)
        )

        harmonized = scaled.multiply(oliETMSlopes).add(oliETMIntercepts)

        return harmonized.copyProperties(image, ["system:time_start"])

    # Calculate ag-year median kNDVI (Jul 1(year) -> Jun 30(year+1))
    def getAnnualKNDVI(year):
        start = ee.Date.fromYMD(year, 7, 1)
        # end-exclusive filterDate, so advance 1 day past Jun 30 to include it
        end = ee.Date.fromYMD(ee.Number(year).add(1), 6, 30).advance(1, "day")

        l89 = (
            ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
            .merge(ee.ImageCollection("LANDSAT/LC09/C02/T1_L2"))
            .filterDate(start, end)
            .filterBounds(aoi)
            .map(prepL89)
            .map(
                lambda img: img.normalizedDifference(["NIR", "RED"])
                .pow(2)
                .tanh()
                .rename("kndvi")
            )
        )

        l57 = (
            ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
            .merge(ee.ImageCollection("LANDSAT/LE07/C02/T1_L2"))
            .filterDate(start, end)
            .filterBounds(aoi)
            .map(prepL57)
            .map(
                lambda img: img.normalizedDifference(["NIR", "RED"])
                .pow(2)
                .tanh()
                .rename("kndvi")
            )
        )

        return l89.merge(l57).median().set("year", year).rename("kndvi")

    kndviYears = ee.List.sequence(kndviMinYear, kndviMaxYear)
    kndviCol = ee.ImageCollection(kndviYears.map(getAnnualKNDVI))

    # BASELINE kNDVI (Yn_bar) :=
    # Mean kNDVI across non-anomalous ag-years only

    # Actual analysis window — the ag-years to compute resistance/resilience FOR.
    analysisYears = ee.List.sequence(start_year, end_year)

    # Frozen baseline window — the ag-years used to work out Yn_bar. Kept
    # separate from analysisYears so extending the analysis later doesn't
    # shift Yn_bar and quietly rewrite already-published results.
    baselineYears = ee.List.sequence(BASELINE_START_YEAR, BASELINE_END_YEAR)

    def calc_Yn_bar(y):
        year = ee.Number(y)
        kndvi = ee.Image(kndviCol.filter(ee.Filter.eq("year", year)).first())
        zScore = (
            ee.Image(zScoreCol.filter(ee.Filter.eq("year", year)).first())
            .resample("bilinear")
            .reproject(crs=kndvi.projection(), scale=30)
        )

        isNormal = zScore.select("zScore").abs().lt(Z_THRESHOLD)
        isForest = startYearTree.lte(year).And(endYearTree.gte(year))

        return kndvi.updateMask(isNormal.And(isForest)).set("year", year)

    Yn_bar = (
        ee.ImageCollection(baselineYears.map(calc_Yn_bar))
        .mean()
        .rename("kndvi_baseline")
    )

    # SIGNED RESISTANCE & RESILIENCE :=
    # Runs only over analysisYears, using the frozen Yn_bar above.
    # Each event year is stored. Stage C averages and classifies.

    def calc_metrics_cols(y):
        year = ee.Number(y)

        kndviYe = ee.Image(kndviCol.filter(ee.Filter.eq("year", year)).first())
        zScore = (
            ee.Image(zScoreCol.filter(ee.Filter.eq("year", year)).first())
            .resample("bilinear")
            .reproject(crs=kndviYe.projection(), scale=30)
        )

        # Only compute on forest pixels during anomalous rainfall ag-years
        isAnomalous = zScore.select("zScore").gt(Z_THRESHOLD)
        isForest = startYearTree.lte(year).And(endYearTree.gte(year))
        eventMask = isAnomalous.And(isForest)

        diffRaw = kndviYe.subtract(Yn_bar)
        diffAbs = diffRaw.abs().max(1e-6)

        # Resistance: signed, computed for ALL anomalous years (both +ve and -ve)
        resistance = (
            Yn_bar.divide(diffAbs)
            .multiply(diffRaw.signum())
            .rename("resistance")
            .updateMask(eventMask)
        )

        # Resilience for every heavy-rain year. Classification keeps only
        # years with eligible = 1 (Ye < 0.95 * Yn_bar).
        kndviNext = ee.Image(kndviCol.filter(ee.Filter.eq("year", year.add(1))).first())
        diffNext = kndviNext.subtract(Yn_bar)
        diffNextAbs = diffNext.abs().max(1e-6)

        resilience = (
            diffAbs.divide(diffNextAbs)
            .multiply(diffNext.signum())
            .rename("resilience")
            .updateMask(eventMask)
        )

        eligible = eventMask.And(kndviYe.lt(Yn_bar.multiply(0.95))).rename(
            "eligible"
        )

        return ee.Image.cat([resistance, resilience, eligible]).set("year", year)

    metricsCol = ee.ImageCollection(analysisYears.map(calc_metrics_cols))

    # One band per ag-year. Averaging and ±20 classes are stage C.
    yearly_bands = []
    for yy in range(start_year, end_year + 1):
        year_img = ee.Image(metricsCol.filter(ee.Filter.eq("year", yy)).first())
        yearly_bands.append(year_img.select("resistance").rename(f"resistance_{yy}"))
        yearly_bands.append(year_img.select("resilience").rename(f"resilience_{yy}"))
        yearly_bands.append(year_img.select("eligible").rename(f"eligible_{yy}"))

    finalOutput = ee.Image.cat(yearly_bands).clip(aoi)

    task_id = export_raster_asset_to_gee(
        finalOutput, OUTPUT_DESC, OUTPUT_ASSET_ID, scale=30, region=aoi
    )

    return task_id
