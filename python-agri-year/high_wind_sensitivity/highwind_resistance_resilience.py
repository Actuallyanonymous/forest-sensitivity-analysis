import ee
from datetime import date

from utilities.constants import AEZ
from utilities.gee_utils import (
    is_gee_asset_exists,
    export_raster_asset_to_gee,
    ee_initialize,
)


def high_wind_sensitivity(aez, start_year=2004, end_year=None, gee_account_id=None):
    """
    * Forest Sensitivity Analysis Pipeline — Script 5b
    * High Windspeed Resistance & Resilience (Harmonized kNDVI + Signed Formulas)
    *
    * AGRICULTURAL YEAR CONVENTION:
    * Year `y` = Jul 1 of `y` -> Jun 30 of `y+1`, matching Scripts 1/3a/3b/
    * fire_index/forest_fire_sensitivity/5a. kNDVI compositing windows and
    * the baseline (Yn_bar) window are all on this convention. Tree-cover
    * asset (Script 1) is also ag-year, so startYear/endYear comparisons
    * against `year` carry no fuzz.
    *
    * For each forest pixel and high-wind year, stores signed resistance
    * and resilience. classify_wind_metrics averages those yearly bands
    * and labels low, medium, and high.
    *
    * Harmonization: Transforms Landsat 8/9 (OLI) to Landsat 5/7 (ETM+)
    * equivalent before computing kNDVI to eliminate sensor-shift bias.
    *
    * Signed resistance (both +ve and -ve events):
    * Resistance = Yn_bar / |Ye - Yn_bar| × sign(Ye - Yn_bar)
    *
    * Resilience = |Ye - Yn_bar| / |Ye+1 - Yn_bar| × sign(Ye+1 - Yn_bar)
    * Stored for every high-wind year. The Ye < 0.95*Yn filter is applied
    * in classify_wind_metrics (eligible_YYYY).
    *
    * Requires:
    * - Forest mask asset from Script 1
    * - Windspeed index asset from Script 5a — must have been exported with
    *   an end_year >= this script's wsMaxYear, or the WSmax_{y} band lookup
    *   below will fail (checked explicitly below rather than left to fail
    *   deep inside the pipeline).
    """
    ee_initialize(gee_account_id)

    if end_year is None:
        raise ValueError(
            "end_year must be specified explicitly — agricultural-year pipelines "
            "cannot infer a safe default."
        )

    TREE_COVER_ASSET = f"projects/corestack-datasets-alpha/assets/datasets/SPEI/Hybrid_Tree_AEZ_{aez}_{str(2004)}_{str(end_year)}"
    WIND_INDEX_ASSET = (
        f"projects/corestack-datasets-alpha/assets/datasets/SPEI/wind_index_AEZ_{aez}"
    )

    OUTPUT_DESC = f"wind_metrics_harmonized_kNDVI_AEZ_{aez}"
    OUTPUT_ASSET_ID = (
        f"projects/corestack-datasets-alpha/assets/datasets/SPEI/{OUTPUT_DESC}"
    )

    if is_gee_asset_exists(OUTPUT_ASSET_ID):
        return None

    WIND_THRESHOLD = 15

    BASELINE_START_YEAR = 2004
    BASELINE_END_YEAR = 2024

    # wind_index (5a) has no baseline concept of its own — it just exports
    # whatever start_year..end_year it was called with. So this script,
    # which does need baseline-year WSmax bands, has to verify the asset
    # actually contains them rather than assume it does.
    wsMinYear = min(start_year, BASELINE_START_YEAR)
    wsMaxYear = max(end_year, BASELINE_END_YEAR)

    kndviMinYear = min(start_year, BASELINE_START_YEAR)
    kndviMaxYear = max(end_year + 1, BASELINE_END_YEAR)

    # kndviMaxYear is the furthest single ag-year fetched anywhere in this
    # script; that ag-year isn't complete until Jun 30 of kndviMaxYear+1.
    # Refuse to silently compute on a partial year. (Landsat, unlike
    # ERA5-Land, is current to near-real-time, so no extra lag buffer
    # needed here — that's handled inside 5a already, for its own data.)
    required_complete_by = date(kndviMaxYear + 1, 6, 30)
    if required_complete_by > date.today():
        raise ValueError(
            f"Agricultural year {kndviMaxYear} (Jul {kndviMaxYear} -> "
            f"Jun {kndviMaxYear + 1}) is not yet complete as of "
            f"{date.today().isoformat()}. Reduce end_year."
        )

    # Tree-mask asset names embed end_year and auto-version; wind-index
    # asset names don't, and silently no-op on re-run unless deleted first
    # (see is_gee_asset_exists guard in max_wind_index). That mismatch
    # means it's easy to expand end_year, regenerate the tree mask, forget
    # to regenerate wind_index, and end up reading a stale wind-index asset
    # against a fresh tree-mask one. Fail loudly instead.
    if not is_gee_asset_exists(TREE_COVER_ASSET):
        raise ValueError(
            f"Tree cover asset not found: {TREE_COVER_ASSET}. Run "
            f"generate_hybrid_tree_mask with end_year={end_year} first."
        )
    if not is_gee_asset_exists(WIND_INDEX_ASSET):
        raise ValueError(f"Wind index asset not found: {WIND_INDEX_ASSET}.")

    requiredWsBands = {f"WSmax_{y}" for y in range(wsMinYear, wsMaxYear + 1)}
    availableWsBands = set(ee.Image(WIND_INDEX_ASSET).bandNames().getInfo())
    missingWsBands = sorted(requiredWsBands - availableWsBands)
    if missingWsBands:
        raise ValueError(
            f"{WIND_INDEX_ASSET} is missing bands: {missingWsBands}. It was "
            f"likely exported with a smaller end_year — re-run max_wind_index "
            f"with end_year >= {wsMaxYear} (delete the existing asset first, "
            "its name doesn't auto-version)."
        )

    aoi = ee.FeatureCollection(AEZ).filter(ee.Filter.eq("ae_regcode", aez)).geometry()
    # aoi = (
    #     ee.FeatureCollection("projects/ext-datasets/assets/datasets/State_pan_india")
    #     .filter(ee.Filter.eq("Name", "Odisha"))
    #     .geometry()
    # )
    # Loading the assets :=

    treeMeta = ee.Image(TREE_COVER_ASSET)
    startYear = treeMeta.select("start_year")
    endYear = treeMeta.select("end_year")

    # Load windspeed index from single multiband asset (Script 5a output)
    windIndex_raw = ee.Image(WIND_INDEX_ASSET)

    wsImages = []
    for y in range(wsMinYear, wsMaxYear + 1):
        wsImages.append(
            windIndex_raw.select(f"WSmax_{y}").rename("windspeed").set("year", y)
        )
    wsCol = ee.ImageCollection(wsImages)

    # LANDSAT HARMONIZATION & kNDVI :=

    # Chastain et al. coefficients (OLI to ETM+)
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

        # Apply mask, select optical bands, and apply Collection 2 scale factors
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

        # Apply mask, select optical bands, and apply Collection 2 scale factors
        scaled = (
            image.updateMask(mask)
            .select(["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"])
            .multiply(0.0000275)
            .add(-0.2)
            .rename(chastainBandNames)
        )

        # Apply Chastain regression model (OLI -> ETM+)
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

    # Load kNDVI for kndviMinYear to kndviMaxYear (need next year for resilience)
    kndviYears = ee.List.sequence(kndviMinYear, kndviMaxYear)
    kndviCol = ee.ImageCollection(kndviYears.map(getAnnualKNDVI))

    # BASELINE kNDVI (Yn_bar):=
    # Mean kNDVI across non-high-wind ag-years only

    analysisYears = ee.List.sequence(start_year, end_year)

    # This is my frozen baseline window — the ag-years I use to WORK OUT
    # what Yn_bar (my "normal" kNDVI reference) is, kept separate from
    # analysisYears on purpose.
    baselineYears = ee.List.sequence(BASELINE_START_YEAR, BASELINE_END_YEAR)

    def calc_kndvi_non_event(y):
        year = ee.Number(y)
        kndvi = kndviCol.filter(ee.Filter.eq("year", year)).first()
        ws = (
            wsCol.filter(ee.Filter.eq("year", year))
            .first()
            .resample("bilinear")
            .reproject(crs=kndvi.projection(), scale=30)
        )
        isNonEvent = ws.lte(WIND_THRESHOLD)
        return kndvi.updateMask(isNonEvent).set("year", year)

    kndviNonEvent = ee.ImageCollection(baselineYears.map(calc_kndvi_non_event))

    Yn_bar = kndviNonEvent.mean().rename("kndvi_baseline")

    # SIGNED RESISTANCE & RESILIENCE :=
    def calc_metrics_col(y):
        year = ee.Number(y)

        kndviYe = kndviCol.filter(ee.Filter.eq("year", year)).first()
        wsYe = (
            wsCol.filter(ee.Filter.eq("year", year))
            .first()
            .resample("bilinear")
            .reproject(crs=kndviYe.projection(), scale=30)
        )

        # Only compute on forest pixels during high-windspeed ag-years
        # Flag ANY year where WSmax crosses the threshold, however briefly
        isForest = startYear.lte(year).And(endYear.gte(year))
        isHighWind = wsYe.gt(WIND_THRESHOLD)
        eventMask = isForest.And(isHighWind)

        diffRaw = kndviYe.subtract(Yn_bar)
        diffAbs = diffRaw.abs().max(1e-6)

        # Resistance: signed, computed for ALL high-wind years
        resistance = (
            Yn_bar.divide(diffAbs)
            .multiply(diffRaw.signum())
            .rename("resistance")
            .updateMask(eventMask)
        )

        kndviNext = kndviCol.filter(ee.Filter.eq("year", year.add(1))).first()
        diffNext = kndviNext.subtract(Yn_bar)
        diffNextAbs = diffNext.abs().max(1e-6)

        resilience = (
            diffAbs.divide(diffNextAbs)
            .multiply(diffNext.signum())
            .rename("resilience")
            .updateMask(eventMask)
        )
        eligible = eventMask.And(kndviYe.lt(Yn_bar.multiply(0.95))).rename("eligible")

        return ee.Image.cat([resistance, resilience, eligible]).set("year", year)

    metricsCol = ee.ImageCollection(analysisYears.map(calc_metrics_col))

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


def _class_of(metric):
    low = metric.lt(-20)
    high = metric.gt(20)
    medium = metric.gte(-20).And(metric.lte(20))
    return (
        low.multiply(1)
        .add(medium.multiply(2))
        .add(high.multiply(3))
        .toInt16()
        .updateMask(metric.mask())
    )


def classify_wind_metrics(aez, start_year=2004, end_year=None, gee_account_id=None):
    """Mean yearly wind bands and label low, medium, and high.

    Resistance uses every high-wind year. Resilience uses only eligible years
    (Ye < 0.95 * Yn). Class codes: 1 low, 2 medium, 3 high.
    """
    ee_initialize(gee_account_id)
    if end_year is None:
        raise ValueError("end_year must be specified explicitly.")

    root = "projects/corestack-datasets-alpha/assets/datasets/SPEI"
    yearly_id = f"{root}/wind_metrics_harmonized_kNDVI_AEZ_{aez}"
    output_desc = f"wind_metrics_classified_AEZ_{aez}"
    output_id = f"{root}/{output_desc}"
    if is_gee_asset_exists(output_id):
        return None
    if not is_gee_asset_exists(yearly_id):
        raise ValueError(f"Yearly metrics asset not found: {yearly_id}")

    img = ee.Image(yearly_id)
    resist, resil = [], []
    for year in range(start_year, end_year + 1):
        resist.append(img.select(f"resistance_{year}"))
        resil.append(
            img.select(f"resilience_{year}").updateMask(img.select(f"eligible_{year}"))
        )
    mean_resistance = ee.ImageCollection(resist).mean().rename("resistance")
    mean_resilience = ee.ImageCollection(resil).mean().rename("resilience")
    aoi = ee.FeatureCollection(AEZ).filter(ee.Filter.eq("ae_regcode", aez)).geometry()
    final_output = (
        mean_resistance.addBands(mean_resilience)
        .addBands(_class_of(mean_resistance).rename("resistance_class"))
        .addBands(_class_of(mean_resilience).rename("resilience_class"))
        .clip(aoi)
    )
    return export_raster_asset_to_gee(
        final_output, output_desc, output_id, scale=30, region=aoi
    )
