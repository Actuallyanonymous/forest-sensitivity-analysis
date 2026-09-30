import ee

from utilities.constants import AEZ
from utilities.gee_utils import (
    ee_initialize,
    is_gee_asset_exists,
    export_raster_asset_to_gee,
)


def classify_wind_resistance(
    aez, start_year=2004, end_year=None, gee_account_id=None
):
    """
    Wind — stage C.

    Mean the per-year bands from high_wind_sensitivity, then class each
    pixel.

    Resistance uses every high-wind year.
    Resilience uses only years with eligible_YYYY = 1 (Ye < 0.95 * Yn_bar).
    Class codes: 1 = low, 2 = medium, 3 = high. Cuts at ±20.
    Exact ±20 is placed in medium.

    start_year and end_year must match the yearly export.
    """

    ee_initialize(gee_account_id)

    if end_year is None:
        raise ValueError(
            "end_year must be specified explicitly — agricultural-year pipelines "
            "cannot infer a safe default."
        )

    yearly_desc = f"wind_metrics_harmonized_kNDVI_AEZ_{aez}"
    yearly_asset = (
        f"projects/corestack-datasets-alpha/assets/datasets/SPEI/{yearly_desc}"
    )
    output_desc = f"{yearly_desc}_Classified"
    output_asset_id = (
        f"projects/corestack-datasets-alpha/assets/datasets/SPEI/{output_desc}"
    )

    if is_gee_asset_exists(output_asset_id):
        return None

    if not is_gee_asset_exists(yearly_asset):
        raise ValueError(
            f"Yearly wind metrics not found: {yearly_asset}. Run "
            f"high_wind_sensitivity with end_year={end_year} first."
        )

    _require_yearly_bands(yearly_asset, start_year, end_year)

    aoi = ee.FeatureCollection(AEZ).filter(ee.Filter.eq("ae_regcode", aez)).geometry()
    img = ee.Image(yearly_asset)
    final_output = _mean_and_classify(img, start_year, end_year).clip(aoi)

    return export_raster_asset_to_gee(
        final_output, output_desc, output_asset_id, scale=30, region=aoi
    )


def _require_yearly_bands(asset_id, start_year, end_year):
    required = []
    for year in range(start_year, end_year + 1):
        required.extend(
            [f"resistance_{year}", f"resilience_{year}", f"eligible_{year}"]
        )
    available = set(ee.Image(asset_id).bandNames().getInfo())
    missing = [band for band in required if band not in available]
    if missing:
        raise ValueError(
            f"{asset_id} is missing bands: {missing}. It is likely an older "
            "2-band mean asset. Delete it and re-run the yearly metrics export."
        )


def _class_of(metric):
    # 1 = low (metric < -20), 2 = medium (exact ±20 included), 3 = high.
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


def _mean_and_classify(img, start_year, end_year):
    resist_list = []
    resil_list = []
    for year in range(start_year, end_year + 1):
        resist_list.append(img.select(f"resistance_{year}"))
        resil_list.append(
            img.select(f"resilience_{year}").updateMask(
                img.select(f"eligible_{year}")
            )
        )

    mean_resistance = ee.ImageCollection(resist_list).mean().rename("resistance")
    mean_resilience = ee.ImageCollection(resil_list).mean().rename("resilience")
    return (
        mean_resistance.addBands(mean_resilience)
        .addBands(_class_of(mean_resistance).rename("resistance_class"))
        .addBands(_class_of(mean_resilience).rename("resilience_class"))
    )
