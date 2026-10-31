# Forest Sensitivity Analysis Pipeline

Google Earth Engine pipeline for forest resistance and resilience to drought,
heavy rainfall, fire, and high wind. An agricultural year runs from **1 July
of year `y` through 30 June of year `y+1`**. `filterDate` is exclusive at the
end, so every script advances one day past 30 June and that day is included.

Checked-in configuration is Odisha, agricultural years **2004–2024**. In FAO
GAUL 2015 the state name is `Orissa`, not `Odisha`. Change `STATE_NAME` and
every asset id together when you run another state. Do not change
`BASELINE_START_YEAR` or `BASELINE_END_YEAR` after results are published:
those years freeze the normal-year kNDVI and the rain/fire z-score.

Built at IIT Delhi under Prof. Aaditeshwar Seth (ICTD Lab / Core Stack).

---

## Two tracks in this repository

| Track | Where | Geometry | Use |
|---|---|---|---|
| State scripts | `scripts/` | FAO GAUL 2015 level 1 (`ADM1_NAME`) | Run in the Earth Engine Code Editor. This is the track behind the state maps. |
| AEZ tasks | `python-agri-year/` | Agro-ecological zone polygons (`ae_regcode`) | Celery tasks. Same formulas, different polygons and asset names. Not a substitute for the GAUL state scripts. |

The sections below describe the state scripts. The Python track is at the end.

---

## What the state scripts export

Defaults below are the Odisha asset ids in the scripts. Scale is metres.

| Step | Script | Asset | Bands |
|---|---|---|---|
| Forest period | `scripts/01_forest_mask/01_hybrid_forest_mask.js` | `Orissa_Hybrid_Tree_Period_2004_2024` | `length`, `start_year`, `end_year` (30 m) |
| Drought, yearly | `scripts/02_drought_sensitivity/02_drought_resistance_resilience.js` | `Orissa_Drought_Metrics_SPEI3_JAS` | `resistance_YYYY`, `resilience_YYYY`, `resistance_class_YYYY`, `resilience_class_YYYY` |
| Drought, final class | `scripts/02_drought_sensitivity/02c_classify_resistance_resilience.js` | `Orissa_Drought_Metrics_ClassMean` | `resistance_class`, `resilience_class` |
| Rain index | `scripts/03_rainfall_sensitivity/03a_export_rainfall_index.js` | `Orissa_Rain_Index_Wet95_Final` | 105 bands: `Hm`, `zScore`, `heavyDays`, `heavyAvg`, `maxDay` for 2004–2024 (5566 m) |
| Rain, yearly | `scripts/03_rainfall_sensitivity/03b_rainfall_resistance_resilience.js` | `Orissa_Rain_Metrics_Harmonized_kNDVI` | same yearly layout as drought |
| Rain, final class | `scripts/03_rainfall_sensitivity/03c_classify_resistance_resilience.js` | `Orissa_Rain_Metrics_ClassMean` | `resistance_class`, `resilience_class` |
| Fire index | `scripts/04_forestfire_sensitivity/04a_export_fire_index.js` | `Orissa_Fire_Index_FRP30` | 105 bands: `FRP`, `zScore`, `fireDays`, `fireAvg`, `maxFRP` for 2004–2024 (1000 m) |
| Fire, yearly | `scripts/04_forestfire_sensitivity/04b_fire_resistance_resilience.js` | `Orissa_Fire_Metrics_Harmonized_kNDVI` | same yearly layout as drought |
| Fire, final class | `scripts/04_forestfire_sensitivity/04c_classify_resistance_resilience.js` | `Orissa_Fire_Metrics_ClassMean` | `resistance_class`, `resilience_class` |
| Wind index | `scripts/05_highwind_sensitivity/05a_export_maxwind_index.js` | `Orissa_Wind_Index` | 63 bands: `WSmax`, `WShoursGT`, `WSmeanGT` for 2004–2024 (11132 m) |
| Wind, yearly | `scripts/05_highwind_sensitivity/05b_highwind_resistance_resilience.js` | `Orissa_Wind_Metrics_Harmonized_kNDVI` | same yearly layout as drought |
| Wind, final class | `scripts/05_highwind_sensitivity/05c_classify_resistance_resilience.js` | `Orissa_Wind_Metrics_ClassMean` | `resistance_class`, `resilience_class` |

Project prefix in the scripts: `projects/sura-496709/assets/`. Drought reads SPEI-3 from `projects/corestack-datasets-alpha/assets/datasets/SPEI_updated/SPEI3`, band `yYYYY_m07_09` (July–September of calendar year `YYYY`, which sits inside agricultural year `YYYY`).

A year is scored only when the next agricultural year's kNDVI exists. With `END_YEAR = 2024`, Script 2 loads kNDVI through agricultural year 2025 (July 2025–June 2026) and can score 2024. Do not raise `END_YEAR` until that next year is complete. ERA5-Land is published about three months behind real time; the Python wind export refuses an unfinished year, and the JavaScript wind export does not.

---

## Forest mask (Script 1)

Hybrid 30 m mask. Agricultural year `y` is forest when enough sources agree.

| Source | Years | Forest test | How it is placed on agricultural year `y` |
|---|---|---|---|
| GLC-FCS30D | through calendar 2022 | classes 51–92 | Mean of calendar years `y` and `y+1`. One band is enough when the other does not exist. Forest if that mean is at least 0.5. |
| Dynamic World | from 2015 | label 1 (trees) | Mode of labels from 1 July `y` through 30 June `y+1`. |
| IndiaSat LULC v3 | from 2017 | `predicted_label` = 6 | Asset `pan_india_lulc_v3_{y}_{y+1}` is already an agricultural year. |

Vote: one source active → 1 vote; two or more active → at least 2 votes. From 2017, IndiaSat is treated as active even when that year's asset is missing, so a missing IndiaSat year still raises the bar to 2. A non-forest year is filled to forest if any other year within ±2 agricultural years is forest. The fill does not turn neighbours into forest.

The export is the most recent unbroken forest run, counted backward from `END_YEAR` until the first non-forest year: `length`, `start_year`, `end_year`. A pixel is inside the forest period in year `y` when `start_year ≤ y ≤ end_year`.

---

## Landsat kNDVI

Used by Scripts 2, 3b, 4b, and 5b.

- Collections: Landsat 5, 7, 8, and 9 Collection 2 Level-2, Tier 1.
- Cloud and cloud shadow: `QA_PIXEL` bits 3 and 4 clear.
- Surface reflectance scale: `0.0000275 * DN - 0.2`.
- Landsat 8/9 are regressed to Landsat 5/7 with the Chastain OLI→ETM+ slopes and intercepts before the index, so the sensor change is not read as a canopy change.
- `kNDVI = tanh(NDVI²)` from NIR and red.
- `Y` for an agricultural year is the **median** of those scene-level kNDVI values from 1 July through 30 June.

`Ȳn` is the **mean**, across non-event agricultural years 2004–2024, of that yearly median. It is not a median of all scenes in the baseline, and it is not recomputed when `START_YEAR` / `END_YEAR` change.

---

## Resistance, resilience, and the final class

Same cuts for every hazard. `Ye` is kNDVI in the event year. `Ye+1` is the next agricultural year, even if that year is also an event.

```
Rt = Ȳn / max(|Ye − Ȳn|, 1e-6) × sgn(Ye − Ȳn)
Re = |Ye − Ȳn| / max(|Ye+1 − Ȳn|, 1e-6) × sgn(Ye+1 − Ȳn)
```

`sgn(0)` is stored as **+1** (`signum().where(eq(0), 1)`). The class is read from the kNDVI ratio, not from the stored Rt or Re. Each band includes its lower edge and excludes its upper edge. `Ye = Ȳn` is resistance class 3.

**Resistance** (every scored event year):

| Class | `Ye / Ȳn` | Colour | Reading |
|---|---|---|---|
| 1 | `< 0.95` | `#8B0000` | More than 5% below normal |
| 2 | `≥ 0.95` and `< 1` | `#A6DBA0` | Below normal by at most 5% |
| 3 | `≥ 1` and `< 1.05` | `#1B7837` | At normal, or less than 5% above |
| 4 | `≥ 1.05` | `#2166AC` | At least 5% above normal |

**Resilience** is stored only when `Ye / Ȳn < 0.95`. Otherwise the class is 0 (no real drop to recover from).

| Class | Next year | Colour | Reading |
|---|---|---|---|
| 0 | drop was under 5%, or the canopy was not below normal | `#BDBDBD` | Not computed |
| 1 | still `< 0.95`, and `Ye+1 ≤ Ye` | `#8B0000` | No recovery |
| 2 | still `< 0.95`, and `Ye+1 > Ye` | `#FFFFFF` | Recovery started, still more than 5% below |
| 3 | `≥ 0.95` and `< 1` | `#A6DBA0` | Within 5% below normal |
| 4 | `≥ 1` and `< 1.05` | `#1B7837` | Back to normal, within 5% above |
| 5 | `≥ 1.05` | `#2166AC` | At least 5% above normal |

Scripts 2c, 3c, 4c, and 5c average the yearly class numbers. Resistance uses every scored year. Resilience leaves class 0 out of the mean and the count; if every year is 0, the final class is 0. A mean of `n.5` rounds up (`add(0.5).floor()`). A pixel with no scored event year is left blank.

---

## Event rules

The hazard grid is bilinearly resampled to 30 m before the rule is applied. A year enters the score only when the pixel is inside its forest period, the event rule is true, and `Ye+1` exists.

| Hazard | Event year | Years in `Ȳn` | Forest mask on `Ȳn` |
|---|---|---|---|
| Drought | SPEI-3 July–September `< −1` | SPEI-3 July–September `≥ −1` | No. A non-drought year counts even if the pixel was not forest. |
| Heavy rain | `Hm` z-score `> +1` | `\|z\| < 1` on `Hm` | Yes |
| Fire | annual fire-day FRP z-score `> +1` | `\|z\| < 1` on that FRP sum | Yes |
| High wind | `WSmax > 15` m/s | `WSmax ≤ 15` m/s | No. Same rule as drought: a non-event year counts even outside the forest period. |

Exactly −1 on SPEI-3, and exactly +1 or −1 on a z-score, are not events. A dry rainfall year (`z < −1`) is not a heavy-rain event and is not removed from the rain baseline for that reason. A low-fire year (`z < −1`) is not a fire event.

### Rain (Script 3a)

CHIRPS daily precipitation.

- A wet day is daily rain `> 1` mm.
- `P95` is the 95th percentile of wet days from **1 July 2000 through 30 June 2025**. It does not move when the analysis years change.
- A heavy day is daily rain `> P95`.
- `Hm` is the sum of rain on heavy days in the agricultural year.
- `heavyDays` is the count of those days (0 if none).
- `heavyAvg` is the mean rain on heavy days, stored as **0** when there was no heavy day. That 0 is a fill, not a storm.
- `maxDay` is the largest daily total in the year.
- The z-score uses the mean and `ee.Reducer.stdDev()` of `Hm` over agricultural years 2004–2024. Where the standard deviation is 0, the z-score is masked and there is no event.

### Fire (Script 4a)

MODIS/061/MOD14A1 `MaxFRP`. The band is multiplied by **0.1** before any threshold, which converts the stored integer to megawatts.

- A fire day is scaled MaxFRP `> 30` MW. The cut is the same at every pixel.
- `FRP` is the sum of daily MaxFRP on fire days. That sum is what is standardised.
- `fireDays` is 0 when the pixel had no fire day.
- `fireAvg` is the mean MaxFRP on fire days, stored as **0** when there was no fire day.
- `maxFRP` is the largest daily value in the year, after the 0.1 scale.
- The z-score uses the 2004–2024 mean and `ee.Reducer.stdDev()` of `FRP`. A zero standard deviation is masked, and that pixel has no fire event.

### Wind (Script 5a)

ERA5-Land hourly 10 m wind:

```
windspeed = sqrt(u_component_of_wind_10m² + v_component_of_wind_10m²)
```

Three values per agricultural year. There is no z-score in this asset.

| Band | Value |
|---|---|
| `WSmax_YYYY` | Maximum hourly speed in the agricultural year, in m/s. No 10 or 15 m/s cut. |
| `WShoursGT_YYYY` | Number of hourly images with speed **strictly above 10 m/s**. |
| `WSmeanGT_YYYY` | Mean speed on those hours only. Empty when no hour exceeded 10 m/s. Not stored as 0. |

The **10 m/s** cut is only for `WShoursGT` and `WSmeanGT`. The forest event, applied in Script 5b, is **`WSmax > 15`**. A year at exactly 15 m/s is a baseline year.

---

## Run order (Code Editor)

Open a script, edit only the configuration block at the top, Run, then submit the task.

```
Script 1 ──┬── Script 2  ── Script 2c
           ├── Script 3b ── Script 3c
           ├── Script 4b ── Script 4c
           └── Script 5b ── Script 5c

Script 3a ── Script 3b
Script 4a ── Script 4b
Script 5a ── Script 5b
SPEI-3 asset ── Script 2
```

Scripts 1, 3a, 4a, and 5a do not depend on each other. Script 2 also needs the SPEI-3 asset. Scripts 3b, 4b, and 5b need Script 1 plus their index. Each `*c` script needs the yearly metrics asset from `*b`, and it selects `resistance_class_YYYY` through `END_YEAR`, so `END_YEAR` must match the years Script `*b` actually wrote.

### 1. Forest mask

`scripts/01_forest_mask/01_hybrid_forest_mask.js`

```javascript
var ADMIN_LEVEL     = 'state';          // or 'district'
var STATE_NAME      = 'Orissa';         // GAUL 2015 name; Odisha is 'Orissa'
var DISTRICT_NAME   = '';
var START_YEAR      = 2004;
var END_YEAR        = 2024;
var OUTPUT_ASSET_ID = 'projects/your-project/assets/Orissa_Hybrid_Tree_Period_2004_2024';
```

Needs `projects/sat-io/open-datasets/GLC-FCS30D` and `projects/corestack-datasets` (IndiaSat). Dynamic World is public.

### 2. Drought

`scripts/02_drought_sensitivity/02_drought_resistance_resilience.js`

```javascript
var TREE_COVER_ASSET = 'projects/your-project/assets/Orissa_Hybrid_Tree_Period_2004_2024';
var SPEI3_ASSET      = 'projects/corestack-datasets-alpha/assets/datasets/SPEI_updated/SPEI3';
var OUTPUT_ASSET_ID  = 'projects/your-project/assets/Orissa_Drought_Metrics_SPEI3_JAS';
var STATE_NAME       = 'Orissa';
var DROUGHT_THRESHOLD = -1.0;   // SPEI-3 July–September strictly below this
```

The script checks that the tree asset exists and that `y2004_m07_09` … `y2024_m07_09` are on the SPEI image. Then run `02c_classify_resistance_resilience.js`.

### 3. Rain, fire, wind

Point `TREE_COVER_ASSET` at Script 1 and the index asset at 3a, 4a, or 5a. Thresholds in the checked-in scripts:

- Rain: `Z_THRESHOLD = 1.0` (event is z `>` 1).
- Fire: `FRP_THRESHOLD = 30` after the 0.1 scale; event is z `>` 1.
- Wind index: `WIND_THRESHOLD = 10` for the hour bands only.
- Wind metrics: `WIND_THRESHOLD = 15` for the event (`WSmax > 15`).

Then run the matching `*c` script.

---

## SPEI

Script 2 does **not** use SPEI-12. The drought event is the 3-month SPEI that ends in September.

`python-agri-year/generate_spei/` builds the inputs:

1. `download_base_datasets.py` — CHIRPS daily precipitation and MODIS/061/MOD16A2GF PET, one GeoTIFF per month. PET is multiplied by the 0.1 scale factor. Download scale in that script is **5500 m**.
2. `generate_ppet_multiband.py` — monthly P − PET from January `start_year` through June `end_year+1`, resampled onto the CHIRPS grid.
3. `compute_spei.R` — CRAN `SPEI` package, log-Logistic, `na.rm = TRUE`. The reference fit is frozen at January 2004–December 2024 (`ref_start`, `ref_end`). Do not move that window when you extend the series; the parameters are refit from those years, not saved as a separate file.
4. `merge_spei.py` — mosaics the per-AEZ GeoTIFFs.

Reported bands:

- SPEI-1: calendar months `y{year}_m{MM}` from January `start_year` through December `end_year`.
- SPEI-3: quarter-end months `y{year}_m01_03`, `m04_06`, `m07_09`, `m10_12`. Script 2 selects `m07_09`.
- SPEI-12: one June-ending value per agricultural year, named `{year}_{year+1}` (July `year` through June `year+1`). This is the annual water-balance layer. It is not the drought rule in Script 2.

`python-agri-year/spei.py` exposes Celery tasks `generate_spei_pipeline`, `run_drought_resistance_resilience`, `run_rainfall_resistance_resilience`, `run_forest_fire_resistance_resilience`, and `run_high_wind_resistance_resilience`. Those tasks take an AEZ id. `generate_spei_pipeline` with no AEZ runs AEZs 0–19 and then merges SPEI-1, SPEI-3, and SPEI-12.

---

## Python AEZ track

`python-agri-year/` repeats the same agricultural-year formulas on AEZ polygons:

- `hybrid_tree_mask.js` — forest period, written to the Core Stack SPEI asset folder as `Hybrid_Tree_AEZ_{aez}_{start}_{end}`.
- `drought_sensitivity/drought_resistance_resilience.py`
- `rainfall_sensitivity/export_rainfall_index.py` and `rainfall_resistance_resilience.py`
- `forestfire_sensitivity/export_fire_index.py` and `forest_fire_resistance_resilience.py`
- `high_wind_sensitivity/export_max_wind_index.py` and `highwind_resistance_resilience.py`

Rain exports at 5566 m, fire at 1000 m, wind at 11132 m, and the kNDVI metrics at 30 m. Python refuses an agricultural year that is not complete (30 June of `end_year+1`). Wind adds the ERA5-Land lag of about three months. MOD16 downloads are refused until 1 March of `end_year+2`, because agricultural year `end_year` needs calendar `end_year+1` PET.

Imports in `spei.py` expect this tree to be installed as `computing.spei` inside the Core Stack app (`utilities.gee_utils`, `utilities.constants.AEZ`, `nrm_app.celery`). Running the files from a bare clone will not resolve those imports.

---

## Prerequisites

- An Earth Engine account that can read:
  - `projects/corestack-datasets` (IndiaSat LULC v3)
  - `projects/sat-io/open-datasets/GLC-FCS30D`
  - `projects/corestack-datasets-alpha` (SPEI-3 asset used by Script 2)
- Public collections: Landsat C2 L2, Dynamic World, CHIRPS daily, MODIS/061/MOD14A1, ERA5-Land hourly. SPEI generation also needs MODIS/061/MOD16A2GF.
- For a local SPEI rebuild: Python with `earthengine-api` and `rasterio`, plus R packages `SPEI`, `raster`, and `terra`.

---

## Acknowledgements

Pipeline developed by Pushkin Mangla (IIT Delhi, 2024CS50081) under the supervision of Prof. Aaditeshwar Seth, Department of CSE, IIT Delhi.

Part of the Core Stack geospatial data framework for climate and forest monitoring.
