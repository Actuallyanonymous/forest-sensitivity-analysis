# Forest Sensitivity Analysis Pipeline

A Google Earth Engine pipeline for computing forest resistance and resilience
to climate disturbances across India.

Built at IIT Delhi under Prof. Aaditeshwar Seth (ICTD Lab / Core Stack).

---

## What this pipeline does

Hybrid 30m forest mask (GLC-FCS30D + Dynamic World + IndiaSat)
↓
Most recent contiguous forest period per pixel (length, start year, end year)
↓
SPEI-12 drought index (from spei-drought-analysis-pipeline)
↓
Drought resistance & resilience (Landsat NDVI vs SPEI baseline)
↓
Heavy rainfall resistance & resilience (CHIRPS anomaly years)
↓
GEE interactive app — all layers visualised

**Output assets per state:**
- `{STATE}_Hybrid_Tree_Period_2003_2022` — 3 bands: `length`, `start_year`, `end_year`
- `{STATE}_Drought_Metrics` — 2 bands: `resistance`, `resilience`
- `{STATE}_Rain_Metrics` — 2 bands: `resistance`, `resilience`

---

## Methodology

### Forest mask
A hybrid 30m annual tree cover mask combining three datasets:
- **GLC-FCS30D** (2003–2022) — classes 51–92 classified as forest
- **Dynamic World** (2015–2022) — label 1 (trees)
- **IndiaSat LULC** (2017–2022) — classes 11, 12, 14 (natural forest types)

Union logic: majority vote among active datasets per year. Years with only
one active dataset (2003–2014) require 1 vote; years with 2–3 active datasets
require at least 2 to agree.

A ±2 year temporal correction fills isolated non-forest pixels surrounded
by forest in neighbouring years.

The most recent *unbroken* forest period per pixel is extracted by iterating
backwards from 2022, stopping at the first non-forest year.

### Resistance & resilience
For both drought and heavy rainfall events, resistance and resilience are
computed using Landsat NDVI time series (30m) relative to a baseline
NDVI (Yn̄) derived from non-event years:

- **Resistance** = Yn̄ / |Ye − Yn̄| (how much NDVI deviated during event)
- **Resilience** = |Ye − Yn̄| / |Y(e+1) − Yn̄| (rate of recovery after event)

**Drought events:** identified using SPEI-12 < -1.0.
SPEI assets come from the companion repository:
https://github.com/Actuallyanonymous/spei-drought-analysis-pipeline

**Rainfall events:** identified using z-score of annual heavy rainfall sum > 1.0,
where heavy days are defined as daily precipitation exceeding the long-term 95th
percentile (CHIRPS).

All metrics are computed only on forest pixels (using the forest mask above)
and averaged across all event years.

---

## Repository structure

scripts/
├── 01_forest_mask/
│   └── 01_hybrid_forest_mask.js          ← Run first. Produces forest period asset.
│
├── 02_drought_sensitivity/
│   └── 02_drought_resistance_resilience.js  ← Requires forest mask + SPEI assets.
│
├── 03_rainfall_sensitivity/
│   └── 03_rainfall_resistance_resilience.js ← Requires forest mask only.
│
└── 04_gee_app/
└── 04_forest_climate_monitor_app.js   ← Visualisation app. Run after all exports.

---

## Prerequisites

- Google Earth Engine account with access to:
  - `projects/corestack-datasets` (IndiaSat LULC)
  - `projects/sat-io/open-datasets/GLC-FCS30D` (GLC forest mask)
- SPEI-12 assets uploaded to your GEE project
  (see: https://github.com/Actuallyanonymous/spei-drought-analysis-pipeline)

---

## Quick start

Run scripts in order. Each script has a `CONFIGURATION` section at the top —
only edit that section.

### Step 1 — Forest mask (GEE Code Editor)

Open `scripts/01_forest_mask/01_hybrid_forest_mask.js`.

Set your config:
```javascript
var ADMIN_LEVEL     = 'state';
var STATE_NAME      = 'Madhya Pradesh';
var OUTPUT_ASSET_ID = 'projects/your-project/assets/MP_Hybrid_Tree_Period_2003_2022';
var OUTPUT_DESC     = 'MP_Hybrid_Tree_Period_2003_2022';
```

Run the script, go to Tasks tab, click Run.
**Expected time:** 2–6 hours at 30m for a full state.

### Step 2 — Drought resistance & resilience (GEE Code Editor)

Open `scripts/02_drought_sensitivity/02_drought_resistance_resilience.js`.

Set your config:
```javascript
var TREE_COVER_ASSET = 'projects/your-project/assets/MP_Hybrid_Tree_Period_2003_2022';
var SPEI_ASSET_PREFIX = 'projects/your-project/assets/SPEI/SPEI12_Madhya_Pradesh';
var OUTPUT_ASSET_ID  = 'projects/your-project/assets/MP_Drought_Metrics';
```

Run and submit task.
**Expected time:** 3–6 hours.

### Step 3 — Rainfall resistance & resilience (GEE Code Editor)

Open `scripts/03_rainfall_sensitivity/03_rainfall_resistance_resilience.js`.

Set your config:
```javascript
var TREE_COVER_ASSET = 'projects/your-project/assets/MP_Hybrid_Tree_Period_2003_2022';
var OUTPUT_ASSET_ID  = 'projects/your-project/assets/MP_Rain_Metrics';
```

Run and submit task.

### Step 4 — GEE app (GEE Code Editor)

Open `scripts/04_gee_app/04_forest_climate_monitor_app.js`.

Set asset paths at the top to match your project, run the script,
then go to Apps → New App to publish.

---

## GEE asset structure

projects/{your-project}/assets/
├── MP_Hybrid_Tree_Period_2003_2022    bands: length, start_year, end_year
├── MP_Drought_Metrics                 bands: resistance, resilience
├── MP_Rain_Metrics                    bands: resistance, resilience
└── SPEI/
└── SPEI12_Madhya_Pradesh          bands: y2004, y2005, ..., y2023

---

## Companion repository

SPEI computation pipeline (Step 2 dependency):
https://github.com/Actuallyanonymous/spei-drought-analysis-pipeline

---

## Acknowledgements

Pipeline developed by Pushkin Mangla (IIT Delhi, 2024CS50081) under the
supervision of Prof. Aaditeshwar Seth, Department of CSE, IIT Delhi.

Part of the Core Stack geospatial data framework for climate and forest monitoring.
