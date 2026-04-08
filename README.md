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
Heavy rainfall index — annual Hm + z-score per pixel (CHIRPS)
↓
Rainfall resistance & resilience (Landsat NDVI vs rainfall baseline)

**Output assets per state:**
- `{STATE}_Hybrid_Tree_Period_2003_2022` — 3 bands: `length`, `start_year`, `end_year`
- `{STATE}_Drought_Metrics` — 2 bands: `resistance`, `resilience`
- `{STATE}_Rain_Index` — 38 bands: `Hm_2004`...`Hm_2022`, `zScore_2004`...`zScore_2022`
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
Both drought and heavy rainfall pipelines follow the same two-stage structure
for methodological consistency:

**Stage 1 — Climate index export**
Compute the climate disturbance index and export as a reusable GEE asset.
- Drought: SPEI-12 (from companion repository below)
- Rainfall: annual heavy rain sum (Hm) and z-score per pixel (Script 3a)

**Stage 2 — Sensitivity metrics**
Load the climate index asset, apply threshold to identify event years,
compute resistance and resilience for forest pixels only.

Metrics are computed using Landsat NDVI time series (30m) relative to a
baseline NDVI (Yn̄) derived from non-event years:

- **Resistance** = Yn̄ / |Ye − Yn̄|
  How strongly the forest maintained its productivity during the event.
- **Resilience** = |Ye − Yn̄| / |Y(e+1) − Yn̄|
  Rate of recovery in the year following the event.

For rainfall metrics both are signed (multiplied by the sign of the NDVI
deviation), so positive = growth/recovery, negative = decline/damage.

**Drought events:** SPEI-12 < -1.0 (moderate to severe drought).

**Rainfall events:** z-score of annual heavy rain sum > 1.0, where heavy
days are defined as daily precipitation exceeding the long-term 95th
percentile (CHIRPS 2000–2023 baseline).

All metrics are computed only on forest pixels and averaged across all
identified event years (2004–2022).

---

## Repository structure

scripts/
├── 01_forest_mask/
│   └── 01_hybrid_forest_mask.js
│       Produces: {STATE}_Hybrid_Tree_Period_2003_2022
│
├── 02_drought_sensitivity/
│   └── 02_drought_resistance_resilience.js
│       Requires: forest mask asset + SPEI-12 assets
│       Produces: {STATE}_Drought_Metrics
│
├── 03_rainfall_sensitivity/
│   ├── 03a_export_rainfall_index.js
│   │   Requires: nothing (public CHIRPS only)
│   │   Produces: {STATE}_Rain_Index
│   │
│   └── 03b_rainfall_resistance_resilience.js
│       Requires: forest mask asset + Rain_Index asset
│       Produces: {STATE}_Rain_Metrics

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
only edit that section. Do not modify anything below it.

### Step 1 — Forest mask (GEE Code Editor)

Open `scripts/01_forest_mask/01_hybrid_forest_mask.js`.

Set your config:
```javascript
var ADMIN_LEVEL     = 'state';            // 'state' or 'district'
var STATE_NAME      = 'Madhya Pradesh';
var DISTRICT_NAME   = '';                 // only if ADMIN_LEVEL = 'district'
var OUTPUT_ASSET_ID = 'projects/your-project/assets/MP_Hybrid_Tree_Period_2003_2022';
var OUTPUT_DESC     = 'MP_Hybrid_Tree_Period_2003_2022';
```

Run the script, go to Tasks tab, click Run.
**Expected time:** 2–6 hours at 30m for a full state.

### Step 2 — Drought resistance & resilience (GEE Code Editor)

Requires: Step 1 asset + SPEI-12 assets from companion repository.

Open `scripts/02_drought_sensitivity/02_drought_resistance_resilience.js`.

Set your config:
```javascript
var TREE_COVER_ASSET  = 'projects/your-project/assets/MP_Hybrid_Tree_Period_2003_2022';
var SPEI_ASSET_PREFIX = 'projects/your-project/assets/SPEI/SPEI12_';
var OUTPUT_ASSET_ID   = 'projects/your-project/assets/MP_Drought_Metrics';
var OUTPUT_DESC       = 'MP_Drought_Metrics';
var STATE_NAME        = 'Madhya Pradesh';
```

Run and submit task.
**Expected time:** 3–6 hours.

### Step 3a — Rainfall index export (GEE Code Editor)

Requires: nothing — uses public CHIRPS dataset only.
Run this in parallel with Steps 1 and 2 if needed.

Open `scripts/03_rainfall_sensitivity/03a_export_rainfall_index.js`.

Set your config:
```javascript
var STATE_NAME      = 'Madhya Pradesh';
var OUTPUT_ASSET_ID = 'projects/your-project/assets/MP_Rain_Index';
var OUTPUT_DESC     = 'MP_Rain_Index';
```

Run and submit task.
**Expected time:** 1–2 hours.
**Output:** 38-band asset — `Hm_2004`...`Hm_2022` and `zScore_2004`...`zScore_2022`.

### Step 3b — Rainfall resistance & resilience (GEE Code Editor)

Requires: Step 1 asset + Step 3a asset.

Open `scripts/03_rainfall_sensitivity/03b_rainfall_resistance_resilience.js`.

Set your config:
```javascript
var TREE_COVER_ASSET = 'projects/your-project/assets/MP_Hybrid_Tree_Period_2003_2022';
var RAIN_INDEX_ASSET = 'projects/your-project/assets/MP_Rain_Index';
var OUTPUT_ASSET_ID  = 'projects/your-project/assets/MP_Rain_Metrics';
var OUTPUT_DESC      = 'MP_Rain_Metrics';
var STATE_NAME       = 'Madhya Pradesh';
```

Run and submit task.
**Expected time:** 3–6 hours.

---

## Execution order and dependencies

Step 3a ─────────────────────────────────┐
▼
Step 1 ──────────┬──► Step 2         Step 3b
│
└──────────────────► Step 3b

Steps 1 and 3a have no dependencies — submit both on the first night.
Steps 2 and 3b both depend on Step 1 — submit after Step 1 completes.

---

## GEE asset structure after all exports

projects/{your-project}/assets/
├── MP_Hybrid_Tree_Period_2003_2022    bands: length, start_year, end_year
├── MP_Drought_Metrics                 bands: resistance, resilience
├── MP_Rain_Index                      bands: Hm_2004...Hm_2022,
│                                             zScore_2004...zScore_2022
├── MP_Rain_Metrics                    bands: resistance, resilience
└── SPEI/
└── SPEI12_Madhya_Pradesh          bands: y2004, y2005, ..., y2023

---

## Companion repository

SPEI-12 computation pipeline (required for Step 2):
https://github.com/Actuallyanonymous/spei-drought-analysis-pipeline

---

## Acknowledgements

Pipeline developed by Pushkin Mangla (IIT Delhi, 2024CS50081) under the
supervision of Prof. Aaditeshwar Seth, Department of CSE, IIT Delhi.

Part of the Core Stack geospatial data framework for climate and forest monitoring.
