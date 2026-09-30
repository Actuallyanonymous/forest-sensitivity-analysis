"""Checks the agricultural-year resistance/resilience scripts follow the
stage-b / stage-c split used by the GEE JavaScript.

Stage B stores signed resistance and resilience for every event year, plus
eligible_YYYY when the pixel is in an event and Ye < 0.95 * Yn_bar.
Stage C averages resistance over event years, averages resilience only where
eligible, and classes the means at ±20 (exact ±20 is medium).
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]

METRIC_FILES = [
    ROOT / "drought_sensitivity" / "drought_resistance_resilience.py",
    ROOT / "rainfall_sensitivity" / "rainfall_resistance_resilience.py",
    ROOT / "forestfire_sensitivity" / "forest_fire_resistance_resilience.py",
    ROOT / "high_wind_sensitivity" / "highwind_resistance_resilience.py",
]

CLASSIFY_FILES = [
    ROOT / "drought_sensitivity" / "classify_resistance_resilience.py",
    ROOT / "rainfall_sensitivity" / "classify_resistance_resilience.py",
    ROOT / "forestfire_sensitivity" / "classify_resistance_resilience.py",
    ROOT / "high_wind_sensitivity" / "classify_resistance_resilience.py",
]


def resistance(yn_bar, ye, eps=1e-6):
    diff = ye - yn_bar
    magnitude = max(abs(diff), eps)
    sign = 0.0 if diff == 0 else (1.0 if diff > 0 else -1.0)
    return yn_bar / magnitude * sign


def resilience(yn_bar, ye, ye_next, eps=1e-6):
    diff = ye - yn_bar
    diff_next = ye_next - yn_bar
    magnitude = max(abs(diff), eps)
    magnitude_next = max(abs(diff_next), eps)
    sign = 0.0 if diff_next == 0 else (1.0 if diff_next > 0 else -1.0)
    return magnitude / magnitude_next * sign


def is_eligible(in_event, ye, yn_bar):
    return bool(in_event) and ye < 0.95 * yn_bar


def class_code(metric):
    if metric < -20:
        return 1
    if metric > 20:
        return 3
    return 2


def aggregate(year_rows):
    """year_rows are (resistance or None, resilience or None, eligible).

    None resistance means the pixel was not an event year (masked), matching
    ImageCollection.mean, which skips masked pixels.
    """
    resist_vals = [row[0] for row in year_rows if row[0] is not None]
    resil_vals = [row[1] for row in year_rows if row[2] and row[1] is not None]
    mean_r = sum(resist_vals) / len(resist_vals)
    mean_d = sum(resil_vals) / len(resil_vals)
    return mean_r, mean_d, class_code(mean_r), class_code(mean_d)


class SignedFormulaTest(unittest.TestCase):
    def test_boundary_resistance_matches_readme(self):
        yn_bar = 0.40
        self.assertAlmostEqual(resistance(yn_bar, 0.95 * yn_bar), -20.0)
        self.assertAlmostEqual(resistance(yn_bar, 1.05 * yn_bar), 20.0)

    def test_resilience_sign_follows_next_year(self):
        yn_bar = 0.40
        ye = 0.30
        self.assertGreater(resilience(yn_bar, ye, 0.41), 0)
        self.assertLess(resilience(yn_bar, ye, 0.39), 0)

    def test_eligibility_is_event_and_strictly_below_95_percent(self):
        yn_bar = 0.40
        self.assertTrue(is_eligible(True, 0.949 * yn_bar, yn_bar))
        self.assertFalse(is_eligible(True, 0.95 * yn_bar, yn_bar))
        self.assertFalse(is_eligible(False, 0.90 * yn_bar, yn_bar))

    def test_exact_bounds_are_medium(self):
        self.assertEqual(class_code(-20), 2)
        self.assertEqual(class_code(20), 2)
        self.assertEqual(class_code(-20.0001), 1)
        self.assertEqual(class_code(20.0001), 3)

    def test_resilience_mean_skips_ineligible_event_years(self):
        yn_bar = 0.40
        # Event, below the 5% line: eligible. Next year just above normal.
        ye_hit = 0.30
        # Event, but only slightly below normal: resistance is stored,
        # resilience is not eligible.
        ye_mild = 0.39
        rows = [
            (
                resistance(yn_bar, ye_hit),
                resilience(yn_bar, ye_hit, yn_bar + 0.001),
                is_eligible(True, ye_hit, yn_bar),
            ),
            (
                resistance(yn_bar, ye_mild),
                resilience(yn_bar, ye_mild, 0.10),
                is_eligible(True, ye_mild, yn_bar),
            ),
            (None, None, False),
        ]
        mean_r, mean_d, r_class, d_class = aggregate(rows)
        self.assertAlmostEqual(
            mean_r, (rows[0][0] + rows[1][0]) / 2
        )
        self.assertAlmostEqual(mean_d, rows[0][1])
        self.assertIn(r_class, (1, 2, 3))
        self.assertIn(d_class, (1, 2, 3))


class StageSplitSourceTest(unittest.TestCase):
    def test_metric_scripts_keep_yearly_bands(self):
        for path in METRIC_FILES:
            source = path.read_text()
            self.assertIn("eligible", source, path.name)
            self.assertIn('rename(f"resistance_{yy}")', source, path.name)
            self.assertIn('rename(f"resilience_{yy}")', source, path.name)
            self.assertIn('rename(f"eligible_{yy}")', source, path.name)
            self.assertIn(".median()", source, path.name)
            self.assertRegex(
                source,
                r"\.mean\(\)\s*\.rename\(\"kndvi_baseline\"\)",
                path.name,
            )
            self.assertIn("Yn_bar.divide(diffAbs)", source, path.name)
            self.assertIn("diffAbs.divide(diffNextAbs)", source, path.name)
            self.assertIn("Yn_bar.multiply(0.95)", source, path.name)
            self.assertIn(".updateMask(eventMask)", source, path.name)
            self.assertNotIn('select("resistance").mean()', source, path.name)
            self.assertNotIn("isNegativeEffect", source, path.name)
            self.assertNotIn("ONLY when Ye < Yn_bar", source, path.name)

    def test_classify_scripts_average_then_cut_at_20(self):
        for path in CLASSIFY_FILES:
            source = path.read_text()
            self.assertIn("metric.lt(-20)", source, path.name)
            self.assertIn("metric.gt(20)", source, path.name)
            self.assertIn("metric.gte(-20)", source, path.name)
            self.assertIn("metric.lte(20)", source, path.name)
            self.assertIn(".updateMask(", source, path.name)
            self.assertIn('f"eligible_{year}"', source, path.name)
            self.assertIn('.rename("resistance_class")', source, path.name)
            self.assertIn('.rename("resilience_class")', source, path.name)
            self.assertIn("ee.ImageCollection(resist_list).mean()", source, path.name)
            self.assertIn("ee.ImageCollection(resil_list).mean()", source, path.name)

    def test_runners_classify_after_yearly_export(self):
        source = (ROOT / "spei.py").read_text()
        self.assertIn("classify_drought_resistance", source)
        self.assertIn("classify_rainfall_resistance", source)
        self.assertIn("classify_fire_resistance", source)
        self.assertIn("classify_wind_resistance", source)


if __name__ == "__main__":
    unittest.main()
