import unittest

import numpy as np
import pandas as pd

from chlorophyll_indices.indices import compute_indices, simulate_bands
from chlorophyll_indices.models import evaluate, fit, leave_group_out_rmse, loo_rmse


class IndexTests(unittest.TestCase):
    def test_index_formulas(self):
        b = {k: np.array([v]) for k, v in {"B3": 0.10, "B4": 0.05, "B5": 0.15, "B6": 0.30, "B7": 0.45}.items()}
        out = {k: float(v[0]) for k, v in compute_indices(b).items()}
        self.assertAlmostEqual(out["NDVI"], 0.40 / 0.50)
        self.assertAlmostEqual(out["GNDVI"], 0.35 / 0.55)
        self.assertAlmostEqual(out["MTCI"], 0.15 / 0.10)
        self.assertAlmostEqual(out["IRECI"], 0.40 / 0.5)
        self.assertAlmostEqual(out["PSSRa"], 9.0)

    def test_band_simulation_point_and_window(self):
        wl = np.arange(400, 901, dtype=float)
        spectrum = (wl / 1000)[:, None]  # reflectance rises linearly with wavelength
        point = simulate_bands(wl, spectrum, "point")
        window = simulate_bands(wl, spectrum, "window")
        self.assertAlmostEqual(float(point["B4"][0]), 0.665)
        self.assertAlmostEqual(float(window["B4"][0]), 0.665, places=3)  # symmetric window keeps the mean
        with self.assertRaises(ValueError):
            simulate_bands(np.arange(400, 600, dtype=float), spectrum[:200], "point")


class ModelTests(unittest.TestCase):
    def test_perfect_fit_and_exact_loo(self):
        x = np.linspace(0, 1, 20)
        stats = fit(x, 3 * x + 2)
        self.assertAlmostEqual(stats["slope"], 3)
        self.assertAlmostEqual(stats["r2"], 1)
        rng = np.random.default_rng(0)
        y = 3 * x + rng.normal(0, 0.1, x.size)
        # Closed-form LOO equals brute-force refitting
        brute = []
        for i in range(x.size):
            keep = np.arange(x.size) != i
            s, c = np.polyfit(x[keep], y[keep], 1)
            brute.append(y[i] - (s * x[i] + c))
        self.assertAlmostEqual(loo_rmse(x, y), float(np.sqrt(np.mean(np.square(brute)))), places=10)

    def test_species_shift_is_penalised_out_of_group(self):
        rng = np.random.default_rng(1)
        x = rng.uniform(0, 1, 60)
        groups = np.repeat(["a", "b", "c"], 20)
        offset = np.select([groups == "a", groups == "b"], [0, 2], 4)  # same slope, different intercepts
        y = x + offset
        self.assertGreater(leave_group_out_rmse(x, y, groups), loo_rmse(x, y))

    def test_evaluate_ranks_by_r2(self):
        x = np.linspace(0, 1, 30)
        df = pd.DataFrame({"good": x, "bad": np.sin(40 * x), "y": 2 * x, "g": np.repeat(["a", "b", "c"], 10)})
        table = evaluate(df, ["bad", "good"], "y", group="g")
        self.assertEqual(table.iloc[0]["index"], "good")


if __name__ == "__main__":
    unittest.main()
