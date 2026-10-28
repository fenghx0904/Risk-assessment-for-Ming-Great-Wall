from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from great_wall_risk.data import load_analysis_input
from great_wall_risk.geodetector import q_statistic
from great_wall_risk.risk import classify_hazard, compute_risk_indices


class AnalysisTests(unittest.TestCase):
    def test_public_input_and_risk_equation(self):
        frame = load_analysis_input(ROOT / "data" / "processed" / "analysis_input.xlsx")
        self.assertEqual(frame.shape, (8656, 19))
        hazard, clusters, thresholds = classify_hazard(frame)
        self.assertTrue(set(hazard.unique()).issubset({1, 2, 3, 4, 5}))
        self.assertEqual(len(clusters), 25)
        self.assertEqual(len(thresholds), 4)

        result = compute_risk_indices(frame, hazard)
        expected = result["H"] * result["V"] / result["P"]
        self.assertTrue(np.allclose(result["R"], expected))

    def test_q_statistic_bounds(self):
        target = np.array([1.0, 1.0, 4.0, 4.0])
        strata = np.array([0, 0, 1, 1])
        self.assertTrue(np.isclose(q_statistic(target, strata), 1.0))
        mixed = q_statistic(target, np.array([0, 1, 0, 1]))
        self.assertGreaterEqual(mixed, 0.0)
        self.assertLessEqual(mixed, 1.0)


if __name__ == "__main__":
    unittest.main()
