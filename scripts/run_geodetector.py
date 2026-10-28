"""Run GeoDetector factor, interaction and discretization-sensitivity analyses."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from great_wall_risk.data import load_analysis_input
from great_wall_risk.geodetector import run_geodetector
from great_wall_risk.risk import classify_hazard


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=REPO_ROOT / "data" / "processed" / "analysis_input.xlsx",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "results" / "generated" / "geodetector",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=9999,
        help="Use 9999 for the manuscript analysis; use 99 for a quick smoke test.",
    )
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    frame = load_analysis_input(args.input)
    frame["H"] = classify_hazard(frame, seed=42)[0]
    tables = run_geodetector(
        frame,
        target="H",
        permutations=args.permutations,
        seed=args.seed,
    )
    for name, table in tables.items():
        table.to_csv(args.output / f"{name}.csv", index=False, encoding="utf-8-sig")
    settings = {
        "target": "H",
        "permutations": args.permutations,
        "seed": args.seed,
        "p_value": "(exceedances + 1) / (permutations + 1)",
        "multiple_testing": "Benjamini-Hochberg across 55 baseline region-factor tests",
        "baseline_discretization": "Jenks natural breaks, five strata; categorical exceptions retained",
        "sensitivity_methods": ["Jenks", "Quantile", "Equal interval"],
        "sensitivity_strata": [3, 4, 5, 6, 7],
    }
    (args.output / "analysis_settings.json").write_text(
        json.dumps(settings, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote GeoDetector results to {args.output}")


if __name__ == "__main__":
    main()
