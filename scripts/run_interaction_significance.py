"""Run permutation significance tests for GeoDetector factor pairs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from great_wall_risk.data import load_analysis_input
from great_wall_risk.geodetector import run_interaction_significance
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
        default=REPO_ROOT / "results" / "generated" / "interaction_significance",
    )
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20261016)
    parser.add_argument(
        "--top-factors",
        type=int,
        default=0,
        help="0 tests all 55 pairs per elevation region; 3 tests the three leading factors.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    frame = load_analysis_input(args.input)
    frame["H"] = classify_hazard(frame, seed=42)[0]
    table = run_interaction_significance(
        frame,
        target="H",
        permutations=args.permutations,
        seed=args.seed,
        top_factors=args.top_factors,
    )
    table.to_csv(
        args.output / "interaction_significance.csv",
        index=False,
        encoding="utf-8-sig",
    )
    settings = {
        "target": "H",
        "permutations": args.permutations,
        "seed": args.seed,
        "top_factors": args.top_factors,
        "tests": len(table),
        "multiple_testing": "Benjamini-Hochberg across all reported pair tests",
    }
    (args.output / "analysis_settings.json").write_text(
        json.dumps(settings, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote {len(table)} interaction tests to {args.output}")


if __name__ == "__main__":
    main()
