"""Run the public risk and GeoDetector workflow in sequence."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from great_wall_risk.data import load_analysis_input
from great_wall_risk.geodetector import run_geodetector
from great_wall_risk.risk import classify_hazard, compute_risk_indices


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=REPO_ROOT / "data" / "processed" / "analysis_input.xlsx",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "results" / "generated" / "complete_run",
    )
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    frame = load_analysis_input(args.input)
    hazard, clusters, _ = classify_hazard(frame, seed=args.seed)
    risk = compute_risk_indices(frame, hazard)
    risk.to_csv(args.output / "analysis_results.csv", index=False, encoding="utf-8-sig")
    clusters.to_csv(args.output / "cluster_summary.csv", index=False, encoding="utf-8-sig")

    tables = run_geodetector(
        risk,
        target="H",
        permutations=args.permutations,
        seed=args.seed,
    )
    for name, table in tables.items():
        table.to_csv(args.output / f"{name}.csv", index=False, encoding="utf-8-sig")
    print(f"Complete workflow finished: {args.output}")


if __name__ == "__main__":
    main()
