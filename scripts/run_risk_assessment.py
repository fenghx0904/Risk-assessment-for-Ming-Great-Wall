"""Run hazard classification and integrated risk calculation."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from great_wall_risk.data import load_analysis_input
from great_wall_risk.risk import classify_hazard, compute_risk_indices


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=REPO_ROOT / "data" / "processed" / "analysis_input.xlsx",
        help="Input Excel table containing the 11 factors and heritage attributes.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "results" / "generated" / "risk",
        help="Directory for generated CSV and JSON files.",
    )
    parser.add_argument("--clusters", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--beta", type=float, default=1.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    frame = load_analysis_input(args.input)
    hazard, cluster_summary, thresholds = classify_hazard(
        frame,
        n_clusters=args.clusters,
        seed=args.seed,
    )
    result = compute_risk_indices(frame, hazard, alpha=args.alpha, beta=args.beta)

    result.to_csv(args.output / "analysis_results.csv", index=False, encoding="utf-8-sig")
    cluster_summary.to_csv(
        args.output / "cluster_summary.csv", index=False, encoding="utf-8-sig"
    )
    metadata = {
        "input": str(args.input.resolve()),
        "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "records": len(result),
        "factors": ["EQ", "LS", "WE", "WAE", "FT", "SA", "NDVI", "HRF", "CL", "RND", "GI"],
        "clusters_per_elevation_region": args.clusters,
        "kmeans_seed": args.seed,
        "hazard_class_thresholds": thresholds,
        "risk_equation": "R = H * V / P",
        "linear_alpha": args.alpha,
        "exponential_beta": args.beta,
    }
    (args.output / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(result)} records to {args.output}")


if __name__ == "__main__":
    main()
