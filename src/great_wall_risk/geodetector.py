"""GeoDetector factor, interaction and robustness analyses."""

from __future__ import annotations

import itertools

import jenkspy
import numpy as np
import pandas as pd

from .constants import CATEGORICAL_FACTORS, ELEVATION_REGIONS, FACTORS


def _as_codes(values: pd.Series) -> np.ndarray:
    """Convert strata labels to contiguous integer codes."""
    codes, _ = pd.factorize(values, sort=True)
    if (codes < 0).any():
        raise ValueError("GeoDetector strata contain missing values.")
    return codes.astype(int)


def discretize(
    series: pd.Series,
    factor: str,
    method: str = "jenks",
    strata: int = 5,
) -> tuple[np.ndarray, int, str]:
    """Discretize a factor while retaining categorical and degenerate inputs."""
    values = pd.to_numeric(series, errors="raise")
    unique = int(values.nunique())
    if unique < 2:
        return _as_codes(values), unique, "constant"
    if factor in CATEGORICAL_FACTORS:
        codes = _as_codes(values)
        return codes, int(np.unique(codes).size), "original categorical"
    if float(values.value_counts(normalize=True).iloc[0]) > 0.95:
        codes = _as_codes(values)
        return codes, int(np.unique(codes).size), "original >95% dominant-value rule"

    requested = min(int(strata), unique)
    method_key = method.lower().replace("_", " ").strip()
    if method_key == "jenks":
        try:
            edges = np.asarray(
                jenkspy.jenks_breaks(values.to_numpy(), n_classes=requested),
                dtype=float,
            )
            if np.unique(edges).size != edges.size:
                raise ValueError("Jenks returned duplicate edges.")
            binned = pd.cut(values, bins=edges, labels=False, include_lowest=True)
            note = "Jenks natural breaks"
        except Exception:
            binned = pd.qcut(
                values.rank(method="first"),
                q=requested,
                labels=False,
                duplicates="drop",
            )
            note = "Jenks fallback: rank-based quantiles"
    elif method_key == "quantile":
        binned = pd.qcut(values, q=requested, labels=False, duplicates="drop")
        note = "quantile; ties retained"
    elif method_key in {"equal interval", "equal-interval"}:
        binned = pd.cut(values, bins=requested, labels=False, include_lowest=True)
        note = "equal interval"
    else:
        raise ValueError(f"Unknown discretization method: {method}")

    codes = _as_codes(binned)
    return codes, int(np.unique(codes).size), note


def q_statistic(target: np.ndarray | pd.Series, codes: np.ndarray) -> float:
    """Return GeoDetector q = 1 - within-stratum SS / total SS."""
    y = np.asarray(target, dtype=float)
    if y.ndim != 1 or y.size != len(codes):
        raise ValueError("Target and strata must be one-dimensional and have equal length.")
    centered = y - y.mean()
    total = float(centered @ centered)
    if total <= 0 or np.unique(codes).size < 2:
        return 0.0
    counts = np.bincount(codes)
    sums = np.bincount(codes, weights=centered)
    between = np.sum((sums[counts > 0] ** 2) / counts[counts > 0])
    return float(np.clip(between / total, 0.0, 1.0))


def permutation_p_value(
    target: np.ndarray | pd.Series,
    codes: np.ndarray,
    permutations: int,
    seed: int,
) -> tuple[float, float]:
    """Upper-tail random-label test using p=(exceedances+1)/(B+1)."""
    observed = q_statistic(target, codes)
    if permutations <= 0:
        return observed, np.nan
    y = np.asarray(target, dtype=float)
    if np.ptp(y) == 0 or np.unique(codes).size < 2:
        return observed, 1.0

    rng = np.random.default_rng(seed)
    centered = y - y.mean()
    total = float(centered @ centered)
    counts = np.bincount(codes)
    present = counts > 0
    exceedances = 0
    for _ in range(permutations):
        sums = np.bincount(codes, weights=rng.permutation(centered))
        permuted_q = float(np.sum((sums[present] ** 2) / counts[present]) / total)
        exceedances += permuted_q >= observed - 1e-12
    return observed, (exceedances + 1.0) / (permutations + 1.0)


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjusted p values, preserving missing entries."""
    result = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.dropna().astype(float)
    if valid.empty:
        return result
    order = np.argsort(valid.to_numpy())
    ranked = valid.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    result.loc[valid.index] = restored
    return result


def interaction_type(q_a: float, q_b: float, q_ab: float) -> str:
    """Classify a GeoDetector pair interaction."""
    tolerance = 1e-12
    if np.isclose(q_ab, q_a + q_b, rtol=0, atol=tolerance):
        return "Independent"
    if q_ab > q_a + q_b:
        return "Nonlinear enhancement"
    if q_ab > max(q_a, q_b) + tolerance:
        return "Bivariate enhancement"
    if q_ab < min(q_a, q_b) - tolerance:
        return "Nonlinear weakening"
    if q_ab < max(q_a, q_b) - tolerance:
        return "Univariate weakening"
    return "Boundary"


def run_geodetector(
    frame: pd.DataFrame,
    target: str = "H",
    factors: list[str] | None = None,
    permutations: int = 9999,
    seed: int = 42,
) -> dict[str, pd.DataFrame]:
    """Run baseline factor/interaction detectors and discretization sensitivity."""
    factors = factors or FACTORS
    factor_rows: list[dict] = []
    interaction_rows: list[dict] = []
    sensitivity_rows: list[dict] = []

    for elevation_code, elevation_name in ELEVATION_REGIONS.items():
        subset = frame.loc[frame["Altitude"].astype(int).eq(elevation_code)].copy()
        y = subset[target].to_numpy(dtype=float)
        baseline_codes: dict[str, np.ndarray] = {}
        baseline_q: dict[str, float] = {}

        for factor_index, factor in enumerate(factors):
            codes, effective, note = discretize(subset[factor], factor, "jenks", 5)
            q_value, raw_p = permutation_p_value(
                y,
                codes,
                permutations=permutations,
                seed=seed + elevation_code * 100 + factor_index,
            )
            baseline_codes[factor] = codes
            baseline_q[factor] = q_value
            factor_rows.append(
                {
                    "Elevation_code": elevation_code,
                    "Elevation": elevation_name,
                    "Factor": factor,
                    "q": q_value,
                    "p_raw": raw_p,
                    "Requested_strata": 5,
                    "Effective_strata": effective,
                    "Discretization": note,
                    "Sample_count": len(subset),
                }
            )

        for factor_a, factor_b in itertools.combinations(factors, 2):
            codes_a = baseline_codes[factor_a]
            codes_b = baseline_codes[factor_b]
            joint_codes = _as_codes(pd.Series(list(zip(codes_a, codes_b))))
            joint_q = q_statistic(y, joint_codes)
            interaction_rows.append(
                {
                    "Elevation_code": elevation_code,
                    "Elevation": elevation_name,
                    "Factor_A": factor_a,
                    "Factor_B": factor_b,
                    "q_A": baseline_q[factor_a],
                    "q_B": baseline_q[factor_b],
                    "q_interaction": joint_q,
                    "Type": interaction_type(
                        baseline_q[factor_a], baseline_q[factor_b], joint_q
                    ),
                }
            )

        for method in ("jenks", "quantile", "equal interval"):
            for requested in range(3, 8):
                for factor in factors:
                    codes, effective, note = discretize(
                        subset[factor], factor, method, requested
                    )
                    sensitivity_rows.append(
                        {
                            "Elevation_code": elevation_code,
                            "Elevation": elevation_name,
                            "Factor": factor,
                            "Method": method,
                            "Requested_strata": requested,
                            "Effective_strata": effective,
                            "q": q_statistic(y, codes),
                            "Note": note,
                        }
                    )

    factor_table = pd.DataFrame(factor_rows)
    factor_table["p_BH"] = benjamini_hochberg(factor_table["p_raw"])
    factor_table["Significant_0.05"] = factor_table["p_BH"].lt(0.05)
    return {
        "factor_detector": factor_table,
        "interaction_detector": pd.DataFrame(interaction_rows),
        "discretization_sensitivity": pd.DataFrame(sensitivity_rows),
    }


def run_interaction_significance(
    frame: pd.DataFrame,
    target: str = "H",
    factors: list[str] | None = None,
    permutations: int = 9999,
    seed: int = 20261016,
    top_factors: int = 0,
) -> pd.DataFrame:
    """Permutation tests for pairwise interaction q values.

    Set ``top_factors`` to 0 to test all 55 pairs in each elevation region.
    A positive value tests combinations among the top-ranked baseline factors.
    """
    factors = factors or FACTORS
    rows: list[dict] = []
    test_index = 0
    for elevation_code, elevation_name in ELEVATION_REGIONS.items():
        subset = frame.loc[frame["Altitude"].astype(int).eq(elevation_code)].copy()
        y = subset[target].to_numpy(dtype=float)
        codes = {
            factor: discretize(subset[factor], factor, "jenks", 5)[0]
            for factor in factors
        }
        q_values = {factor: q_statistic(y, code) for factor, code in codes.items()}
        selected = factors
        if top_factors > 0:
            selected = sorted(factors, key=q_values.get, reverse=True)[:top_factors]
        for factor_a, factor_b in itertools.combinations(selected, 2):
            joint = _as_codes(pd.Series(list(zip(codes[factor_a], codes[factor_b]))))
            joint_q, raw_p = permutation_p_value(
                y,
                joint,
                permutations=permutations,
                seed=seed + test_index,
            )
            test_index += 1
            rows.append(
                {
                    "Elevation_code": elevation_code,
                    "Elevation": elevation_name,
                    "Factor_A": factor_a,
                    "Factor_B": factor_b,
                    "q_A": q_values[factor_a],
                    "q_B": q_values[factor_b],
                    "q_interaction": joint_q,
                    "p_raw": raw_p,
                    "Type": interaction_type(q_values[factor_a], q_values[factor_b], joint_q),
                }
            )
    table = pd.DataFrame(rows)
    table["p_BH"] = benjamini_hochberg(table["p_raw"])
    table["Significant_0.05"] = table["p_BH"].lt(0.05)
    return table
