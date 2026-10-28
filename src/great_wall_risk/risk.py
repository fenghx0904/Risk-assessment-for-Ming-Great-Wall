"""Hazard classification and integrated risk equations."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

from .constants import ELEVATION_REGIONS, FACTORS


def _classification_thresholds(cluster_scores: np.ndarray, classes: int = 5) -> np.ndarray:
    """Build common class thresholds from equally weighted cluster means."""
    ordered = np.sort(np.asarray(cluster_scores, dtype=float))
    if ordered.size == 0:
        raise ValueError("No cluster scores were produced.")
    groups = np.array_split(ordered, min(classes, ordered.size))
    return np.asarray(
        [
            (lower[-1] + upper[0]) / 2
            for lower, upper in zip(groups[:-1], groups[1:])
            if lower[-1] < upper[0]
        ],
        dtype=float,
    )


def classify_hazard(
    frame: pd.DataFrame,
    factors: list[str] | None = None,
    elevation_regions: dict[int, str] | None = None,
    n_clusters: int = 5,
    seed: int = 42,
) -> tuple[pd.Series, pd.DataFrame, list[float]]:
    """Classify hazard separately within each elevation region.

    All observations are standardized together. K-means is then fitted within
    each elevation region. Cluster hazard scores are the means of the globally
    min-max-normalized factor values. The cluster means share one set of five
    class thresholds across the entire corridor.
    """
    factors = factors or FACTORS
    elevation_regions = elevation_regions or ELEVATION_REGIONS

    scaled = pd.DataFrame(
        StandardScaler().fit_transform(frame[factors].fillna(0)),
        columns=factors,
        index=frame.index,
    )

    factor_min = frame[factors].min()
    factor_range = (frame[factors].max() - factor_min).replace(0, np.nan)
    normalized = ((frame[factors] - factor_min) / factor_range).fillna(0)
    hazard_score = normalized.mean(axis=1)

    cluster_labels = pd.Series(pd.NA, index=frame.index, dtype="Int64")
    summary_rows: list[dict[str, float | int | str]] = []
    for elevation_code, elevation_name in elevation_regions.items():
        mask = frame["Altitude"].astype(int).eq(elevation_code)
        if int(mask.sum()) < max(10, n_clusters):
            raise ValueError(
                f"Elevation region {elevation_name} has only {int(mask.sum())} records; "
                f"at least {max(10, n_clusters)} are required."
            )
        labels = KMeans(
            n_clusters=n_clusters,
            random_state=seed,
            n_init=10,
        ).fit_predict(scaled.loc[mask])
        cluster_labels.loc[mask] = labels
        for label in range(n_clusters):
            member_index = scaled.loc[mask].index[labels == label]
            summary_rows.append(
                {
                    "Elevation_code": elevation_code,
                    "Elevation": elevation_name,
                    "Cluster": label,
                    "Sample_count": len(member_index),
                    "Hazard_score": float(hazard_score.loc[member_index].mean()),
                }
            )

    summary = pd.DataFrame(summary_rows)
    thresholds = _classification_thresholds(summary["Hazard_score"].to_numpy(), classes=5)
    summary["H"] = np.searchsorted(
        thresholds,
        summary["Hazard_score"].to_numpy(),
        side="right",
    ) + 1

    class_lookup = {
        (int(row.Elevation_code), int(row.Cluster)): int(row.H)
        for row in summary.itertuples()
    }
    hazard_class = pd.Series(
        [
            class_lookup[(int(elevation), int(cluster))]
            for elevation, cluster in zip(frame["Altitude"], cluster_labels)
        ],
        index=frame.index,
        name="H",
        dtype=int,
    )
    summary = summary.sort_values(["Elevation_code", "Cluster"]).reset_index(drop=True)
    return hazard_class, summary, thresholds.tolist()


def compute_risk_indices(
    frame: pd.DataFrame,
    hazard_class: pd.Series,
    alpha: float = 0.5,
    beta: float = 1.0,
) -> pd.DataFrame:
    """Append the hazard, vulnerability, protection and risk indices.

    The equations match the manuscript workflow:
    V = (material vulnerability + preservation status / 6 * 5) / 2
    P = protection level / 4 * 5
    R = H * V / P

    Three alternative formulations are retained for sensitivity analysis.
    """
    result = frame.copy()
    result["H"] = hazard_class.astype(int)
    result["V"] = (
        result["Caozhizhi"].astype(float)
        + result["Baohuzhuangtai"].astype(float) / 6.0 * 5.0
    ) / 2.0
    result["P"] = result["Baohudengjizhi"].astype(float) / 4.0 * 5.0
    result["R"] = result["H"] * result["V"] / result["P"]

    normalized_protection = result["P"] / 5.0
    result["R_linear"] = result["H"] * result["V"] * (1.0 - alpha * normalized_protection)
    result["R_exponential"] = result["H"] * result["V"] * np.exp(-beta * normalized_protection)
    result["R_reference"] = result["H"] * result["V"]
    return result
