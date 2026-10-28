"""Input loading and validation."""

from pathlib import Path

import pandas as pd

from .constants import ELEVATION_REGIONS, FACTORS, INPUT_COLUMNS


def load_analysis_input(path: str | Path) -> pd.DataFrame:
    """Load the model table and validate the fields used by the public code."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    frame = pd.read_excel(path)
    missing = [column for column in INPUT_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")

    if frame["FID"].isna().any() or frame["FID"].duplicated().any():
        raise ValueError("FID must be non-null and unique.")

    numeric_columns = [
        "FID",
        *FACTORS,
        "Altitude",
        "Baohudengjizhi",
        "Baohuzhuangtai",
        "Caizhiyuanshi",
        "Caozhizhi",
    ]
    converted = frame[numeric_columns].apply(pd.to_numeric, errors="coerce")
    if converted.isna().any().any():
        bad = converted.columns[converted.isna().any()].tolist()
        raise ValueError(f"Non-numeric or missing values found in: {bad}")
    frame[numeric_columns] = converted

    unknown_regions = sorted(set(frame["Altitude"].astype(int)) - set(ELEVATION_REGIONS))
    if unknown_regions:
        raise ValueError(f"Unexpected Altitude codes: {unknown_regions}")

    if (frame["Baohudengjizhi"] <= 0).any():
        raise ValueError("Baohudengjizhi must be positive because it is used in the risk denominator.")

    return frame
