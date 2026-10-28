"""Shared column names and study settings."""

FACTORS = ["EQ", "LS", "WE", "WAE", "FT", "SA", "NDVI", "HRF", "CL", "RND", "GI"]

ELEVATION_REGIONS = {
    1: "ELER",  # 0-500 m
    2: "LER",   # 500-1000 m
    3: "MER",   # 1000-1500 m
    4: "HER",   # 1500-2000 m
    5: "EHER",  # >2000 m
}

CATEGORICAL_FACTORS = {"SA", "CL"}

INPUT_COLUMNS = [
    "FID",
    *FACTORS,
    "Altitude",
    "Baohudengji",
    "Baohudengjizhi",
    "Baohuzhuangtai",
    "Caizhiyuanshi",
    "Caizhi",
    "Caozhizhi",
]
