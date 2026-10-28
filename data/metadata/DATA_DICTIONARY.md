# Data dictionary

Each row in `data/processed/analysis_input.xlsx` represents one assessment unit. The public file contains 8,656 rows and 19 columns.

| Field | Type | Meaning |
| --- | --- | --- |
| FID | integer | Unique assessment-unit identifier |
| EQ | numeric | Earthquake density value assigned to the unit |
| LS | numeric | Landslide density value assigned to the unit |
| WE | numeric | Wind-erosion value assigned to the unit |
| WAE | numeric | Water-erosion value assigned to the unit |
| FT | numeric | Freeze–thaw value assigned to the unit |
| SA | integer category | Salinization class (0–4 in this table) |
| NDVI | numeric | Normalized difference vegetation index |
| HRF | integer | Heavy-rainfall frequency value (0–10 in this table) |
| CL | binary integer | Cropland indicator (0 or 1) |
| RND | numeric | Road-network density value |
| GI | numeric | Grazing-intensity value |
| Altitude | integer category | Elevation region: 1 ELER, 2 LER, 3 MER, 4 HER, 5 EHER |
| Baohudengji | text | Protection-level label from the heritage attributes |
| Baohudengjizhi | integer | Protection-level score used to calculate P (1–4) |
| Baohuzhuangtai | integer | Preservation-status score used to calculate V (1–6) |
| Caizhiyuanshi | integer | Original wall-material category code |
| Caizhi | text | Wall-material label |
| Caozhizhi | integer | Material vulnerability score used to calculate V (1–4 in this table) |

The public code adds these derived fields:

| Field | Meaning |
| --- | --- |
| H | Hazard class (1 very low to 5 very high) |
| V | Vulnerability index |
| P | Protection-capacity index |
| R | Main integrated risk, H × V / P |
| R_linear | Linear protection-effect sensitivity formulation |
| R_exponential | Exponential protection-effect sensitivity formulation |
| R_reference | Reference formulation H × V |

Values are published as model inputs. They should not be interpreted as the native resolution or precision of the original source product.
