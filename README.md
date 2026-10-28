# GIS-driven environmental risk assessment for the Ming Great Wall

This repository contains the processed model input and reproducible Python workflow used to assess environmental risk along the Ming Great Wall. The analysis combines 11 environmental and human-activity factors with wall material, preservation status and protection level. It produces elevation-stratified hazard classes, vulnerability and protection indices, integrated risk, GeoDetector factor effects, pair interactions and discretization-sensitivity results.

[Data](#data) · [Methods](#methods) · [Installation](#installation) · [Run the analysis](#run-the-analysis) · [Outputs](#outputs) · [Citation](#citation)

## Repository structure

```text
data/
  processed/analysis_input.xlsx   Final model table used by the public scripts
  metadata/                       Field definitions, provenance and checksum
docs/                             Methods and raw-raster download guidance
scripts/                          Command-line entry points
src/great_wall_risk/              Reusable analysis functions
results/generated/                Recreated outputs (ignored by Git)
legacy/                           Earlier scripts preserved for traceability
tests/                            Reproducibility checks
```

## Data

`data/processed/analysis_input.xlsx` contains 8,656 assessment records. Each row is one assessment unit identified by `FID`. The table has no missing values in the fields required by the public workflow and includes the following 11 factors:

| Code | Factor |
| --- | --- |
| EQ | Earthquake density |
| LS | Landslide density |
| WE | Wind erosion |
| WAE | Water erosion |
| FT | Freeze–thaw |
| SA | Salinization |
| NDVI | Normalized difference vegetation index |
| HRF | Heavy-rainfall frequency |
| CL | Cropland distribution |
| RND | Road-network density |
| GI | Grazing intensity |

The large factor rasters are not stored in GitHub. The processed table is sufficient to reproduce the risk and GeoDetector calculations. See [`data/metadata/DATA_DICTIONARY.md`](data/metadata/DATA_DICTIONARY.md) for every field and [`data/metadata/DATA_PROVENANCE.md`](data/metadata/DATA_PROVENANCE.md) for source and period notes. Raw-raster download information will be added to [`docs/BAIDU_NETDISK_UPLOAD_CHECKLIST.md`](docs/BAIDU_NETDISK_UPLOAD_CHECKLIST.md) after the archive is uploaded.

## Methods

The workflow standardizes all 11 factors across the full dataset and runs K-means separately in five elevation regions: ELER (0–500 m), LER (500–1000 m), MER (1000–1500 m), HER (1500–2000 m) and EHER (>2000 m). Five clusters are fitted in each region with `random_state=42` and `n_init=10`. Cluster hazard scores are based on the mean of the globally min–max-normalized factors and are converted to five hazard classes using common thresholds across the corridor.

The main risk formulation is:

```text
V = (material vulnerability + preservation status / 6 × 5) / 2
P = protection level / 4 × 5
R = H × V / P
```

The GeoDetector workflow reports factor q-statistics, all 55 pair interactions within each elevation region, permutation p-values with Benjamini–Hochberg correction, and sensitivity to Jenks, quantile and equal-interval discretization using 3–7 strata. Full details are in [`docs/METHODS.md`](docs/METHODS.md).

## Installation

Python 3.11 is recommended.

```bash
git clone https://github.com/fenghx0904/Risk-assessment-for-Ming-Great-Wall.git
cd Risk-assessment-for-Ming-Great-Wall
python -m venv .venv
```

Activate the environment and install the pinned dependencies:

```bash
# Windows
.venv\Scripts\activate

# macOS or Linux
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Run the analysis

Run the complete public workflow:

```bash
python scripts/run_all.py --permutations 9999
```

For a quick installation check, use fewer permutations:

```bash
python scripts/run_all.py --permutations 99 --output results/generated/smoke_test
```

Individual stages can also be run separately:

```bash
python scripts/run_risk_assessment.py
python scripts/run_geodetector.py --permutations 9999
python scripts/run_interaction_significance.py --permutations 9999
```

The interaction-significance command tests all 55 factor pairs in each elevation region by default (275 tests). Use `--top-factors 3` to test only combinations among the three leading factors in each region.

All scripts accept `--help`. Paths are resolved from the repository root, so the commands work without editing local absolute paths.

## Outputs

Generated files are written under `results/generated/` and are ignored by Git. The main files are:

| File | Contents |
| --- | --- |
| `analysis_results.csv` | Original input plus H, V, P, R and three alternative risk formulations |
| `cluster_summary.csv` | Cluster size, score and hazard class by elevation region |
| `factor_detector.csv` | Factor q, raw permutation p and BH-adjusted p |
| `interaction_detector.csv` | Joint q and interaction type for every factor pair |
| `discretization_sensitivity.csv` | q under three methods and 3–7 requested strata |
| `interaction_significance.csv` | Permutation tests for pair interactions |
| `analysis_settings.json` | Parameters needed to interpret and repeat the run |

## Reproducibility notes

- The public input checksum is recorded in `data/metadata/SHA256SUMS.txt`.
- The official manuscript settings use 9,999 permutations. A smaller value is only a smoke test.
- Permutation tests shuffle the constructed hazard labels within each elevation region. They are descriptive robustness checks and do not preserve spatial autocorrelation or establish causality.
- `legacy/` preserves earlier scripts and tables for provenance. New users should run the entry points in `scripts/`.
- Large rasters, temporary GIS files, credentials, caches and generated results are excluded by `.gitignore`.

## Citation

The paper citation and DOI will be added after publication. Until then, please cite this repository using [`CITATION.cff`](CITATION.cff) and include the repository URL and access date.

## License

Unless a file states otherwise, the repository is released under the [Creative Commons Attribution 4.0 International License](LICENSE). Third-party source datasets retain their original licenses and terms. Users are responsible for obtaining source data under those terms and providing the required attribution.

## Related research from our group

Our group also develops methods for automatically mapping linear cultural heritage from remote-sensing imagery. [SGGWSeg](https://github.com/2022jiangjiazheng/SGGWSeg) is a skeleton-guided dual-task framework that jointly learns semantic segmentation and skeleton representations to improve the continuity and completeness of Great Wall remains extraction, especially for narrow or fragmented structures in complex backgrounds.
