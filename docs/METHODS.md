# Analysis methods

## Hazard classification

The 11 factors are standardized once across all assessment units with `StandardScaler`. K-means is then fitted separately within each of the five elevation regions using five clusters, `random_state=42` and `n_init=10`.

Each factor is also min–max normalized across the full dataset. A cluster hazard score is the mean normalized factor value across all members and all 11 factors. The 25 cluster scores are sorted, split into five equally weighted groups and separated at the midpoint between adjacent groups. Those common thresholds assign H from 1 to 5.

## Vulnerability, protection and risk

```text
V = (Caozhizhi + Baohuzhuangtai / 6 × 5) / 2
P = Baohudengjizhi / 4 × 5
R = H × V / P
```

Sensitivity formulations are:

```text
R_linear      = H × V × (1 - 0.5 × P/5)
R_exponential = H × V × exp(-1.0 × P/5)
R_reference   = H × V
```

## GeoDetector

For target H and factor X divided into strata h, the factor detector is:

```text
q = 1 - sum_h(N_h × variance_h) / (N × variance)
```

The baseline uses Jenks natural breaks with five requested strata. SA and CL retain their original categories. Constant factors and factors dominated by one value (>95%) are not artificially divided.

Pair interactions use the Cartesian combination of two factor-stratum labels. Interaction type is determined by comparing joint q with the two single-factor q values and their sum.

## Significance and sensitivity

The permutation test shuffles H within each elevation region and uses the upper-tail p-value `(exceedances + 1) / (B + 1)`. The manuscript setting is B = 9,999. Baseline factor tests are adjusted together with Benjamini–Hochberg correction.

Discretization sensitivity compares Jenks natural breaks, quantiles and equal intervals with 3–7 requested strata. Quantile ties may reduce the effective number of strata, which is recorded in the output.

These tests describe the decomposition of a constructed hazard index. The unrestricted permutations do not preserve spatial autocorrelation and should not be interpreted as causal evidence or independent model validation.
