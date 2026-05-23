# Mouse epilepsy connectome experiment (Daria / Anderson cohort)

Apply the §4.5 workflow from Liao's thesis (originally implemented for human data in `epileptor-experiment/runEZ1_section45.m`) to mouse coarse-grained connectomes (66 brain-atlas regions per animal).

## Prerequisites

- **MATLAB Optimization Toolbox** — `fsolve` for Epileptor fixed points (`tvb` / default schemes).
- **Brain Connectivity Toolbox (BCT)** — optional for stability-only runs; **required** for betweenness/closeness in `computeNetworkCentralities` and for the betweenness panel in `compareAndersonVsArnold`. BCT is **not** included in git; install locally under `../2019_03_03_BCT/` — see [docs/BCT.md](../docs/BCT.md).

## Layout

| Path | Role |
|------|------|
| `setupMousePaths.m` | Adds linsync root, optional `../2019_03_03_BCT/` (BCT), and this folder to the MATLAB path |
| `loadMouseConnectome.m` | Loads `data/<mouseId>/fine_family_labelled_coarse.csv` and returns `K`, region labels, and a `info.isTrivial` flag for empty rows |
| `compareMouseHeatmap.m` | Sanity check: render `imagesc(K)` for every mouse so you can compare visually with `data/<mouseId>/coarse_connectome_<mouseid>.png` |
| `runMouseSection45.m` | Per-mouse driver -- mirrors `runEZ1_section45.m`. Healthy 1-D Epileptor on K, computes \(D(\to i)\), \(D(k \to)\), per-node critical excitability \(x^{c}_{0,i}\) |
| `runAllMiceSection45.m` | Loops over all 6 mice, ranks nodes within each mouse and aggregated across mice, saves CSV summary + cross-mouse overview |
| `compareAndersonVsArnold.m` | For each Anderson mouse, plots per-node \(D(\to i)\) and \(D(k \to)\) against the Arnold cohort mean ± SD; flags nodes with \|z\| ≥ threshold |
| `data/<mouseId>/` | Coarse connectome CSV, reference PNG, and the underlying compact connectome (unused here) |
| `results/` | Per-mouse `.mat` results, scatter figures, summary CSV, overview figure |

## Inputs

For each mouse, the relevant file is:

```
data/<mouseId>/fine_family_labelled_coarse.csv
```

This is a 67-line CSV: header row + 66 data rows. Column 1 is the row label (`L-CORTEX_VISUAL`, ..., `R-BACKGROUND`); the remaining 66 columns hold the symmetric weighted adjacency.

Several "rows" are placeholders that have no real connectivity (`L-BACKGROUND`, `R-BACKGROUND`, and a few `*_MASK` rows). `loadMouseConnectome` flags them in `info.isTrivial` and `runMouseSection45` skips them in the per-node sweep.

## Workflow per mouse

1. Load the 66×66 raw weighted matrix `K_raw`.
2. Normalise it. Schemes available via `applyConnectomeNormalisation` in the linsync root:
    - `tvb` (default) — `normal()`: zero diagonal, truncate at the 95th percentile of off-diagonal entries, rescale to `[0,1]`.
    - `none` — pass `K_raw` through unchanged.
    - `parkes` — `nctpy.utils.matrix_normalization` from Parkes et al. (Nat Protoc 2024, [doi:10.1038/s41596-024-01023-w](https://doi.org/10.1038/s41596-024-01023-w)): `A_norm = A / (|λ(A)|_max + c)`, with `c = ParkesC` (default 1). The continuous-time `-I` subtraction described in the paper is applied implicitly inside linsync's `con2cov` (which solves `dX = -X(I − C) dt + dW`), so we deliberately stop at the rescaling step. With this scheme the script **skips the Epileptor fixed-point solve and the per-node `x_0^c` sweep**: `A_norm` is fed straight in as the coupling matrix `C` and only `D(→i)` and `D(k→)` are computed (Liao & Lizier-style stability centralities on the Parkes-normalised connectome).
3. Solve the healthy-state 1-D Epileptor fixed point with `x0_i = -2.3` for all i.
4. Form the effective coupling matrix `C` (`CouplingMatrix.m` in linsync root).
5. Compute \(\Omega\) via `covariancesGaussianNet(C, false, MaxK, false, 1)` and read \(D(\to i) = \Omega_{ii}\). Repeat with \(C^{\top}\) for \(D(k \to)\).
6. For each non-trivial node, sweep \(x_{0,i}\) upward from −2.3 (others held at −2.3) until \(\rho(C) \geq 1\); coarse step 0.01 + bisection to 1e-4 in \(x_0\).
7. Render the §4.5 two-panel scatter and persist `_results.mat` / `_figure.{fig,png}`.

## Cross-mouse standout nodes

`runAllMiceSection45.m` collects the per-mouse vectors and produces:

- `results/section45_summary_topnodes_<scheme>.csv` — for every region: mean rank across mice (by \(x^{c}_{0,i}\), \(D(\to i)\), \(D(k \to)\)) plus the raw per-mouse values. Sorted by aggregate \(x^{c}_{0,i}\) rank (most epileptogenic first).
- `results/section45_summary_overview_<scheme>.png` — region × mouse heatmaps for \(x^{c}_{0,i}\), \(D(\to i)\), \(D(k \to)\). Useful for spotting nodes that are consistently low-\(x^c\) across animals (i.e. structurally susceptible regardless of mouse).

## Anderson vs Arnold per-node comparison

`compareAndersonVsArnold.m` consumes the per-mouse `_results.mat` files written by `runAllMiceSection45` and, for each Anderson mouse, produces a 2-panel figure (one panel per centrality) overlaying:

- shaded band — Arnold mean ± SD per node
- blue line — Arnold mean per node
- grey dots — individual Arnold values per node
- red markers — the Anderson mouse's per-node value
- black ring + label — nodes whose \|z\| = \|(Anderson − Arnold mean) / Arnold SD\| exceeds `ZThreshold` (default 2)

Outputs (under `results/`):

- `compare_anderson_vs_arnold_<mouseId>_<scheme>.{fig,png}` — per-Anderson figures
- `compare_anderson_vs_arnold_<mouseId>_<scheme>_outliers.csv` — list of flagged nodes with both centrality values, Arnold mean/SD, and z-scores, sorted by largest \|z\|
- `compare_anderson_vs_arnold_<scheme>.mat` — packed matrices for downstream analysis

Run after the per-mouse pipeline:

```matlab
compareAndersonVsArnold;                      % uses tvb scheme, z=2
compareAndersonVsArnold('ZThreshold', 1.5);   % stricter outlier net
```

## Quick start

```matlab
cd mouse-epilepsy-daria-anderson-experiment
setupMousePaths();

compareMouseHeatmap                 % visual sanity check vs reference PNGs
results = runMouseSection45('Anderson_1');   % single mouse
runAllMiceSection45                 % all 6 mice + cross-mouse summary
compareAndersonVsArnold             % Anderson-vs-Arnold per-node comparison

% Parkes normalisation (no Epileptor solve / no x0 sweep)
runMouseSection45('Anderson_1', 'Normalisation', 'parkes');
runAllMiceSection45;                                % edit `normalisation = 'parkes'` first
compareAndersonVsArnold('Normalisation', 'parkes');
```

Per-mouse runtime is similar to the human Epileptor case (~10 min / mouse); all six mice take ~1 hour with default settings.
