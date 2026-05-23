# Mouse epilepsy connectome experiment (Daria / Anderson cohort)

Apply the §4.5 workflow from Liao's thesis (originally implemented for human data in `epileptor-experiment/runEZ1_section45.m`) to mouse coarse-grained connectomes (66 brain-atlas regions per animal).

## Prerequisites

- **MATLAB Optimization Toolbox** — `fsolve` for Epileptor fixed points (`tvb` / `none` schemes).
- **Brain Connectivity Toolbox (BCT)** — optional for stability-only runs; **required** for betweenness/closeness in `computeNetworkCentralities` and for the betweenness panel in `compareAndersonVsArnold`. BCT is **not** included in git; install locally under `../2019_03_03_BCT/` — see [docs/BCT.md](../docs/BCT.md).

## Cohort

Five mice have the coarse CSV required by the pipeline:

`Anderson_1`, `Anderson_2`, `Arnold_3`, `Arnold_4`, `Arnold_5`

`Arnold_2` has other data files but **no** `fine_family_labelled_coarse.csv` and is skipped automatically.

---

## Managing experiment runs (recommended workflow)

Use a **`.properties` config file** per experiment. Each experiment has:

- **One normalisation scheme** (`column`, `parkes`, or `tvb`)
- **One parameter set** (Epileptor bounds, Parkes constant, comparison settings, etc.)
- **Its own results folder** under `results/<experiment.name>/`, so runs never overwrite each other

### Quick start

```matlab
cd mouse-epilepsy-daria-anderson-experiment
setupMousePaths();

% Run the three shipped baseline experiments (column, parkes, tvb)
% Warning: tvb is slow (~10 min/mouse × 5 mice)
runAllMouseExperiments

% Or a single experiment:
runMouseExperiment('configs/initial_column.properties')

% Create a new timestamped config, edit it, then run:
newMouseExperiment('parkes')   % -> configs/yyyy-MM-dd_HH-mm-ss_parkes.properties
runMouseExperiment('configs/yyyy-MM-dd_HH-mm-ss_parkes.properties')
```

### Shipped configs

| Config file | `experiment.name` | Normalisation |
|-------------|-------------------|---------------|
| `configs/initial_column.properties` | `initial_column` | `column` |
| `configs/initial_parkes.properties` | `initial_parkes` | `parkes` |
| `configs/initial_tvb.properties` | `initial_tvb` | `tvb` |

Copy `configs/experiment.template.properties` or any `initial_*.properties` file to start a new study with different `parkes.c`, `col.scale`, Epileptor search bounds, or pipeline toggles.

### What `runMouseExperiment` runs

When pipeline flags are `true` (defaults), the orchestrator runs these steps **in order** for the scheme in the config:

| Step | Script | Output (examples) |
|------|--------|-------------------|
| Optional QC | `compareMouseHeatmap` | `mouse_heatmaps_overview_*.{fig,png}` |
| Per-mouse §4.5 | `runAllMiceSection45` | `section45_<mouse>_<scheme>_results.mat`, figures, summary CSV/overview |
| Centrality correlations | `compareCentralityMeasures` | `centrality_corr_*` |
| Anderson vs Arnold | `compareAndersonVsArnold` | `compare_anderson_vs_arnold_*` |
| Network \(D_{\mathrm{st}}\) | `compareDstAcrossCohort` | `D_st_cohort_*` |
| Console reports | `reportTopNodes`, `reportAndersonOutliers` | `reports.log` |

After §4.5, the orchestrator checks that every available mouse produced a `section45_*_<scheme>_results.mat` file. Step success and timing are recorded in `run_manifest.mat`.

Disable steps with `pipeline.run*` keys in the properties file (e.g. `pipeline.runHeatmap=true` for connectome QC only).

---

## Results folder layout

Each experiment writes to **`results/<experiment.name>/`** (not the legacy flat `results/` root, unless you call scripts manually without an experiment name).

Example after `runMouseExperiment('configs/initial_column.properties')`:

```
results/initial_column/
  experiment.properties          # exact config copied at run start
  experiment_parameters.mat      # full parameter snapshot (MATLAB)
  experiment_parameters.json     # same snapshot (human-readable)
  run_manifest.mat               # pipeline steps + runParams + cfg
  reports.log                    # reportTopNodes / reportAndersonOutliers output

  section45_Anderson_1_column_results.mat
  section45_Anderson_1_column_figure.{fig,png}
  ...                            # one result set per mouse
  section45_summary_topnodes_column.csv
  section45_summary_overview_column.{fig,png}
  section45_summary_column.mat

  centrality_corr_*.mat / .fig / .png
  compare_anderson_vs_arnold_*
  D_st_cohort_column.csv / .mat / .fig / .png
```

Legacy workflows that call `runAllMiceSection45` without `ResultsDir` still write into the flat `results/` directory. Prefer the orchestrator for new work.

---

## Parameter provenance (what gets recorded)

Every orchestrated run records **all parameters** used to generate the results, at the experiment level and inside individual `.mat` artefacts.

### Experiment-level files

| File | Contents |
|------|----------|
| `experiment.properties` | Raw Java-style config (`key=value`) as run |
| `experiment_parameters.mat` | Struct `runParams` — see below |
| `experiment_parameters.json` | Same struct as JSON (easy diff/review in git or editors) |
| `run_manifest.mat` | `manifest` (step names, success, duration), `cfg` (parsed config), `runParams` (updated with `finishedAt` and `pipelineSteps`) |

`runParams` is built by `mouseExperimentRunParameters.m` and includes:

- **Metadata:** `experimentName`, `experimentDescription`, `configFile`, `recordedAt`, `resultsDir`
- **Environment:** MATLAB version, computer arch, hostname, user, toolkit paths
- **Cohort:** list of mouse IDs and count
- **`normalisation`:** scheme for this experiment
- **`section45`:** `parkesC`, `colScale`, `x0Base` / `x0Upper` / `x0Step`, `bisectTol`, `maxK`, `tau0`, `topK`, `discreteTime` (= `false`), `fsolve` tolerances, classical centrality defaults (`alphaPR`, `alphaKZ`, …)
- **`comparison`:** `corrType`, `alpha` (Bonferroni), `zThreshold`, centrality-plot options
- **`pipeline`:** which steps were enabled (`runSection45`, `runCentralityCorr`, …)
- **`output`:** `saveResults`, `plot`, `verbose`
- **`config`:** full parsed struct from `loadMouseExperimentConfig`

Parameters are saved **at the start** of the run (`experiment_parameters.*`) and **again at the end** (with `finishedAt` and pipeline step outcomes).

### Embedded in each output `.mat`

So a single per-mouse file remains self-describing if copied elsewhere:

| File pattern | Field |
|--------------|--------|
| `section45_<mouse>_<scheme>_results.mat` | `results.runParameters` |
| `section45_summary_<scheme>.mat` | `runParameters` |
| `compare_anderson_vs_arnold_<scheme>.mat` | `runParameters` |
| `centrality_corr_<scheme>.mat` | `runParameters` |
| `D_st_cohort_<scheme>.mat` | `runParameters` |

Manual runs (without the orchestrator) still attach `runParameters` built from each script’s `inputParser` options, but omit experiment metadata unless you pass `'RunParameters', runParams` yourself.

**Inspect parameters in MATLAB:**

```matlab
load('results/initial_column/experiment_parameters.mat', 'runParams')
disp(runParams.section45)
disp(runParams.comparison)

% Or open experiment_parameters.json in any text editor
```

---

## Configuration file reference (`.properties`)

Files live in `configs/`. Syntax: `key=value`, `#` comments, one key per line.

### Experiment identity

| Key | Description |
|-----|-------------|
| `experiment.name` | Folder name under `results/` (e.g. `initial_column` or `2026-05-23_14-30-00_column`) |
| `experiment.description` | Free-text note stored in `runParams` |

### Normalisation and §4.5

| Key | Default | Description |
|-----|---------|-------------|
| `normalisation` | `column` | `column`, `parkes`, `tvb`, or `none` |
| `parkes.c` | `1.0` | Parkes rescaling constant \(c\) |
| `col.scale` | `0.95` | Target column sum for column normalisation |
| `x0.base` | `-2.3` | Healthy excitability (Epileptor schemes) |
| `x0.upper` | `-1.0` | Upper bound of per-node \(x_0^c\) search |
| `x0.step` | `0.01` | Coarse sweep step |
| `bisect.tol` | `1e-4` | Bisection tolerance in \(x_0\) |
| `max.k` | `1e8` | `covariancesGaussianNet` iteration cap |
| `tau0` | `6667` | 1-D Epileptor slow timescale |
| `top.k` | `10` | Top nodes in summary CSV / reports |

### Comparisons and reports

| Key | Default | Description |
|-----|---------|-------------|
| `corr.type` | `Spearman` | Correlation type in `compareCentralityMeasures` |
| `alpha` | `0.05` | Family-wise level for Bonferroni outlier flagging in `compareAndersonVsArnold` |
| `z.threshold` | `2` | Minimum \|z\| shown in `reportAndersonOutliers` |

### Pipeline and output

| Key | Default | Description |
|-----|---------|-------------|
| `pipeline.runHeatmap` | `false` | Run `compareMouseHeatmap` (scheme-independent QC) |
| `pipeline.runSection45` | `true` | Per-mouse + cohort summary |
| `pipeline.runCentralityCorr` | `true` | Stability vs classical centrality correlations |
| `pipeline.runAndersonVsArnold` | `true` | Strain comparison figures and outlier CSVs |
| `pipeline.runDstCohort` | `true` | Network-level \(D_{\mathrm{st}}\) bar chart |
| `pipeline.runReports` | `true` | `reportTopNodes` + `reportAndersonOutliers` → `reports.log` |
| `pipeline.stopOnError` | `false` | If `true`, abort the experiment on first failed step |
| `save.results` | `true` | Write `.mat` / figures / CSVs |
| `plot` | `true` | Generate figures during §4.5 |
| `verbose` | `true` | Per-node console output in §4.5 |

---

## Script layout

| Path | Role |
|------|------|
| `setupMousePaths.m` | Add linsync root, BCT, and this folder to the path; optional `ExperimentName` → `results/<name>/` |
| `resolveMouseResultsDir.m` | Resolve `ResultsDir` name-value for all writers |
| `loadMouseExperimentConfig.m` | Parse `.properties` → MATLAB struct |
| `mouseExperimentRunParameters.m` | Build / merge / save `runParams` provenance |
| `newMouseExperiment.m` | Create timestamped config from template |
| `runMouseExperiment.m` | **Orchestrator** — full pipeline from one config file |
| `runAllMouseExperiments.m` | Run every `configs/*.properties` except the template |
| `loadMouseConnectome.m` | Load `data/<mouseId>/fine_family_labelled_coarse.csv` |
| `runMouseSection45.m` | Per-mouse §4.5 driver |
| `runAllMiceSection45.m` | Loop cohort + cross-mouse summary |
| `compareMouseHeatmap.m` | Visual QC vs reference PNGs |
| `compareCentralityMeasures.m` | Correlation heatmaps / bar charts |
| `compareAndersonVsArnold.m` | Anderson vs Arnold per-node plots + outliers |
| `compareDstAcrossCohort.m` | Cohort \(D_{\mathrm{st}}\) table and figure |
| `reportTopNodes.m` | Print top regions from summary CSV |
| `reportAndersonOutliers.m` | Print outlier tables from comparison CSVs |
| `configs/` | Experiment configs |
| `data/<mouseId>/` | Connectome CSV and reference images |
| `results/<experimentName>/` | Outputs for one experiment run |

---

## Inputs

For each mouse, the relevant file is:

```
data/<mouseId>/fine_family_labelled_coarse.csv
```

This is a 67-line CSV: header row + 66 data rows. Column 1 is the row label (`L-CORTEX_VISUAL`, ..., `R-BACKGROUND`); the remaining 66 columns hold the symmetric weighted adjacency.

Placeholder rows (`L-BACKGROUND`, `R-BACKGROUND`, `*_MASK`) are flagged in `info.isTrivial` and skipped in the per-node \(x_0^c\) sweep.

## Workflow per mouse (§4.5)

1. Load the 66×66 raw weighted matrix `K_raw`.
2. Normalise via `applyConnectomeNormalisation` (see scheme table below).
3. **Epileptor schemes (`tvb`, `none`):** solve healthy fixed point, form `C`, compute \(\Omega\), sweep \(x_{0,i}^c\) per non-trivial node.
4. **Linear schemes (`parkes`, `column`):** use normalised `K` as `C` directly (no Epileptor solve / no \(x_0^c\) sweep).
5. Compute \(D(\to i)\), \(D(k \to)\) from `covariancesGaussianNet` (continuous-time: `discreteTime = false`).
6. Save results, figures, and `runParameters`.

### Normalisation schemes

| Scheme | Epileptor + \(x_0^c\) sweep? | Notes |
|--------|------------------------------|-------|
| `tvb` | Yes | TVB-style 95th-percentile truncate + rescale to `[0,1]` |
| `none` | Yes | Raw `K` (no rescaling) |
| `parkes` | No | \(A / (\|\lambda\|_{\max} + c)\); symmetric \(K\) ⇒ \(D(k\to) \equiv D(\to i)\) |
| `column` | No | Column-sum scaling; \(D(k\to) \neq D(\to i)\) in general |

## Cross-mouse summary

`runAllMiceSection45` produces (under the experiment results folder):

- `section45_summary_topnodes_<scheme>.csv` — mean ranks and per-mouse values
- `section45_summary_overview_<scheme>.{fig,png}` — region × mouse heatmaps (1–3 panels depending on scheme)
- `section45_summary_<scheme>.mat` — packed matrices + `runParameters`

## Anderson vs Arnold comparison

`compareAndersonVsArnold` loads per-mouse `section45_*_results.mat` and, for each Anderson mouse, plots **three panels**: \(D(\to i)\), \(D(k \to)\), and betweenness centrality vs the Arnold cohort mean ± SD.

Outlier **flagging** uses **Bonferroni correction** at `alpha` from the config (default `0.05`), not a simple \|z\| > `z.threshold`. The `z.threshold` key controls what `reportAndersonOutliers` **displays** in `reports.log`.

Outputs:

- `compare_anderson_vs_arnold_<mouseId>_<scheme>.{fig,png}`
- `compare_anderson_vs_arnold_<mouseId>_<scheme>_outliers.csv`
- `compare_anderson_vs_arnold_<scheme>.mat` (includes `runParameters`)

---

## Manual / legacy workflow

You can still call scripts directly. Outputs go to flat `results/` unless you pass an explicit directory:

```matlab
setupMousePaths();
resultsDir = setupMousePaths('ExperimentName', 'my_manual_run');

compareMouseHeatmap('ResultsDir', resultsDir);
runAllMiceSection45('Normalisation', 'column', 'ResultsDir', resultsDir);
compareAndersonVsArnold('Normalisation', 'column', 'ResultsDir', resultsDir);
```

Per-mouse §4.5 files will contain `results.runParameters` built from that script’s options only (no full experiment config unless you pass `'RunParameters', ...`).

**Runtime:** ~10 min per mouse for `tvb`; ~1 hour for all five mice at default settings.
