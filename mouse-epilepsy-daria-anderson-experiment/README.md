# Mouse epilepsy connectome experiment (Daria / Anderson cohort)

Apply the **stability centralities** workflow (per-node \(D(\to i)\), \(D(k \to)\), and critical \(x_0^c\) where applicable) to mouse coarse-grained connectomes (66 brain-atlas regions per animal). The human reference implementation is `epileptor-experiment/runEZ1_stabilityCentralities.m` (formerly thesis §4.5 / `runEZ1_section45`).

## Prerequisites

- **MATLAB Optimization Toolbox** — `fsolve` for Epileptor fixed points (`tvb` / `none` schemes).
- **Brain Connectivity Toolbox (BCT)** — optional for stability-only runs; **required** for betweenness/closeness in `computeNetworkCentralities` and for the betweenness panel in `compareAndersonVsArnold`. BCT is **not** included in git; install locally under `../2019_03_03_BCT/` — see [docs/BCT.md](../docs/BCT.md).

## Cohort

The reference cohort has five coarse connectomes:

`Anderson_1`, `Anderson_2`, `Arnold_3`, `Arnold_4`, `Arnold_5`

`Arnold_2` has other data files but **no** `fine_family_labelled_coarse.csv` and is skipped automatically. An extra folder such as `Arnold_1` is included if its CSV is present (`listAvailableMice` discovers all `data/<mouseId>/fine_family_labelled_coarse.csv` files).

**Git checkout:** this repository currently ships connectome CSVs for **`Arnold_5` only** (`data/Arnold_5/fine_family_labelled_coarse.csv`). Add the other four mice under `data/<mouseId>/` locally (same filename) for the full cohort; `listAvailableMice` discovers whatever is present. With a partial cohort, per-mouse stability and cohort \(D_{st}\) plots still run; Anderson-vs-Arnold comparison and related reports are skipped until at least one Anderson and one Arnold mouse are available.

---

## Managing experiment runs (recommended workflow)

Use a **`.properties` config file** per experiment. Each experiment has:

- **One normalisation scheme** (`column`, `parkes`, or `tvb`)
- **One parameter set** (Epileptor bounds, Parkes constant, comparison settings, etc.)
- **Its own results folder** under `results/<experiment.name>_<yyyy-mm-dd_HHMM>/` (config name + run timestamp), so runs never overwrite each other

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
| Per-mouse stability centralities | `runAllMiceStabilityCentralities` | `stabilityCentralities_<mouse>_<scheme>_results.mat`, `_figure.{fig,png}`, summary CSV/overview |
| Centrality correlations | `compareCentralityMeasures` | `centrality_corr_*` (only if cohort step succeeded) |
| Anderson vs Arnold | `compareAndersonVsArnold` | `compare_anderson_vs_arnold_*` |
| Left–right asymmetry | `compareLeftRightAsymmetry` | `compare_LR_asymmetry_*`, `LR_asymmetry_*` |
| Comparison summary figures (`column` only) | `renderComparisonFigures` | `comparison_figures/region_susc_outliers`, `laterality_influence_outliers` (`.fig`, `.png` each) |
| Network \(D_{\mathrm{st}}\) | `compareDstAcrossCohort` | `D_st_cohort_*` |
| Console reports | `reportTopNodes`, `reportAndersonOutliers` | `reports_<run>.log` |

After the per-mouse step, the orchestrator checks that **every** available mouse produced a `stabilityCentralities_*_<scheme>_results.mat` file. If any mouse fails, `runAllMiceStabilityCentralities` errors out (no partial summary), comparison/report steps are **skipped**, and `run_manifest.mat` records which step failed.

A complete single-scheme experiment has the same **types** of artefact as [`results/mouse_experiment_reference/`](results/mouse_experiment_reference/). That folder is a **flat, multi-scheme snapshot** (column + parkes + tvb together) saved with the legacy `section45_*` prefix; new orchestrated runs use **`stabilityCentralities_*`** for the same per-mouse and summary files.

At the end of `runMouseExperiment`, `assertMouseExperimentOutputs` checks that every required file for the enabled pipeline steps exists (see table below). Anderson outlier CSVs are optional (only written when Bonferroni outliers exist).

Disable steps with `pipeline.run*` keys in the properties file (e.g. `pipeline.runHeatmap=true` for connectome QC only).

### Expected outputs per scheme (checked by `assertMouseExperimentOutputs`)

For each of the five cohort mice, with `save.results=true` and all pipeline steps enabled:

| Category | Files (prefix `stabilityCentralities` unless noted) |
|----------|---------------------------------------------------|
| Per mouse (×5) | `<prefix>_<mouse>_<scheme>_results.mat`, `_figure.{fig,png}` |
| Cohort summary | `<prefix>_summary_topnodes_<scheme>.csv`, `_summary_overview_<scheme>.{fig,png}`, `_summary_<scheme>.mat` |
| Centrality comparison | `centrality_corr_<mouse>_<scheme>.{fig,png}` (×5), `centrality_corr_mean_<scheme>.{fig,png}`, `centrality_corr_bars_<scheme>.{fig,png}`, `centrality_corr_<scheme>.mat` |
| Anderson vs Arnold | Three fig/png per **Anderson** mouse (`Anderson_1`, `Anderson_2`): `*_D_to_i`, `*_D_k_to`, `*_BC` (BC when BCT available) vs Arnold cohort mean±SD; plus `compare_anderson_vs_arnold_<scheme>.mat`, optional `*_outliers.csv` per Anderson mouse |
| Left–right asymmetry | `compare_LR_asymmetry_<anderson>_<scheme>.{fig,png}`, optional `*_outliers.csv`; `LR_asymmetry_groupTest_<scheme>.{csv,fig,png}`; `LR_asymmetry_systematic_<scheme>.{csv,fig,png}`; `LR_asymmetry_<scheme>.mat` |
| Comparison summary figures (`column` only) | `comparison_figures/region_susc_outliers`, `laterality_influence_outliers` (`.fig`, `.png` at 1500×988 px) |
| Network \(D_{\mathrm{st}}\) | `D_st_cohort_<scheme>.{csv,fig,png,mat}` |
| QC (optional) | `mouse_heatmaps_overview_reference.{fig,png}` |
| Provenance | `experiment_<run>.properties`, `experiment_parameters_<run>.{mat,json}`, `run_manifest.mat`, `reports_<run>.log` (`<run>` = results subfolder name) |

---

## Results folder layout

Each experiment writes to **`results/<experiment.name>/`** (not the legacy flat `results/` root, unless you call scripts manually without an experiment name).

Example after `runMouseExperiment('configs/initial_column.properties')` (if run at 2026-05-23 14:30):

```
results/initial_column_2026-05-23_1430/
  experiment_initial_column_2026-05-23_1430.properties
  experiment_parameters_initial_column_2026-05-23_1430.mat
  experiment_parameters_initial_column_2026-05-23_1430.json
  run_manifest.mat               # pipeline steps + runParams + cfg
  reports_initial_column_2026-05-23_1430.log

  stabilityCentralities_Anderson_1_column_results.mat
  stabilityCentralities_Anderson_1_column_figure.{fig,png}
  ...                            # one result set per mouse (5 mice)
  stabilityCentralities_summary_topnodes_column.csv
  stabilityCentralities_summary_overview_column.{fig,png}
  stabilityCentralities_summary_column.mat

  compare_anderson_vs_arnold_Anderson_*_column_D_to_i.{fig,png}
  compare_anderson_vs_arnold_Anderson_*_column_D_k_to.{fig,png}
  compare_anderson_vs_arnold_Anderson_*_column_BC.{fig,png}
  compare_anderson_vs_arnold_Anderson_*_column_outliers.csv
  compare_anderson_vs_arnold_column.mat

  compare_LR_asymmetry_Anderson_*_column.{fig,png}
  compare_LR_asymmetry_Anderson_*_column_outliers.csv
  LR_asymmetry_groupTest_column.{csv,fig,png}
  LR_asymmetry_systematic_column.{csv,fig,png}
  LR_asymmetry_column.mat

  comparison_figures/            # column scheme only
    region_susc_outliers.{fig,png}
    laterality_influence_outliers.{fig,png}

  centrality_corr_<mouse>_column.{fig,png}
  centrality_corr_mean_column.{fig,png}
  centrality_corr_bars_column.{fig,png}
  centrality_corr_column.mat

  D_st_cohort_column.csv / .mat / .fig / .png
  mouse_heatmaps_overview_reference.{fig,png}   # if pipeline.runHeatmap=true
```

Legacy workflows that call `runAllMiceStabilityCentralities` without `ResultsDir` still write into the flat `results/` directory. Prefer the orchestrator for new work.

---

## Parameter provenance (what gets recorded)

Every orchestrated run records **all parameters** used to generate the results, at the experiment level and inside individual `.mat` artefacts.

### Experiment-level files

| File | Contents |
|------|----------|
| `experiment_<run>.properties` | Raw Java-style config (`key=value`) as run |
| `experiment_parameters_<run>.mat` | Struct `runParams` — see below |
| `experiment_parameters_<run>.json` | Same struct as JSON (easy diff/review in git or editors) |
| `run_manifest.mat` | `manifest` (step names, success, duration), `cfg` (parsed config), `runParams` (updated with `finishedAt` and `pipelineSteps`) |

`runParams` is built by `mouseExperimentRunParameters.m` and includes:

- **Metadata:** `experimentName`, `experimentDescription`, `configFile`, `recordedAt`, `resultsDir`
- **Environment:** MATLAB version, computer arch, hostname, user, toolkit paths
- **Cohort:** list of mouse IDs and count
- **`normalisation`:** scheme for this experiment
- **`stabilityCentralities`:** `parkesC`, `colScale`, `x0Base` / `x0Upper` / `x0Step`, `bisectTol`, `maxK`, `tau0`, `topK`, `discreteTime` (= `false`), `fsolve` tolerances, classical centrality defaults (`alphaPR`, `alphaKZ`, …)
- **`comparison`:** `corrType`, `alpha` (Bonferroni), `zThreshold`, centrality-plot options
- **`pipeline`:** which steps were enabled (`runStabilityCentralities`, `runCentralityCorr`, …)
- **`output`:** `saveResults`, `plot`, `verbose`
- **`config`:** full parsed struct from `loadMouseExperimentConfig`

Parameters are saved **at the start** of the run (`experiment_parameters_<run>.*`) and **again at the end** (with `finishedAt` and pipeline step outcomes).

### Embedded in each output `.mat`

So a single per-mouse file remains self-describing if copied elsewhere:

| File pattern | Field |
|--------------|--------|
| `stabilityCentralities_<mouse>_<scheme>_results.mat` | `results.runParameters` |
| `stabilityCentralities_summary_<scheme>.mat` | `runParameters` |
| `compare_anderson_vs_arnold_<scheme>.mat` | `runParameters` |
| `LR_asymmetry_<scheme>.mat` | `runParameters` |
| `centrality_corr_<scheme>.mat` | `runParameters` |
| `D_st_cohort_<scheme>.mat` | `runParameters` |

Manual runs (without the orchestrator) still attach `runParameters` built from each script’s `inputParser` options, but omit experiment metadata unless you pass `'RunParameters', runParams` yourself.

**Inspect parameters in MATLAB:**

```matlab
load('results/initial_column_2026-05-23_1430/experiment_parameters_initial_column_2026-05-23_1430.mat', 'runParams')
disp(runParams.stabilityCentralities)
disp(runParams.comparison)

% Or open experiment_parameters_<run>.json in any text editor
```

---

## Configuration file reference (`.properties`)

Files live in `configs/`. Syntax: `key=value`, `#` comments, one key per line.

### Experiment identity

| Key | Description |
|-----|-------------|
| `experiment.name` | Base name for the results folder (e.g. `initial_column`); at run time the folder becomes `initial_column_2026-05-23_1430` |
| `experiment.description` | Free-text note stored in `runParams` |

### Normalisation and stability centralities

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
| `laterality.basis` | `signed` | Basis for the per-Anderson L–R z-score and Bonferroni flagging in `compareLeftRightAsymmetry`: `signed` uses \(\mathrm{LI}_{\mathrm{signed}} = L - R\) (raw difference, matches the manuscript spec); `norm` uses \(\mathrm{LI}_{\mathrm{norm}} = (L-R)/(L+R)\). Group-level Welch and systematic sign-rank are always on \(\mathrm{LI}_{\mathrm{norm}}\). |

### Pipeline and output

| Key | Default | Description |
|-----|---------|-------------|
| `pipeline.runHeatmap` | `false` | Run `compareMouseHeatmap` (scheme-independent QC) |
| `pipeline.runStabilityCentralities` | `true` | Per-mouse + cohort summary (`pipeline.runSection45` is an alias) |
| `pipeline.runCentralityCorr` | `true` | Stability vs classical centrality correlations |
| `pipeline.runAndersonVsArnold` | `true` | Strain comparison figures and outlier CSVs |
| `pipeline.runLRAsymmetry` | `true` | Left–right laterality (Anderson vs Arnold per pair + group tests) |
| `pipeline.runDstCohort` | `true` | Network-level \(D_{\mathrm{st}}\) bar chart |
| `pipeline.runReports` | `true` | `reportTopNodes` + `reportAndersonOutliers` → `reports_<run>.log` |
| `pipeline.stopOnError` | `false` | If `true`, abort the experiment on first failed step |
| `save.results` | `true` | Write `.mat` / figures / CSVs |
| `plot` | `true` | Generate figures during stability-centralities step |
| `verbose` | `true` | Per-node console output in stability-centralities step |

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
| `runMouseStabilityCentralities.m` | Per-mouse stability centralities driver |
| `runAllMiceStabilityCentralities.m` | Loop cohort + cross-mouse summary (`stabilityCentralities_*` filenames) |
| `mouseExperimentResultPrefix.m` | Returns `stabilityCentralities` (result artefact prefix) |
| `mouseExperimentFolderName.m` | Build `configName_yyyy-mm-dd_HHMM` results subfolder |
| `mouseExperimentProvenanceFile.m` | Suffixed provenance paths (`experiment_<run>.properties`, etc.) |
| `listPerMouseResultFiles.m` | Find per-mouse `*_results.mat` (stabilityCentralities or legacy section45) |
| `assertMouseExperimentOutputs.m` | Post-run checklist vs reference artefact types |
| `compareMouseHeatmap.m` | Visual QC vs reference PNGs |
| `compareCentralityMeasures.m` | Correlation heatmaps / bar charts |
| `compareAndersonVsArnold.m` | Anderson vs Arnold per-node plots + outliers |
| `compareLeftRightAsymmetry.m` | L–R laterality indices; Anderson vs Arnold per pair + group/systematic views |
| `renderComparisonFigures.m` | Cohort comparison summary figures; uses in-memory comparison summaries only |
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

## Workflow per mouse (stability centralities)

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

`runAllMiceStabilityCentralities` produces (under the experiment results folder):

- `stabilityCentralities_summary_topnodes_<scheme>.csv` — mean ranks and per-mouse values. Includes an `is_trivial` column for placeholder rows (`L-/R-BACKGROUND`, `*_MASK`, `*_MASKS`); `reportTopNodes` filters these out by default (pass `'IncludeTrivial', true` to keep them).
- `stabilityCentralities_summary_overview_<scheme>.{fig,png}` — region × mouse heatmaps (1–3 panels depending on scheme)
- `stabilityCentralities_summary_<scheme>.mat` — packed matrices + `runParameters`

## Anderson vs Arnold comparison

`compareAndersonVsArnold` loads per-mouse `stabilityCentralities_*_results.mat` and, for each Anderson mouse, writes **three separate figures** (one per metric): \(D(\to i)\), \(D(k \to)\), and betweenness centrality vs the Arnold cohort mean ± SD.

## Left–right hemisphere asymmetry

`compareLeftRightAsymmetry` pairs each `L-<region>` node with `R-<region>` (same suffix) and computes, per mouse and metric:

- **\( \mathrm{LI}_{\mathrm{norm}} = (L - R) / (L + R) \)** — normalised laterality index in \([-1, +1]\) (positive ⇒ higher on the L-labelled side)
- **\( \mathrm{LI}_{\mathrm{signed}} = L - R \)** — signed difference in the metric's own units

Metrics: \(D(\to i)\), \(D(k \to)\), and BC (when BCT was used). Three comparison views:

1. **Per Anderson mouse** — z-score of \(\mathrm{LI}_{\langle\mathrm{basis}\rangle}\) vs Arnold mean ± SD per region pair; Bonferroni over (pairs × metrics); figures `compare_LR_asymmetry_<mouse>_<scheme>.*`. The basis is set by `laterality.basis` in the config (default `signed`, i.e. the raw difference \(L - R\); set `laterality.basis = norm` to revert to the bounded index). Outlier-CSV columns adapt to the basis: `z_<metric>` / `p_bonf_<metric>` for `norm`, `z_signed_<metric>` / `p_bonf_signed_<metric>` for `signed`; both `LI_norm_*` and `LI_signed_*` raw values are always emitted.
2. **Group-level** — Welch t-test of \(\mathrm{LI}_{\mathrm{norm}}\) per pair × metric (basis-independent); `LR_asymmetry_groupTest_<scheme>.*`
3. **Systematic direction** — cohort-wide mean \(\mathrm{LI}_{\mathrm{norm}}\) and sign-rank vs zero (basis-independent; answers whether epileptic mice show a consistent L vs R bias); `LR_asymmetry_systematic_<scheme>.*`

Full packed results: `LR_asymmetry_<scheme>.mat`. Disable with `pipeline.runLRAsymmetry=false` in the config.

## Comparison summary figures (`column` only)

After `compareAndersonVsArnold` and `compareLeftRightAsymmetry`, `runMouseExperiment` calls `renderComparisonFigures` for the **`column`** scheme. Two summary figures are written under `comparison_figures/` (`.fig` plus 1500×988 px PNG at 200 dpi):

| File | Content |
|------|---------|
| `region_susc_outliers` | Bonferroni-significant \(D(\to i)\) excess over Arnold mean |
| `laterality_influence_outliers` | Significant \(D(k\to)\) L–R pairs (\(\mathrm{LI}_{\mathrm{signed}}\)) |

Data come from the same in-memory quantities as the comparison outlier CSVs (no `readtable` of those files). Titles, legends, and cohort counts are derived from the summary struct (`cohortGroupInfo`), not hard-coded strain names.

For a different cohort, set mouse ID prefixes in the properties file (defaults match the shipped Anderson / Arnold study):

```properties
cohort.case.prefix=Anderson
cohort.control.prefix=Arnold
```

Example for `Case_1`, `Case_2` vs `Ctrl_1`…`Ctrl_4`: use `cohort.case.prefix=Case` and `cohort.control.prefix=Ctrl`.

To regenerate manually:

```matlab
ava = compareAndersonVsArnold('Normalisation', 'column', 'ResultsDir', resultsDir, ...
    'CasePrefix', 'Anderson', 'ControlPrefix', 'Arnold');
lr  = compareLeftRightAsymmetry('Normalisation', 'column', 'ResultsDir', resultsDir, ...
    'CasePrefix', 'Anderson', 'ControlPrefix', 'Arnold');
renderComparisonFigures(ava, lr, 'OutputDir', fullfile(resultsDir, 'comparison_figures'));
```

Outlier **flagging** uses **Bonferroni correction** at `alpha` from the config (default `0.05`), not a simple \|z\| > `z.threshold`. The `z.threshold` key controls what `reportAndersonOutliers` **displays** in `reports_<run>.log`.

Outputs:

- `compare_anderson_vs_arnold_<mouseId>_<scheme>_D_to_i.{fig,png}`
- `compare_anderson_vs_arnold_<mouseId>_<scheme>_D_k_to.{fig,png}`
- `compare_anderson_vs_arnold_<mouseId>_<scheme>_BC.{fig,png}` (when BC is available)
- `compare_anderson_vs_arnold_<mouseId>_<scheme>_outliers.csv`
- `compare_anderson_vs_arnold_<scheme>.mat` (includes `runParameters`)

---

## Manual / legacy workflow

You can still call scripts directly. Outputs go to flat `results/` unless you pass an explicit directory:

```matlab
setupMousePaths();
resultsDir = setupMousePaths('ExperimentName', 'my_manual_run');

compareMouseHeatmap('ResultsDir', resultsDir);
runAllMiceStabilityCentralities('Normalisation', 'column', 'ResultsDir', resultsDir);
ava = compareAndersonVsArnold('Normalisation', 'column', 'ResultsDir', resultsDir);
lr  = compareLeftRightAsymmetry('Normalisation', 'column', 'ResultsDir', resultsDir);
renderComparisonFigures(ava, lr, 'OutputDir', fullfile(resultsDir, 'comparison_figures'));
```

Per-mouse stability-centralities files will contain `results.runParameters` built from that script’s options only (no full experiment config unless you pass `'RunParameters', ...`).

**Runtime:** ~10 min per mouse for `tvb`; ~1 hour for all five mice at default settings.
