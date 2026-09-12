# Connectome stability pipeline — runnable demo

A complete, self-contained example of the connectome cohort pipeline, on a
**synthetic** cohort. No real patient or animal data is needed, and the whole
run takes a few seconds.

Use it to check your installation, to see the shape of the outputs before
pointing the pipeline at your own data, or as the template for a new study.

## Run it

```matlab
cd connectome-demo
setupConnectomePaths();          % puts the linsync root on the path
makeSyntheticCohort('data');     % writes data/Ctrl_1 ... data/Case_2
runConnectomeExperiment('configs/demo.properties');
```

That is the whole demo. Outputs land in
`results/demo_column_<yyyy-mm-dd_HHMM>/`.

## What the synthetic cohort contains

`makeSyntheticCohort` builds six subjects on a 26-node bilateral atlas
(12 regions × 2 hemispheres, plus one disconnected `BACKGROUND` placeholder
per hemisphere so the trivial-node handling is exercised):

| Group | Subjects | Connectome |
|-------|----------|------------|
| Control | `Ctrl_1` … `Ctrl_4` | dense within-hemisphere edges, sparse cross-hemisphere edges, strong homotopic links |
| Case | `Case_1`, `Case_2` | the same, with **`L-HIPPOCAMPUS`** and **`L-CORTEX_TEMPORAL`** up-weighted by 1.8× |

The planted effect is left-sided and confined to two regions, so a working
pipeline should recover it twice over: as elevated `D(→i)` for those nodes in
the case-vs-control comparison, and as a left–right asymmetry in the
laterality analysis. If you change `EffectSize` to `1.0`, both should stop
finding anything — a useful null check.

Everything is seeded (`'Seed', 42`), so the cohort is identical on every
machine.

## What the run produces

```
results/demo_column_<stamp>/
  provenance/      config snapshot, parameters (.mat + .json), run_manifest, reports log
  qc/              connectome heatmaps for every subject
  stability_centralities/
                   per-subject results + figures, cohort summary CSV, overview heatmap
  centrality_correlation/
                   stability vs classical centralities, per subject and averaged
  cohort_comparison/
                   per-node Case-vs-Ctrl figures and Bonferroni outlier CSVs
  hemispheric_asymmetry/
                   per-subject laterality figures, group test, systematic bias
  dst_cohort/      network-level D_st per subject, with a Welch test
  comparison_figures/
                   the two publication-style summary figures
```

Start with `provenance/reports_<run>.log` — it holds the top-ranked regions
and the flagged outliers in plain text.

## Adapting it to your own study

1. Copy this folder's layout: a project root holding `configs/`, `data/` and
   `results/`.
2. Put one folder per subject under `data/`, each with a connectome file
   (labelled CSV, or a plain numeric matrix). Name the folders so the group
   is in the prefix, e.g. `Patient_01` and `Control_01`.
3. Create a config from the annotated template and edit it:

   ```matlab
   newConnectomeExperiment('Name', 'pilot', 'Normalisation', 'column');
   ```

   The keys you are most likely to change are `data.file`,
   `cohort.case.prefix` / `cohort.control.prefix`,
   `hemisphere.left.prefix` / `hemisphere.right.prefix` and
   `trivial.label.suffixes`.
4. Run it:

   ```matlab
   runConnectomeExperiment('configs/pilot.properties');
   ```

For the mathematics behind each step, see
[`docs/connectome-stability.md`](../docs/connectome-stability.md).
