# Connectome stability analysis — methods reference

This document explains what the connectome cohort pipeline in `linsync`
computes and why each choice was made. It is written to stand on its own: you
should be able to read it without opening the source.

For how to *run* the pipeline, see [`connectome-demo/README.md`](../connectome-demo/README.md)
and the annotated [`configs/experiment.template.properties`](../configs/experiment.template.properties).

---

## Contents

1. [What the pipeline answers](#1-what-the-pipeline-answers)
2. [Input data](#2-input-data)
3. [Normalisation](#3-normalisation)
4. [From a connectome to a coupling matrix](#4-from-a-connectome-to-a-coupling-matrix)
5. [Stability centralities and $D_{\mathrm{st}}$](#5-stability-centralities-and-d_st)
6. [Critical excitability $x^{c}_{0,i}$](#6-critical-excitability-xc_0i)
7. [Classical centralities for comparison](#7-classical-centralities-for-comparison)
8. [Case-vs-control statistics](#8-case-vs-control-statistics)
9. [Hemispheric asymmetry](#9-hemispheric-asymmetry)
10. [Network-level comparison](#10-network-level-comparison)
11. [Provenance](#11-provenance)
12. [Limitations](#12-limitations)
13. [Function index](#13-function-index)

---

## 1. What the pipeline answers

Given a cohort of weighted structural connectomes — one per subject, split
into a *case* group and a *control* group — the pipeline asks:

* **Per node**: which regions are dynamically fragile, in the sense that
  network-wide fluctuations concentrate on them (or radiate from them)?
* **Per subject pair of hemispheres**: is that fragility lateralised?
* **Per network**: how far from stability is the whole system?
* **Between groups**: which of those quantities differ between cases and
  controls by more than the control cohort's own spread?

Everything is derived from the connectome plus a dynamical model — no
functional recordings are required.

---

## 2. Input data

One folder per subject:

```
<data root>/<subjectId>/<connectome file>
```

The connectome file is either a **labelled CSV** (header row of node names,
first column repeating those names as row labels, then an $N \times N$ block
of weights) or a **plain numeric matrix** with no labels, in which case nodes
are named `node_001 …`. Both are read by `loadConnectome`.

Group membership comes from the subject ID prefix (`cohort.case.prefix`,
`cohort.control.prefix`), so a study needs no separate group table and a
cohort can grow by dropping in a folder.

### Trivial nodes

A node is marked **trivial** when either

* its off-diagonal row and column are exactly zero (it is disconnected, so
  per-node measures are meaningless), or
* its label ends with one of `trivial.label.suffixes` — atlas placeholder
  rows such as background or mask entries.

The label rule exists because placeholder rows often carry tiny non-zero
weights. After a normalisation that rescales by the maximum, those tiny
weights can become comparable to real ones and the placeholder marches to the
top of the ranking with no anatomical meaning behind it.

Trivial nodes are **never deleted** — indices stay aligned with the source
file — but they are excluded from rankings, from the excitability sweep, and
from the multiple-comparison families.

---

## 3. Normalisation

A raw connectome carries arbitrary units (streamline counts, fibre densities).
Before it can act as the coupling matrix of a stable linear system it has to
be rescaled. `applyConnectomeNormalisation` offers four schemes, and
`normalisationSchemeInfo` reports what each one implies downstream.

| Scheme | Formula | Epileptor solve | Symmetric? |
|--------|---------|-----------------|------------|
| `column` | each column rescaled to sum to `col.scale` | no | no |
| `parkes` | $A / (\lvert\lambda\rvert_{\max} + c)$ | no | yes |
| `tvb` | truncate at the `tvb.percentile` percentile of the upper triangle, zero the diagonal, rescale to $[0,1]$ | yes | yes |
| `none` | raw weights | yes | yes |

### The symmetry trade-off

Structural connectomes are usually symmetric. A symmetry-preserving
normalisation therefore produces a symmetric coupling matrix, and then
$D(k \to) \equiv D(\to i)$ identically — the influence centrality carries no
information the susceptibility centrality did not already have.

**Column normalisation breaks that degeneracy.** Dividing each column by its
own sum means a node with many weak inputs and a node with few strong inputs
end up scaled differently, so the matrix becomes asymmetric and the two
centralities separate. The interpretation is that each node has a fixed total
input gain, which it distributes over its afferents — a standard assumption in
network-control work. The cost is that the weights no longer mean "fibre
count" in any direct sense.

Column normalisation also guarantees $\rho(C) \le$ `col.scale`, because for a
non-negative matrix the spectral radius is bounded by the maximum column sum.
That is why `col.scale = 0.95` is the default: stability is structural, not
something to be tuned per subject.

### The Parkes constant

`parkes.c` shifts the whole spectrum: $A/(\rho(A)+c)$ has spectral radius
$\rho(A)/(\rho(A)+c) < 1$ for any $c > 0$. It is a stability margin, not a
free scientific parameter — changing it rescales every node's $D(\to i)$ in
much the same way and barely moves the *ranking*, which is what the
comparisons use.

---

## 4. From a connectome to a coupling matrix

### Linear schemes (`column`, `parkes`)

The normalised connectome $K$ is used directly as the coupling matrix
$C = K$. No neural model is assumed; the weights are read as linear gains.
The run errors out if $\rho(C) \ge 1$, since the stationary covariance would
not exist.

### Epileptor schemes (`tvb`, `none`)

Each region $i$ carries a single slow Epileptor variable $z_i(t)$, coupled
through $K$:

$$\tau_0 \dot{z}_i = G(z_i) - 4 x_{0,i} - \sum_{j} K_{ij}\left[F(z_j) - F(z_i)\right],$$

with

$$F(z) = \tfrac{1}{4}\left(-\tfrac{16}{3} - \sqrt{8z - \tfrac{629.6}{27}}\right), \qquad G(z) = 4F(z) - z .$$

This is the standard one-dimensional reduction of the Epileptor used for
seizure-onset analyses; the square-root term is what produces its
characteristic bifurcation as the excitability $x_0$ rises.

`healthyEpileptorCoupling` holds every node at the healthy excitability
$x_{0,i} = x_0^{\mathrm{base}}$ (default $-2.3$), solves $\dot z = 0$ with
`fsolve` (Optimization Toolbox, tolerances $10^{-14}$, warm start
$z = 3\cdot\mathbf{1}$), and linearises around the solution $z^{*}$. The
resulting effective coupling matrix is

$$C_{ii} = 1 - \tfrac{1}{\tau_0} + \tfrac{1}{\tau_0}\left(4 + \sum_j K_{ij} - K_{ii}\right) F'(z_i^{*}),$$

$$C_{ij} = -\tfrac{1}{\tau_0} K_{ji} F'(z_j^{*}) \quad (i \neq j), \qquad F'(z) = \frac{-1}{\sqrt{8z - 629.6/27}} .$$

If the all-healthy network is already unstable ($\rho(C) \ge 1$), the run
stops: everything downstream would be undefined.

Note the index transpose in the off-diagonal term. `linsync`'s convention is
that $C_{ij}$ is the weight of the edge **from** $i$ **to** $j$, i.e. matrices
act on row vectors.

### The dynamics both routes feed

Whatever produced $C$, the fluctuations $X(t)$ around the operating point are
modelled as the continuous-time Ornstein–Uhlenbeck process

$$\mathrm{d}X(t) = -X(t)(I - C)\theta\,\mathrm{d}t + \zeta\,\mathrm{d}w(t),$$

with $\theta = \zeta = 1$ (they enter only as constant multipliers) and $w$ a
Wiener process of identity covariance. This is the same model the rest of the
toolkit uses. The pipeline always runs continuous time
(`discreteTime = false`).

---

## 5. Stability centralities and $D_{\mathrm{st}}$

The stationary covariance $\Omega$ of that process is the unique
positive-definite solution of the Lyapunov equation, computed by
`covariancesGaussianNet` via the power series

$$\Omega = \frac{\zeta^{2}}{2\theta} \sum_{m=0}^{\infty} 2^{-m} \sum_{u=0}^{m} \binom{m}{u} (C^{u})^{\top} C^{\,m-u},$$

with a closed form in the symmetric case. Convergence requires $\rho(C) < 1$.

`computeStabilityCentralities(C, MaxK)` then returns

$$D(\to i) = \Omega_{ii}, \qquad D(k \to) = \left(\Omega_{C^{\top}}\right)_{kk}, \qquad D_{\mathrm{st}} = \frac{1}{N}\,\mathrm{trace}\,\Omega .$$

**Interpretation.**

* $D(\to i)$ — **stability susceptibility centrality**. The steady-state
  variance the network deposits *into* node $i$ along all upstream paths.
  Large means: fluctuations anywhere in the network arrive here amplified.
* $D(k \to)$ — **stability influence centrality**, obtained from the
  covariance of the *transposed* coupling. Large means: fluctuations starting
  here spread far.
* $D_{\mathrm{st}}$ — the scalar **deviation from stability** of the whole
  network: the per-node average variance amplification.

$D_{\mathrm{st}}$ is worth reporting alongside $\rho(C)$ because two networks
pinned at the same spectral radius can still amplify very differently — the
leading eigenvalue is only one mode, and the trace of $\Omega$ integrates over
all of them.

These are centralities in the same generalised sense as PageRank or
eigenvector centrality, but explicitly **dynamical**: they encode what the
linearised dynamics on the graph does, not just who is connected to whom.

The convergence error code of each power series is recorded in the results
(`err_fwd`, `err_trans`) so a partially-converged run is detectable after the
fact.

---

## 6. Critical excitability $x^{c}_{0,i}$

For the Epileptor schemes the pipeline also computes a model-based ground
truth for "how epileptogenic is this node", independent of the covariance.

For each non-trivial node $i$: hold every other node at $x_0^{\mathrm{base}}$,
raise $x_{0,i}$, and record the smallest value at which the linearised network
loses stability, $\rho(C) \ge 1$. That value is $x^{c}_{0,i}$. **Low means
fragile**: only a small excitability increase is needed to destabilise the
whole network through this node.

`sweepCriticalExcitability` searches with a coarse forward sweep at step
`x0.step` followed by bisection to `bisect.tol` — cheap where the threshold is
far away, precise where it matters. A single isolated `fsolve` failure is
tolerated and skipped rather than being read as instability, so one numerical
glitch cannot truncate the sweep; a sustained failure *is* treated as seizure
onset, because losing the fixed point is itself the bifurcation.

Nodes that stay stable all the way to `x0.upper` get `NaN`.

$x^{c}_{0,i}$ is what validates the centralities: if $D(\to i)$ really
measures dynamical fragility, it should correlate negatively with
$x^{c}_{0,i}$ across nodes. The per-subject scatter figure plots exactly that.

---

## 7. Classical centralities for comparison

`computeNetworkCentralities` computes, on the same $C$: weighted in/out
degree, left and right eigenvector centrality, PageRank, Katz, self-
communicability $\mathrm{diag}(e^{C})$, betweenness, and in/out closeness.
The last three come from the Brain Connectivity Toolbox; without BCT they are
`NaN` and everything else still runs.

`compareCentralityMeasures` correlates every pair of measures across nodes
(Spearman by default — these measures live on different scales and are
heavy-tailed, so Pearson would be dominated by a few hubs), per subject and
averaged over the cohort.

The practical finding on normalised structural connectomes is that most
classical centralities collapse onto $D(\to i)$, while **betweenness** stays
partly independent, because it reflects path structure rather than
accumulated weight. That is why betweenness is the one classical measure
carried into the group comparisons by default (`defaultCohortMetrics`).

---

## 8. Case-vs-control statistics

`compareCohortGroups` compares each case subject, node by node, against the
control cohort:

$$z_{i,m} = \frac{x_{i,m} - \mu^{\mathrm{ctrl}}_{i,m}}{\sigma^{\mathrm{ctrl}}_{i,m}}$$

for node $i$ and metric $m$, with two-tailed normal p-values.

### Why Bonferroni, and not $|z| \ge 2$

With $N$ nodes and several metrics there are hundreds of implicit tests per
subject. A bare $|z| \ge 2$ rule would be expected to "find" a dozen regions
in pure noise. Flagging therefore uses **Bonferroni correction at `alpha`**
over the family of (non-trivial nodes × metrics) tested for that subject:

$$p^{\mathrm{bonf}} = \min(1, m \cdot p), \qquad m = \#\{\text{eligible tests}\} .$$

`z.threshold` in the config controls only what the console reports *display*.
It cannot promote a node the correction rejected.

A node that is trivial in **any control subject** is excluded from the family
for every case subject, because a control mean built on a placeholder row is
not a meaningful reference.

The correction is deliberately conservative. With a handful of controls the
per-node SD is itself badly estimated, so the honest reading of a flagged
region is "worth looking at", not "established". The summary figure
(`plotCohortOutlierSummary`) separates regions flagged in **all** case
subjects from those flagged in only one, for the same reason.

Implementation: `zScoreOutliers` uses `erfc`/`erfinv` rather than
`normcdf`/`norminv`, so no Statistics Toolbox licence is needed.

---

## 9. Hemispheric asymmetry

`pairHemisphereNodes` matches each `<left prefix><NAME>` node to
`<right prefix><NAME>`. `compareHemisphericAsymmetry` then computes, per pair
and metric, two laterality indices:

$$\mathrm{LI}_{\mathrm{norm}} = \frac{L - R}{L + R}, \qquad \mathrm{LI}_{\mathrm{signed}} = L - R .$$

Positive means the left-labelled node carries more of the metric.

The normalised index is bounded in $[-1,1]$ and scale-free, which is what
makes it comparable across metrics and subjects. The signed index keeps the
effect in the metric's own units, which is what makes a specific finding
interpretable. Which one drives the per-subject z-scores is set by
`laterality.basis`; the group-level tests always use the normalised index,
because pooling across subjects requires the scale-free version.

Three views, answering different questions:

| View | Test | Question |
|------|------|----------|
| Per case subject | z of $\mathrm{LI}$ vs control mean ± SD, Bonferroni over (pairs × metrics) | which regions are unusual in *this* subject |
| Group level | Welch t-test of $\mathrm{LI}_{\mathrm{norm}}$ per pair, Bonferroni | which regions differ *between the groups* |
| Systematic direction | signed-rank of pooled $\mathrm{LI}_{\mathrm{norm}}$ vs zero, per group | does one hemisphere dominate *overall* |

They can disagree, and the disagreement is informative: a region can be
extreme in one subject without the groups differing, and the groups can differ
without any systematic side bias.

Welch (`welchTTest`) rather than Student's t because equal variances are
untestable at these cohort sizes. Signed-rank (`signRankTest`) rather than a
t-test because laterality indices are bounded and typically not normal. Both
are implemented in base MATLAB.

Pairs whose left or right member is trivial in any control subject are dropped
before any testing.

---

## 10. Network-level comparison

`compareDstAcrossCohort` collects $D_{\mathrm{st}}$ and $\rho(C)$ for every
subject, tabulates them, and runs a Welch t-test between the groups.

With `column` or `parkes` normalisation, $\rho(C)$ is pinned by construction
and is identical across subjects — which is exactly what makes
$D_{\mathrm{st}}$ the informative quantity: any difference between subjects is
about the shape of the spectrum, not its leading edge.

---

## 11. Provenance

Every orchestrated run records what produced it, in three places:

* `provenance/experiment_<run>.properties` — the config exactly as run, with
  the results folder appended.
* `provenance/experiment_parameters_<run>.{mat,json}` — the full parameter
  record (`experimentRunParameters`): parameters, environment, MATLAB version,
  cohort, pipeline flags. Written at run start *and* at run end, the second
  time with the per-step outcomes.
* `provenance/run_manifest.mat` — per-step name, success, duration and error.

The same record is embedded inside every analysis `.mat`, so a result file
copied elsewhere still says how it was made.

Results folders are named `<experiment.name>_<yyyy-mm-dd_HHMM>`, so re-running
a config never overwrites an earlier run.

---

## 12. Limitations

* **Linear.** Everything follows from a linearisation around one operating
  point. It describes the onset of instability, not the seizure itself.
* **Structural only.** No functional data enters. The model supplies the
  dynamics; the connectome supplies the structure.
* **Small cohorts.** The per-node control SD is estimated from a handful of
  subjects. Bonferroni correction makes the flagging conservative but cannot
  repair a badly estimated denominator.
* **Normalisation is a modelling choice**, not preprocessing. `column` and
  `parkes` will not generally agree on rankings, and neither is "the" right
  answer — running both and reporting the overlap is the honest approach, and
  is why the batch runner exists.
* **Independence is assumed** by the Bonferroni family, while node metrics on
  a connectome are strongly correlated. The correction is therefore
  conservative rather than exact.
* **Normal null.** The per-subject z-scores assume normality of the control
  distribution per node, which is not testable at these cohort sizes.

---

## 13. Function index

**Data and normalisation**

| Function | Role |
|----------|------|
| `loadConnectome` | read one subject's connectome; trivial-node and symmetry checks |
| `listConnectomeSubjects` | discover subjects under a data root |
| `applyConnectomeNormalisation` | dispatch over normalisation schemes |
| `normalisationSchemeInfo` | canonical scheme name and its implications |
| `matrixNormalizationColumn`, `matrixNormalizationParkes`, `matrixNormalizationTvb` | the individual schemes |
| `makeSyntheticCohort` | generate a synthetic cohort for demos and tests |

**Model and centralities**

| Function | Role |
|----------|------|
| `oneDepileptor`, `CouplingMatrix` | 1-D Epileptor field and its linearisation |
| `healthyEpileptorCoupling` | healthy fixed point and effective $C$ |
| `evalEpileptorStability`, `findCriticalX0` | per-node stability test and threshold search |
| `sweepCriticalExcitability` | the $x^{c}_{0,i}$ sweep over all nodes |
| `covariancesGaussianNet` | stationary covariance $\Omega$ |
| `computeStabilityCentralities` | $D(\to i)$, $D(k \to)$, $D_{\mathrm{st}}$ |
| `computeNetworkCentralities` | classical graph centralities |

**Cohort analyses**

| Function | Role |
|----------|------|
| `runStabilityCentralities` | one subject, end to end |
| `runCohortStabilityCentralities` | the whole cohort plus the cross-subject summary |
| `loadCohortResults`, `packCohortMetric` | load a run and gather one metric across it |
| `compareCohortGroups` | per-node case vs control |
| `compareHemisphericAsymmetry` | laterality, three views |
| `compareCentralityMeasures` | stability vs classical centralities |
| `compareDstAcrossCohort` | network-level $D_{\mathrm{st}}$ |
| `reportTopNodes`, `reportCohortOutliers` | console reports |

**Statistics**

| Function | Role |
|----------|------|
| `zScoreOutliers` | z-scores plus Bonferroni flagging |
| `welchTTest`, `signRankTest` | toolbox-free two-sample and signed-rank tests |
| `pairwiseCorrelation` | NaN-tolerant Pearson/Spearman matrix |
| `rankWithNaN` | ranking that sends NaN to the bottom |
| `pairHemisphereNodes`, `lateralityIndex`, `stripHemispherePrefix` | hemisphere pairing and laterality |
| `cohortGroupMask`, `cohortGroupInfo` | case/control split and its metadata |

**Experiment framework**

| Function | Role |
|----------|------|
| `runConnectomeExperiment` | the orchestrator |
| `runAllConnectomeExperiments` | batch over a configs folder |
| `newConnectomeExperiment` | scaffold a config from the template |
| `loadExperimentConfig` | parse a `.properties` file |
| `experimentRunParameters` | build, merge and save provenance |
| `experimentResultsLayout`, `experimentResultsDir`, `locateExperimentResultsFile` | results folder structure |
| `experimentFolderName`, `experimentProvenanceFile`, `experimentResultPrefix` | artefact naming |
| `listSubjectResultFiles` | discover a run's per-subject outputs |
| `setupConnectomePaths` | put the toolkit (and BCT) on the path |

**Figures**

| Function | Role |
|----------|------|
| `plotConnectomeHeatmaps` | input quality control |
| `plotStabilityCentralitiesFigure`, `plotStabilityScatterFigure` | per-subject figures |
| `plotCohortComparisonFigure` | one case subject against the control band |
| `plotCohortOutlierSummary`, `renderCohortComparisonFigures` | publication-style summaries |
| `collectCohortOutlierRows`, `collectLateralityOutlierRows` | findings as flat rows |
| `cohortFigStyle`, `applyCohortAxesStyle`, `cohortSubjectColors` | shared styling |
| `saveFigureBoth` | write `.fig` and `.png` together |
