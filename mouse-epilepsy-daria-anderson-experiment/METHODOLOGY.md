# Mouse epilepsy connectome experiment — Methodology, theory, and design decisions

This document is the self-contained reference for the
`mouse-epilepsy-daria-anderson-experiment/` pipeline. It is written so that
the methodology, scientific goals, mathematical model, normalisation choices,
statistical decisions, and output schema can be explained to readers (human
or LLM) **without giving them access to the source code**.

For the *operational* manual (how to install, run, and read pipeline outputs)
see [`README.md`](README.md). This document complements it with the
*"why"* behind every design choice.

---

## Table of contents

1. [Scientific goal](#1-scientific-goal)
2. [Input data: mouse coarse connectomes](#2-input-data-mouse-coarse-connectomes)
3. [Dynamical model and linearisation](#3-dynamical-model-and-linearisation)
4. [Per-node stability centralities and \(D_{\mathrm{st}}\)](#4-per-node-stability-centralities-and-d_st)
5. [Critical excitability \(x_{0,i}^{c}\)](#5-critical-excitability-x_0ic)
6. [Connectome normalisation schemes](#6-connectome-normalisation-schemes)
7. [Classical comparison centralities](#7-classical-comparison-centralities)
8. [Cross-cohort statistics: Anderson vs Arnold](#8-cross-cohort-statistics-anderson-vs-arnold)
9. [Left–right hemisphere asymmetry](#9-leftright-hemisphere-asymmetry)
10. [Network-level \(D_{\mathrm{st}}\) cohort comparison](#10-network-level-d_st-cohort-comparison)
11. [Pipeline architecture and provenance model](#11-pipeline-architecture-and-provenance-model)
12. [Output artefact schema](#12-output-artefact-schema)
13. [Key empirical findings on this cohort](#13-key-empirical-findings-on-this-cohort)
14. [Tradeoffs, limitations, and design decisions](#14-tradeoffs-limitations-and-design-decisions)
15. [Glossary of symbols](#15-glossary-of-symbols)

---

## 1. Scientific goal

The experiment asks: **given mouse coarse-grained structural brain
connectomes from two cohorts (Anderson and Arnold) — where the Anderson
strain is suspected to have an epileptic phenotype — can a per-region
"how close to seizure" centrality, computed entirely from the connectome
plus a simple network-dynamical model, identify regions whose role differs
between the two cohorts in a biologically plausible way?**

This is operationalised by:

- Treating each connectome as the coupling structure of a network of
  identically-tuned 1-D Epileptor (or purely linear) units.
- Computing the **stationary covariance** of small linear fluctuations
  around the network's healthy fixed point.
- Reading off per-node "stability centralities" \(D(\to i)\) and
  \(D(k \to)\) from the diagonal of that covariance and its transpose.
- Reading off a scalar **network-level \(D_{\mathrm{st}}\)** as the trace
  of the same covariance.
- For Epileptor-based schemes, also sweeping the per-node excitability
  \(x_{0,i}\) until the network loses linear stability, yielding a
  **critical excitability** \(x_{0,i}^{c}\) per region.

The output supports two kinds of downstream comparison:

1. **Per-mouse / cross-mouse rankings** of the most "susceptible",
   "influential", or "epileptogenic" regions.
2. **Anderson vs Arnold** strain contrasts — per-node z-scores with
   Bonferroni correction, and a scalar \(D_{\mathrm{st}}\) cohort plot.

The pipeline applies the **same** stability-centralities recipe used in
the human Epileptor experiment (`epileptor-experiment/`) and tied back to
the toolkit's core formulas in the [linsync README §1.6](../README.md#16-deviation-from-stability-d_st).

---

## 2. Input data: mouse coarse connectomes

### 2.1 File layout

For each mouse, the relevant input is a single CSV at:

```
data/<mouseId>/fine_family_labelled_coarse.csv
```

It is a **67-line file**: one header row plus 66 data rows. Column 1
holds the row label and columns 2–67 hold a 66×66 weighted adjacency.
Row labels (`L-CORTEX_VISUAL`, …, `R-BACKGROUND`) must equal the column
labels in the same order; this is checked at load time and a warning is
issued on mismatch.

Each connectome comes from DTI-derived structural tractography that has
been coarse-grained from the fine atlas into 66 standard brain-atlas
regions. The matrix is structurally symmetric: across the shipped cohort
\(\lVert K - K^{\top}\rVert_{\max}\) sits in \([1.8, 3.6]\times 10^{-12}\)
absolute — float-roundoff-level noise from however the upstream pipeline
emitted the CSV, and \(\sim 10^{-16}\) relative to
\(\max_{ij} K_{ij} \sim 2\)–\(3\times 10^{4}\).

After the `tvb` normalisation (`normal()`,
[§6](#6-connectome-normalisation-schemes)) this collapses to
\(\lVert K_{\text{tvb}} - K_{\text{tvb}}^{\top}\rVert_{\max} \sim 2\times 10^{-16}\)
on every shipped mouse — machine epsilon at the \(K_{\text{tvb}}\) scale
\(\max_{ij}K_{\text{tvb},ij} = 1\). That symmetry then propagates through
the Epileptor coupling and Lyapunov solves to give
\(\max_{i} |D(\to i) - D(k\to)| \leq 3\times 10^{-14}\) on every shipped
mouse (cf. [§6.2](#62-the-symmetryasymmetry-tradeoff-and-why-it-matters)).
The mechanism is *not* active symmetrisation: `normal()` zeroes the
diagonal, clips at the 95th-percentile threshold \(uu\), and divides by
\(\max_{ij}\); each step preserves symmetry on a symmetric input but
does not impose it. The collapse from \(\sim 10^{-12}\) (raw) to
\(\sim 10^{-16}\) (post-`normal`) is driven by the divide-by-\(\max_{ij}\)
scaling.

> **Truncation-straddle edge case.** The clip step can in principle
> assign \(uu\) to \(K_{ij}\) but not to \(K_{ji}\) when one entry of a
> near-symmetric pair sits just above \(uu\) and the other just below.
> On the shipped cohort this happens on at most 2 entries (only
> `Arnold_3`), the introduced asymmetry is bounded by the input's own
> asymmetry, and the subsequent divide-by-\(\max_{ij}\) keeps the
> post-`normal` symmetry at machine epsilon. So the machine-precision
> identity \(D(\to i) \equiv D(k\to)\) is preserved, but the symmetry
> guarantee is "up to input asymmetry / \(\max_{ij} K_{ij}\)", not
> "exactly symmetric on the nose".

### 2.2 Reference cohort

The canonical cohort is five animals:

`Anderson_1`, `Anderson_2`, `Arnold_3`, `Arnold_4`, `Arnold_5`

with `Arnold_2` excluded automatically (other files present but the
coarse CSV is missing). `Arnold_1` is an optional extra animal that is
included whenever its CSV is present (so the operational cohort is often
6 when the optional CSV ships locally).

The `Anderson_*` mice are the candidate epileptic strain; the `Arnold_*`
mice are the control cohort. The strain split is detected by name prefix
in the comparison scripts.

### 2.3 Trivial nodes

Some atlas regions are placeholders (`L-BACKGROUND`, `R-BACKGROUND`, and
the `*_MASK` rows). These have all-zero connectivity (ignoring the
diagonal) and are flagged as **trivial**. Trivial nodes:

- are kept in the matrix (so node indices stay 1..66 and labels line up
  across mice),
- are **skipped** in the per-node \(x_{0,i}^{c}\) sweep (no Epileptor
  fixed point to solve there),
- are **omitted** from `compareAndersonVsArnold` figures (x-axis shows
  connected brain regions only), and
- are **excluded** from the test family used for Bonferroni correction
  in `compareAndersonVsArnold`.

### 2.4 Provenance check on load

`loadMouseConnectome` enforces three sanity checks before returning:

1. matrix is square,
2. row labels equal column labels in the same order,
3. symmetry error \(\lVert K - K^{\top}\rVert_{\max}\) is recorded (and
   reported in verbose logs).

These are intentional: any silent reshape of the CSV would otherwise
desynchronise per-node results from region names, which is the most
common failure mode for atlas-based analyses.

---

## 3. Dynamical model and linearisation

The same model is used in the human Epileptor experiment; this section
re-states it in the mouse context.

### 3.1 1-D Epileptor + structural coupling

Each region \(i \in \{1, \dots, N\}\) (with \(N=66\)) carries a single
slow Epileptor variable \(z_i(t)\). Coupling between regions is mediated
through the (normalised) connectome matrix \(K\):

\[
\tau_0 \dot{z}_i \;=\; G(z_i) \;-\; 4 x_{0,i} \;-\; \sum_{j=1}^{N} K_{ij}\!\left[F(z_j) - F(z_i)\right],
\]

with

\[
F(z) = \tfrac{1}{4}\!\left(-\tfrac{16}{3} - \sqrt{8z - \tfrac{629.6}{27}}\right), \qquad G(z) = 4 F(z) - z.
\]

Per the linsync convention, \(K_{ij}\) is the directed edge weight from
node \(i\) (source) to node \(j\) (target). \(\tau_0\) is the slow
timescale (default `tau0 = 6667`). The default healthy excitability is
\(x_{0,i} = x_0^{\text{base}} = -2.3\) for every node.

This is the standard 1-D reduction of the Epileptor model used for
seizure-onset analyses. The square-root term in \(F\) is what produces
the characteristic Epileptor bifurcation as \(x_0\) is raised.

### 3.2 Healthy fixed point

The pipeline holds every node at \(x_{0,i} = x_0^{\text{base}}\) and
solves \(\dot{z}=0\) numerically with `fsolve` (Optimization Toolbox
required), warm-starting from \(Z_0 = 3 \cdot \mathbf{1}\). The solver
tolerances are tight (`TolFun = TolX = 1e-14`) and the solution is
rejected if `fsolve` reports a non-positive exit flag.

This gives a per-node fixed point \(z_i^{*}\). If even the
all-healthy network cannot reach a fixed point with \(\rho(C) < 1\),
the pipeline errors out with `healthyEpileptorCoupling:Unstable` — the
input network is already linearly unstable and the rest of the analysis
is undefined.

### 3.3 Jacobian / effective coupling matrix \(C\)

Linearising the network ODE around \(z^{*}\) gives a Jacobian whose
discretised one-step form (with \(\Delta t = 1\)) the toolkit calls the
**effective coupling matrix** \(C\):

\[
C_{ii} \;=\; 1 - \tfrac{1}{\tau_0} + \tfrac{1}{\tau_0}\!\left(4 + \textstyle\sum_{j} K_{ij} - K_{ii}\right) F'(z_i^{*}),
\]

\[
C_{ij} \;=\; -\tfrac{1}{\tau_0}\, K_{ji}\, F'(z_j^{*}) \quad (i \neq j),
\]

with \(F'(z) = -1/\sqrt{8z - 629.6/27}\).

\(C\) is the matrix that drives all downstream stability and covariance
calculations. The spectral radius \(\rho(C) = \max_k |\lambda_k(C)|\) is
the linear-stability indicator: \(\rho(C) < 1\) ⇔ the fixed point is
linearly stable under the canonical OU dynamics described below.

For the **linear schemes** (`parkes`, `column`), the pipeline **bypasses
the Epileptor solve entirely** and uses the normalised connectome \(K\)
*as* \(C\). This is a substantive modelling choice, justified by the
fact that those normalisations encode "synaptic gain" directly; see
[§6](#6-connectome-normalisation-schemes) for the tradeoffs.

### 3.4 Canonical continuous-time OU dynamics

Whatever the source of \(C\), the toolkit assumes the following
continuous-time Ornstein–Uhlenbeck dynamics for the linearised
fluctuations \(X(t) = z(t) - z^{*}\):

\[
\mathrm{d}X(t) = -X(t)\,(I - C)\,\theta\,\mathrm{d}t + \zeta\,\mathrm{d}w(t),
\]

with \(\theta = \zeta = 1\) (fixed in the toolkit; they appear only as
constant multipliers) and \(w(t)\) a multivariate Wiener process with
identity covariance. This is the same model documented in the linsync
README and the underlying synchronizability paper.

The pipeline always uses `discreteTime = false`; the discrete-time VAR
analogue is supported by the toolkit but not exercised here.

---

## 4. Per-node stability centralities and \(D_{\mathrm{st}}\)

### 4.1 Full stationary covariance \(\Omega\)

For the OU process of §3.4, the stationary covariance is the unique
positive-definite solution to the Lyapunov equation; the toolkit computes
it via the power series

\[
\Omega \;=\; \tfrac{\zeta^{2}}{2\theta} \sum_{m=0}^{\infty} 2^{-m} \sum_{u=0}^{m} \binom{m}{u} (C^{u})^{\!\top} C^{\,m-u},
\]

implemented in `covariancesGaussianNet` (with `MaxK = 1e8` iteration cap
by default). The same routine handles the symmetric special case in
closed form. Convergence requires \(\rho(C) < 1\); this is checked
explicitly and the pipeline errors out otherwise.

`computeStabilityCentralities(C, MaxK)` then returns:

\[
\boxed{\;\; D(\to i) \;=\; \Omega_{ii}, \qquad D(k \to) \;=\; \Omega^{\!\top}_{kk} = (\Omega_{C^{\top}})_{kk}, \qquad D_{\mathrm{st}} \;=\; \tfrac{1}{N} \mathrm{trace}(\Omega). \;\;}
\]

The second equality for \(D(k \to)\) requires computing the covariance
of the **transposed** coupling \(C^{\top}\) (i.e. flipping source and
target). For symmetric \(C\) the two diagonals coincide and
\(D(k \to) \equiv D(\to i)\); see §6 on why this matters for the choice
of normalisation.

The error code (`err_fwd`, `err_trans`) from each call is recorded so
that downstream consumers can detect a partially-converged power series.

### 4.2 Interpretation

- **\(D(\to i) = \Omega_{ii}\)** — the **stability susceptibility
  centrality**. It is the steady-state variance the network injects
  *into* node \(i\) via *all* upstream paths under the OU dynamics. A
  large \(D(\to i)\) means small fluctuations propagating through the
  network arrive at node \(i\) amplified; node \(i\) is "downstream
  of a lot".

- **\(D(k \to) = (\Omega_{C^{\top}})_{kk}\)** — the **stability
  influence centrality**. It is the steady-state variance that
  fluctuations originating *at* node \(k\) deposit through the network.
  Large \(D(k\to)\) ⇔ node \(k\) is "upstream of a lot".

- **\(D_{\mathrm{st}} = \tfrac{1}{N}\mathrm{trace}(\Omega)\)** — the
  scalar **Deviation from Stability** for the whole network. It is the
  per-node-averaged variance amplification under the same dynamics; the
  same quantity studied in the 2026 stability paper.

The choice of trace-based \(D_{\mathrm{st}}\) (as opposed to a synchronisation
projection) is deliberate: stability and synchronizability are two
different toolkit "modes" and this pipeline is unambiguously in
stability mode. The CLAUDE.md file at the repo root makes this explicit:
do **not** use the old ncomp `Stability()` helper, which mixes the two.

### 4.3 Why "centrality"?

\(D(\to i)\) and \(D(k\to)\) are network-derived per-node scalars; they
are centrality measures in the same generalised sense as PageRank or
eigenvector centrality (and reduce to those in special regimes; see §12).
They differ in being explicitly **dynamical**: they encode "what the
linearised dynamics on this graph does" rather than just topology.

---

## 5. Critical excitability \(x_{0,i}^{c}\)

For the Epileptor-based schemes (`tvb`, `none`), the pipeline goes beyond
the healthy-state covariance and asks, **for each non-trivial node
\(i\), how much does the excitability of that single node have to rise
before the network as a whole loses linear stability?**

This is the per-node **critical excitability** \(x_{0,i}^{c}\): the
smallest \(x_{0,i}\), with every other \(x_{0,j} = x_0^{\text{base}}\),
such that the new healthy fixed point yields \(\rho(C) \geq 1\).

### 5.1 Search algorithm

1. **Coarse forward sweep.** Step \(x_{0,i}\) from \(x_0^{\text{base}}\)
   upward in increments of \(\Delta x_0 = \text{`x0Step`}\)
   (default `0.01`) until either `fsolve` fails or
   \(\rho(C) \geq 1\). This bracket gives the *last stable* and *first
   unstable* \(x_{0,i}\) values.
2. **Bisection refinement.** Bisect the bracket to a tolerance
   `BisectTol` in \(x_0\) (default `1e-4`). At each test point the
   fixed-point solve is warm-started from the previous solution to keep
   `fsolve` close to the right branch of \(F\).
3. **Skip trivial nodes.** Nodes with all-zero connectivity are recorded
   as `NaN` and the sweep is skipped.
4. **One-failure tolerance.** A single isolated `fsolve` failure inside
   the coarse sweep is treated as a transient numerical glitch and
   skipped without declaring instability; two failures in a row are
   treated as unstable. This prevents a single missed-convergence step
   from truncating the sweep early.

The upper bound `x0Upper = -1.0` is high enough that every connected
mouse node reaches instability before the search runs out of room; if a
node never destabilises within \([x_0^{\text{base}}, x_0^{\text{upper}}]\),
\(x_{0,i}^{c}\) is returned as `NaN`.

### 5.2 Interpretation

A **smaller** \(x_{0,i}^{c}\) means node \(i\) is **closer to its
seizure threshold** when held at the healthy baseline. The summary
ranking
`stabilityCentralities_summary_topnodes_<scheme>.csv` sorts by
**mean rank of \(x_{0,i}^{c}\) across mice (ascending)**, so the top of
the table is the most "epileptogenic" regions in the consensus across
mice.

For the linear schemes (`parkes`, `column`), \(x_{0,i}^{c}\) is **not
defined** (no Epileptor model in play) and the corresponding columns of
the summary CSV are `NaN`. The summary code falls back to ranking by
\(D(\to i)\) in that case (and `reportTopNodes` warns and falls back to
`susc` if the user asks for `x0` sorting under a linear scheme).

---

## 6. Connectome normalisation schemes

`applyConnectomeNormalisation` dispatches over four schemes. Choice of
scheme is the single most consequential decision in the pipeline; the
shipped configs deliberately exercise three of the four side-by-side.

| Scheme   | Symmetric \(C\)? | Epileptor + \(x_{0}^{c}\) sweep? | Default extra parameter | What it does |
|----------|------------------|---------------------------------|------------------------|--------------|
| `tvb`    | Yes              | Yes                             | (none)                 | Zero the diagonal, truncate every entry above the 95th percentile of the off-diagonal upper-triangle, divide by \(\max_{ij} K_{ij}\). Same `normal()` recipe used in the TVB human Epileptor pipeline. |
| `none`   | Yes              | Yes                             | (none)                 | Use \(K_{\text{raw}}\) unchanged. |
| `parkes` | Yes              | No                              | `parkes.c = 1.0`       | \(K = A / (\rho(A) + c)\) — the Parkes et al. (Nature Protocols 2024) prescription. Symmetric input gives symmetric \(K\); the continuous-time \(-I\) subtraction is supplied implicitly by `con2cov`. |
| `column` | No (in general)  | No                              | `col.scale = 0.95`     | Rescale each column to sum to `col.scale`, breaking the input's symmetry. The prescription Liao & Lizier (2026) use in Appendix H Fig 13. |

### 6.1 Why pick a scheme?

Each scheme is a *modelling claim* about how the structural connectome
should map to coupling gains in a dynamical model. They give qualitatively
different answers:

- **`tvb`** mimics the truncation+normalise pipeline used by The Virtual
  Brain; numerical results match the human Epileptor experiment by
  construction. It is the slowest scheme because the per-node
  \(x_{0,i}^{c}\) sweep is run (≈10 min/mouse).

- **`none`** is the no-rescale baseline; useful as a control and for
  diagnosing whether downstream differences are driven by the
  normalisation versus the connectome itself.

- **`parkes`** pins the spectral radius just below 1 by construction
  (\(\rho(K) \to \rho(A)/(\rho(A)+c) < 1\)), so the healthy
  configuration is *guaranteed* linearly stable for any positive `c`.
  It preserves the input symmetry. This is the scheme that interfaces
  most cleanly with the classical centrality literature (much of which
  assumes undirected graphs).

- **`column`** pins every column sum to `col.scale = 0.95`, so each node
  receives the same *total* upstream weight regardless of how that weight
  is distributed. This **breaks the input's symmetry**, making \(C\)
  genuinely directed. \(D(\to i)\) and \(D(k\to)\) then carry different
  information, and the in/out classical centralities decouple from each
  other and from \(D\).

### 6.2 The symmetry/asymmetry tradeoff (and why it matters)

DTI-derived structural connectomes are symmetric by construction. Under
a symmetry-preserving normalisation (`tvb`, `none`, `parkes`):

- \(D(\to i) \equiv D(k\to)\) at machine precision — the in/out
  distinction the stability centralities are designed to expose is
  trivially absent.
- All power-series-based centralities (in/out-DC, EC^L/R, PR, KZ, SelfC,
  \(D(\to i)\)) collapse onto a leading-eigenvector-driven ranking with
  pairwise Spearman correlation \(\rho \geq 0.96\) on this dataset
  (see [§12](#12-key-empirical-findings-on-this-cohort)).
- The only classical centrality that carries genuinely independent
  per-node information is **Betweenness Centrality (BC)** at
  \(\rho \approx 0.78\).

Under the asymmetry-imposing normalisation (`column`):

- \(D(\to i)\) and \(D(k\to)\) decouple (\(\rho \approx 0.72\)).
- In-style measures (EC^L, in-DC, KZ, CC^in) sign-flip to
  \(\rho \approx -0.78\) with \(D(\to i)\) — the Liao & Lizier (2026)
  Fig 12 → Fig 13 transition reproduced on the mouse data.
- The two strains' \(D_{\mathrm{st}}\) values **separate**: the highest
  Arnold value is below the lowest Anderson value. (Caveat: \(n=2\) vs
  \(n=3\); see §13.)

This is the load-bearing methodological choice. The shipped reference
analysis (`results/mouse_experiment_reference/centrality_comparison_analysis.md`)
documents the full Parkes-vs-column comparison in detail; the headline is
that **conclusions about whether \(D\) "adds information" depend on the
normalisation choice**, and the column normalisation is the regime in
which the stability centralities can carry independent signal on this
data.

### 6.3 The Parkes-`c` non-knob

An earlier version of the analysis hoped that lowering Parkes `c` to
push \(\rho \to 1\) would decouple \(D\) from EC. **It does not.** A
sweep across seven orders of magnitude of `c` on Anderson_1 yields
Spearman \(\rho(D, \mathrm{EC})\) in `[0.93, 1.00]`. The argument is
simple: for symmetric \(C\), \(\Omega = (\zeta^{2}/2)(I-C)^{-1}\), so
\(D(\to i) = (\zeta^{2}/2)\sum_{k} v_{k}(i)^{2}/(1-\lambda_{k})\). Parkes
`c` rescales the spectrum uniformly; the eigenvectors and spectral-gap
ratios are invariant. Conclusion: `c` is the wrong knob; **breaking
symmetry** is the right knob. That is what `column` does.

---

## 7. Classical comparison centralities

`computeNetworkCentralities(C_healthy)` is called inside
`runMouseStabilityCentralities` after the stability centralities, and
attaches the following per-node measures to `results.centralities`:

| Field name           | Symbol        | What it is                                                                  |
|----------------------|---------------|-----------------------------------------------------------------------------|
| `in_strength`        | in-DC         | Weighted in-degree \(\sum_{j} C_{j,i}\)                                     |
| `out_strength`       | out-DC        | Weighted out-degree \(\sum_{j} C_{i,j}\)                                    |
| `eigenvector_right`  | \(\mathrm{EC}^{R}\) | Leading right eigenvector of \(C\) (out-style)                              |
| `eigenvector_left`   | \(\mathrm{EC}^{L}\) | Leading right eigenvector of \(C^{\top}\) (in-style)                        |
| `pagerank`           | PR            | PageRank with damping `AlphaPR = 0.85`, dangling-node uniform teleportation |
| `katz`               | KZ            | Katz centrality with attenuation `AlphaKZ/rho(C) = 0.5/\rho(C)`             |
| `self_comm`          | SelfC         | Self-communicability \(\mathrm{diag}(\exp(C))\)                             |
| `betweenness`        | BC            | Weighted betweenness on the length matrix \(L_{ij} = 1/C_{ij}\) (BCT)       |
| `closeness_in`       | \(\mathrm{CC}^{in}\)  | \((N-1) / \sum_{j} d(j \to i)\) (BCT)                                       |
| `closeness_out`      | \(\mathrm{CC}^{out}\) | \((N-1) / \sum_{j} d(i \to j)\) (BCT)                                       |

Conventions:

- `C(i,j)` is the directed edge weight from source \(i\) to target \(j\)
  (linsync convention).
- All measures are computed on `C_healthy` (Epileptor schemes) or the
  normalised `K` (linear schemes), with the diagonal zeroed first
  (`IgnoreSelfLoops = true`).
- Betweenness and closeness require the
  [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/);
  if BCT is not on the path the corresponding fields are returned as
  `NaN` and a single warning is issued. The remaining centralities still
  run.

These are the same set used by Liao & Lizier (2026) for their
centrality-comparison analysis (their Appendix G). The point of running
them is to test whether the stability centralities reduce to a classical
centrality on this dataset (the answer depends on normalisation; see
§§6.2 and 12).

`compareCentralityMeasures` then computes the M×M Spearman correlation
matrix (M = 12 measures) per mouse, the mean and SD across mice, and
both a heatmap (with hierarchical-clustering reorder, matching the
paper's Figs 12–13 convention) and a bar chart comparing \(D(\to i)\)
and \(D(k\to)\) against each classical centrality (paper's Fig 4
analogue, aggregated across mice instead of across a randomisation
parameter).

---

## 8. Cross-cohort statistics: Anderson vs Arnold

`compareAndersonVsArnold` is the per-node strain contrast. For a given
normalisation scheme, and each Anderson mouse \(a\):

1. Load all per-mouse `stabilityCentralities_*_<scheme>_results.mat`
   files in the chosen results folder.
2. Verify that all mice share the same parcellation (66 labels in the
   same order); warn on mismatch.
3. Build per-node Arnold cohort statistics: mean and SD across the
   \(n_{\mathrm{Arn}}\) Arnold animals for each of three metrics:
   - \(D(\to i)\)  (always present)
   - \(D(k \to)\)  (always present, but ≡ \(D(\to i)\) under symmetric \(C\))
   - BC           (only if BCT was available; otherwise the panel is
     silently skipped and the comparison family shrinks)
4. For Anderson mouse \(a\), compute per-node z-scores
   \(z_{a,i} = (x_{a,i} - \mu_{\text{Arn},i}) / \sigma_{\text{Arn},i}\),
   skipping nodes where the Arnold SD is zero or `NaN` and skipping
   nodes that are trivial in *any* Arnold mouse.
5. Convert to two-sided p-values via the normal approximation:
   \(p_{a,i} = 2(1 - \Phi(|z_{a,i}|))\).

### 8.1 Bonferroni correction (the actual flagging rule)

The family for the correction is **per Anderson mouse**: it is the union
of all non-trivial nodes \(\times\) all included metrics. So with
\(n_{\mathrm{nontrivial}} \approx 64\) and three metrics (when BC is
available) the family size is \(m \approx 192\); with BC missing the
family shrinks to \(m \approx 128\) and the corrected threshold tightens
accordingly.

A node is flagged as an **outlier** for mouse \(a\) and metric \(M\) if
the Bonferroni-corrected p-value satisfies \(m \cdot p_{a,i}^{(M)} < \alpha\),
with `alpha = 0.05` by default. Equivalently, in z-units this is
\(|z_{a,i}^{(M)}| \geq z_{\alpha / 2m} = \Phi^{-1}(1 - \alpha/(2m))\); the
plot title displays this **corrected** threshold rather than the
nominal `ZThreshold` from the config.

### 8.2 Why Bonferroni rather than \(|z| \geq 2\)?

The plot would otherwise flag dozens of nodes per Anderson mouse purely
from multiple testing. Bonferroni is conservative but transparent:
the entire family of node × metric tests for that Anderson mouse is the
denominator, and the corrected threshold scales with how many metrics
are present.

The `z.threshold` config key is **decoupled** from the flagging rule.
It controls only what `reportAndersonOutliers` chooses to **display** in
the human-readable `reports.log`. This separation is intentional: the
Bonferroni rule is the statistical decision (calibrated by `alpha`); the
z-threshold is just a verbosity knob.

### 8.3 What gets saved

For each Anderson mouse \(a\), one figure per metric (up to three:
\(D(\to i)\), \(D(k\to)\), BC):

- **shaded band** — Arnold mean ± SD per node
- **blue line** — Arnold mean
- **grey dots** — individual Arnold values at each node
- **red markers** — the Anderson mouse's per-node value
- **black ring + label** — Bonferroni-flagged outlier nodes (z-score
  above the marker; x-tick label **bold** and coloured red)
- trivial atlas placeholders (`*BACKGROUND`, `*_MASK`) are not shown on
  the x-axis

The Bonferroni-flagged nodes are also written to
`compare_anderson_vs_arnold_<mouseId>_<scheme>_outliers.csv`, sorted by
the smallest of the three Bonferroni p-values for that node. The CSV
contains every test statistic and Bonferroni p-value for the flagged
nodes (so a reader can verify which metric drove the flag).

### 8.4 Multi-mouse summarisation

`reportAndersonOutliers` reads the per-Anderson CSVs and produces a
plain-text report that:

1. lists each mouse's outliers (filtered by `ZThreshold` for *display*,
   sorted by \(|z|\)), and
2. lists **shared outliers** — regions flagged across multiple Anderson
   mice. The shared list is sorted by how many mice flagged each region
   (descending), then alphabetically. This is the most useful column for
   biological interpretation: a node that lights up in both Anderson
   mice for the same metric is more likely to reflect strain biology
   than experimental noise.

---

## 9. Left–right hemisphere asymmetry

`compareLeftRightAsymmetry` addresses a complementary question to §8:
**after** per-node strain contrasts, do epileptic (Anderson) mice show
*different left–right balance* on corresponding homotopic regions than
controls (Arnold)? This is motivated by Daria's hemispheric intervention
(one hemisphere edited — labelled `L-` in the CSV though possibly
anatomical right), affecting hippocampus and other seizure-onset regions.

### 9.1 Pairing homotopic regions

The coarse atlas uses `L-<SUFFIX>` and `R-<SUFFIX>` labels (e.g.
`L-HIPPOCAMPAL FORMATION` / `R-HIPPOCAMPAL FORMATION`). The pipeline
matches pairs by shared suffix. Pairs where either side is trivial in
any Arnold mouse (same rule as §8) are excluded from testing.

### 9.2 Laterality indices

For each mouse, region pair \(p\), and metric \(M \in \{D(\to i), D(k\to), \mathrm{BC}\}\):

\[
\mathrm{LI}_{\mathrm{norm}}(p, M) = \frac{L - R}{L + R}, \qquad
\mathrm{LI}_{\mathrm{signed}}(p, M) = L - R,
\]

with \(L,R\) the metric values at the left- and right-labelled nodes.
\(\mathrm{LI}_{\mathrm{norm}} \in [-1, +1]\): positive ⇒ higher on the
`L-` side; zero ⇒ symmetric. `NaN` when \(L+R \le 0\) or either side is
missing. **Primary score for statistics:** \(\mathrm{LI}_{\mathrm{norm}}\);
\(\mathrm{LI}_{\mathrm{signed}}\) is stored for interpretation in metric
units.

### 9.3 Per-Anderson-mouse view

For each Anderson mouse \(a\), per pair and metric:

\[
z_{a,p}^{(M)} = \frac{\mathrm{LI}_{\mathrm{norm},a,p}^{(M)} - \mu_{\mathrm{Arn},p}^{(M)}}{\sigma_{\mathrm{Arn},p}^{(M)}},
\]

using Arnold cohort mean and SD of \(\mathrm{LI}_{\mathrm{norm}}\) at
that pair. Bonferroni correction over the family **(all non-trivial pairs
× all included metrics)** for that mouse — same convention as
§8.1. Outliers are written to
`compare_LR_asymmetry_<mouse>_<scheme>_outliers.csv` (optional).

Figures (`compare_LR_asymmetry_<mouse>_<scheme>.{fig,png}`) plot
\(\mathrm{LI}_{\mathrm{norm}}\) per pair with Arnold mean ± SD, a dashed
\(\mathrm{LI}=0\) reference, and Bonferroni-flagged pairs labelled.

### 9.4 Group-level view

Per pair × metric, a **Welch t-test** compares \(\mathrm{LI}_{\mathrm{norm}}\)
between Anderson and Arnold mice (`ttest2`, unequal variance). Bonferroni
correction over all pair × metric tests in the run. Results:
`LR_asymmetry_groupTest_<scheme>.csv` and a scatter figure of Anderson
vs Arnold pair means (highlighted when Bonferroni-significant).

### 9.5 Systematic direction (which hemisphere dominates?)

To answer whether epileptic mice show a **cohort-wide** bias toward one
side (relevant to the L/R labelling question), the pipeline pools all
\(\mathrm{LI}_{\mathrm{norm}}\) values per cohort × metric (all pairs ×
all mice in that cohort) and reports mean ± SD plus a **sign-rank test**
against zero (`signrank`). If mean \(> 0\) and \(p < \alpha\), the
inferred dominant side is **`L`** (CSV label); if mean \(< 0\), **`R`**;
otherwise **`none`**. Outputs: `LR_asymmetry_systematic_<scheme>.{csv,fig,png}`.

### 9.6 Packed summary

`LR_asymmetry_<scheme>.mat` holds pair indices, all per-mouse
\(\mathrm{LI}_{\mathrm{norm}}\) / \(\mathrm{LI}_{\mathrm{signed}}\)
matrices, Arnold statistics, per-Anderson z/p tables, group-test results,
systematic summaries, and `runParameters`.

Pipeline flag: `pipeline.runLRAsymmetry` (default `true`), run after
`compareAndersonVsArnold` when both strains are present.

---

## 10. Network-level \(D_{\mathrm{st}}\) cohort comparison

`compareDstAcrossCohort` is the simplest piece of the pipeline and the
one most directly answering the original scientific question. It loads
`results.D_st_healthy` from every per-mouse file, classifies each mouse
as Anderson or Arnold by name prefix, and:

1. Writes a CSV with `mouse`, `D_st`, `rho_C` (each computed from
   `eig(C_healthy)`).
2. Prints cohort means, SDs, and ranges.
3. Runs a **Welch t-test** (`ttest2` with `Vartype='unequal'`) when
   both cohorts have \(n \geq 2\). The Welch flavour is used because
   the two cohort SDs are different; equal-variance pooling would be
   wrong here.
4. Plots a per-mouse bar chart with horizontal dashed lines at each
   cohort mean, value labels above each bar, and bar colour by strain.
   Title carries the Welch \(t\) and \(p\) when present.

The key sanity check in the figure is the relative position of the two
horizontal lines and the *overlap* between the bars of each strain.
Under `column` normalisation on the reference cohort the bars do not
overlap (smallest Anderson > largest Arnold); under `parkes` they do.

This is a scalar/network-level analog of the per-node Anderson-vs-Arnold
analysis and serves as a sanity check: if \(D_{\mathrm{st}}\) does not
discriminate cohorts, no per-node test that aggregates to \(D_{\mathrm{st}}\)
is going to discriminate them meaningfully either.

---

## 11. Pipeline architecture and provenance model

### 10.1 The orchestrator pattern

The pipeline is driven by **one config file per experiment**:

- A `.properties` file under `configs/` pins **one normalisation scheme**
  and a complete set of Epileptor / comparison / pipeline parameters.
- `runMouseExperiment('configs/<name>.properties')` orchestrates the
  full pipeline for that config.
- `runAllMouseExperiments()` runs every non-template `.properties` file
  in `configs/` in sequence (typically the three shipped baselines:
  `initial_column`, `initial_parkes`, `initial_tvb`).

This was a deliberate departure from earlier scripted runs that mutated
a flat `results/` folder. The key properties of the new design are:

1. **One scheme per experiment.** Schemes are sufficiently different
   methodologically that mixing their outputs in the same folder is
   confusing; separating them at the config level makes that explicit.
2. **Timestamped results folders.** The orchestrator writes to
   `results/<experiment.name>_<yyyy-mm-dd_HHMM>/`, so two runs of the
   same config never overwrite each other.
3. **Per-step manifest.** Each pipeline step is wrapped in
   `runPipelineStep`, which records `name`, `success`, `durationSec`,
   and `errorMessage`. By default a failure does **not** abort the
   pipeline (`pipeline.stopOnError=false`) — subsequent steps that need
   the failed output skip themselves, and the manifest records why.
4. **Cohort post-flight check.** After per-mouse runs, the orchestrator
   verifies that every available mouse produced a results file. If any
   mouse failed, the per-mouse step is recorded as a failure and
   downstream comparison/report steps are *skipped* rather than
   silently producing a misleading partial summary.
5. **Output assertion.** At the end of a run with all expected steps
   succeeding, `assertMouseExperimentOutputs` checks that every required
   filename exists, and errors if any are missing. (Anderson outlier
   CSVs are exempt because they are only written when outliers exist.)
6. **Skip-when-impossible.** If the cohort has no Anderson or no Arnold
   mice (a common case when a developer drops in only a partial CSV),
   `compareAndersonVsArnold` is marked **skipped** (counts as success)
   and the pipeline continues.

### 10.2 Provenance: every parameter, every run

A guiding principle of the pipeline is that **every output should be
self-describing**. The provenance model has four layers:

| Layer                      | Where it lives                                                                                                                            |
|----------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| Raw config                 | `experiment_<run>.properties` (verbatim config copy; `<run>` = results subfolder name; actual folder appended at run start)               |
| Structured snapshot        | `experiment_parameters_<run>.{mat,json}` — the full `runParams` struct (saved at start *and* end of the run)                                |
| Run manifest               | `run_manifest.mat` — `{ manifest, cfg, runParams }` with per-step success/duration/error                                                  |
| Embedded in result `.mat`  | Every per-mouse / summary / comparison `.mat` file carries a `runParameters` field with the full provenance struct                        |

`runParams` contains (at minimum):

- **Metadata** — `experimentName`, `experimentDescription`, `configFile`,
  `recordedAt`, `resultsDir`.
- **Environment** — MATLAB version, computer arch, hostname, user,
  toolkit / experiment root paths.
- **Cohort** — full list of mouse IDs and count.
- **`normalisation`** — scheme name for this experiment.
- **`stabilityCentralities`** — `parkesC`, `colScale`, `x0Base`,
  `x0Upper`, `x0Step`, `bisectTol`, `maxK`, `tau0`, `topK`,
  `discreteTime` (always `false`), fsolve tolerances, classical
  centrality defaults (`alphaPR`, `alphaKZ`, `ignoreSelfLoops`).
- **`comparison`** — `corrType`, `alpha`, `zThreshold`, plus
  comparison-script options like `compareCentralityReorder`.
- **`pipeline`** — boolean for each pipeline step (stability centralities,
  centrality correlations, Anderson vs Arnold, L–R asymmetry, \(D_{\mathrm{st}}\),
  reports, optional heatmap) plus
  `stopOnError`.
- **`output`** — `saveResults`, `plot`, `verbose`.
- **`config`** — the full parsed config struct from
  `loadMouseExperimentConfig`.
- **`finishedAt`**, **`pipelineSteps`** — appended at run end.

The JSON copy is included specifically so that diffs across experiments
work in plain text editors and version control without needing MATLAB
to read them.

### 10.3 Failure modes the architecture handles

- **Incomplete cohort.** `runAllMiceStabilityCentralities` errors out
  (`IncompleteCohort`) if any mouse failed its per-mouse run.
  Comparison steps are skipped rather than producing a misleading
  partial summary.
- **Missing BCT.** Classical centralities that depend on BCT
  (`betweenness`, `closeness_in/out`) are returned as `NaN` with a
  single warning. `compareAndersonVsArnold` silently skips the BC figure.
- **Linear scheme \(\rho(C) \geq 1\).** `runMouseStabilityCentralities`
  errors with `LinearSchemeUnstable` and a hint to raise `parkes.c`
  (for `parkes`) or lower `col.scale` (for `column`).
- **Power series didn't converge** — `err_fwd != 0` or `err_trans != 0`
  surface as warnings and are recorded in the result struct so consumers
  can decide whether to discard.
- **Strain not represented** — `compareAndersonVsArnold` is skipped
  (with the skip *counting as success* in the manifest) so partial
  cohorts can still run the per-mouse and \(D_{\mathrm{st}}\) steps
  without failing the experiment.

---

## 12. Output artefact schema

A complete `runMouseExperiment` run produces the following files under
`results/<experiment.name>_<yyyy-mm-dd_HHMM>/`, grouped into category
subfolders (`provenance/`, `stability_centralities/`, `anderson_vs_arnold/`,
`lr_asymmetry/`, `centrality_correlation/`, `dst_cohort/`, `qc/`,
`comparison_figures/`). Paths below are relative to that experiment root.
Loaders also accept the legacy flat layout (same filenames in the root).
All filenames are parameterised by the `<scheme>` (`column`, `parkes`, `tvb`,
or `none`).

### 12.1 Provenance (`provenance/`)

| File                            | Contents                                              |
|---------------------------------|-------------------------------------------------------|
| `experiment_<run>.properties`         | Verbatim config + appended actual results folder      |
| `experiment_parameters_<run>.mat`     | `runParams` struct                                    |
| `experiment_parameters_<run>.json`    | `runParams` struct as JSON                            |
| `run_manifest.mat`                    | `{ manifest, cfg, runParams }`                        |
| `reports_<run>.log`                   | `reportTopNodes` + `reportAndersonOutliers` console   |

### 12.2 Per-mouse stability centralities (`stability_centralities/`, ×5–6)

For each mouse `<mouseId>`:

| File                                                                  | Contents                                                                            |
|-----------------------------------------------------------------------|-------------------------------------------------------------------------------------|
| `stabilityCentralities_<mouseId>_<scheme>_results.mat`                | `results` struct (see schema below)                                                 |
| `stabilityCentralities_<mouseId>_<scheme>_figure.{fig,png}`           | Per-mouse bar / scatter figure                                                      |

The `results` struct contains:

- `mouseId`, `labels`, `isTrivial`, `symmetryError`, `N`
- `tau0`, `x0_base`, `x0_upper`, `x0_step`, `bisectTol`
- `normalisation`, `parkesC`, `colScale`, `isParkes`, `isColumn`,
  `isLinearScheme`
- `K_raw`, `K` (raw and normalised connectome)
- `z_fixed_healthy`, `C_healthy`, `rho_healthy` (Epileptor; the first is
  `NaN` for linear schemes)
- `Omega`, `OmegaTranspose` (full stationary covariance)
- `D_susceptibility`, `D_influence`, `D_st_healthy`
- `centralities` (struct of classical centralities; see §7)
- `x0_crit`, `rho_at_crit` (Epileptor schemes only; `NaN` otherwise)
- `err_fwd`, `err_trans`, `sweepTime`
- `runParameters` (full provenance struct)

### 12.3 Cross-mouse summary

| File                                                            | Contents                                                                  |
|-----------------------------------------------------------------|---------------------------------------------------------------------------|
| `stabilityCentralities_summary_topnodes_<scheme>.csv`           | Per-mouse + mean rank of `x0_crit`, `D_susc`, `D_infl` for every region. Carries an `is_trivial` column (`1` for `BACKGROUND` / `*_MASK` placeholder rows) so consumers like `reportTopNodes` can filter them out of the top-K |
| `stabilityCentralities_summary_overview_<scheme>.{fig,png}`     | Region × mouse heatmaps (1–3 panels; see below)                           |
| `stabilityCentralities_summary_<scheme>.mat`                    | Packed matrices `D_susc_all`, `D_infl_all`, `x0_crit_all`, mean ranks, `runParameters` |

Overview heatmap panel counts:

- `parkes`: 1 panel — \(D(\to i)\) only (since \(D(k\to) \equiv D(\to i)\))
- `column`: 2 panels — \(D(\to i)\) and \(D(k\to)\)
- `tvb` / `none`: 3 panels — \(x_{0,i}^{c}\), \(D(\to i)\), \(D(k\to)\)

### 12.4 Centrality correlations (`centrality_correlation/`)

| File                                                       | Contents                                                                          |
|------------------------------------------------------------|-----------------------------------------------------------------------------------|
| `centrality_corr_<mouseId>_<scheme>.{fig,png}` (×N)        | Per-mouse Spearman correlation heatmap                                            |
| `centrality_corr_mean_<scheme>.{fig,png}`                  | Mean correlation heatmap across mice                                              |
| `centrality_corr_bars_<scheme>.{fig,png}`                  | Bar chart of `corr(D(→i), <classical>)` and `corr(D(k→), <classical>)`            |
| `centrality_corr_<scheme>.mat`                             | Packed summary: `corrAll`, `meanCorr`, `sdCorr`, ordering, labels, `runParameters` |

### 12.5 Anderson vs Arnold (`anderson_vs_arnold/`)

| File                                                                       | Contents                                                                  |
|----------------------------------------------------------------------------|---------------------------------------------------------------------------|
| `compare_anderson_vs_arnold_<andersonMouse>_<scheme>_D_to_i.{fig,png}` (per Anderson) | \(D(\to i)\) vs Arnold cohort mean ± SD                          |
| `compare_anderson_vs_arnold_<andersonMouse>_<scheme>_D_k_to.{fig,png}` (per Anderson) | \(D(k\to)\) vs Arnold cohort mean ± SD                             |
| `compare_anderson_vs_arnold_<andersonMouse>_<scheme>_BC.{fig,png}` (per Anderson, when BC present) | BC vs Arnold cohort mean ± SD                         |
| `compare_anderson_vs_arnold_<andersonMouse>_<scheme>_outliers.csv` (optional) | Bonferroni-flagged outlier rows                                            |
| `compare_anderson_vs_arnold_<scheme>.mat`                                  | Cohort matrices, mean/SD, per-mouse z/p tables, `runParameters`           |

### 12.6 Left–right asymmetry (`lr_asymmetry/`)

| File                                                                       | Contents                                                                  |
|----------------------------------------------------------------------------|---------------------------------------------------------------------------|
| `compare_LR_asymmetry_<andersonMouse>_<scheme>.{fig,png}` (per Anderson) | Per-pair \(\mathrm{LI}_{\mathrm{norm}}\) vs Arnold mean ± SD             |
| `compare_LR_asymmetry_<andersonMouse>_<scheme>_outliers.csv` (optional)  | Bonferroni-flagged pair × metric rows                                     |
| `LR_asymmetry_groupTest_<scheme>.{csv,fig,png}`                            | Welch t-tests per pair × metric; scatter of cohort means                  |
| `LR_asymmetry_systematic_<scheme>.{csv,fig,png}`                           | Cohort-wide mean \(\mathrm{LI}_{\mathrm{norm}}\), sign-rank vs zero       |
| `LR_asymmetry_<scheme>.mat`                                               | Full packed summary + `runParameters`                                     |

### 12.7 \(D_{\mathrm{st}}\) cohort (`dst_cohort/`)

| File                                                       | Contents                                                                          |
|------------------------------------------------------------|-----------------------------------------------------------------------------------|
| `D_st_cohort_<scheme>.csv`                                 | Table: `mouse`, `D_st`, `rho_C`                                                   |
| `D_st_cohort_<scheme>.{fig,png}`                           | Per-mouse bar chart with strain mean lines, Welch t-test in title (if applicable) |
| `D_st_cohort_<scheme>.mat`                                 | Full summary struct + `runParameters`                                             |

### 12.8 Optional QC

| File                                                       | Contents                                                                          |
|------------------------------------------------------------|-----------------------------------------------------------------------------------|
| `mouse_heatmaps_overview_reference.{fig,png}`              | Multi-panel `imagesc` of every mouse's raw \(K\), permuted to match the reference PNGs in `data/<mouseId>/coarse_connectome_*.png` |

Only written if `pipeline.runHeatmap=true`.

---

## 13. Key empirical findings on this cohort

(These are summarised from `results/mouse_experiment_reference/centrality_comparison_analysis.md`,
which is the load-bearing analysis document for the experiment. Numbers
below are from the n=5 reference cohort: Anderson_1, Anderson_2,
Arnold_1, Arnold_3, Arnold_4. Add Arnold_5 and the qualitative picture
is unchanged.)

### 12.1 Under `parkes`: \(D\) collapses onto eigenvector centrality

On Parkes-normalised symmetric mouse connectomes, the pairwise Spearman
correlations between every classical power-series-based centrality
(in/out-DC, \(\mathrm{EC}^{L/R}\), PR, KZ, SelfC) and the stability
centralities \(D(\to i), D(k\to)\) are **all ≥ 0.96**. The lone
exception is **BC at \(\rho \approx 0.78\)**.

Why: three properties of the data + regime combine to collapse the
power series \(\Omega = (\zeta^2/2) \sum_{m,u} 2^{-m} \binom{m}{u} (C^u)^{\top} C^{m-u}\)
onto the leading-eigenvector outer product.

1. \(C\) is exactly symmetric (DTI is undirected; Parkes preserves
   symmetry). \(D(\to i) \equiv D(k\to)\) is forced.
2. Significant spectral gap: \(\rho_2 / \rho_1 \in [0.49, 0.68]\) across
   the cohort.
3. Operating well below the stability boundary (Parkes pins
   \(\rho(C) = \rho(A)/(\rho(A)+c)\) close to but below 1, but higher-order
   terms are still negligible).

So all measures asymptote to \(C^u \approx \lambda_1^u v v^{\top}\) and
\(\mathrm{diag}(\Omega)_i \approx \text{const}\cdot v_i^2\), which under
Spearman ranking is just eigenvector centrality.

**This is not a bug.** It is the textbook regime in which the stability
centralities should reduce to EC, and BC is the only classical
centrality that adds anything.

### 12.2 Under `column`: \(D\) decouples and the strains separate

Re-running the entire pipeline with `column` normalisation flips the
qualitative picture:

- The spectrum is no longer dominated by the leading eigenvalue (every
  column sum is pinned to 0.95, so \(\rho(C) = 0.95\) by construction);
  higher-order terms in the series matter.
- \(D(\to i)\) and \(D(k\to)\) decouple: Spearman \(\rho \approx 0.72\).
- In-style classical measures **sign-flip** to negative correlation
  with \(D(\to i)\): \(\mathrm{EC}^{L}\) at \(-0.78\), in-DC at
  \(-0.75\), KZ at \(-0.81\). Out-style measures keep moderate positive
  correlation.
- \(D_{\mathrm{st}}\) separates the strains: Anderson cohort at
  \(0.792\pm 0.010\), Arnold cohort at \(0.768\pm 0.007\); the two
  ranges do not overlap.
- Per-Anderson Bonferroni outliers triple in count and converge on
  epileptologically interpretable regions: **hippocampus, amygdala,
  hypothalamus** light up under \(D(\to i)\) and BC across both Anderson
  mice. Under Parkes the hippocampus is a BC-only finding on a single
  Anderson mouse.

This is the Liao & Lizier (2026) Fig 12 → Fig 13 transition reproduced
on the mouse data, and it is the empirical justification for treating
`column` as the primary scheme of interest.

### 12.3 \(D(k\to)\) does something \(D(\to i)\) does not

Under `column`, a few standout per-node patterns are visible *only* in
\(D(k\to)\):

- `R-CORTEX_SOMATOSENSORY` on Anderson_1 has \(z_{D(\to i)} = +3.9\)
  but \(z_{D(k\to)} = +22.5\) — massively elevated influence capacity
  with only mildly elevated susceptibility.
- The same region on Anderson_2 has \(z_{D(k\to)} = -7.9\) — dramatically
  *depressed* influence capacity. Same region, opposite direction in
  the two animals.

Neither pattern is visible to \(D(\to i)\) alone, and neither is visible
at all under `parkes` (where \(D(k\to) \equiv D(\to i)\) by construction).

### 12.4 Under `tvb`: trivial placeholder rows inflate \(\Omega_{ii}\)

The TVB normalisation (zero diagonal, 95th-percentile truncate,
rescale to \([0,1]\)) leaves the placeholder rows (`L-BACKGROUND`,
`R-BACKGROUND`, `*_MASK`) at exactly zero connectivity, so their
column of \(C\) collapses to a single self-loop term. This makes
their \(\Omega_{ii}\) numerically identical across mice
(\(D(\to i) \approx 589.6\) for both `*BACKGROUND` rows) and large
enough to dominate the \(D(\to i)\) ranking. Under `parkes` and
`column` the trivial rows are still trivial but their \(D(\to i)\)
sits below the genuinely connected regions, so the issue is
TVB-specific.

This is a presentation artefact: the per-node `runMouseStabilityCentralities`
results still mark `isTrivial(i) = true` on those rows, and the
comparison / outlier scripts (`compareAndersonVsArnold`,
`compareLeftRightAsymmetry`) already exclude them from the test
family. The cohort summary CSV
(`stabilityCentralities_summary_topnodes_<scheme>.csv`) therefore
carries an `is_trivial` column and `reportTopNodes` filters those
rows out of the printed top-K by default (pass
`'IncludeTrivial', true` to keep them).

### 12.5 Under `tvb`: critical excitability \(x_{0,i}^{c}\) is compressed

In the n = 6 reference cohort all non-trivial nodes destabilise
within roughly \([-2.06, -1.0]\); many of the high-rank regions are
clustered within \(\approx 0.02\) of \(-2.057\), well below the
upper-bound `x0.upper = -1.0`. The thesis (Patient 1, Figure 4.2)
shows a comparatively wider per-node \(x_{0}^{c}\) spread under the
TVB pipeline. The compression observed here is consistent with the
TVB 95th-percentile truncation flattening the in-strength
distribution \(\sum_{j} K_{ji}\) — which by §4.6 of the thesis (and
Eq. 4.13) is the dominant input to \(C_{ii}\) and therefore to the
critical excitability. Practical implication: the rank order of
\(x_{0,i}^{c}\) under `tvb` is meaningful, but the absolute spread
should not be interpreted as biologically meaningful resilience
margins; for comparisons across cohorts, prefer rank-based statistics
(as the summary CSV's `mean_rank_x0` already does).

### 12.6 Honest caveats

- **Column normalisation is an imposed asymmetry.** DTI is undirected;
  the in/out distinction in the column scheme is a mathematical
  convention ("each node receives the same total upstream weight"), not
  biological directionality. The non-tautological part of the §12.2
  result is that the *specific nodes* \(D(\to i)\) flags under
  `column` overlap with the BC-flagged nodes under `parkes`
  (hippocampus, amygdala, etc.). That convergence on the same biology
  under two methodologically independent lenses is the actual evidence
  that \(D\) is detecting something real.
- **Tiny n.** Anderson \(n=2\), Arnold \(n=3\)–4. The strain
  \(D_{\mathrm{st}}\) separation needs replication on a larger cohort
  before being treated as a result rather than a hypothesis.
- **Welch t-test under tiny n.** The cohort \(D_{\mathrm{st}}\) chart
  prints a Welch \(t\) and \(p\) value automatically when both cohorts
  have \(n \geq 2\), but these numbers are largely cosmetic given the
  sample sizes involved.

---

## 14. Tradeoffs, limitations, and design decisions

### 13.1 Why a per-config experiment isolation rather than a single mega-config

Earlier iterations of the pipeline ran every scheme in a single flat
`results/` directory (this is what the legacy `section45_*` prefix in
`results/mouse_experiment_reference/` reflects). That layout made it
impossible to (a) reproduce a single experiment by re-running, and (b)
record which parameters belonged to which output. The
one-config-per-experiment, timestamped-folder design is a deliberate
reaction to that.

### 13.2 Why three normalisation schemes are shipped

Each of `column`, `parkes`, `tvb` answers a different methodological
question:

- `column` is the **primary scientific scheme** — it is the only one
  in which \(D\) can carry independent per-node signal on this data,
  and the only one in which \(D_{\mathrm{st}}\) separates the strains
  on the reference cohort.
- `parkes` is the **literature-comparison scheme** — it cleanly
  interfaces with classical (undirected) centrality work and is the
  regime where \(D\) reduces to EC. Running `parkes` alongside
  `column` is what surfaces the asymmetry effect.
- `tvb` is the **Epileptor reference scheme** — it matches the TVB
  pipeline used in the human Epileptor experiment, with the full
  per-node \(x_{0,i}^{c}\) sweep. It's slow (~10 min/mouse) but is the
  most directly model-driven of the three.

Shipping all three together as `initial_column.properties`,
`initial_parkes.properties`, and `initial_tvb.properties` lets the
qualitative comparison from §12 be reproduced with one
`runAllMouseExperiments` call.

### 13.3 Why `discreteTime = false` is locked

Both modes are mathematically supported by `covariancesGaussianNet`,
but the toolkit standard for stability analyses on physiological
networks is continuous-time OU, matching the linsync README and the
human Epileptor experiment. Mixing them in the same pipeline would
require careful annotation of every output; keeping `discreteTime`
fixed in `runParams` removes a category of footguns.

### 13.4 Why Bonferroni rather than FDR

Bonferroni is conservative and transparent. With \(n_{\text{Arn}} = 3\)
or 4, the Arnold cohort's per-node SD estimate is noisy enough that a
less conservative procedure (e.g. Benjamini–Hochberg FDR) would let in
nodes whose flagged status is sensitive to which single Arnold mouse is
included. The cost is missing some real effects; that's an acceptable
trade-off given the cohort size. The raw and Bonferroni-corrected
p-values are *both* saved in the outliers CSV so a user who wants to
re-do the multiple-testing correction can do so without re-running the
analysis.

### 13.5 Why \(D(k\to)\) is computed even when ≡ \(D(\to i)\)

For symmetric \(C\), \(\Omega = \Omega^{\top}\) and the two diagonals
coincide at machine precision. The pipeline still computes both, for
three reasons:

1. **Consistency check.** A mismatch above machine precision indicates a
   bug somewhere upstream (asymmetric numerical artefacts, mismatched
   `MaxK`, etc.).
2. **Schema stability.** Every per-mouse `.mat` always has both fields,
   regardless of scheme; downstream consumers don't need to special-case
   symmetric cases.
3. **Cheap insurance.** The second `covariancesGaussianNet` call is
   the same cost as the first; the savings would be negligible.

### 13.6 Why the per-node \(x_{0,i}^{c}\) sweep is skipped for linear schemes

The Epileptor solve is the source of the per-node Jacobian \(C\); under
`parkes` / `column`, \(C = K\) directly. There is no "single-node
seizure" picture because there is no Epileptor to perturb. The sweep is
returned as a vector of `NaN`s rather than an arbitrary stand-in so
that downstream consumers (CSV writers, plot panels, summary rankings)
visibly degrade rather than silently inventing values.

### 13.7 Why Anderson and Arnold are detected by name prefix

The directory layout and the mouse IDs themselves are the only
machine-readable signal in the input data that tells us which strain a
mouse belongs to. Detecting by `startsWith(string(mouseId), 'Anderson')`
keeps that mapping out of any config file, so adding a new Anderson or
Arnold animal is a matter of dropping the CSV into `data/<mouseId>/`
without touching any other configuration. Costs: any mouse whose name
doesn't start with `Anderson` or `Arnold` is silently *not classified*,
and the strain comparison is skipped if either prefix is absent.

### 13.8 Why config files use Java properties syntax

The MATLAB code parses configs with `java.util.Properties`, which
matches the way the parent linsync toolkit handles its own configs.
This is intentional cross-toolkit consistency; the format is also
trivially diffable and ships with both `#` comments and the
`key.with.dots = value` convention that maps cleanly to the
structured `runParams` field names.

### 13.9 What is NOT in scope for this experiment

- **No fitting of Epileptor parameters.** \(\tau_0\), \(x_0^{\text{base}}\),
  and the \(F\) coefficients are fixed to their literature defaults.
- **No subject-specific simulation timeseries.** The covariance is
  obtained analytically from \(C\); no time-domain integration is run.
- **No multi-scale (fine atlas) analysis.** The pipeline operates on the
  pre-coarsened 66-region matrix. The raw fine atlas `compact_idx_to_name.csv`
  and `connectome_compact_labelled.csv` are present in `data/<mouseId>/`
  but unused.
- **No effective / functional connectivity.** Inputs are structural
  (DTI). Whether \(D\)'s in/out asymmetry would carry biological
  meaning under a directed effective-connectivity input (transfer
  entropy, DCM, etc.) is the natural follow-up but is out of scope here.

---

## 15. Glossary of symbols

| Symbol                 | Meaning                                                                                                                  |
|------------------------|--------------------------------------------------------------------------------------------------------------------------|
| \(N\)                  | Number of regions (66 for the coarse mouse atlas)                                                                        |
| \(K\)                  | Normalised connectome matrix; \(K_{ij}\) is the directed weight from source \(i\) to target \(j\)                        |
| \(K_{\text{raw}}\)     | Raw connectome from the input CSV (before normalisation)                                                                 |
| \(C\)                  | Effective coupling matrix used for stability / covariance calculations (\(=K\) for linear schemes, Jacobian otherwise)   |
| \(\rho(C)\)            | Spectral radius of \(C\); \(\rho(C) < 1\) ⇔ healthy fixed point is linearly stable                                       |
| \(\Omega\)             | Stationary covariance of the linearised OU dynamics around the healthy fixed point                                       |
| \(D(\to i)\)           | Stability susceptibility centrality of node \(i\); \(= \Omega_{ii}\)                                                     |
| \(D(k \to)\)           | Stability influence centrality of node \(k\); \(= (\Omega_{C^{\top}})_{kk}\)                                              |
| \(D_{\mathrm{st}}\)    | Deviation from Stability; \(= \mathrm{trace}(\Omega) / N\)                                                               |
| \(x_{0,i}\)            | Per-node Epileptor excitability                                                                                          |
| \(x_0^{\text{base}}\)  | Healthy baseline excitability (default \(-2.3\))                                                                         |
| \(x_{0,i}^{c}\)        | Critical excitability — smallest \(x_{0,i}\) making \(\rho(C) \geq 1\) with all others held at \(x_0^{\text{base}}\)     |
| \(\tau_0\)             | 1-D Epileptor slow timescale (default \(6667\))                                                                          |
| \(F(z), G(z)\)         | Epileptor coupling and self-dynamics functions; see §3.1                                                                 |
| `parkes.c`, \(c\)      | Parkes normalisation constant; \(K = A / (\rho(A) + c)\)                                                                 |
| `col.scale`            | Target column sum for column normalisation (default \(0.95\))                                                            |
| `alpha`                | Family-wise level for Bonferroni outlier flagging (default \(0.05\))                                                     |
| `z.threshold`          | Display-only threshold in `reportAndersonOutliers`; does **not** control flagging                                        |
| \(\mathrm{LI}_{\mathrm{norm}}\) | Laterality index \((L-R)/(L+R)\) for homotopic pairs; see §9                                                      |
| \(\mathrm{LI}_{\mathrm{signed}}\) | Signed L–R difference \(L-R\) on a metric; see §9                                                              |
| BC, EC, KZ, PR, etc.   | Classical centralities; see §7                                                                                           |

---

## Where to go next

- The operational manual (how to install, run, list expected outputs):
  [`README.md`](README.md).
- The full empirical comparison of normalisation regimes, with numerical
  tables and the per-region biological interpretation:
  [`results/mouse_experiment_reference/centrality_comparison_analysis.md`](results/mouse_experiment_reference/centrality_comparison_analysis.md).
- The mathematical reference for \(\Omega\), \(D_{\mathrm{st}}\), and
  the OU / VAR dynamics: [linsync README §1.6](../README.md#16-deviation-from-stability-d_st).
- The human-cohort analogue of this pipeline (using patient connectomes
  and Epileptor-only): [`epileptor-experiment/README.md`](../epileptor-experiment/README.md).
