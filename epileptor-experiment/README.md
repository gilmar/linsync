# Epileptor network experiment (1D epileptor + structural connectivity)

This folder is self-contained except for the parent **linsync** toolkit. Run scripts from MATLAB with the current folder set anywhere; each driver calls `setupEpileptorPaths()` to add:

- the **linsync** root (parent of this directory), which provides covariance / \(D_{\mathrm{st}}\) machinery plus shared Epileptor helpers (`oneDepileptor`, `CouplingMatrix`, `normal`, `findCriticalX0`, `computeStabilityCentralities`, etc.), and  
- this **epileptor-experiment** directory (data loaders and drivers).

## Prerequisites

- **MATLAB Optimization Toolbox** — required. The drivers solve the epileptor fixed point with **`fsolve`**, which is part of this toolbox.

## Layout

| Path | Role |
|------|------|
| `runEZ1.m` | Parameter sweep: fixed point → coupling matrix **C** → \(D_{\mathrm{st}}\), eigenvalues, optional transpose-Ω diagnostics |
| `runEZ1_stabilityCentralities.m` | Stability centralities for Patient 1 (thesis Fig. 4.2): per-node \(D(\to i)\) and \(D(k \to)\) from the healthy P1 network, per-node critical \(x^{c}_{0,i}\) sweep (\(\rho(C) \ge 1\)), two-panel scatter plot |
| `JL_playing.m` | Single baseline + virtual resections, \(D_{\mathrm{st}}\) and per-node Ω diagonals |
| `main1K12_JLadjusted.m` | Resection scenarios and critical \(x_0\) for EZ (fsolve failure boundary) |
| `loadPatientWeights.m`, `setupEpileptorPaths.m` | Data loading and path setup |
| `results/` | `.mat` outputs from the drivers (folder created by `setupEpileptorPaths`) |
| `data/connectivity_<patient>/weights.txt` | Structural weights (sample: `P1`) |

## Data

Place patient structural matrices under `data/connectivity_<ID>/weights.txt`. Optionally add `weights.txt` in this folder as a fallback for quick tests.

## See also

Repository [README.md](../README.md) §1.6 for \(D_{\mathrm{st}}\) definition.
