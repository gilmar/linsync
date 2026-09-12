# Brain Connectivity Toolbox (external dependency)

The [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/) (BCT) is **not** bundled with
this repository. Install it locally if you want the classical path-based graph centralities
(`betweenness_wei`, `distance_wei`) in the connectome pipeline.

BCT is optional. Without it, the stability centralities, $D_{\mathrm{st}}$, the critical-excitability
sweep and every cohort comparison still run; only the betweenness and closeness fields come back as
`NaN`, and `computeNetworkCentralities` issues a one-time warning.

## What needs it

| Function | BCT required? |
|----------|---------------|
| `computeNetworkCentralities` | for the `betweenness` and `closeness_*` fields only |
| `runStabilityCentralities` | no — it calls the above and tolerates `NaN` |
| `compareCentralityMeasures` | only to include betweenness/closeness in the correlation matrix |
| `compareCohortGroups` | only for the betweenness metric; it is dropped automatically when absent |
| everything else | no |

## Install (one-time)

1. Download BCT from the official site:
   - [bctnet — Getting started](https://sites.google.com/site/bctnet/getting-started), or
   - [NITRC BCT releases](https://www.nitrc.org/projects/bct) (e.g. the **2019-03-03** `BCT.zip`).

2. Extract it so the `.m` files sit **directly** inside a folder whose name contains `BCT`, beside the
   `linsync` root or beside your study folder:

   ```
   linsync/2019_03_03_BCT/betweenness_wei.m
   linsync/2019_03_03_BCT/distance_wei.m
   ...
   ```

   `setupConnectomePaths` looks for a `*BCT*` folder in, and next to, both the toolkit root and the
   project root. Any other location works too — pass it explicitly:

   ```matlab
   setupConnectomePaths('BctRoot', '/path/to/BCT');
   ```

3. Check it resolved:

   ```matlab
   info = setupConnectomePaths();
   assert(info.hasBct, 'BCT is not on the MATLAB path');
   ```

If BCT is missing, `setupConnectomePaths` issues a **warning**, not an error. Re-run the stability
step after installing it if you need the path-based centralities in the saved results — they are
computed once, when each subject is analysed, and stored in the per-subject `.mat`.

## Citation

When publishing work that uses BCT, cite:

> Rubinov M, Sporns O (2010). Complex network measures of brain connectivity: Uses and
> interpretations. *NeuroImage* 52:1059–1069.

See [bctnet](https://sites.google.com/site/bctnet/) for updates and per-measure documentation
(`doc <function name>` in MATLAB).
