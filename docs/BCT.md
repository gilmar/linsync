# Brain Connectivity Toolbox (external dependency)

The [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/) (BCT) is **not** vendored in this repository. Install it locally if you run the mouse epilepsy experiment scripts that need classical graph centralities (`betweenness_wei`, `distance_wei`).

## When you need it

| Script | BCT required? |
|--------|----------------|
| `runMouseStabilityCentralities` | Optional — betweenness and closeness are `NaN` without BCT; stability centralities still run |
| `computeNetworkCentralities` | Yes, for `betweenness` and `closeness_*` fields |
| `compareCentralityMeasures` | Only if results were computed with BCT available |
| `compareAndersonVsArnold` | Only for the betweenness panel (third subplot) |
| `epileptor-experiment/` | No |

## Install (one-time)

1. Download **BCT** from the official site:
   - [bctnet — Getting started](https://sites.google.com/site/bctnet/getting-started), or
   - [NITRC BCT releases](https://www.nitrc.org/projects/bct) (e.g. **2019-03-03** `BCT.zip`).

   The [MathWorks File Exchange entry](https://www.mathworks.com/matlabcentral/fileexchange/61173-brain-connectivity-toolbox) is an optional per-machine install via the Add-On Manager; for this repo, extracting the zip to the path below is simpler and matches `setupMousePaths`.

2. Extract the archive so the `.m` files sit **directly** in this folder (not in a nested `BCT/` subfolder unless you adjust the path):

   ```
   linsync/2019_03_03_BCT/betweenness_wei.m
   linsync/2019_03_03_BCT/distance_wei.m
   ...
   ```

3. In MATLAB, from the mouse experiment:

   ```matlab
   cd('mouse-epilepsy-daria-anderson-experiment')
   setupMousePaths();
   assert(exist('betweenness_wei', 'file') == 2, 'BCT not on path');
   ```

If BCT is missing, `setupMousePaths` issues a **warning** (not an error). Re-run `runMouseStabilityCentralities` after installing BCT if you need full classical centrality fields in the saved results.

## Citation

When publishing work that uses BCT, cite:

> Rubinov M, Sporns O (2010). Complex network measures of brain connectivity: Uses and interpretations. *NeuroImage* 52:1059–1069.

See [bctnet](https://sites.google.com/site/bctnet/) for updates and measure documentation (`doc function_name` in MATLAB).
