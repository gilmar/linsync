# Figures from Lizier et al. 2023

This folder contains scripts / parameters files to recreate the results of:

J.T. Lizier, F.M. Atay, and J. Jost,
_"TBA"_,
2023

## Experiment 1 - convergence of numerical results from simulations - discrete time, no delay

The parameter settings for recreating experiment 1 are contained in `exp1Parameters.m`.

The script to run the experiments is `exp1ConvergenceOfSimulations.m`.
Run that script in Matlab from this folder:
```matlab
figure1ConvergenceOfSimulations
```

An `.eps` and `.mat` file is saved for the figure.

*Crucially* - note that the script adjusts the number of repeat runs 
downwards from 2000 in the parameters file (matching the paper) to 10 -- this
is done to ensure that the results run quickly for you. 
Also, the script removes the largest number of time series samples (L=1,000,000)
for the same reason.
These shorter runs take about 3 minutes 
on my machine.

Obviously you can remove those restrictions -- running for 2000 repeats and with more time series samples
will take about 100 hours (the longest number of time series samples really 
blows this out)
so you would be advised to run for more repeats in parallel fashion on a cluster (instructions below).

In any case, the shorter runs give a good idea of the main trends already,
though there are some fluctuations / std error differences compared to the longer runs.

To run this on a cluster, first set up your cluster environment as per the
[cluster instructions](/cluster).

Then following the basic steps outlined there to get an experiment running:
1. Copy `exp1Parameters.m` here over the top of `cluster/parameters.m`,
1. From your shell environment run `runManyProcesses.sh 1 1 200`
and then follow the other instructions at the [cluster instructions](/cluster)
page to combine the results together etc.

## Experiment 2 - convergence of numerical results from simulations - discrete time, fixed delays

The parameter settings for recreating experiment 2 are contained in `exp2Parameters.m`.

The script to run the experiments is `exp2ConvergenceOfSimulations.m`.
Run that script in Matlab from this folder:
```matlab
figure2ConvergenceOfSimulations
```

An `.eps` and `.mat` file is saved for the figure.

*Crucially* - note that the script adjusts the number of repeat runs 
downwards from 2000 in the parameters file (matching the paper) to 10 -- this
is done to ensure that the results run faster for you.
Also, the script removes the largest number of time series samples (L=1,000,000)
for the same reason.
These shorter runs take about 45 minutes on my machine.

Obviously you can remove those restrictions -- running for 2000 repeats and with more time series samples
will take about 800 hours (the longest number of time series samples really 
blows this out)
so you would be advised to run for more repeats in parallel fashion on a cluster (instructions below).

In any case, the shorter runs give a good idea of the main trends already,
though there are some fluctuations / std error differences compared to the longer runs.

To run this on a cluster, first set up your cluster environment as per the
[cluster instructions](/cluster).

Then following the basic steps outlined there to get an experiment running:
1. Copy `exp2Parameters.m` here over the top of `cluster/parameters.m`,
1. From your shell environment run `runManyProcesses.sh 1 1 200`
and then follow the other instructions at the [cluster instructions](/cluster)
page to combine the results together etc.

For the saved plot here I altered the use of `jet` palette to `hsv` in .

