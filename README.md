# odin-dust-plots
Code to recreate plots in 'Reproducible parallel inference and simulation of stochastic state space models using odin, dust, and mcstate'

## Figure 1
Runs a speedup test of four models across a number of cores.

Run `figure1.R` to reproduce.

This file contains all the code needed to produce timings, but also
includes the data we observed on a test machine used to reproduce the
exact plot.

## Figure 2
Drawn by hand, no data.

## Figure 3
Runs the volatility model with 10 particles, highlights a single
trajectory, and draws the average with a loess curve.

Run `figure3.R` to reproduce.

## Figures 4 and 5
Runs the SIR model with beta = 0.2, gamma = 0.1 (and therefore R0 = 2)
against case incidence data simulated from the same model. This code is
identical to that in the mcstate SIR models vignette, with some changes
to the plotting functions. Figure 4 plots are simulation runs, a run of the
particle filter, fit to incidence data, and a fit plus forecasting.
Figure 5 is a plot of pMCMC inference of beta, gamma and R0 over four
chains.

Run `figure4-figure5.R` to reproduce. Figure 4 comes as four panels
which were collated and labeled manually (adding a forecast/nowcast
divide in panel D at t = 50). The legend of figure 5 was enlarged and
panels labelled manually.

## Table 1
These were timed using `/usr/bin/time`. The equivalent pomp and libBi
models we used can be found in the corresponding subdirectories.
