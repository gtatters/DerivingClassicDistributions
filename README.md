# Deriving Classical Distributions from Normal Sampling

An interactive Shiny app for BIOL 3P96 (Biostatistics) at Brock University.

## What this app does

The t, F, and chi-square distributions are not arbitrary — they arise naturally
from sampling normal populations. This app simulates that process and lets you
watch each distribution emerge from first principles.

- **t distribution** — standardised sample means from a normal population;
  approaches the standard normal as sample size grows
- **F distribution** — ratio of two independent sample variances;
  shaped by the degrees of freedom from each group
- **Chi-square distribution** — sum of squared standard normal values;
  connects to variance estimation and goodness-of-fit tests

Press **Resample** to draw a new set of experiments. A previous-sample overlay
helps you see how much the histogram shifts by chance.

## How to use

1. Choose a distribution from the dropdown menu.
2. Adjust the sample size(s) with the sliders.
3. Press **Resample** to generate a new batch of simulated values.
4. Compare the histogram to the red theoretical curve and note how closely
   they match — and how the shape changes with degrees of freedom.

## Learning goals

- Understand that these distributions are not "made up" — they follow
  mathematically from normal random sampling
- See how degrees of freedom control the shape of each distribution
- Observe the t distribution converging to normal as n increases

## Course context

Developed for BIOL 3P96 — Biostatistics, Brock University.
Built with R and Shiny (base R graphics only).