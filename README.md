# Deriving Classical Distributions from Normal Sampling

An interactive Shiny app for BIOL 3P96 (Biostatistics) at Brock University.

## What this app does

The test statistics used throughout the course do not follow arbitrary
distributions: each one arises naturally from random sampling. This app
simulates that process and lets you watch each distribution emerge from first
principles, then compares it with the theoretical (calculated) distribution.

- **t distribution**: standardised sample means from a normal population;
  approaches the standard normal as sample size grows
- **F distribution**: ratio of two independent sample variances;
  shaped by the degrees of freedom from each group
- **Chi-square distribution**: sum of squared standard normal values;
  connects to variance estimation and goodness-of-fit tests
- **Binomial distribution**: number of successes in *n* independent trials,
  each with probability *p*; compared with the exact probabilities from
  `dbinom()`
- **Pearson's correlation (r)**: correlation between unrelated *X* and *Y*,
  converted to the *t* statistic that `cor.test()` uses and compared with the
  *t* distribution on *n* − 2 degrees of freedom
- **Spearman's rank correlation (r<sub>S</sub>)**: the correlation of the ranks
  of unrelated *X* and *Y*; shows how far *r<sub>S</sub>* wanders from zero by
  chance alone
- **Mann-Whitney U**: calculated as in Whitlock & Schluter (rank both groups
  together, find U<sub>1</sub> and U<sub>2</sub>, and keep the larger);
  compared with the exact distribution of *U*

In every simulation the null hypothesis is true, so each histogram is a null
distribution.

## Critical values

Each plot shows the critical value(s) at α = 0.05 in two ways:

- **Grey dashed lines**: the calculated critical values (e.g. from `qt()`,
  `qf()`, `qchisq()`, or the exact binomial and *U* distributions). These change
  only when the sliders change the degrees of freedom or sample sizes.
- **Orange dotted lines**: the same cut-offs found directly from the simulated
  values (the values that leave 5% of the simulations beyond them). These
  bounce around the grey lines each time you resample, and settle closer to
  them as the number of simulations increases.

Because the null hypothesis is true in every simulation, about 5% of simulated
values fall beyond the critical values: these are the false positives. For the
binomial and *U*, whose values are whole numbers, the proportion is at most 5%.

## How to use

1. Choose a distribution from the dropdown menu.
2. Adjust the sample size(s), probability or number of simulations with the
   sliders.
3. Press **Resample** to generate a new batch of simulated values. The previous
   batch is drawn underneath in grey, so you can see how much the histogram
   shifts by chance.
4. Compare the histogram with the red theoretical curve or points, and the
   orange simulated critical values with the grey calculated ones.
5. Read the R code beneath the plot: it shows how the values were simulated and
   how to draw the histogram and critical values yourself.
6. Choose higher number of simulations to observe how the empirically drawn
   better approximate the theoretical values.

## Learning goals

- Understand that these distributions are not "made up": they follow
  mathematically from random sampling
- See how degrees of freedom and sample size control the shape of each
  distribution
- Observe the *t* distribution converging to normal as *n* increases
- Understand where critical values come from, and that about 5% of results
  exceed them by chance when the null hypothesis is true
- Connect each distribution to the test that uses it: *t*-tests, ANOVA,
  chi-square tests, the binomial test, correlation tests and the Mann-Whitney
  U test

## Course context

Developed for BIOL 3P96 (Biostatistics), Brock University.
Built with R and Shiny (base R graphics only).