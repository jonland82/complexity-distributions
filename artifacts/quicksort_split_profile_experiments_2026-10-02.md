# Finite-split and profile experiments

Date: October 2, 2026.

## One-step finite-size split model

The previous experiment found that interpolating recursive generations improved the central CDF but understated finite-size skewness and the right tail. Here the correction is the exact finite-size split rule, with each child law replaced by a moment-matched approximation to the terminal Quicksort law $Q$:

\[
\widetilde C_n=n-1+\mu_I+\sqrt{v_I}Q_1+\mu_{n-1-I}+\sqrt{v_{n-1-I}}Q_2,
\qquad I\sim\operatorname{Uniform}\{0,\ldots,n-1\}.
\]

For child sizes below 3, the law is deterministic. A second variant uses exact child PMFs through size 8 and $Q$ above size 8. Both models preserve the exact finite-size mean and variance. They are one-step approximations, not exact Quicksort laws or diffusion processes.

Finite-size targets are PMFs from the exact split recurrence through $n=256$, calculated by floating-point FFT convolution. Their means and variances agree with the closed formulas to numerical precision. The direct recursive baseline uses the previously fitted fractional times, evaluated here with fresh recursive populations. Two independent seeds each generated a one-million-value terminal population and two million model draws per size. The table averages the two measured distances; the target columns are exact up to FFT roundoff.

| $n$ | Direct fractional KS | Split with $Q$ children KS | Split with exact children through 8 KS | Exact skew | $Q$-child skew | Exact $P(Z>3)$ | $Q$-child $P(Z>3)$ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 0.03940 | 0.03985 | **0.03290** | 1.077 | 1.004 | 0.01176 | 0.00989 |
| 32 | 0.02199 | **0.01862** | 0.01881 | 1.041 | 0.959 | 0.01071 | 0.00984 |
| 64 | 0.01304 | 0.01051 | **0.01022** | 0.989 | 0.921 | 0.01047 | 0.00957 |
| 128 | 0.00846 | **0.00602** | 0.00616 | 0.945 | 0.896 | 0.00994 | 0.00949 |
| 256 | 0.00497 | **0.00354** | 0.00397 | 0.913 | 0.879 | 0.00966 | 0.00934 |

The split rule improves the CDF at sizes 32 through 256; exact small-child laws give the largest gain at size 16. The $Q$-child model also moves skew and the $Z>3$ tail toward the exact values, but still understates both. Its excess kurtosis versus the target is about 1.23 versus 1.42 at size 16 and 1.21 versus 1.29 at size 256. Pure $Q$ children do not improve KS at size 16. The size-256 KS difference is small relative to Monte Carlo resolution, although its direction repeated across seeds. Approximate 95% Dvoretzky–Kiefer–Wolfowitz uncertainty from two million draws alone is 0.00096 per sampled model CDF; terminal-population error adds uncertainty.

## Profile dependence of next-step variance

In 100,000 independent uniform external-slot insertion histories, the exact next-step martingale increment variance was recorded at sizes 16, 64, and 256. Half the histories fitted a quadratic function of standardized cost; the other half evaluated how much of the across-tree variance this scalar predictor explained.

| $n$ | Mean next-step variance | Coefficient of variation across trees | Held-out $R^2$ from cost |
|---:|---:|---:|---:|
| 16 | 0.007230 | 0.589 | 0.862 |
| 64 | 0.001077 | 0.440 | 0.767 |
| 256 | 0.000110 | 0.332 | 0.716 |

The simulated means agree with the exact innovation-variance identity: 0.007230, 0.001074, and 0.000110, respectively. Cost explains much, but leaves 14% to 28% of the across-tree variance unexplained in this test. This supports assessing a scalar variance approximation empirically; it does not make the coupled cost process Markov or validate a diffusion approximation.

## Interpretation and next experiment

The finite-size split is a useful structural correction. The remaining tail and moment gap is likely influenced by replacing the child laws with $Q$, especially for moderate child sizes. The next controlled comparison should use size-dependent child approximations and test a mapping from child size to recursive generation on held-out sizes. Separately, fit a scalar variance function to profile data and evaluate its one-time CDFs **and joint increments** against the full-profile jump process before treating it as a diffusion model.

Reproduce with `python quicksort_split_experiment.py`, `python quicksort_split_experiment.py --seed 20261042 --draws 2000000`, and `python quicksort_profile_experiment.py`. Full numeric outputs are the adjacent JSON files.
