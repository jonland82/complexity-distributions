# Fractional Quicksort distribution evolution

Date: October 1, 2026

## Question and method

Can interpolation between adjacent laws in the recursive family improve its approximation to finite-size Quicksort?

For integer $k\geq1$ and $0\leq\alpha\leq1$, define a mixture of the *unstandardized* recursive laws,

$$
\nu_{k+\alpha}^{\rm mix}=(1-\alpha)\nu_k+\alpha\nu_{k+1},
\qquad
v_{k+\alpha}=(1-\alpha)\left[1-(2/3)^k\right]
 +\alpha\left[1-(2/3)^{k+1}\right].
$$

The comparison model is $\nu_{k+\alpha}^{\rm mix}$ divided by $\sqrt{v_{k+\alpha}}$. This is a continuous interpolation of one-time distributions. It is not a continuous-time Quicksort process or a diffusion.

The experiment used the exact finite-size Quicksort split recurrence, evaluated by floating-point FFT convolution, for $n=16,32,64,128,256,512$. The $n=1024$ target used 400,000 training and 400,000 independent evaluation runs. Independent recursive populations of one million values per generation were generated with seeds `20261020` and `20261021`; the first was used to choose the time and the second to evaluate it. The terminal-$Q$ comparator was generation 35. The CDF-fit time minimizes training Kolmogorov distance on a grid with 0.05-generation spacing from 1 to 14. The best integer generation is also selected on training data. All displayed comparisons use the independent evaluation population. For exact finite-size targets, the Kolmogorov calculation includes both sides of each probability jump. The reported Wasserstein-1 values integrate the absolute CDF gap numerically over $[-6,12]$.

Reproduce with `python quicksort_fractional_experiment.py`. The script also writes the complete numerical results to `quicksort_fractional_results.json`. The local run took 116 seconds, including 96 seconds for exact PMFs through $n=512$.

## CDF results

All distances below are Kolmogorov distances. The integer comparator is the best generation chosen on the training population. Rows through $n=512$ use exact finite-size target PMFs up to floating-point roundoff; only $n=1024$ uses a simulated target.

| $n$ | CDF-fit time | Best integer | Fitted mixture | Best integer | Terminal $Q$ | Variance-aligned mixture |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 1.90 | 2 | 0.03954 | 0.04909 | 0.06060 | 0.08340 |
| 32 | 2.75 | 3 | 0.02207 | 0.02312 | 0.03110 | 0.02312 |
| 64 | 2.90 | 3 | 0.01302 | 0.01327 | 0.01767 | 0.01300 |
| 128 | 4.15 | 4 | 0.00774 | 0.00767 | 0.01030 | 0.00821 |
| 256 | 4.70 | 5 | 0.00483 | 0.00500 | 0.00607 | 0.00511 |
| 512 | 6.55 | 7 | 0.00294 | 0.00290 | 0.00357 | 0.00286 |
| 1,024 (simulated target) | 11.00 | 11 | 0.00309 | 0.00309 | 0.00299 | 0.00275 |

At $n=16$, interpolation reduces the held-out CDF gap by 0.00955, or 19%, relative to the best integer generation. The much larger variance-aligned gap shows that variance is a poor rule for choosing time at this size. At $n=32$, the fractional gain is 0.00105. At $n\geq64$, the fitted mixture and best integer are separated by at most 0.00025 in this run, well within the recursive population's Monte Carlo uncertainty. At $n=1024$, the fitted time is unstable: its CDF gap is comparable with the roughly 0.0025 two-sample 95% Kolmogorov noise scale for 400,000 versus one million observations.

For Wasserstein-1, fitted mixture versus best integer is 0.04971 versus 0.05055 at $n=16$, and 0.03289 versus 0.03402 at $n=32$. At $n=256$, terminal $Q$ has the smaller Wasserstein-1 value (0.00994 versus 0.01070 for the fitted mixture). The CDF improvement therefore does not consistently improve this second distance.

## Shape and tail results

| $n$ | Target skew | Mixture skew | Target excess kurtosis | Mixture excess kurtosis | Target $P(Z>3)$ | Mixture $P(Z>3)$ |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 1.0773 | 0.8304 | 1.4206 | 0.4170 | 0.01176 | 0.00791 |
| 32 | 1.0414 | 0.8234 | 1.4981 | 0.6704 | 0.01071 | 0.00808 |
| 64 | 0.9894 | 0.8215 | 1.4321 | 0.6836 | 0.01047 | 0.00819 |
| 128 | 0.9451 | 0.8261 | 1.3526 | 0.8635 | 0.00994 | 0.00866 |
| 256 | 0.9126 | 0.8274 | 1.2901 | 0.9084 | 0.00966 | 0.00872 |
| 512 | 0.8906 | 0.8415 | 1.2471 | 1.0599 | 0.00944 | 0.00892 |

The fitted interpolation understates skewness, kurtosis, and the probability of costs more than three standard deviations above the mean. At $n=16$, it predicts 0.79% above that threshold versus the exact 1.18%; at $n=256$, 0.87% versus 0.97%. The mismatch shrinks with $n$, but tuning the interpolation time alone does not reproduce the early finite-size shape.

## Interpretation

Fractional mixing is useful mainly at the smallest sizes. It smooths the coarse jump between generations and gives a clear CDF gain at $n=16$. It does not resolve the missing finite-size skewness and far-right tail. The next model change should account for finite-size split effects or use the exact size-indexed tree/recurrence evolution; finer tuning of the same mixture is unlikely to repair the shape discrepancy. The advantage of individual intermediate generations over terminal $Q$ at larger sizes requires more precise recursive populations before drawing quantitative conclusions.
