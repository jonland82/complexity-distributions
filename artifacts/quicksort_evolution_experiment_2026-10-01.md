# Quicksort distribution-evolution experiment

Date: October 1, 2026

## Question

Does the recursively generated family

$$
\nu_0=\delta_0,
\qquad
\nu_{k+1}=\mathcal T(\nu_k),
$$

where

$$
\mathcal T(\nu)
=\operatorname{Law}\left(
U X_1+(1-U)X_2+
\frac{1+2U\log U+2(1-U)\log(1-U)}{\sigma}
\right),
$$

provide useful intermediate approximations to the finite-size Quicksort distributions on its way to the Quicksort limit law $Q$?

Here $U$ is uniform on $(0,1)$, $X_1,X_2$ are independent with law $\nu$, and

$$
\sigma^2=7-\frac{2\pi^2}{3}.
$$

The objects compared in this experiment must be distinguished:

- $C_n$ is the finite-size Quicksort comparison count.
- $\widehat C_n=(C_n-\mu_n)/\sqrt{v_n}$ is its standardized finite-size shape.
- $Q$ is the continuous limiting law of $\widehat C_n$ as $n\to\infty$.
- $\nu_k$ is the recursive-generation family converging to $Q$; it is not assumed to equal a finite-size Quicksort law.

## Method

The experiment imported `quicksort_samples` and `exact_mean_variance` from `make_figures.py`.

- Recursive trajectory: 1,000,000 population values for each $k=1,\ldots,25$.
- Recursive seed: `20261001`.
- Independent terminal-$Q$ approximation: 1,000,000 values after 35 population iterations.
- Terminal-$Q$ seed: `20261002`.
- Finite-size Quicksort seed at size $n$: `20261010 + n`.
- Finite-size runs: 100,000 to 300,000, depending on $n$.
- Python 3.13.5, NumPy 2.1.3, and SciPy 1.15.3.

The theoretical variance of the recursive family is

$$
\operatorname{Var}(\nu_k)=1-\left(\frac23\right)^k.
$$

Each $\nu_k$ was divided by this theoretical standard deviation before its shape was compared with $\widehat C_n$. The variance-aligned generation was

$$
k_{\mathrm{var}}(n)
=\frac{\log(1-a_n/\sigma^2)}{\log(2/3)},
\qquad
a_n=\frac{v_n}{(n+1)^2}.
$$

The principal metric was two-distribution Kolmogorov distance. Wasserstein-1 distance, skewness, and excess kurtosis were also recorded. Kolmogorov and Wasserstein distances were evaluated on a dense union of 20,000 equal-probability quantiles. The displayed Monte Carlo threshold is the approximate 95% two-sample Kolmogorov threshold for the finite-size sample and the one-million-value recursive population.

The main run took 30.10 seconds after compilation. A separate exact-PMF check through $n=32$ took 8.6 seconds.

## Main results

`k_near` is the integer nearest to $k_var`. `k_best` minimizes the measured Kolmogorov distance among generations 1 through 25. At large $n$, several generations are closer than Monte Carlo resolution, so `k_best` is unstable and should not be interpreted literally.

| $n$ | Runs | $k_{\rm var}$ | `k_near` | `k_best` | KS near | W1 near | KS best | W1 at best | KS to $Q$ | W1 to $Q$ | KS to normal | 95% KS noise |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 300,000 | 1.434 | 1 | 2 | 0.125513 | 0.145009 | 0.048890 | 0.050224 | 0.060828 | 0.060696 | 0.109577 | 0.002831 |
| 32 | 300,000 | 2.402 | 2 | 3 | 0.027108 | 0.037750 | 0.025364 | 0.034814 | 0.032651 | 0.040090 | 0.087348 | 0.002831 |
| 64 | 300,000 | 3.522 | 4 | 3 | 0.015638 | 0.025175 | 0.013942 | 0.024243 | 0.018554 | 0.026023 | 0.074722 | 0.002831 |
| 128 | 300,000 | 4.751 | 5 | 4 | 0.009290 | 0.017277 | 0.008390 | 0.017192 | 0.011245 | 0.017611 | 0.067210 | 0.002831 |
| 256 | 300,000 | 6.057 | 6 | 4 | 0.006262 | 0.010637 | 0.005944 | 0.012483 | 0.007407 | 0.011376 | 0.064236 | 0.002831 |
| 512 | 250,000 | 7.420 | 7 | 6 | 0.004714 | 0.008089 | 0.003801 | 0.007312 | 0.004930 | 0.007366 | 0.061466 | 0.003041 |
| 1,024 | 200,000 | 8.826 | 9 | 25 | 0.001841 | 0.002746 | 0.001629 | 0.002307 | 0.002120 | 0.003383 | 0.058954 | 0.003331 |
| 2,048 | 150,000 | 10.266 | 10 | 7 | 0.001718 | 0.003329 | 0.001618 | 0.002706 | 0.002335 | 0.003878 | 0.057213 | 0.003766 |
| 4,096 | 100,000 | 11.732 | 12 | 12 | 0.001425 | 0.001960 | 0.001425 | 0.001960 | 0.002225 | 0.002564 | 0.057614 | 0.004511 |

The independent terminal-$Q$ population had estimated skewness `0.858617` and excess kurtosis `1.209753`. The theoretical skewness of $Q$ is approximately `0.854882`.

## Moment comparison at the variance-aligned generation

| $n$ | Skewness of $\widehat C_n$ | Skewness of $\nu_{k_{near}}$ | Excess kurtosis of $\widehat C_n$ | Excess kurtosis of $\nu_{k_{near}}$ |
|---:|---:|---:|---:|---:|
| 16 | 1.071337 | 0.883263 | 1.401910 | -0.328380 |
| 32 | 1.037097 | 0.824728 | 1.473793 | 0.401794 |
| 64 | 0.984413 | 0.823571 | 1.427045 | 0.848482 |
| 128 | 0.952007 | 0.826952 | 1.381502 | 0.932715 |
| 256 | 0.912736 | 0.836767 | 1.284831 | 1.017381 |
| 512 | 0.891978 | 0.836397 | 1.243603 | 1.058607 |
| 1,024 | 0.871119 | 0.853257 | 1.187311 | 1.150599 |
| 2,048 | 0.854755 | 0.841355 | 1.134714 | 1.102227 |
| 4,096 | 0.859928 | 0.848226 | 1.189250 | 1.139616 |

These moment estimates show that the recursive family captures the main CDF shape faster than it captures finite-size skewness and kurtosis. The discrepancy is material at small and moderate sizes and becomes small at larger sizes.

## Exact finite-size checks

Finite-size PMFs were generated through the exact Quicksort recurrence using floating-point convolution. The target PMF was therefore exact up to floating-point roundoff; the recursive populations remained Monte Carlo approximations.

| $n$ | Nonzero PMF support | $k_{\rm var}$ | `k_near` | KS near | W1 near | `k_best` | KS best | W1 at best | Exact skewness | Exact excess kurtosis |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 83 | 1.434 | 1 | 0.125630 | 0.144534 | 2 | 0.048770 | 0.050179 | 1.077339 | 1.420605 |
| 32 | 394 | 2.402 | 2 | 0.027972 | 0.037284 | 3 | 0.023332 | 0.033866 | 1.041394 | 1.498077 |

These exact checks confirm that the early-size findings are not artifacts of finite-size Quicksort simulation.

## Interpretation

The experiment does **not** show that an intermediate $\nu_k$ approximates $Q$ better than $Q$ approximates itself. The comparison is

$$
d(\widehat C_n,\nu_k)
\quad\text{versus}\quad
d(\widehat C_n,Q).
$$

An intermediate recursive distribution can be closer to a finite-size Quicksort distribution because the latter has not yet reached its limiting shape.

The evidence supports the following conclusions:

1. The recursive family is a useful model of the evolution toward $Q$.
2. At sizes 16 through roughly 512, a suitable intermediate generation usually improves the Kolmogorov fit relative to using terminal $Q$ prematurely. The Wasserstein improvement is smaller and is essentially gone by size 512.
3. The recursive model is dramatically closer than a normal distribution. At $n=4096$, its variance-aligned KS distance was approximately `0.0014`, versus approximately `0.0576` for the normal.
4. Integer recursive generations are too coarse at small sizes. At $n=16$, variance alignment selects generation 1, while generation 2 is much closer. A fractional-generation or continuous-time interpolation is therefore worth investigating.
5. At $n\ge1024$, differences among nearby recursive generations and terminal $Q$ fall below the Monte Carlo resolution of this run. Reported best-generation values in this regime are not statistically identifiable.
6. The recursive family does not exactly reproduce finite-size Quicksort. It under-represents finite-size skewness and kurtosis at small and moderate sizes, even when its central CDF is close.

## Status and next step

This is positive empirical evidence for a Quicksort-specific evolving family, not a proof that $\nu_{k(n)}$ equals or asymptotically gives the best approximation to the finite-size law.

The next focused step is to construct a continuous or fractional interpolation between $\nu_k$ and $\nu_{k+1}$, align its time parameter with $n$, and determine whether it improves both the CDF and the finite-size skewness and kurtosis. A diffusion model should be evaluated only after this interpolation is defined, using the exact finite-size recurrence or coupled tree process as the reference.
