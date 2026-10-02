"""Compare a continuous mixture of Quicksort recursion generations with C_n.

Run: python quicksort_fractional_experiment.py
The exact finite-size PMFs are computed by the Quicksort split recurrence.
Recursive laws are Monte Carlo populations; train and evaluation populations
use independent random seeds.
"""

from __future__ import annotations

import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.signal import fftconvolve

from make_figures import exact_mean_variance, quicksort_samples


OUT = Path(__file__).resolve().parent / "artifacts" / "quicksort_fractional_results.json"
SIGMA = np.sqrt(7 - 2 * np.pi**2 / 3)
SIZES = (16, 32, 64, 128, 256, 512, 1024)
EXACT_MAX = 512
POPULATION = 1_000_000
FINITE_RUNS = 400_000
LAST_GENERATION = 14
Q_GENERATION = 35


def exact_pmfs(max_n: int):
    """Floating-point PMFs from the exact finite-size split recurrence."""
    laws = [(0, np.array([1.0])), (0, np.array([1.0]))]
    for n in range(2, max_n + 1):
        low = n - 1 + min(laws[j][0] + laws[n - 1 - j][0] for j in range(n))
        high = n * (n - 1) // 2
        pmf = np.zeros(high - low + 1)
        for j in range((n + 1) // 2):
            left_low, left = laws[j]
            right_low, right = laws[n - 1 - j]
            convolution = fftconvolve(left, right)
            offset = n - 1 + left_low + right_low - low
            multiplier = 1 if j == n - 1 - j else 2
            pmf[offset : offset + convolution.size] += multiplier * convolution / n
        if np.min(pmf) < -1e-11:
            raise ValueError(f"Large negative recurrence probability at n={n}")
        pmf = np.maximum(pmf, 0)
        pmf /= pmf.sum()
        laws.append((low, pmf))
    return laws


def population_trajectory(seed: int):
    rng = np.random.default_rng(seed)
    pool = np.zeros(POPULATION)
    laws = {}
    for k in range(1, Q_GENERATION + 1):
        u = rng.random(POPULATION)
        left = pool[rng.integers(POPULATION, size=POPULATION)]
        right = pool[rng.integers(POPULATION, size=POPULATION)]
        toll = 1 + 2 * (u * np.log(u) + (1 - u) * np.log1p(-u))
        pool = u * left + (1 - u) * right + toll / SIGMA
        pool -= pool.mean()
        if k <= LAST_GENERATION or k == Q_GENERATION:
            # The population noise in its variance is removed; the recursion
            # itself always uses the unadjusted pool from the previous step.
            adjusted = pool * np.sqrt((1 - (2 / 3) ** k) / pool.var())
            laws[k] = np.sort(adjusted)
    return laws


def standard_target(values, probabilities, n):
    mean, variance = exact_mean_variance(n)
    z = (np.asarray(values, dtype=np.float64) - mean) / np.sqrt(variance)
    p = np.asarray(probabilities, dtype=np.float64)
    return z, p / p.sum()


def simulated_target(n, repetitions, seed):
    counts = quicksort_samples(n, repetitions, seed)
    values, frequencies = np.unique(counts, return_counts=True)
    return standard_target(values, frequencies, n)


def theoretical_variance(k):
    return 1 - (2 / 3) ** k


def split_time(t):
    k = min(int(np.floor(t)), LAST_GENERATION - 1)
    return k, float(t - k)


def model_cdf(x, t, laws):
    k, alpha = split_time(t)
    variance = (1 - alpha) * theoretical_variance(k) + alpha * theoretical_variance(k + 1)
    scaled = np.asarray(x) * np.sqrt(variance)
    first = np.searchsorted(laws[k], scaled, side="right") / POPULATION
    second = np.searchsorted(laws[k + 1], scaled, side="right") / POPULATION
    return (1 - alpha) * first + alpha * second


def integer_cdf(x, k, laws):
    return np.searchsorted(laws[k], np.asarray(x) * np.sqrt(theoretical_variance(k)), side="right") / POPULATION


def q_cdf(x, laws):
    return np.searchsorted(laws[Q_GENERATION], x, side="right") / POPULATION


def ks(z, p, model):
    observed = np.cumsum(p)
    predicted = model(z)
    return float(max(np.max(np.abs(observed - predicted)),
                     np.max(np.abs(observed - p - predicted))))


def w1(z, p, model):
    grid = np.linspace(-6, 12, 9001)
    target_cdf = np.concatenate(([0.0], np.cumsum(p)))[np.searchsorted(z, grid, side="right")]
    return float(np.trapezoid(np.abs(target_cdf - model(grid)), grid))


def shape(z, p):
    mean = np.dot(p, z)
    centered = z - mean
    variance = np.dot(p, centered**2)
    return (float(np.dot(p, centered**3) / variance**1.5),
            float(np.dot(p, centered**4) / variance**2 - 3),
            float(p[z > 2].sum()),
            float(p[z > 3].sum()))


def model_shape(t, laws):
    k, alpha = split_time(t)
    variance = (1 - alpha) * theoretical_variance(k) + alpha * theoretical_variance(k + 1)
    first = laws[k] / np.sqrt(variance)
    second = laws[k + 1] / np.sqrt(variance)
    # Both populations have been centered and assigned their theoretical variance.
    m3 = (1 - alpha) * np.mean(first**3) + alpha * np.mean(second**3)
    m4 = (1 - alpha) * np.mean(first**4) + alpha * np.mean(second**4)
    tails = 1 - model_cdf(np.array([2.0, 3.0]), t, laws)
    return float(m3), float(m4 - 3), float(tails[0]), float(tails[1])


def fit_cdf_time(z, p, laws):
    times = np.round(np.arange(1.0, LAST_GENERATION + 0.0001, 0.05), 2)
    scores = np.array([ks(z, p, lambda x, t=t: model_cdf(x, t, laws)) for t in times])
    winner = int(np.argmin(scores))
    return float(times[winner]), float(scores[winner])


def evaluate(z, p, t, integer_k, q_laws, shape_laws):
    mixture = lambda x: model_cdf(x, t, q_laws)
    integer = lambda x: integer_cdf(x, integer_k, q_laws)
    terminal = lambda x: q_cdf(x, q_laws)
    return {
        "ks_mixture": ks(z, p, mixture),
        "ks_integer": ks(z, p, integer),
        "ks_q": ks(z, p, terminal),
        "w1_mixture": w1(z, p, mixture),
        "w1_integer": w1(z, p, integer),
        "w1_q": w1(z, p, terminal),
        "target_shape": shape(z, p),
        "mixture_shape": model_shape(t, shape_laws),
        "integer_shape": model_shape(float(integer_k), shape_laws),
    }


def main():
    started = perf_counter()
    pmfs = exact_pmfs(EXACT_MAX)
    print(f"Exact PMFs through n={EXACT_MAX}: {perf_counter() - started:.1f}s", flush=True)
    train_laws = population_trajectory(20261020)
    test_laws = population_trajectory(20261021)
    print(f"Independent recursive populations: {perf_counter() - started:.1f}s", flush=True)
    rows = []
    for n in SIZES:
        if n <= EXACT_MAX:
            low, p = pmfs[n]
            train_target = standard_target(np.arange(low, low + p.size), p, n)
            test_target = train_target
            target_type = "exact PMF"
        else:
            train_target = simulated_target(n, FINITE_RUNS, 20261030 + n)
            test_target = simulated_target(n, FINITE_RUNS, 20261040 + n)
            target_type = "independent simulations"
        mean, variance = exact_mean_variance(n)
        ratio = variance / ((n + 1) ** 2 * SIGMA**2)
        k_var = np.log1p(-ratio) / np.log(2 / 3)
        k = max(1, min(int(np.floor(k_var)), LAST_GENERATION - 1))
        alpha_var = (ratio - theoretical_variance(k)) / (theoretical_variance(k + 1) - theoretical_variance(k))
        t_var = k + alpha_var
        t_best, train_score = fit_cdf_time(*train_target, train_laws)
        integer_scores = [ks(*train_target, lambda x, k=k: integer_cdf(x, k, train_laws))
                          for k in range(1, LAST_GENERATION + 1)]
        k_best = int(np.argmin(integer_scores)) + 1
        z, p = test_target
        row = {
            "n": n,
            "target_type": target_type,
            "k_var": float(k_var),
            "t_var": float(t_var),
            "t_fit": t_best,
            "fit_train_ks": train_score,
            "k_best_integer": k_best,
            "fitted": evaluate(z, p, t_best, k_best, test_laws, test_laws),
            "variance_aligned": evaluate(z, p, t_var, int(round(k_var)), test_laws, test_laws),
        }
        rows.append(row)
        print(f"n={n:4d} {target_type:23s} t_fit={t_best:5.2f} "
              f"KS fitted={row['fitted']['ks_mixture']:.5f} "
              f"integer={row['fitted']['ks_integer']:.5f} "
              f"Q={row['fitted']['ks_q']:.5f} "
              f"skew={row['fitted']['target_shape'][0]:.3f}/"
              f"{row['fitted']['mixture_shape'][0]:.3f}", flush=True)
    payload = {
        "method": "Adjacent raw recursive-law mixture, then standardized to unit theoretical variance; train and test recursive populations are independent.",
        "population_per_generation": POPULATION,
        "finite_runs_per_split": FINITE_RUNS,
        "exact_max": EXACT_MAX,
        "q_generation": Q_GENERATION,
        "elapsed_seconds": round(perf_counter() - started, 2),
        "rows": rows,
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT} in {payload['elapsed_seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
