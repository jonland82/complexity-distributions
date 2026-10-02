"""One-step finite-split corrections to the Quicksort limit approximation.

Run: python quicksort_split_experiment.py
The finite-size targets are exact PMFs (up to FFT roundoff). Model CDFs use
independent Monte Carlo samples, while reported skewness and kurtosis use
the one-step split identity and child-law moments to reduce simulation noise.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from make_figures import exact_mean_variance
from quicksort_fractional_experiment import exact_pmfs, population_trajectory, shape


ROOT = Path(__file__).resolve().parent
SIZES = (16, 32, 64, 128, 256)
DIRECT_TIMES = {16: 1.90, 32: 2.75, 64: 2.90, 128: 4.15, 256: 4.70}
DRAW_COUNT = 2_000_000
EXACT_CHILD_MAX = 8
SEED = 20261002


def moments_of_pool(pool):
    centered = pool - pool.mean()
    variance = np.mean(centered**2)
    return np.mean(centered**3) / variance**1.5, np.mean(centered**4) / variance**2


def model_moments(n, means, variances, child_m3, child_m4):
    """Third/fourth moments of the standardized one-step split law."""
    m = n - 1 - np.arange(n)
    i = np.arange(n)
    a = (n - 1 + means[i] + means[m] - means[n]) / np.sqrt(variances[n])
    b = np.sqrt(np.maximum(variances[i], 0) / variances[n])
    c = np.sqrt(np.maximum(variances[m], 0) / variances[n])
    v = b**2 + c**2
    third = np.mean(a**3 + 3 * a * v + child_m3[i] * b**3 + child_m3[m] * c**3)
    fourth = np.mean(
        a**4 + 6 * a**2 * v
        + 4 * a * (child_m3[i] * b**3 + child_m3[m] * c**3)
        + child_m4[i] * b**4 + child_m4[m] * c**4 + 6 * b**2 * c**2
    )
    return float(third), float(fourth - 3)


def sample_split(n, means, variances, q_pool, pmfs, rng, exact_children, draws):
    i = rng.integers(0, n, size=draws)
    j = n - 1 - i
    left = means[i] + np.sqrt(np.maximum(variances[i], 0)) * q_pool[rng.integers(q_pool.size, size=draws)]
    right = means[j] + np.sqrt(np.maximum(variances[j], 0)) * q_pool[rng.integers(q_pool.size, size=draws)]
    if exact_children:
        for m in range(2, EXACT_CHILD_MAX + 1):
            low, p = pmfs[m]
            left_mask = i == m
            right_mask = j == m
            left[left_mask] = low + np.searchsorted(np.cumsum(p), rng.random(left_mask.sum()))
            right[right_mask] = low + np.searchsorted(np.cumsum(p), rng.random(right_mask.sum()))
    return np.sort((n - 1 + left + right - means[n]) / np.sqrt(variances[n]))


def sample_direct(t, laws, rng, draws):
    k = int(np.floor(t))
    alpha = t - k
    which = rng.random(draws) < alpha
    indices = rng.integers(laws[k].size, size=draws)
    values = laws[k][indices].copy()
    values[which] = laws[k + 1][indices[which]]
    variance = (1 - alpha) * (1 - (2 / 3) ** k) + alpha * (1 - (2 / 3) ** (k + 1))
    return np.sort(values / np.sqrt(variance))


def compare(z, p, samples):
    target_right = np.cumsum(p)
    target_left = target_right - p
    model_left = np.searchsorted(samples, z, side="left") / samples.size
    model_right = np.searchsorted(samples, z, side="right") / samples.size
    return {
        "ks": float(max(np.max(np.abs(target_left - model_left)),
                        np.max(np.abs(target_right - model_right)))),
        "tail_2": float(1 - np.searchsorted(samples, 2, side="right") / samples.size),
        "tail_3": float(1 - np.searchsorted(samples, 3, side="right") / samples.size),
    }


def main(seed=SEED, draws=DRAW_COUNT):
    rng = np.random.default_rng(seed)
    pmfs = exact_pmfs(max(SIZES))
    laws = population_trajectory(seed + 1)
    q_pool = laws[35]
    q_m3, q_m4 = moments_of_pool(q_pool)
    means = np.array([exact_mean_variance(m)[0] for m in range(max(SIZES) + 1)])
    variances = np.array([exact_mean_variance(m)[1] for m in range(max(SIZES) + 1)])
    variances[:3] = 0  # The exact law has no randomness until n=3.
    q3 = np.full(max(SIZES) + 1, q_m3)
    q4 = np.full(max(SIZES) + 1, q_m4)
    hybrid3, hybrid4 = q3.copy(), q4.copy()
    for m in range(3, EXACT_CHILD_MAX + 1):
        low, p = pmfs[m]
        z = (np.arange(low, low + p.size) - means[m]) / np.sqrt(variances[m])
        hybrid3[m], excess, _, _ = shape(z, p)
        hybrid4[m] = excess + 3

    rows = []
    for n in SIZES:
        low, p = pmfs[n]
        z = (np.arange(low, low + p.size) - means[n]) / np.sqrt(variances[n])
        target_shape = shape(z, p)
        direct = sample_direct(DIRECT_TIMES[n], laws, rng, draws)
        pure = sample_split(n, means, variances, q_pool, pmfs, rng, False, draws)
        hybrid = sample_split(n, means, variances, q_pool, pmfs, rng, True, draws)
        row = {
            "n": n,
            "target": {"skew": target_shape[0], "excess_kurtosis": target_shape[1],
                       "tail_2": target_shape[2], "tail_3": target_shape[3]},
            "direct_recursive": compare(z, p, direct),
            "split_q_children": compare(z, p, pure),
            "split_exact_small_children": compare(z, p, hybrid),
        }
        for label, third, fourth in (("split_q_children", q3, q4),
                                     ("split_exact_small_children", hybrid3, hybrid4)):
            row[label]["skew"], row[label]["excess_kurtosis"] = model_moments(
                n, means, variances, third, fourth)
        rows.append(row)
        print(f"n={n:3d} KS direct/pure/hybrid="
              f"{row['direct_recursive']['ks']:.5f}/"
              f"{row['split_q_children']['ks']:.5f}/"
              f"{row['split_exact_small_children']['ks']:.5f} "
              f"skew target/pure/hybrid={target_shape[0]:.3f}/"
              f"{row['split_q_children']['skew']:.3f}/"
              f"{row['split_exact_small_children']['skew']:.3f}", flush=True)

    result = {
        "seed": seed, "draws_per_model_and_size": draws,
        "q_population": int(q_pool.size), "q_empirical_skew": q_m3,
        "q_empirical_excess_kurtosis": q_m4 - 3,
        "exact_child_max": EXACT_CHILD_MAX,
        "direct_times_from_prior_fitted_experiment": DIRECT_TIMES,
        "rows": rows,
    }
    suffix = "" if seed == SEED else f"_seed{seed}"
    out = ROOT / "artifacts" / f"quicksort_split_results_2026-10-02{suffix}.json"
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {out}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--draws", type=int, default=DRAW_COUNT)
    options = parser.parse_args()
    main(options.seed, options.draws)
