"""Count-binned Jensen-Shannon distance between exact and model laws.

Run: python quicksort_interpolation_information_experiment.py
The model's continuous standardized cost is converted back to a count and
rounded to the nearest integer before comparing distributions.
"""

import numpy as np

from make_figures import exact_mean_variance
from quicksort_fractional_experiment import exact_pmfs, population_trajectory
from quicksort_split_experiment import DIRECT_TIMES, sample_direct, sample_split


SIZES = (16, 32, 64)
DRAWS = 2_000_000
SEEDS = (20261002, 20261042)


def js_bits(target_low, target_p, simulated_counts):
    low = min(target_low, int(simulated_counts.min()))
    high = max(target_low + len(target_p) - 1, int(simulated_counts.max()))
    p = np.zeros(high - low + 1)
    p[target_low - low : target_low - low + len(target_p)] = target_p
    q = np.bincount(simulated_counts - low, minlength=len(p)).astype(float) / len(simulated_counts)
    mixture = (p + q) / 2
    p_mask, q_mask = p > 0, q > 0
    return float((np.sum(p[p_mask] * np.log2(p[p_mask] / mixture[p_mask]))
                  + np.sum(q[q_mask] * np.log2(q[q_mask] / mixture[q_mask]))) / 2)


def main():
    pmfs = exact_pmfs(max(SIZES))
    means = np.array([exact_mean_variance(n)[0] for n in range(max(SIZES) + 1)])
    variances = np.array([exact_mean_variance(n)[1] for n in range(max(SIZES) + 1)])
    variances[:3] = 0
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        laws = population_trajectory(seed + 1)
        print(f"seed={seed}  n  JS direct  JS split-Q  JS split-small-exact (bits)", flush=True)
        for n in SIZES:
            low, p = pmfs[n]
            models = (
                sample_direct(DIRECT_TIMES[n], laws, rng, DRAWS),
                sample_split(n, means, variances, laws[35], pmfs, rng, False, DRAWS),
                sample_split(n, means, variances, laws[35], pmfs, rng, True, DRAWS),
            )
            values = []
            for z in models:
                counts = np.rint(means[n] + np.sqrt(variances[n]) * z).astype(np.int64)
                values.append(js_bits(low, p, counts))
            print(f"{n:3d}  " + "  ".join(f"{v:.6f}" for v in values), flush=True)


if __name__ == "__main__":
    main()
