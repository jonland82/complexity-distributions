"""Measure whether current cost predicts the next insertion variance.

Run: python quicksort_profile_experiment.py
"""

import json
from pathlib import Path

import numpy as np
from numba import njit

from make_figures import exact_mean_variance


CHECKPOINTS = (16, 64, 256)
RUNS = 100_000
SEED = 20261002
OUT = Path(__file__).resolve().parent / "artifacts" / "quicksort_profile_results_2026-10-02.json"


@njit
def simulate_profiles(runs, checkpoints, seed):
    np.random.seed(seed)
    costs = np.empty((len(checkpoints), runs), dtype=np.int64)
    innovations = np.empty((len(checkpoints), runs), dtype=np.float64)
    for run in range(runs):
        depths = np.empty(checkpoints[-1] + 1, dtype=np.int64)
        depths[0] = 0
        cost = 0
        checkpoint = 0
        for n in range(1, checkpoints[-1] + 1):
            slot = np.random.randint(n)
            depth = depths[slot]
            cost += depth
            depths[slot] = depth + 1
            depths[n] = depth + 1
            if n == checkpoints[checkpoint]:
                second = 0.0
                for j in range(n + 1):
                    second += depths[j] * depths[j]
                mean_depth = (cost + 2 * n) / (n + 1)
                innovations[checkpoint, run] = (
                    second / (n + 1) - mean_depth * mean_depth
                ) / (n + 2) ** 2
                costs[checkpoint, run] = cost
                checkpoint += 1
                if checkpoint == len(checkpoints):
                    break
    return costs, innovations


def main():
    costs, innovations = simulate_profiles(RUNS, np.array(CHECKPOINTS), SEED)
    rng = np.random.default_rng(SEED + 1)
    order = rng.permutation(RUNS)
    train, test = order[:RUNS // 2], order[RUNS // 2:]
    rows = []
    for j, n in enumerate(CHECKPOINTS):
        mean, variance = exact_mean_variance(n)
        z = (costs[j] - mean) / np.sqrt(variance)
        y = innovations[j]
        x = np.stack((np.ones(RUNS), z, z**2), axis=1)
        coefficients = np.linalg.lstsq(x[train], y[train], rcond=None)[0]
        fitted = x[test] @ coefficients
        constant = np.mean(y[train])
        baseline_mse = np.mean((y[test] - constant) ** 2)
        r2 = 1 - np.mean((y[test] - fitted) ** 2) / baseline_mse
        # Exact empirical conditional means by comparison count are a flexible
        # upper benchmark, evaluated only where the training count was seen.
        unique, inverse = np.unique(costs[j, train], return_inverse=True)
        sums = np.bincount(inverse, weights=y[train])
        number = np.bincount(inverse)
        locations = np.searchsorted(unique, costs[j, test])
        observed = locations < unique.size
        observed[observed] &= unique[locations[observed]] == costs[j, test][observed]
        predictions = np.full(test.size, constant)
        predictions[observed] = sums[locations[observed]] / number[locations[observed]]
        empirical_r2 = 1 - np.mean((y[test] - predictions) ** 2) / baseline_mse
        rows.append({
            "n": n,
            "mean_innovation_variance": float(np.mean(y)),
            "innovation_variance_cv": float(np.std(y) / np.mean(y)),
            "quadratic_cost_r2": float(r2),
            "empirical_cost_r2": float(empirical_r2),
            "same_cost_coverage": float(np.mean(observed)),
        })
        print(f"n={n:3d} CV(V)={rows[-1]['innovation_variance_cv']:.3f} "
              f"cost quadratic R2={r2:.3f} empirical cost R2={empirical_r2:.3f}", flush=True)
    OUT.write_text(json.dumps({"runs": RUNS, "seed": SEED, "rows": rows}, indent=2) + "\n")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
