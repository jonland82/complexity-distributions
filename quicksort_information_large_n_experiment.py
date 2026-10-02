"""Estimate predictive profile information for larger Quicksort trees.

The exact depth profile is simulated under uniform external-slot insertion.
Two independent halves estimate the conditional depth law given cost. A pooled
plug-in entropy and held-out cross entropy expose finite-sample bias.

Run: python quicksort_information_large_n_experiment.py --runs 500000
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numba import njit


SIZES = (16, 32, 64, 128, 256, 512)
MAX_DEPTH = 64
SEED = 20261002
OUT = Path(__file__).resolve().parent / "artifacts" / "quicksort_information_large_n_2026-10-02.json"


@njit
def collect(n, runs_per_half, seed):
    max_cost = n * (n - 1) // 2
    joint = np.zeros((2, max_cost + 1, MAX_DEPTH), dtype=np.float64)
    counts = np.zeros((2, max_cost + 1), dtype=np.int64)
    profile_entropy = np.zeros(2)
    deepest = 0
    for half in range(2):
        np.random.seed(seed + half)
        for _ in range(runs_per_half):
            depths = np.empty(n + 1, dtype=np.int64)
            depths[0] = 0
            cost = 0
            for size in range(1, n + 1):
                slot = np.random.randint(size)
                depth = depths[slot]
                cost += depth
                depths[slot] = depth + 1
                depths[size] = depth + 1
            profile = np.zeros(MAX_DEPTH, dtype=np.int64)
            for slot in range(n + 1):
                depth = depths[slot]
                if depth >= MAX_DEPTH:
                    raise ValueError("Depth exceeded MAX_DEPTH")
                profile[depth] += 1
                if depth > deepest:
                    deepest = depth
            counts[half, cost] += 1
            for depth in range(deepest + 1):
                slots = profile[depth]
                if slots:
                    probability = slots / (n + 1)
                    joint[half, cost, depth] += probability
                    profile_entropy[half] -= probability * np.log2(probability)
    return joint, counts, profile_entropy, deepest


def entropy_given_cost(joint, counts):
    totals = counts.sum(axis=0)
    mask = totals > 0
    conditional = joint[:, mask, :].sum(axis=0) / totals[mask, None]
    terms = np.zeros_like(conditional)
    np.log2(conditional, out=terms, where=conditional > 0)
    entropy = -np.sum(conditional * terms, axis=1)
    return float(np.dot(totals[mask], entropy) / totals.sum())


def heldout_cross_entropy(joint, counts, alpha=10.0):
    values = []
    for train in (0, 1):
        test = 1 - train
        global_depth = joint[train].sum(axis=0)
        prior = (global_depth + 1e-12) / (global_depth.sum() + 1e-12 * MAX_DEPTH)
        q = (joint[train] + alpha * prior[None, :]) / (counts[train, :, None] + alpha)
        values.append(float(-np.sum(joint[test] * np.log2(q)) / counts[test].sum()))
    return float(np.mean(values)), values


def main(runs, seed=SEED):
    rows = []
    for n in SIZES:
        joint, counts, profile_entropy, deepest = collect(n, runs // 2, seed + n * 10)
        prof_h = float(profile_entropy.sum() / runs)
        pooled_h = entropy_given_cost(joint, counts)
        cross_h, directions = heldout_cross_entropy(joint, counts)
        row = {
            "n": n,
            "runs": runs,
            "deepest_external_slot_seen": deepest,
            "profile_entropy": prof_h,
            "conditional_depth_entropy_plugin": pooled_h,
            "conditional_depth_entropy_crossfit": cross_h,
            "predictive_bits_plugin": pooled_h - prof_h,
            "predictive_bits_crossfit": cross_h - prof_h,
            "crossfit_directions": directions,
        }
        rows.append(row)
        print(f"n={n:3d} predictive bits plug-in/cross-fit="
              f"{row['predictive_bits_plugin']:.5f}/"
              f"{row['predictive_bits_crossfit']:.5f} "
              f"deepest={deepest}", flush=True)
    out = OUT if seed == SEED else OUT.with_name(f"{OUT.stem}_seed{seed}.json")
    out.write_text(json.dumps({"seed": seed, "rows": rows}, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {out}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=int, default=500_000)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    if args.runs < 2 or args.runs % 2:
        parser.error("--runs must be an even integer at least 2")
    main(args.runs, args.seed)
