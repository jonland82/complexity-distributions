"""Reproduce the two figures in quicksort_complexity_distribution.tex.

Requires numpy, scipy, matplotlib, and numba. All random seeds are fixed.
The finite-size samples use the paper's uniform-pivot split recurrence.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from numba import njit
from scipy.stats import lognorm, norm


OUT = Path(__file__).resolve().parent
BLUE = "#21618c"
ORANGE = "#d35400"
GREEN = "#287c5d"


@njit
def quicksort_samples(n, repetitions, seed):
    """Independent comparison counts from uniform pivot choices at each node."""
    np.random.seed(seed)
    counts = np.empty(repetitions, dtype=np.int64)
    stack = np.empty(n, dtype=np.int64)
    for run in range(repetitions):
        top = 1
        stack[0] = n
        total = 0
        while top:
            top -= 1
            size = stack[top]
            if size < 2:
                continue
            total += size - 1
            left = np.random.randint(size)
            right = size - 1 - left
            if left > 1:
                stack[top] = left
                top += 1
            if right > 1:
                stack[top] = right
                top += 1
        counts[run] = total
    return counts


def exact_mean_variance(n):
    j = np.arange(1, n + 1, dtype=np.float64)
    harmonic = np.sum(1 / j)
    harmonic_2 = np.sum(1 / j**2)
    mean = 2 * (n + 1) * harmonic - 4 * n
    variance = 7 * n**2 + 13 * n - 2 * (n + 1) * harmonic - 4 * (n + 1) ** 2 * harmonic_2
    return mean, variance


def matched_lognormal(n):
    mean, variance = exact_mean_variance(n)
    log_variance = np.log1p(variance / mean**2)
    log_mean = np.log(mean) - log_variance / 2
    return lognorm(s=np.sqrt(log_variance), scale=np.exp(log_mean))


def quicksort_limit_samples(repetitions=400_000, iterations=25, seed=20260929):
    """Population iteration for Y = U Y1 + (1-U) Y2 + g(U)."""
    rng = np.random.default_rng(seed)
    pool = np.zeros(repetitions)
    for _ in range(iterations):
        u = rng.random(repetitions)
        a = pool[rng.integers(repetitions, size=repetitions)]
        b = pool[rng.integers(repetitions, size=repetitions)]
        toll = 1 + 2 * (u * np.log(u) + (1 - u) * np.log1p(-u))
        pool = u * a + (1 - u) * b + toll
        pool -= pool.mean()
    sigma = np.sqrt(7 - 2 * np.pi**2 / 3)
    return pool / sigma


def empirical_cdf(samples, grid):
    return np.searchsorted(np.sort(samples), grid, side="right") / samples.size


def plot_histograms(samples):
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.9))
    scores = {}
    for ax, (n, counts) in zip(axes, samples.items()):
        values = counts.astype(float)
        mean = values.mean()
        std = values.std()
        log_values = np.log(values)
        log_mean = log_values.mean()
        log_std = log_values.std()
        normal = norm(loc=mean, scale=std)
        lognormal = lognorm(s=log_std, scale=np.exp(log_mean))
        grid = np.linspace(values.min(), values.max(), 900)
        ax.hist(values, bins=35 if n == 64 else 48, density=True,
                color="#c7d7df", edgecolor="white", linewidth=0.4,
                label=f"{len(values):,} runs")
        ax.plot(grid, normal.pdf(grid), color=BLUE, lw=1.8, label="Normal")
        ax.plot(grid, lognormal.pdf(grid), color=ORANGE, lw=1.8, label="Lognormal")
        ax.set(xlabel="Comparisons", title=f"n = {n}")
        ax.grid(axis="y", color="#e1e5e8", linewidth=0.6)
        ax.set_axisbelow(True)
        scores[n] = (normal.logpdf(values).mean(), lognormal.logpdf(values).mean())
    axes[0].set_ylabel("Density")
    axes[1].legend(loc="upper right", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "quicksort_fits.pdf", bbox_inches="tight")
    fig.savefig(OUT / "quicksort_fits.svg", bbox_inches="tight")
    plt.close(fig)
    return scores


def plot_limit(counts_256, counts_4096, limit):
    grid = np.linspace(-2.7, 3.7, 700)
    panels = []
    for n, counts in ((256, counts_256), (4096, counts_4096)):
        mean, variance = exact_mean_variance(n)
        standardized = (counts - mean) / np.sqrt(variance)
        model = matched_lognormal(n)
        model_cdf = model.cdf(mean + grid * np.sqrt(variance))
        panels.append((f"n = {n:,}", empirical_cdf(standardized, grid), model_cdf))
    panels.append(("Limit", empirical_cdf(limit, grid), norm.cdf(grid)))

    fig, axes = plt.subplots(2, 3, figsize=(7.1, 3.8), sharex="col",
                             gridspec_kw={"height_ratios": [2, 1], "hspace": 0.09})
    for column, (title, quicksort_cdf, model_cdf) in enumerate(panels):
        ax = axes[0, column]
        ax.plot(grid, quicksort_cdf, color=GREEN, lw=1.8, label="Quicksort")
        ax.plot(grid, model_cdf, color=ORANGE, lw=1.8, ls="--", label="Lognormal")
        ax.set_title(title, fontsize=10)
        ax.set_xlim(-2.7, 3.7)
        ax.set_ylim(0, 1)
        ax.set_xticks([-2, 0, 2])
        ax.grid(color="#e1e5e8", linewidth=0.6)
        ax.set_axisbelow(True)
        difference = quicksort_cdf - model_cdf
        lower = axes[1, column]
        lower.plot(grid, difference, color=GREEN, lw=1.6)
        lower.axhline(0, color="#89969c", lw=0.8, ls=":")
        lower.set_ylim(-0.06, 0.08)
        lower.set_xlim(-2.7, 3.7)
        lower.set_xticks([-2, 0, 2])
        lower.set_xlabel("Standardized cost", fontsize=8)
        lower.grid(axis="y", color="#e1e5e8", linewidth=0.6)
        lower.set_axisbelow(True)
        if column > 0:
            ax.tick_params(labelleft=False)
            lower.tick_params(labelleft=False)
    axes[0, 0].set_ylabel("Cumulative probability", fontsize=8)
    axes[1, 0].set_ylabel("CDF gap", fontsize=8)
    axes[0, 2].legend(["Quicksort limit Q", "Normal limit"], loc="lower right",
                      frameon=False, fontsize=6.8)
    limit_gap = panels[-1][1] - panels[-1][2]
    maximum = np.argmax(np.abs(limit_gap))
    axes[1, 2].plot(grid[maximum], limit_gap[maximum], marker="o",
                    color=ORANGE, markersize=3.5)
    fig.subplots_adjust(left=0.095, right=0.985, top=0.92, bottom=0.13,
                        wspace=0.15, hspace=0.09)
    fig.savefig(OUT / "quicksort_limit.pdf", bbox_inches="tight")
    fig.savefig(OUT / "quicksort_limit.svg", bbox_inches="tight")
    plt.close(fig)
    return np.max(np.abs(limit_gap))


def main():
    plt.rcParams.update({"font.family": "serif", "font.size": 9,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    counts_64 = quicksort_samples(64, 30_000, 20260927)
    counts_256 = quicksort_samples(256, 30_000, 20260928)
    counts_4096 = quicksort_samples(4096, 30_000, 20260930)
    scores = plot_histograms({64: counts_64, 256: counts_256})
    limit = quicksort_limit_samples()
    estimated_gap = plot_limit(counts_256, counts_4096, limit)
    for n, (normal_score, lognormal_score) in scores.items():
        print(f"n={n} simulated mean log density: normal={normal_score:.6f}, "
              f"lognormal={lognormal_score:.6f}, difference={lognormal_score-normal_score:.6f}")
    print(f"Estimated limiting CDF gap (Monte Carlo): {estimated_gap:.4f}")


if __name__ == "__main__":
    main()
