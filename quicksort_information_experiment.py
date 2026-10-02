"""Exact bits lost when tree depth profile is compressed to Quicksort cost.

Run: python quicksort_information_experiment.py --max-n 16
The insertion histories are counted exactly with integer multiplicities.
"""

import argparse
from collections import defaultdict
from math import factorial, log2


def entropy(weights):
    total = sum(weights)
    return -sum((w / total) * log2(w / total) for w in weights if w)


def cost(profile, n):
    return sum(depth * count for depth, count in enumerate(profile)) - 2 * n


def metrics(states, n):
    histories = factorial(n)
    by_cost = defaultdict(int)
    by_cost_and_next_depth = defaultdict(lambda: defaultdict(int))
    depth_given_profile = 0.0
    for profile, multiplicity in states.items():
        c = cost(profile, n)
        by_cost[c] += multiplicity
        depth_given_profile += (multiplicity / histories) * entropy(profile)
        for depth, slots in enumerate(profile):
            if slots:
                by_cost_and_next_depth[c][depth] += multiplicity * slots
    depth_given_cost = sum(
        (count / histories) * entropy(by_cost_and_next_depth[c].values())
        for c, count in by_cost.items()
    )
    profile_given_cost = entropy(states.values()) - entropy(by_cost.values())
    return {
        "n": n,
        "profiles": len(states),
        "cost_values": len(by_cost),
        "profile_bits_given_cost": profile_given_cost,
        "next_depth_bits_given_cost": depth_given_cost,
        "next_depth_bits_given_profile": depth_given_profile,
        "predictive_bits_lost": depth_given_cost - depth_given_profile,
    }


def evolve(states, n):
    following = defaultdict(int)
    for profile, multiplicity in states.items():
        for depth, slots in enumerate(profile):
            if slots == 0:
                continue
            next_profile = list(profile)
            next_profile[depth] -= 1
            if depth + 1 == len(next_profile):
                next_profile.append(0)
            next_profile[depth + 1] += 2
            while next_profile[-1] == 0:
                next_profile.pop()
            following[tuple(next_profile)] += multiplicity * slots
    assert sum(following.values()) == factorial(n + 1)
    return following


def main(max_n):
    states = {(1,): 1}
    print("n profiles costs H(profile|cost) H(next depth|cost) "
          "H(next depth|profile) predictive loss (bits)", flush=True)
    for n in range(max_n + 1):
        row = metrics(states, n)
        print(f"{n:2d} {row['profiles']:8d} {row['cost_values']:5d} "
              f"{row['profile_bits_given_cost']:15.8f} "
              f"{row['next_depth_bits_given_cost']:18.8f} "
              f"{row['next_depth_bits_given_profile']:21.8f} "
              f"{row['predictive_bits_lost']:18.8f}", flush=True)
        if n < max_n:
            states = evolve(states, n)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=16)
    main(parser.parse_args().max_n)
