# How predictive information in the Quicksort tree changes with size

Date: October 2, 2026.

## Question

Does the one-step information omitted by a comparison-count-only description persist or shrink as the binary search tree grows? The earlier exact calculation stopped at $n=16$, where the gap was 0.31545 bits. This report extends exact profile enumeration through $n=26$ and estimates the gap by simulation through $n=512$.

## Definition and exact computation

Let $K_n(d)$ count external slots at depth $d$, $C_n$ be the internal path length (the Quicksort comparison count), and $D_{n+1}$ be the depth of the next uniformly selected external slot. Then

\[
C_n=\sum_d dK_n(d)-2n,
\qquad
P(D_{n+1}=d\mid K_n)=\frac{K_n(d)}{n+1}.
\]

The predictive information discarded by retaining only $C_n$ is

\[
I_n:=I(D_{n+1};K_n\mid C_n)
=H(D_{n+1}\mid C_n)-H(D_{n+1}\mid K_n).
\]

For each depth profile $K$, the script counts the number $w_n(K)$ of insertion histories leading to it. Its exact integer recurrence is

\[
w_{n+1}(K-e_d+2e_{d+1})
\mathrel{+}=w_n(K)K(d).
\]

Here $e_d$ is the unit vector for depth $d$. The weights sum to $n!$ at every size. Probabilities are obtained by dividing by $n!$; logarithms are then evaluated numerically in base 2. Merging ordered trees by depth profile is valid here because that profile determines the distribution of the next insertion depth. Run `python quicksort_information_experiment.py --max-n 26` to reproduce the table.

## Results

| $n$ | Distinct profiles | $H(K_n\mid C_n)$, bits | $H(D_{n+1}\mid C_n)$, bits | $H(D_{n+1}\mid K_n)$, bits | $I_n$, bits |
|---:|---:|---:|---:|---:|---:|
| 6 | 9 | 0.222222 | 1.640600 | 1.519402 | 0.121198 |
| 7 | 16 | 0.272415 | 1.766221 | 1.637474 | 0.128747 |
| 8 | 28 | 0.613249 | 1.982253 | 1.731670 | 0.250584 |
| 9 | 50 | 0.827930 | 2.068808 | 1.815447 | 0.253362 |
| 10 | 89 | 1.140614 | 2.183977 | 1.890087 | 0.293890 |
| 11 | 159 | 1.428930 | 2.253310 | 1.955573 | 0.297737 |
| 12 | 285 | 1.754409 | 2.328558 | 2.013754 | 0.314804 |
| 13 | 510 | 2.060027 | 2.380021 | 2.066478 | 0.313544 |
| 14 | 914 | 2.378592 | 2.433414 | 2.114751 | **0.318662** |
| 15 | 1,639 | 2.687144 | 2.475009 | 2.159102 | 0.315907 |
| 16 | 2,938 | 2.998236 | 2.515405 | 2.199951 | 0.315453 |
| 17 | 5,269 | 3.302137 | 2.549483 | 2.237718 | 0.311765 |
| 18 | 9,451 | 3.604695 | 2.582143 | 2.272794 | 0.309349 |
| 19 | 16,952 | 3.900877 | 2.610768 | 2.305513 | 0.305254 |
| 20 | 30,410 | 4.193784 | 2.637979 | 2.336144 | 0.301835 |
| 21 | 54,555 | 4.481216 | 2.662586 | 2.364903 | 0.297683 |
| 22 | 97,871 | 4.764633 | 2.685837 | 2.391973 | 0.293864 |
| 23 | 175,586 | 5.043126 | 2.707287 | 2.417511 | 0.289775 |
| 24 | 315,016 | 5.317407 | 2.727554 | 2.441656 | 0.285899 |
| 25 | 565,168 | 5.587095 | 2.746481 | 2.464531 | 0.281950 |
| 26 | 1,013,976 | 5.852575 | 2.764400 | 2.486245 | 0.278155 |

The first nonzero gap occurs at $n=6$. Over the tested range it reaches its largest value, 0.318662 bits, at $n=14$, then decreases at every size through $n=26$. The decline from $n=14$ to $n=26$ is 0.040508 bits, or 12.7%. During the same interval, $H(K_n\mid C_n)$ rises from 2.378592 to 5.852575 bits.

## Larger trees: sampled estimates

Exact profile enumeration grows rapidly. To test the trend farther out, two independent runs each simulated 500,000 uniform-slot insertion histories at every listed size. For each sampled tree, its depth profile gives the exact next-depth probabilities and entropy. Averaging profiles with the same cost estimates $P(D_{n+1}\mid C_n)$. The pooled plug-in entropy tends to underestimate this conditional entropy. A separate held-out cross-entropy calculation, fitting on one half of the histories and evaluating on the other, tends to overestimate it. The pair is a bias diagnostic, **not a confidence interval**. Values below average the two independent runs.

| $n$ | Exact $I_n$ where available | Plug-in estimate, bits | Held-out estimate, bits |
|---:|---:|---:|---:|
| 16 | 0.31545 | 0.31540 | 0.31585 |
| 32 | — | 0.25661 | 0.25756 |
| 64 | — | 0.18583 | 0.18736 |
| 128 | — | 0.12756 | 0.12983 |
| 256 | — | 0.08515 | 0.08837 |
| 512 | — | 0.05619 | 0.06076 |

The $n=16$ estimate closely matches the exact value, and both estimators decline at every larger tested size. At $n=512$, the gap is roughly 0.06 bits per next insertion, about 80% below the exact peak near $n=14$. The two seeds agree closely; the growing distance between plug-in and held-out estimates shows why the larger-size values should be treated as approximate. Reproduce with `python quicksort_information_large_n_experiment.py --runs 500000` and the same command with `--seed 20261042`. Full outputs are in the adjacent JSON files.

## Interpretation and limit

The total cost conceals increasingly many bits about the external-depth profile through $n=26$, yet the portion useful for predicting the *next insertion depth* shrinks after $n=14$ in the exact range and continues to shrink in the sampled range. The evidence is consistent with $I_n\to0$. It does **not** prove that limit or a decay rate.

This metric is for the natural tree insertion coupling and next-step prediction. It is distinct from the Jensen–Shannon divergence between a continuous one-time interpolation and the finite-size count law. The exact profile state count exceeds one million by $n=26$, making direct enumeration increasingly costly.

The next theoretical question is the asymptotic order of $I_n$. A proof of decay or another long-run behavior would determine the appropriate paper claim about scalar state compression.
