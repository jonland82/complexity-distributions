# Information in Quicksort's discrete profile

Date: October 2, 2026.

Let $K_n(d)$ be the number of external tree slots at depth $d$, $C_n=\sum_d dK_n(d)-2n$, and $D_{n+1}$ the next insertion depth. Conditional on the profile,

\[
P(D_{n+1}=d\mid K_n)=K_n(d)/(n+1).
\]

The predictive information discarded by retaining only comparison count is

\[
I(D_{n+1};K_n\mid C_n)
=H(D_{n+1}\mid C_n)-H(D_{n+1}\mid K_n)
=\mathbb E\,D_{\rm KL,2}
\left(P(D_{n+1}\mid K_n)\Vert P(D_{n+1}\mid C_n)\right).
\]

Exact integer enumeration of all external-slot insertion histories, merged by depth profile, gives:

| $n$ | $H(K_n\mid C_n)$, bits | $H(D_{n+1}\mid C_n)$, bits | $H(D_{n+1}\mid K_n)$, bits | Predictive bits lost |
|---:|---:|---:|---:|---:|
| 6 | 0.222222 | 1.640600 | 1.519402 | 0.121198 |
| 8 | 0.613249 | 1.982253 | 1.731670 | 0.250584 |
| 16 | 2.998236 | 2.515405 | 2.199951 | 0.315453 |

This is a loss caused by compressing the discrete profile to cost. It is present even if the cost distribution is represented exactly. The fractional recursive interpolation specifies only one-time laws, so a predictive information loss for that interpolation is undefined until a joint process or coupling is specified.

For a separate measure of one-time distribution mismatch, map each continuous model draw back to a comparison count and round to the nearest integer. With an equal prior over “exact Quicksort” and “model,” the Jensen–Shannon divergence in bits is the information one rounded count gives about which source generated it. Two independent two-million-draw runs gave:

| $n$ | Fractional interpolation, bits | Finite-size split correction, bits |
|---:|---:|---:|
| 16 | 0.00509 | 0.00186 (exact child laws through 8) |
| 32 | 0.00235 | 0.00135 (terminal-law children) |
| 64 | 0.00085 | 0.00070 (terminal-law children) |

These small divergence values describe marginal fit after the stated rounding convention. They are not the predictive bits above and cannot be added to or subtracted from them. The $n=64$ improvement is small relative to Monte Carlo resolution. Without rounding, the exact finite-size law is discrete and the interpolation is continuous, so ordinary KL divergence is infinite and does not provide a useful finite bit count.

Reproduce with `python quicksort_information_experiment.py --max-n 16` and `python quicksort_interpolation_information_experiment.py`.
