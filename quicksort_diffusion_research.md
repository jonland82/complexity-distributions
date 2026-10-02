# Quicksort size evolution and diffusion: research notes

Date: October 1, 2026.

## Findings

The natural coupling across input sizes gives a convergent martingale. Its limit is the non-normal Quicksort law, rather than an equilibrium maintained by persistent noise. The exact next-step noise depends on the external depth profile of a growing binary search tree. Cost alone is not Markov under this coupling; an exact counterexample appears below.

A second, complementary evolution is already implicit in Quicksort's limiting fixed-point equation. Iterating its distributional operator from a point mass produces a nonparametric family that converges to the Quicksort limit. This is a more natural candidate for studying shape evolution than searching for another fixed named family such as the normal, lognormal, or beta. The iteration is indexed by recursive generation, not directly by input size, so its relationship to the exact finite-size laws must be measured rather than assumed.

A diffusion approximation remains plausible for large-size corrections. Published results already establish a Gaussian limit for the shrinking residual between the normalized finite-size cost and its non-Gaussian terminal limit. This distinguishes Gaussian late noise from a Gaussian limiting cost distribution.

It is also possible to construct diffusions having a point mass at zero as their initial law and $Q$ as their terminal law. Such endpoint constructions are not unique and, by themselves, do not identify Quicksort's finite-size evolution. No scalar diffusion reproducing the finite-size Quicksort laws has been derived here. The conditional variance below identifies a candidate coefficient, and the profile process gives an exact alternative with additional state.

## Context and scope

The discussion concerns the comparison count of uniform-pivot Quicksort on distinct keys, as specified in [the existing paper](quicksort_complexity_distribution.tex). The evolution variable is input size, not elapsed execution time within one sort.

The paper states two results:

- For every integer $n\ge4$, the optimal ordinary lognormal has higher expected log-density on the finite Quicksort law than the optimal normal. This is a relative score ranking, not a theorem about monotonic approximation error as size changes.
- After centering and scaling, Quicksort approaches a non-normal law $Q$, whereas the moment-matched lognormal approaches a standard normal. The paper also gives the normal standardized limit for the log-score-optimal lognormal, standardized using its own moments. A finite-size advantage can therefore coexist with a positive limiting distributional gap.

An earlier exploratory beta comparison used theoretical count bounds padded by half a comparison. Independent training/evaluation samples favored beta over lognormal over normal at sizes 64, 256, and 4096. Beta-minus-lognormal score differences were approximately 0.0351, 0.0259, and 0.0190, with paired evaluation standard errors 0.00082, 0.00070, and 0.00047. Full-PMF floating-point calculations also favored beta at sizes 4 through 32. These are discussion results, not an all-size theorem; their former scripts and outputs were removed at the user's request. The moment-matched beta with those bounds also loses its skewness in the limit. More fundamentally, normal, lognormal, and beta are static approximation families, not starting laws selected by Quicksort's dynamics. A family can fit particular finite sizes reasonably well without containing the limiting law $Q$.

The user then proposed studying the evolution itself: find a distributional family that develops toward $Q$, perhaps with a diffusion-like continuous interpolation. A Langevin process targeting the density of $Q$ is possible, but specifying the correct equilibrium does not derive Quicksort's size evolution. The appropriate order is to identify the exact or recursively generated marginal family first and only then ask whether a diffusion accurately interpolates it. The theory below therefore separates the exact size coupling, a fixed-point distributional evolution, and endpoint-only diffusion constructions.

This task authorizes a new Markdown research note only. The paper, README, code, figures, and existing artifacts are not updated. Exact verification code was run through standard input and is included below for reproducibility.

## 1. Notation and a natural coupling

Let

$$
\mu_n=2(n+1)H_n-4n,
\qquad
v_n=7n^2+13n-2(n+1)H_n-4(n+1)^2H_n^{(2)}.
$$

Here $H_0=H_0^{(2)}=0$, and the paper's comparison counts have $C_0=C_1=0$.

Construct an ordered binary search tree by starting with one empty root slot. At size $n$, choose uniformly among its $n+1$ external slots and insert a node there. A slot at depth $d$ becomes an internal node with two external children at depth $d+1$. Root depth is zero. The internal path length has the same marginal law as Quicksort's comparison count. This construction is described in [Grübel and Kabluchko, Section 5.5.1](https://arxiv.org/pdf/1410.0469).

Write $\mathcal F_n$ for the information in the tree insertion history through size $n$. This is an abstract tree/rank coupling. It does not condition on numerical key spacings; if actual numerical key values were also revealed, insertion probabilities into their gaps would not generally be uniform.

Let $D_{n+1}$ be the depth selected for the next insertion. Under this coupling,

$$
C_{n+1}=C_n+D_{n+1}.
$$

This couples the marginal laws from the paper. It does not prescribe shared randomness for independently executed randomized-pivot sorts.

Define the external depth profile and its moments:

$$
K_n(d)=\#\{\text{external slots at depth }d\},
\qquad S_{r,n}=\sum_d d^r K_n(d).
$$

Counting slots and summing depths gives

$$
S_{0,n}=n+1,
\qquad S_{1,n}=C_n+2n,
\qquad
\Pr(D_{n+1}=d\mid \mathcal F_n)=\frac{K_n(d)}{n+1}.
$$

Consequently, by direct calculation,

$$
\mathbb E[D_{n+1}\mid \mathcal F_n]=\frac{C_n+2n}{n+1},
$$

$$
\operatorname{Var}(D_{n+1}\mid \mathcal F_n)
=\frac{S_{2,n}}{n+1}-\left(\frac{C_n+2n}{n+1}\right)^2.
$$

The conditional mean closes in cost and size. The conditional variance needs additional profile information.

## 2. Exact normalized increments

Use the martingale normalization

$$
M_n=\frac{C_n-\mu_n}{n+1}.
$$

It differs slightly from the paper's unit-variance normalization. Write

$$
a_n=\operatorname{Var}(M_n)=\frac{v_n}{(n+1)^2},
\qquad \widehat C_n=\frac{M_n}{\sqrt{a_n}}\quad(n\ge3).
$$

The exact mean identity

$$
\mu_{n+1}=\frac{n+2}{n+1}\mu_n+\frac{2n}{n+1}
$$

and the cost increment give

$$
\boxed{
M_{n+1}-M_n
=\frac{D_{n+1}-(C_n+2n)/(n+1)}{n+2}.
}
$$

Thus

$$
\boxed{\mathbb E[M_{n+1}-M_n\mid \mathcal F_n]=0,}
$$

$$
\boxed{
\mathbb E[(M_{n+1}-M_n)^2\mid \mathcal F_n]
=\frac{S_{2,n}/(n+1)-[(C_n+2n)/(n+1)]^2}{(n+2)^2}.
}
$$

These formulas derive the martingale property directly. The classical martingale limit was established by [Régnier](https://www.numdam.org/item/ITA_1989__23_3_335_0/); the normalized coupling and convergence are also reviewed in [Fill and Janson, Section 7](https://arxiv.org/pdf/math/0105248).

Let

$$
\sigma^2=7-\frac{2\pi^2}{3}=0.420263732607\ldots.
$$

The martingale is bounded in second moment, so it converges almost surely and in mean square to $Y$, with variance $\sigma^2$. The paper's standardized limit is $Q=Y/\sigma$.

The exact starting point is $M_0=0$. The normalized process remains zero through size 2; genuine randomness begins at size 3. Nothing has to be inferred backward to obtain this starting condition once the coupling is specified.

### Standardization introduces a drift, but not an equilibrium mechanism

For $n\ge3$,

$$
\widehat C_{n+1}
=\sqrt{\frac{a_n}{a_{n+1}}}\,\widehat C_n
+\frac{M_{n+1}-M_n}{\sqrt{a_{n+1}}}.
$$

The conditional mean increment is therefore

$$
\mathbb E[\widehat C_{n+1}-\widehat C_n\mid \mathcal F_n]
=\left(\sqrt{\frac{a_n}{a_{n+1}}}-1\right)\widehat C_n.
$$

This apparent restoring term comes from changing the normalization. Both its strength and the innovation variance diminish at large sizes. It is not evidence that the natural process is a constant-noise Langevin process equilibrating to $Q$.

## 3. An exact scalar Markov counterexample

At size 6, consider the two positive-probability cost histories

$$
H_A:(C_1,\ldots,C_6)=(0,1,3,6,10,11),
$$

$$
H_B:(C_1,\ldots,C_6)=(0,1,3,6,9,11).
$$

Both have the same current state $(n,C_n)=(6,11)$. Exact enumeration of all 720 equally likely insertion histories gives probabilities $1/45$ and $1/90$ for these two cost histories. Their conditional next insertion depth laws differ:

| Next depth $d$ | Given $H_A$ | Given $H_B$ |
|---:|---:|---:|
| 1 | 0 | 1/7 |
| 2 | 3/7 | 0 |
| 3 | 1/7 | 2/7 |
| 4 | 1/7 | 4/7 |
| 5 | 2/7 | 0 |

For example, $\Pr(C_7=12\mid H_A)=0$, whereas $\Pr(C_7=12\mid H_B)=1/7$. Hence the scalar cost process is not Markov under this coupling, even when size is included as deterministic time. The same conclusion applies to its invertible size-dependent normalizations.

The corresponding external depths are, up to ordering,

$$
A:(2,2,2,3,4,5,5),\qquad B:(1,3,3,4,4,4,4).
$$

Both sum to 23, but their sums of squares are 87 and 83. Their insertion-depth variances are $80/49$ and $52/49$; their conditional normalized increment variances are respectively $5/196$ and $13/784$.

This rules out an exact scalar Markov description of these coupled paths. It does not rule out a diffusion approximation, a different coupling, or a scalar model that reproduces only the one-time marginal distributions.

## 4. The total remaining noise is finite and shrinks

Martingale increments are orthogonal in second moment. Algebra using the paper's exact variance gives

$$
a_n=7-4H_n^{(2)}-\frac{2H_n+1}{n+1}-\frac{6}{(n+1)^2}.
$$

The unconditional innovation variance is exactly

$$
\boxed{
\eta_n:=\mathbb E[(M_{n+1}-M_n)^2]
=a_{n+1}-a_n
=\frac{2H_n-1}{(n+1)(n+2)}
+\frac{2}{(n+1)^2}-\frac{6}{(n+2)^2}.
}
$$

In particular,

$$
\eta_n\sim\frac{2\log n}{n^2}.
$$

The remaining mean-square error is

$$
\boxed{
\mathbb E[(Y-M_n)^2]
=\sigma^2-a_n
=\frac{2H_n+1}{n+1}+\frac{6}{(n+1)^2}
-4\left(\frac{\pi^2}{6}-H_n^{(2)}\right)
\sim\frac{2\log n}{n}.
}
$$

This agrees with the exact coupled error formula in [Bindjeme and Fill, Theorem 1.4](https://arxiv.org/pdf/1201.6445). It is a coupling-specific mean-square error, not the optimal Wasserstein distance between the two marginal laws.

| $n$ | $\operatorname{Var}(M_n)$ | Remaining variance $\mathbb E[(Y-M_n)^2]$ |
|---:|---:|---:|
| 4 | 0.032222222 | 0.388041510 |
| 64 | 0.319507387 | 0.100756346 |
| 256 | 0.384216090 | 0.036047642 |
| 4096 | 0.416653485 | 0.003610248 |
| Limit | 0.420263733 | 0 |

So $Q$ is the standardized terminal random value of this process. Individual normalized paths converge; they do not keep wandering with nonzero constant noise at the limit.

## 5. Gaussian residuals are established theory

[Neininger, Theorem 1.1](https://arxiv.org/pdf/1207.4556) proves

$$
\sqrt{\frac{n}{2\log n}}(M_n-Y)
\Rightarrow\mathcal N(0,1).
$$

Thus the small remaining error has an asymptotically Gaussian shape, even though $Y$ itself is non-Gaussian. This supports investigating diffusion descriptions of late corrections.

[Grübel and Kabluchko, Section 5.5.1](https://arxiv.org/pdf/1410.0469) give a stronger conditional version in their external-path-length normalization, and [Sulzbach](https://arxiv.org/pdf/1412.3508) studies martingale tails and higher moments. These strengthen the case for Gaussian late residuals; they do not identify the finite-size residual as exactly independent Gaussian noise.

### Independence matters when working backward

Under this coupling,

$$
M_n=\mathbb E[Y\mid \mathcal F_n],
\qquad \mathbb E[Y-M_n\mid \mathcal F_n]=0.
$$

This is an exact interpretation of working backward from a terminal variable: finite-size values are conditional expectations as tree information is revealed. The terminal marginal law alone does not specify the filtration or its relation to $Y$.

The residual is uncorrelated with $M_n$, but not with $Y$:

$$
\operatorname{Cov}(M_n,Y-M_n)=0,
\qquad \operatorname{Cov}(Y,Y-M_n)=\sigma^2-a_n.
$$

Replacing $M_n=Y-(Y-M_n)$ by $Y$ minus an independent Gaussian with the same residual variance would incorrectly give variance $\sigma^2+(\sigma^2-a_n)$, instead of $a_n$. A residual CLT alone is not a finite-size independent-noise reconstruction formula.

## 6. What kind of evolution model remains viable?

### The exact size-indexed family

There is already an exact family that converges to $Q$; it does not belong to a standard named parametric class. Define

$$
P_n=\operatorname{Law}\left(\frac{M_n}{\sigma}\right)
=\operatorname{Law}\left(\frac{C_n-\mu_n}{(n+1)\sigma}\right).
$$

Then

$$
P_0=P_1=P_2=\delta_0,
\qquad P_n\Rightarrow\operatorname{Law}(Q),
$$

and

$$
\operatorname{Var}(P_n)=\frac{a_n}{\sigma^2}\uparrow1.
$$

This normalization is preferable for evolution questions because it retains the deterministic starting point and lets variance accumulate toward its terminal value. By contrast, the paper's unit-variance normalization $\widehat C_n=M_n/\sqrt{a_n}$ is undefined while $a_n=0$ and forces every later marginal to have variance one, hiding this part of the evolution.

The distributions $P_n$ are generated exactly by the finite-size Quicksort recurrence

$$
C_n\overset d=n-1+C_{I_n}+C'_{n-1-I_n}.
$$

They are discrete at every finite $n$. Thus an ordinary nondegenerate Brownian diffusion cannot reproduce all of them exactly at deterministic positive times. The exact size evolution is the tree/profile jump process; a scalar diffusion can only be a smoothed or asymptotic interpolation.

### A recursive distribution family converging to $Q$

Quicksort's limiting fixed-point equation supplies a second useful evolution. Put

$$
g(u)=1+2u\log u+2(1-u)\log(1-u),
$$

and, for a centered probability law $\nu$, define

$$
\mathcal T(\nu)
=\operatorname{Law}\left(
U X_1+(1-U)X_2+\frac{g(U)}{\sigma}
\right),
$$

where $U$ is uniform on $(0,1)$ and $X_1,X_2$ are independent with law $\nu$, independently of $U$. The standardized Quicksort law is a fixed point:

$$
\mathcal T(\operatorname{Law}(Q))=\operatorname{Law}(Q).
$$

Starting from a point mass gives the recursively generated family

$$
\nu_0=\delta_0,
\qquad \nu_{k+1}=\mathcal T(\nu_k).
$$

On centered laws with finite second moment, the transform contracts squared Wasserstein distance by at most the factor

$$
\mathbb E[U^2+(1-U)^2]=\frac23.
$$

Consequently $\nu_k$ converges to the unique centered finite-variance fixed point $\operatorname{Law}(Q)$. Since $\operatorname{Var}(g(U)/\sigma)=1/3$, its variance evolves explicitly as

$$
v_{k+1}=\frac23v_k+\frac13,
\qquad v_0=0,
\qquad
\boxed{v_k=1-\left(\frac23\right)^k.}
$$

This family develops the skewness and tail structure of $Q$ through the same split-and-recombine operation that defines the Quicksort limit. It is therefore a better shape-evolution candidate than a transition among normal, lognormal, and beta families.

The distinction between the indices matters. The exact family $P_n$ is indexed by input size; $\nu_k$ is indexed by recursive population iteration. There is no established identity $k=k(n)$. A variance-matched or empirically optimized mapping can be tested, but it must not be presented as exact without evidence.

The existing `quicksort_limit_samples` function in `make_figures.py` already simulates $\mathcal T^k(\delta_0)$ for 25 iterations, but retains only the last population. Saving each intermediate population would expose this evolution directly.

### A diffusion with endpoint $Q$ exists, but is not unique

An endpoint-only construction confirms that diffusion is possible in principle. Let $F_Q$ be the CDF of $Q$, let $\Phi$ be the standard normal CDF, and define the quantile map

$$
r(z)=F_Q^{-1}(\Phi(z)).
$$

For Brownian motion $B_t$, set

$$
X_t=\mathbb E[r(B_1)\mid B_t],
\qquad 0\le t\le1.
$$

Then $X_0=\mathbb E Q=0$ and $X_1=r(B_1)\sim Q$. Writing

$$
u(t,x)=\mathbb E\left[r\left(x+\sqrt{1-t}\,Z\right)\right],
\qquad Z\sim\mathcal N(0,1),
$$

gives $X_t=u(t,B_t)$. The backward heat equation and Itô's formula yield a driftless, time-inhomogeneous diffusion. When $u(t,\cdot)$ is inverted, it can be written

$$
dX_t=s(t,X_t)\,dW_t,
\qquad
s(t,x)=\partial_xu\left(t,u^{-1}(t,x)\right).
$$

This [Bass-martingale construction](https://www.numdam.org/item/SPS_1983__17__221_0/) gives a concrete continuous family

$$
\delta_0\longrightarrow\operatorname{Law}(Q).
$$

It is selected by a Brownian quantile construction, not by Quicksort's finite-size recurrence. Infinitely many other processes share the same endpoint laws, so agreement at $t=0$ and $t=1$ is insufficient. Its intermediate marginals should be treated as a candidate and compared with $P_n$, not identified with them in advance.

### Exact profile process and continuous-time embedding

The external profile is a sufficient Markov state. At each insertion,

$$
K_{n+1}=K_n-e_D+2e_{D+1},
\qquad \Pr(D=d\mid K_n)=K_n(d)/(n+1),
$$

where $e_d$ is the profile with one slot at depth $d$ and zero elsewhere.

Give each external slot an independent rate-one splitting clock. The resulting continuous-time branching process has generator

$$
\mathcal L f(K)=\sum_d K(d)
\bigl[f(K-e_d+2e_{d+1})-f(K)\bigr].
$$

Its embedded insertion chain is exactly the profile process above. This is an exact jump process, not a Brownian diffusion. The branching/Yule embedding is described in [Sulzbach, Section 2.3](https://arxiv.org/pdf/1412.3508).

Its total event rate at size $n$ is $n+1$. Consequently log-size is a natural approximate continuous clock; the embedding's actual clock is random, so it must not simply be identified with deterministic $\log n$.

Tracking only cost and $S_2$ does not automatically close all fluctuation equations: an insertion changes $S_2$ by $D^2+4D+2$, whose variance depends on third and fourth profile moments.

### A scalar marginal approximation

Averaging over hidden profiles at a given current cost gives an exact discrete-time transition kernel for one-time marginals:

$$
\Pr_{\mathrm{projected}}(C_{n+1}=c+d\mid C_n=c)
=\mathbb E\left[\frac{K_n(d)}{n+1}\,\middle|\,C_n=c\right].
$$

A Markov chain using these time-dependent kernels has the correct marginal laws by total probability, but does not reproduce the natural process's joint history laws.

Define the corresponding candidate local variance for normalized cost:

$$
A(n,m)=\mathbb E\left[
\frac{S_{2,n}/(n+1)-[(C_n+2n)/(n+1)]^2}{(n+2)^2}
\,\middle|\,M_n=m\right].
$$

A formal diffusion approximation would be

$$
dZ_n=\sqrt{A(n,Z_n)}\,dW_n,
\qquad
\partial_n p(m,n)=\frac12\partial_m^2[A(n,m)p(m,n)].
$$

Here $n$ is interpolated continuously. This is a candidate, not a derived exact PDE: neglecting higher conditional jump moments and temporal discretization must be justified. The coefficient also requires the hidden-profile conditional distribution; it is not yet an explicit closed function of $n$ and $m$.

With deterministic log-size $t=\log n$, the candidate noise variance per unit $t$ is $nA(n,m)$. Its mean scales like $2\log n/n=2te^{-t}$, so the late effective noise dies away. This asymptotic statement concerns the mean coefficient, not a proof of a uniform pointwise approximation in $m$.

### Why a deterministic Gaussian variance fit is insufficient

Replacing the conditional variance by its unconditional value $\eta_n$ and using independent Gaussian innovations can match second moments while missing the shape.

For example, start with the exact $M_4$ and replace all later innovations by independent centered Gaussians with total variance $\sigma^2-a_4$. The terminal variance is correct, but its skewness is only about $0.007069$, versus $0.854882$ for the actual $Q$. The calculation uses

$$
a_4=29/900,\qquad \mathbb E[M_4^3]=13/6750,
$$

and the fact that adding an independent centered Gaussian preserves the third central moment. This is a test of that particular surrogate, not an impossibility result for state-dependent diffusions.

The exact third-moment increment makes the missing dependence visible:

$$
\mathbb E[M_{n+1}^3-M_n^3\mid \mathcal F_n]
=3M_n V_n+T_n,
$$

where $V_n=\mathbb E[(M_{n+1}-M_n)^2\mid\mathcal F_n]$ and $T_n=\mathbb E[(M_{n+1}-M_n)^3\mid\mathcal F_n]$. Correlation between current cost and future variance, as well as finite jump asymmetry, matters for the evolving shape.

### The earlier equilibrium construction remains a different model

If $q$ is the density of $Q$, the artificial Langevin construction

$$
dX_t=D\,\partial_x\log q(X_t)\,dt+\sqrt{2D}\,dW_t
$$

has $q$ as a stationary density. Its density equation is

$$
\partial_t p=D\,\partial_x\left[p\,\partial_x\log(p/q)\right].
$$

These equations specify an equilibrium target, not Quicksort's natural size coupling. They also require knowledge or estimation of $q$; writing them does not produce elementary closed coefficients. Starting at equilibrium cannot identify a unique finite-size history.

This stationary construction should not be confused with either $\nu_{k+1}=\mathcal T(\nu_k)$ or the Bass martingale. The fixed-point iteration approaches $Q$ through Quicksort's recursive smoothing transform. The Bass martingale reaches $Q$ at a finite terminal time from $\delta_0$. The stationary Langevin process, if initialized with $Q$, remains distributed as $Q$ and reveals no earlier shape.

## 7. Status and the next focused experiment

Established here by derivation and exact finite verification:

- The conditional drift and variance of normalized increments.
- The scalar Markov counterexample at size 6.
- Exact unconditional innovation and remaining-variance formulas, consistent with the cited theory.
- Failure of the independent-Gaussian-future surrogate started at size 4 to preserve limiting skewness.
- Identification of the exact evolving marginal family $P_n$, beginning at $\delta_0$ and converging to $Q$.
- A Quicksort-specific recursive family $\nu_k=\mathcal T^k(\delta_0)$ converging to $Q$, with variance $1-(2/3)^k$.
- An endpoint-only Bass diffusion from $\delta_0$ to $Q$, distinguished from the finite-size Quicksort process.

Established in the cited literature:

- Martingale convergence to the Quicksort limit.
- The Gaussian limit for the rescaled martingale residual.
- Stronger conditional and higher-moment versions of residual limit results.

Still unresolved here:

- An explicit accurate scalar variance function $A(n,m)$.
- Error bounds for a scalar Fokker–Planck approximation.
- How far toward small sizes such an approximation can be trusted.
- Whether the recursive-generation marginals $\nu_k$ approximate the exact size marginals $P_n$ under a useful mapping $k=k(n)$.
- Whether the intermediate Bass marginals approximate $P_n$ after a variance or time alignment.
- A finite-size reconstruction from $Q$ alone without selecting additional dynamics; the Bass construction is one selection, not an identification theorem.

The next experiment should first expose the distributional evolution already present in the repository. Modify `quicksort_limit_samples` to retain the population after every fixed-point iteration. At each iteration record variance, skewness, kurtosis, and distance to a high-accuracy terminal approximation of $Q$. Separately simulate the natural coupled insertion/profile process and record $P_n$, $M_n$, its increment, and $S_2$ at logarithmically spaced size checkpoints. Compare $\nu_k$ with $P_n$ using both variance matching and direct held-out CDF optimization; report the resulting $k(n)$ rather than assuming one.

Only after that comparison should a continuous diffusion be fitted. Compare three baselines: the exact full-profile jump process, the recursive family $\nu_k$, and a scalar diffusion using a fitted version of $A(n,m)$. An endpoint-only Bass interpolation can be included as a fourth baseline. Evaluate one-time CDFs, skewness, kurtosis, tail probabilities, and joint increments. For late-size comparisons, check the normalized residual against a Gaussian, accounting for the finite endpoint used to approximate $Y$.

No such evolution comparison, large coupled simulation, or diffusion fit was run for this note. The existing figure code runs only the terminal fixed-point population, and the checks below are finite exact calculations rather than the proposed experiment.

## Appendix: reproducible exact checks

Run this block with Python. It writes no files and requires only the standard library. It enumerates insertion choices through size 6, merges equivalent states with multiplicities, and uses exact rational arithmetic.

```python
from collections import Counter
from fractions import Fraction as F
from math import factorial

states = {((0,), ()): 1}  # (ordered external depths, cost history) -> count
for n in range(1, 7):
    following = Counter()
    for (depths, history), count in states.items():
        cost = history[-1] if history else 0
        assert sum(depths) == cost + 2 * (n - 1)
        for i, depth in enumerate(depths):
            updated = depths[:i] + (depth + 1, depth + 1) + depths[i + 1:]
            following[updated, history + (cost + depth,)] += count
    states = following
    assert sum(states.values()) == factorial(n)

    harmonic = sum((F(1, k) for k in range(1, n + 1)), F(0))
    harmonic2 = sum((F(1, k * k) for k in range(1, n + 1)), F(0))
    mean = 2 * (n + 1) * harmonic - 4 * n
    variance = (7 * n * n + 13 * n - 2 * (n + 1) * harmonic
                - 4 * (n + 1)**2 * harmonic2)
    observed_variance = F(0)
    observed_innovation = F(0)
    for (depths, history), count in states.items():
        weight = F(count, factorial(n))
        observed_variance += weight * ((history[-1] - mean) / (n + 1))**2
        depth_mean = F(sum(depths), n + 1)
        depth_variance = F(sum(d * d for d in depths), n + 1) - depth_mean**2
        observed_innovation += weight * depth_variance / (n + 2)**2
    expected_innovation = ((2 * harmonic - 1) / ((n + 1) * (n + 2))
                           + F(2, (n + 1)**2) - F(6, (n + 2)**2))
    assert observed_variance == variance / (n + 1)**2
    assert observed_innovation == expected_innovation

for history in ((0, 1, 3, 6, 10, 11), (0, 1, 3, 6, 9, 11)):
    count = sum(v for (depths, h), v in states.items() if h == history)
    depth_counts = Counter()
    for (depths, h), multiplicity in states.items():
        if h == history:
            for d in depths:
                depth_counts[d] += multiplicity
    print(history, "probability", F(count, factorial(6)),
          "next-depth law",
          {d: F(v, 7 * count) for d, v in sorted(depth_counts.items())})

print("Exact recurrence, variance, innovation, and history checks passed.")
```

These checks passed for every size 1 through 6. The counterexample laws were exactly those displayed in Section 3. The variance table and Gaussian-surrogate skewness were separately evaluated from the displayed formulas.
