# Quicksort's Lognormal Fit and Non-Normal Limit

[Read the paper online](https://jonland82.github.io/complexity-distributions/) · [Download the PDF](quicksort_complexity_distribution.pdf)

Best case, worst case, and average case give three numbers for an algorithm's cost. They do not say how often the costs *between* those numbers occur. This paper follows that missing shape for the number of comparisons made by randomized-pivot Quicksort.

## The story

**Random splits create a distribution.** For $n$ distinct keys, Quicksort compares the pivot with the other $n-1$ keys, then recurses on two smaller groups. If $C_n$ is the total comparison count,

$$
C_n \overset{d}{=} n-1+C_{I_n}+C'_{n-1-I_n}.
$$

Here $I_n$ is uniformly distributed over the integers from $0$ to $n-1$, and the prime marks an independent recursive copy. Starting with $C_0=C_1=0$, this recurrence describes the *whole distribution*, not only its average.

**A lognormal curve wins the finite-size comparison.** Comparison counts have a longer right side: very uneven splits can make a run costly. The paper compares an ordinary normal fit with an ordinary lognormal fit using the expected log-density assigned to the true counts. For a positive random variable $X$, the difference between the best scores is

$$
\Delta(X)
=\frac12\log\left(\frac{\mathrm{Var}(X)}{\mathrm{Var}(\log X)}\right)
-\mathbb E[\log X].
$$

The main finite-size result is $\Delta(C_n)>0$ for **every integer $n\ge4$**. Thus the lognormal beats the normal under this scoring rule at every such size. The first figure shows what both fitted curves look like over experimental histograms; the theorem concerns the exact distribution rather than those samples.

**The winning fit still misses the limiting shape.** As inputs grow, the mean comparison count grows roughly like $n\log n$, while its standard deviation grows like $n$. A lognormal with the same mean and variance therefore becomes nearly normal once centered and scaled. Quicksort's centered and scaled count does something else:

$$
\widehat C_n=\frac{C_n-\mu_n}{\sqrt{v_n}}\Rightarrow Q,
\qquad
\widehat L_n=\frac{L_n-\mu_n}{\sqrt{v_n}}\Rightarrow\mathcal N(0,1),
$$

where $\mu_n=\mathbb E[C_n]$, $v_n=\mathrm{Var}(C_n)$, and $L_n$ is the mean-and-variance-matched lognormal. The Quicksort limit $Q$ keeps a rightward skew of about $0.855$; the standardized lognormal's skew tends to zero. Their largest cumulative-probability gap stays positive. The second figure shows that gap at two finite sizes and in the limit.

The point is a useful distinction: **a model can fit better at every finite size and still miss the shape that remains after standardization.** The claims here concern the uniform-pivot comparison-count model described above. Proofs and the exact certificate specification are in the appendices.

## Read and reproduce

- [Live HTML paper](https://jonland82.github.io/complexity-distributions/), with native MathML and section navigation
- [PDF paper](quicksort_complexity_distribution.pdf)
- [LaTeX source](quicksort_complexity_distribution.tex)

The website is built from the LaTeX source. Install [Pandoc](https://pandoc.org/installing.html) and the Python packages in `requirements.txt`, then run:

```sh
python -m pip install -r requirements.txt
python make_figures.py
python build_html.py
```

`make_figures.py` regenerates the plots from fixed random seeds. `build_html.py` writes both `index.html` for GitHub Pages and `quicksort_complexity_distribution.html` for local reading. The site uses local SVG figures and native MathML, so reading it does not require a network connection.
