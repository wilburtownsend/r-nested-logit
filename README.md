## rnestedlogit

This repo contains an R function for drawing `nested extreme value' random variables --- i.e., those with this CDF:
```math
\mathbb{P}
\left[\left(X_{j}\right)_{j \in \mathbf{J}} \leq 
\left(x_{j}\right)_{j \in \mathbf{J}}\right]
= \exp\left(
\sum_{C \in \mathbf{C}} 
\left(
\sum_{j \in C} \exp\left(\frac{-x_{j}}{1-\sigma}\right)
\right)^{1-\sigma}
\right),
```
where $\mathbf{C}$ is a partition on $\mathbf{J}$ and the parameter $σ$ measures within-group correlation

(These random variables are sometimes called 'nested logit' random variables, because they are the latent variables which generate the nested logit model.)

For methods, see Appendix F.2 of my [job market paper](https://wilburtownsend.github.io/papers/visas.pdf). (Please cite that paper if you use this function in research!) In short, the random variable can be decomposed into a nest-level component plus an IID component. The IID component is a scaled Gumbel, so the challenge is drawing the nest-level component. To do so, we approximate its PDF using a Fourier transform.

The function requires the libraries `fourierin`, `pracma` and `extraDistr`.

Typical usage is 
```
X = rnestedlogit(N, σ, nests),
```
where `N` is the number of vectors you want to draw, `σ` is the coefficient measuring within-nest correlation, and `nests` is a vector indicating, for each product, the nest of that product (which must be coded as an integer). The output `X` is an `N`-by-`length(nests)` matrix.

There are additional optional keyword arguments. The PDF of the nest-level component is approximated with a Fourier transform, and the approximate PDF can be negative, imaginary, or it can be positive at the boundary of its domain. The function ensures these failures are not too egregious, with tolerance dictated by the keyword `tol`. If the `stopifnot()` assertions are causing the function to fail (as is likely with a very low value of `σ`), try setting `plot_pdf=TRUE` and/or `plot_imaginary=TRUE`, to plot the approximate PDF. You may then want to play with the following keyword arguments, which are passed directly to the `fourierin` package: `lower_int` (default -50), `upper_int` (default 50), `lower_eval` (default -5), `upper_eval` (default 10), and `resolution` (default 2^15).

If you have any problems, feel free to drop a Github issue.