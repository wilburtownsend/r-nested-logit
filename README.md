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
\right).
```

(These random variables are sometimes called 'nested logit' random variables, because they are the latent variables which generate the nested logit model.)

For methods, see Appendix F.2 of my [job market paper](https://wilburtownsend.github.io/papers/visas.pdf). (Please cite that paper if you use this function in research!) In short, the function uses a Fourier transform to approximate the required PDF, and then uses inverse-transform sampling to draw from that PDF.
