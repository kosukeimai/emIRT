# emIRT [![Build Status](https://travis-ci.org/kosukeimai/emIRT.svg?branch=master)](https://travis-ci.org/kosukeimai/emIRT)
emIRT: EM Algorithms for Estimating Item Response Theory Models


Various Expectation-Maximization (EM) algorithms are implemented
for item response theory (IRT) models. The current implementation includes IRT
models for binary and ordinal responses, along with dynamic and hierarchical IRT
models with binary responses. The latter two models are derived and implemented
using variational EM.

The current release of the R package is available on
[CRAN](https://cran.r-project.org/web/packages/emIRT/).

The github versions of the R package are available with

    library("devtools")
    install_github("kosukeimai/emIRT")
    install_github("kosukeimai/emIRT", ref ="development")
