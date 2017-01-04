# emIRT: EM Algorithms for Estimating Item Response Theory Models [![Build Status](https://travis-ci.org/kosukeimai/emIRT.svg?branch=master)](https://travis-ci.org/kosukeimai/emIRT)

Various Expectation-Maximization (EM) algorithms are implemented for item
response theory (IRT) models. The current implementation includes IRT models for
binary and ordinal responses, along with dynamic and hierarchical IRT models
with binary responses. Variational algorithms for scaling network and text data
are also included.  The details about the methods implemented in this package are 
described in Imai, Lo, and Olmsted. (2016). "[Fast Estimation of Ideal Points with Massive Data.](http://imai.princeton.edu/research/fastideal.html)" *American Political Science Review*, Vol. 110, No. 4 
    (December), pp. 631-656.

The current release of the R package is available on
[CRAN](https://cran.r-project.org/web/packages/emIRT/).

The github versions of the R package are available with

    library("devtools")
    install_github("kosukeimai/emIRT")
    install_github("kosukeimai/emIRT", ref ="development")
