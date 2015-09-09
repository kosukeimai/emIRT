// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat getEx2x2(const arma::mat &Ex,
                   const arma::mat &Vx,
                   const int N,
                   const int D
                   ) ;
