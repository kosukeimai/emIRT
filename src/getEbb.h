// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat getEbb(const arma::mat &Eb,
                 const arma::mat &Vb,
                 const int J,
                 const int D
                 ) ;
