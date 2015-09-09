// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEBA_H
#define GETEBA_H

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat getEba(const arma::mat &Eb,
                 const arma::mat &Ea,
                 const arma::mat &Vb2,
                 const int J,
                 const int D,
                 const bool asEM
                 ) ;

#endif
