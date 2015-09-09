// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEX2X2_ORDIRT_H
#define GETEX2X2_ORDIRT_H

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat getEx2x2_ordIRT(const arma::mat &Ex,
                   const arma::mat &Vx,
                   const int N) ;
#endif
