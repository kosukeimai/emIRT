// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVBETA2_ORDIRT_HPP
#define GETVBETA2_ORDIRT_HPP

#include <RcppArmadillo.h>

arma::cube getVb2_ordIRT(const arma::mat &Ex,
                 const arma::mat &Ex2x2,
                 const arma::mat &sigma,
                 const arma::mat &Edd,
                 const int nJ
                 );

#endif
