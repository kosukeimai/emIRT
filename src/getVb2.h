// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVBETA2_HPP
#define GETVBETA2_HPP

#include <RcppArmadillo.h>

arma::mat getVb2(const arma::mat &Ex,
                 const arma::mat &Ex2x2,
                 const arma::mat &sigma
                 ) ;

#endif
