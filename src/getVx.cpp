// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>
#include "getEbb.h"

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getVx(const arma::mat &Eb,
                const arma::mat &Ebt,
                const arma::mat &Ebb,
                const arma::mat &sigma
                ) {
    arma::mat Vx = inv_sympd(inv_sympd(sigma) + Ebb) ;
    return(Vx) ;
}
