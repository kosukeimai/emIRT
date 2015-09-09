// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getVb2(const arma::mat &Ex,
                 const arma::mat &Ex2x2,
                 const arma::mat &sigma
                 ) {
    arma::mat Vb2 = inv_sympd(inv_sympd(sigma) + Ex2x2) ;
    return(Vb2) ;
}
