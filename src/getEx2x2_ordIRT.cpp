// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEx2x2_ordIRT(const arma::mat &Ex,
                   const arma::mat &Vx,
                   const int N) {

    arma::mat Ex2x2(2, 2, arma::fill::zeros) ;
    arma::mat tmp(1, 1) ;
    tmp.fill(N) ;

    arma::mat S = tmp % Vx + Ex.t() * Ex ;

    Ex2x2(0, 0) = N ;
    Ex2x2(1, 1) = S(0,0) ;

    arma::mat sums = sum(Ex, 0) ;
    Ex2x2(1, 0) = sums(0,0) ;
    Ex2x2(0, 1) = sums(0,0) ;

    return(Ex2x2) ;

}
