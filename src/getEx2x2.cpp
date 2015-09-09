// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEx2x2(const arma::mat &Ex,
                   const arma::mat &Vx,
                   const int N,
                   const int D
                   ) {
    arma::mat Ex2x2(D + 1, D + 1, arma::fill::zeros) ;
    arma::mat tmp(D, D) ;
    tmp.fill(N) ;
    arma::mat S = tmp % Vx + Ex.t() * Ex ;

    // S.print("S") ;
    // (tmp % Vx).print("tmp") ;

    Ex2x2(0, 0) = N ;
    Ex2x2.submat(1, 1, D, D) = S ;

    arma::mat sums = sum(Ex, 0) ;
    for (int d = 0; d < D ; d++) {
        Ex2x2(d + 1, 0) = sums(0,d) ;
        Ex2x2(0, d + 1) = sums(0,d) ;
    }
    return(Ex2x2) ;
}
