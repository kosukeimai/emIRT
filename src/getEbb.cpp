// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEbb(const arma::mat &Eb,
                 const arma::mat &Vb,
                 const int J,
                 const int D
                 ) {
    arma::mat tmp(D, D) ;
    tmp.fill(J) ;
    arma::mat Ebb = tmp % Vb + Eb.t() * Eb ;
    // Ebb.print("Ebb") ;
    // (tmp % Vb).print("vi-spec cor Ebb") ;
    return(Ebb) ;
}
