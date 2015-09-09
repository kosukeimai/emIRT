// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEb2(const arma::mat &Eystar,
                 const arma::mat &Ex,
                 const arma::mat &Vb2,
                 const arma::mat &mubeta,
                 const arma::mat &sigmabeta,
                 const int J,
                 const int D,
                 const bool asEM
                 ) {

    arma::mat ones(Ex.n_rows, 1) ;
    ones.ones() ;

    arma::mat Eb2(J, D + 1) ;

    arma::mat Ex2 = Ex ;
    Ex2.insert_cols(0, ones) ;

    arma::mat B = Vb2 ;

    if (asEM) {
        B = inv_sympd(inv_sympd(sigmabeta) + trans(Ex2) * Ex2) ;
    }

#pragma omp parallel for
    for (int j = 0; j < J; j++) {
        Eb2.row(j) = trans(B * (inv_sympd(sigmabeta) * mubeta + trans(Ex2) * Eystar.col(j))) ;
    }

    return(Eb2) ;
}
