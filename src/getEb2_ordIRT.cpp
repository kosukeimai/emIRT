// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEb2_ordIRT(const arma::mat &Ezstar,
                 const arma::mat &Ex,
                 const arma::cube &Vb2,
                 const arma::mat &mubeta,
                 const arma::mat &sigmabeta,
                 const arma::mat &Edd,
                 const int J
                 ) {

    arma::mat ones(Ex.n_rows, 1) ;
    ones.ones() ;
    arma::mat Ex2 = Ex ;
    Ex2.insert_cols(0, ones) ;

    arma::mat Eb2(J, 2) ;

#pragma omp parallel for
    for (int j = 0; j < J; j++) {
        Eb2.row(j) = trans(Vb2.slice(j) * (inv_sympd(sigmabeta) * mubeta + Edd(j,0) * trans(Ex2) * Ezstar.col(j))) ;
    }

    return(Eb2) ;
}
