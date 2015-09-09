// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEba(const arma::mat &Eb,
                 const arma::mat &Ea,
                 const arma::mat &Vb2,
                 const int J,
                 const int D,
                 const bool asEM
                 ) {
    arma::mat Eba(D, 1) ;

    // Vb2.print("getEba Vb2") ;

    for (int d = 0 ; d < D ; d++) {
        double q = 0.0 ;
        for (int j = 0 ; j < J ; j++) {
            q += Eb(j, d) * Ea(j, 0) ;
            // Just VI
            if (!asEM) {
                // Rcout << Vb2(0, d + 1) << std::endl ;
                q += Vb2(0, d + 1) ;
            }
        }
        Eba(d, 0) = q ;
    }

    // Eba.print("Eba") ;
    return(Eba) ;
}
