// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEx(const arma::mat &Eystar,
                const arma::mat &Eb, // framework dependent
                const arma::mat &Vx,
                const arma::mat &Eba, // framework dependent
                const arma::mat &mu,
                const arma::mat &sigma,
                const int N,
                const int J,
                const int D,
                const bool asEM
                ) {
    arma::mat Ebt = trans(Eb) ;
    arma::mat x(N, D) ;
    x.fill(0.0) ;

    arma::mat A = Vx ;
    // A.print("getEx Avi") ;

    if (asEM) {
        // Ebt.print("ebt") ;
        // (Ebt * Eb).print("ebteb") ;

        arma::mat tmp = inv_sympd(sigma) + Ebt * Eb ;
        //  tmp.print("tmp") ;

// Before v0.0.9
//        arma::mat A = inv_sympd(tmp) ;

		// Following Rothenberg, Jo, Seo's testing, v. 0.0.11
        A = inv_sympd(tmp);
        // A.print("getEx Aem") ;
    }

#pragma omp parallel for
    for (int n = 0 ; n < N; n++) {
        x.row(n) = trans(A * (inv_sympd(sigma) * mu +
                               Ebt * trans(Eystar.row(n)) - Eba
                               )
                         ) ;
    }
    return(x) ;
}
