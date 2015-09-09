// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>
//#include <RcppTN.h>
#include "etn1.h"

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEystar(const arma::mat &alpha,
                    const arma::mat &beta,
                    const arma::mat &x,
                    const arma::mat &y,
                    const int D,
                    const int N,
                    const int J
                    ) {
    arma::mat ystars(N, J, arma::fill::zeros) ;


    // Main Calculation
#pragma omp parallel for
    for (int n = 0; n < N; n++) {
        arma::mat thisx = x.row(n) ;
        arma::mat theseystars = arma::mat(1, J, arma::fill::zeros) ;

        for (int j = 0; j < J; j++) {
            double q1 = 0.0 ;
            q1 += alpha(j, 0) ;
            for (int d = 0; d < D; d++) {
                q1 += thisx(0, d) * beta(j, d) ;
            }
            // defaults to untruncated

            double low = y(n,j) == 1 ? 0.0 : R_NegInf ;
            double high = y(n,j) == -1 ? 0.0 : R_PosInf ;

            // Rcout << n << ':' << j << std::endl ;
            // Rcout << q1 << " " << low << " " << high << std::endl ;
//            theseystars(0, j) = RcppTN::etn1(q1, 1.0, low, high) ;
            theseystars(0, j) = etn1(q1, 1.0, low, high) ;

        }
        ystars.row(n) = theseystars ;
    }

    return(ystars) ;
}
