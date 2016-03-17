// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>
#include "etn1.h"

using namespace Rcpp;

void getEystar_endorseIRT (const arma::mat &alpha,
                           const arma::mat &beta,
                           const arma::mat &w,
                           const arma::mat &theta,
                           const arma::mat &gamma,
                           const arma::mat &y,
                           const int N,
                           const int J,
                           arma::mat &ystars,
                           const arma::mat &theta2,
                           const arma::mat &w2
                           ) {
    // arma::mat ystars(N, J, arma::fill::zeros) ;

    // Main Calculation

#pragma omp parallel for
    for (int n = 0; n < N; n++) {

        arma::mat theseystars = arma::mat(1, J, arma::fill::zeros) ;

        for (int j = 0; j < J; j++) {
            // defaults to untruncated
            double low = y(n,j) == 1 ? 0.0 : R_NegInf ;
            double high = y(n,j) == 0 ? 0.0 : R_PosInf ;

            // mu

            double mu =
                alpha(j, 0) +
                beta(n, 0) -
                gamma(0, 0) * (
                               pow(theta(n, 0), 2) +
                               pow(w(j, 0), 2) -
                               2 * theta(n, 0) * w(j, 0)
                               ) ;

            double exp = etn1(mu, 1.0, low, high) ;

            double exp2 = exp ;

            // if (isinf(exp2) || isnan(exp2)) {
            //     if (low == 0.0) exp2 =   .0001 ;
            //     if (high == 0.0) exp2 = -.0001 ;
            // }


            // if (isnan(exp2)) {
            //     exp2 = 0 ;
            // }

            // if (isinf(exp)) {
            //     Rcout << n << " " << j << " " << low << " " << high << " " << mu << " " << exp << " " << exp2 << std::endl ;
            // }

            // if (isnan(exp)) {
            //     Rcout << n << " " << j << " " << low << " " << high << " " << mu << " " << exp << " " << exp2 << std::endl ;
            // }

            theseystars(0, j) = exp2 ;
        }
        ystars.row(n) = theseystars ;
    }

    arma::uvec check = find_nonfinite(ystars) ;

    if (check.n_elem > 0) {
        ystars.elem(check).print("ystar check") ;
    }

    // return(ystars) ;
}
