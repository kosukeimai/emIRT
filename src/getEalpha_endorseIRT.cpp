// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

void getEalpha_endorseIRT (const arma::mat &ystar,
                           const arma::mat &beta,
                           const arma::mat &theta,
                           const arma::mat &w,
                           const arma::mat &gamma,
                           const arma::mat &mu,
                           const arma::mat &sigma,
                           const int N,
                           const int J,
                           arma::mat &Ealpha,
                           arma::mat &Valpha,
                           const arma::mat &theta2,
                           const arma::mat &w2
                           ) {
    // arma::mat Ealpha(J, 1) ;
    // Ealpha.fill(0.0) ;

    Valpha.fill(pow((N + 1/sigma(0, 0)), -1)) ;

#pragma omp parallel for
    for (int j = 0 ; j < J ; j++) {

        double q1 = mu(0, 0) / sigma(0, 0) ;

        for (int n = 0 ; n < N ; n++) {

            q1 = q1 +
                (ystar(n, j) -
                 beta(n, 0) +
                 gamma(0, 0) * (pow(theta(n, 0), 2)
                                - 2 * theta(n, 0) * w(j, 0)
                                + pow(w(j, 0), 2)
                                )
                 ) ;

        }

        Ealpha(j, 0) = Valpha(j, 0) * (q1) ;
    }

    // return(Ealpha) ;
}
