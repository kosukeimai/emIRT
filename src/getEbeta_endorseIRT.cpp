// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

void getEbeta_endorseIRT (const arma::mat &ystar,
                          const arma::mat &alpha,
                          const arma::mat &theta,
                          const arma::mat &w,
                          const arma::mat &gamma,
                          const arma::mat &mu,
                          const arma::mat &sigma,
                          const int N,
                          const int J,
                          arma::mat &Ebeta,
                          arma::mat &Vbeta,
                          const arma::mat theta2,
                          const arma::mat w2
                          ) {
    // arma::mat Ebeta(N, 1) ;
    // Ebeta.fill(0.0) ;

    // arma::mat Vbeta(1, 1) ;
    // Vbeta(0, 0) = pow((J + 1/sigma(0, 0)), -1) ;

    Vbeta.fill(pow((J + 1/sigma(0, 0)), -1)) ;

#pragma omp parallel for
    for (int n = 0 ; n < N ; n++) {

        double q1 = mu(0, 0) / sigma(0, 0) ;

        for (int j = 0 ; j < J ; j++) {

            q1 = q1 +
                (ystar(n, j) -
                 alpha(j, 0) +
                 gamma(0, 0) * (pow(theta(n, 0), 2)
                                - 2 * theta(n, 0) * w(j, 0)
                                + pow(w(j, 0), 2)
                                )
                 ) ;

        }

        Ebeta(n, 0) = Vbeta(n, 0) * (q1) ;

    }

    // return(Ebeta) ;
}
