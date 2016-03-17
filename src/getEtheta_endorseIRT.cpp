// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

void getEtheta_endorseIRT (const arma::mat &ystar,
                           const arma::mat &alpha,
                           const arma::mat &beta,
                           const arma::mat &w,
                           const arma::mat &gamma,
                           const arma::mat &mu,
                           const arma::mat &sigma,
                           const int N,
                           const int J,
                           const arma::mat &oldEtheta,
                           arma::mat &Etheta,
                           arma::mat &Vtheta,
                           const arma::mat &w2,
                           const arma::mat &w3,
                           const arma::mat &gamma2
                           ) {
    // arma::mat Etheta(N, 1) ;
    // Etheta.fill(0.0) ;

#pragma omp parallel for
    for (int n = 0 ; n < N ; n++) {
        double v1 = 2.0 / sigma(0, 0) ;

        double q1 = (2 / sigma(0, 0)) * (oldEtheta(n, 0) - mu(0, 0)) ;

        for (int j = 0 ; j < J ; j++) {

            // 2nd Deriv Term

            v1 = v1 +
                4 * (
                     (gamma(0, 0) *
                      (ystar(n, j)
                       - alpha(j, 0)
                       - beta(n, 0)
                       )
                      ) +
                     pow(gamma(0, 0), 2) *
                     (
                      3 * pow(oldEtheta(n, 0), 2)
                      - 6 * oldEtheta(n, 0) * w(j, 0)
                      + 3 * pow(w(j, 0), 2)
                      )
                     ) ;


            // 1st Deriv Term

            q1 = q1 +
                4 * (
                     (
                      (gamma(0, 0)) *
                      (oldEtheta(n, 0) - w(j, 0)) *
                      (ystar(n, j) - alpha(j, 0) - beta(n, 0))
                      ) +
                     (pow(gamma(0, 0), 2) *
                      (pow(oldEtheta(n, 0), 3)
                       - 3 * pow(oldEtheta(n, 0), 2) * w(j, 0)
                       + 3 * oldEtheta(n, 0) * pow(w(j, 0), 2)
                       - pow(w(j, 0), 3)
                       )
                      )
                     ) ;

        }
        double T = pow(v1 / 2, -1.0) ;
        double t = T * (oldEtheta(n, 0) * v1 - q1) / 2 ;

        Etheta(n, 0) = t ;
        Vtheta(n, 0) = T ;
    }

    // Returns
    // return(Etheta) ;
}
