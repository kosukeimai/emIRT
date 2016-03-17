// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

void getEw_endorseIRT (const arma::mat &ystar,
                       const arma::mat &alpha,
                       const arma::mat &beta,
                       const arma::mat &theta,
                       const arma::mat &gamma,
                       const arma::mat &mu,
                       const arma::mat &sigma,
                       const int N,
                       const int J,
                       const arma::mat &oldEw,
                       arma::mat &Ew,
                       arma::mat &Vw,
                       const arma::mat &gamma2,
                       const arma::mat &theta2,
                       const arma::mat &theta3
                       ) {
    // arma::mat Ew(J, 1) ;
    // Ew.fill(0.0) ;

#pragma omp parallel for
    for (int j = 0 ; j < J ; j++) {
        double v1 = 2.0 / sigma(0, 0) ;

        double q1 = (2 / sigma(0, 0)) * (oldEw(j, 0) - mu(0, 0)) ;

        for (int n = 0 ; n < N ; n++) {

            // 2nd Deriv Term
            v1 = v1 +
                4 * (
                     (gamma(0, 0) *
                      (ystar(n, j)
                       - alpha(j, 0)
                       - beta(n, 0)
                       )
                      ) +
                     (3 *
                      pow(gamma(0, 0), 2) *
                      (
                       (pow(theta(n, 0), 2)
                        - 2 * oldEw(j, 0) * theta(n, 0)
                        + pow(oldEw(j, 0), 2)
                        )
                       )
                      )
                     )  ;

            // 1st Deriv Term
            q1 = q1 -
                4 * (
                     (
                      (gamma(0, 0)) *
                      (theta(n, 0) - oldEw(j, 0)) *
                      (ystar(n, j) - alpha(j, 0) - beta(n, 0))
                     ) +
                     (pow(gamma(0, 0), 2) *
                      (pow(theta(n, 0), 3)
                       - 3 * pow(theta(n, 0), 2) * oldEw(j, 0)
                       + 3 * theta(n, 0) * pow(oldEw(j, 0), 2)
                       - pow(oldEw(j, 0), 3)
                       )
                      )
                     ) ;

        }
        double T = pow(v1 / 2, -1.0) ;
        double t = T * (oldEw(j, 0) * v1 - q1) / 2 ;

        Ew(j, 0) = t ;
        Vw(j, 0) = T ;
    }

    // Returns
    // return(Ew) ;
}
