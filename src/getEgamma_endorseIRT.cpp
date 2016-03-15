// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat getEgamma_endorseIRT (const arma::mat &ystar,
                                const arma::mat &alpha,
                                const arma::mat &beta,
                                const arma::mat &theta,
                                const arma::mat &w,
                                const arma::mat &mu,
                                const arma::mat &sigma,
                                const int N,
                                const int J,
                                arma::mat &Vgamma,
                                const arma::mat &theta2,
                                const arma::mat &theta3,
                                const arma::mat &theta4,
                                const arma::mat &w2,
                                const arma::mat &w3,
                                const arma::mat &w4
                                ) {
    arma::mat Egamma(1, 1) ;
    Egamma.fill(0.0) ;

    double v1 = 1 / sigma(0, 0) ;
    double q1 = mu(0, 0) / sigma(0, 0) ;

    // #pragma omp parallel for
    for (int n = 0 ; n < N ; n++) {
        for (int j = 0 ; j < J ; j++) {

            v1 = v1 +
                (
                 pow(theta(n, 0), 4) -
                 - 4 * pow(theta(n, 0), 3) * w(j, 0)
                 + 6 * pow(theta(n, 0), 2) * pow(w(j, 0), 2)
                 - 4 * theta(n, 0) * pow(w(j, 0), 3)
                 + pow(w(j, 0), 4)
                 ) ;

            // v1 = v1 +
            //     (
            //      theta4(n, 0)
            //      - 4 * theta3(n, 0) * w(j, 0)
            //      + 6 * theta2(n, 0) * w2(j, 0)
            //      - 4 * theta(n, 0) * w3(j, 0)
            //      + w4(j, 0)
            //      ) ;


            q1 = q1 +
                (pow(theta(n, 0), 2)
                 - 2 * theta(n, 0) * w(j, 0)
                 + pow(w(j, 0), 2)
                 ) * (alpha(j, 0)
                      + beta(n, 0)
                      - ystar(n, j)
                      ) ;
        }
    }
    double G = pow(v1, -1) ;
    double g = G * q1 ;

    Egamma(0, 0) = g ;
    Egamma(0, 0) = 1.0 ;

    Vgamma(0, 0) = G ;
    Vgamma(0, 0) = 0.0 ;

    return(Egamma) ;
}
