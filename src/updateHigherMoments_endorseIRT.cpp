// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

void  updateHigherMoments_endorseIRT (const arma::mat &Etheta,
                                      const arma::mat &Vtheta,
                                      const arma::mat &Ew,
                                      const arma::mat &Vw,
                                      const arma::mat &Egamma,
                                      const arma::mat &Vgamma,
                                      arma::mat &Etheta2,
                                      arma::mat &Etheta3,
                                      arma::mat &Etheta4,
                                      arma::mat &Ew2,
                                      arma::mat &Ew3,
                                      arma::mat &Ew4,
                                      arma::mat &Egamma2
                                      ) {

    // theta
    Etheta2 = pow(Etheta, 2) + Vtheta ;

    Etheta3 = pow(Etheta, 3) +
        3 * Etheta % Vtheta ;

    Etheta4 = pow(Etheta, 4) +
        6 * pow(Etheta, 2) % Vtheta +
        3 * pow(Vtheta, 2) ;

    // w
    Ew2 = pow(Ew, 2) + Vw ;

    Ew3 = pow(Ew, 3) +
        3 * Ew % Vw ;

    Ew4 = pow(Ew, 4) +
        6 * pow(Ew, 2) % Vw +
        3 * pow(Vw, 2) ;

    // gamma
    Egamma2 = pow(Egamma, 2) + Vgamma ;
}
