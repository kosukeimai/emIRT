// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETETHETA_ENDORSEIRT_H
#define GETETHETA_ENDORSEIRT_H

#include <RcppArmadillo.h>

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
                           ) ;


#endif
