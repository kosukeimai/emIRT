// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEBETA_ENDORSEIRT_H
#define GETEBETA_ENDORSEIRT_H

#include <RcppArmadillo.h>

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
                          ) ;


#endif
