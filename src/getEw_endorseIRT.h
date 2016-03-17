// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEW_ENDORSEIRT_H
#define GETEW_ENDORSEIRT_H

#include <RcppArmadillo.h>

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
                       ) ;


#endif
