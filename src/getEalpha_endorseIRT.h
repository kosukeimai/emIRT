// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEALPHA_ENDORSEIRT_H
#define GETEALPHA_ENDORSEIRT_H

#include <RcppArmadillo.h>

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
                           ) ;


#endif
