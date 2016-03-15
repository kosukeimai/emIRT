// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETYESTAR_ENDORSEIRT_H
#define GETYESTAR_ENDORSEIRT_H

#include <RcppArmadillo.h>

void getEystar_endorseIRT (const arma::mat &alpha,
                           const arma::mat &beta,
                           const arma::mat &w,
                           const arma::mat &theta,
                           const arma::mat &gamma,
                           const arma::mat &y,
                           const int N,
                           const int J,
                           arma::mat &ystars,
                           const arma::mat &theta2,
                           const arma::mat &w2
                           ) ;

#endif
