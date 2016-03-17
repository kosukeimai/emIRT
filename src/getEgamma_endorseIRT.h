// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEGAMMA_ENDORSEIRT_H
#define GETEGAMMA_ENDORSEIRT_H

#include <RcppArmadillo.h>

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
                                ) ;

#endif
