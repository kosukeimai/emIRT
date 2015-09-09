// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEX_H
#define GETEX_H

#include <RcppArmadillo.h>

arma::mat getEx(const arma::mat &Eystar,
                const arma::mat &Eb,
                const arma::mat &Vx,
                const arma::mat &Eba,
                const arma::mat &mu,
                const arma::mat &sigma,
                const int N,
                const int J,
                const int D,
                const bool asEM
                ) ;

#endif
