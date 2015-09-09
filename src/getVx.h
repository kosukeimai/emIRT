// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVX_HPP
#define GETVX_HPP

#include <RcppArmadillo.h>

arma::mat getVx(const arma::mat &Eb,
                const arma::mat &Ebt,
                const arma::mat &Ebb,
                const arma::mat &sigma
                ) ;

#endif
