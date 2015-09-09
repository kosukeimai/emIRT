// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEB2_H
#define GETEB2_H

#include <RcppArmadillo.h>

arma::mat getEb2(const arma::mat &Eystar,
                 const arma::mat &Ex,
                 const arma::mat &Vb2,
                 const arma::mat &mubeta,
                 const arma::mat &sigmabeta,
                 const int J,
                 const int D,
                 const bool asEM
                 ) ;

#endif
