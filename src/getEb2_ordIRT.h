// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEB2_ORDIRT_H
#define GETEB2_ORDIRT_H

#include <RcppArmadillo.h>

arma::mat getEb2_ordIRT(const arma::mat &Ezstar,
                 const arma::mat &Ex,
                 const arma::cube &Vb2,
                 const arma::mat &mubeta,
                 const arma::mat &sigmabeta,
                 const arma::mat &Edd,
                 const int J
                 );

#endif
