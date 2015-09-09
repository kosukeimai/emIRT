// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETED_ORDIRT_H
#define GETED_ORDIRT_H

#include <RcppArmadillo.h>

arma::mat getEd_ordIRT(const arma::mat &Ex,
                 const arma::mat &Exx,
                 const arma::mat &Eb,
                 const arma::mat &Ebb,
                 const arma::mat &Etau,
                 const arma::mat &Ett,
                 const arma::mat &Ezstar,
                 const arma::mat &Ezzstar,
                 const int nN,
                 const int nJ
                 );

#endif
