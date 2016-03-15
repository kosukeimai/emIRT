// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVB_POISIRT_H
#define GETVB_POISIRT_H

#include <RcppArmadillo.h>

void getVbeta(arma::mat &Vbeta,
				   const arma::mat &Ex,
				   const arma::mat &Vx,
                   const arma::mat &exi,
                   const arma::mat &i,
				   const double beta_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) ;

#endif
