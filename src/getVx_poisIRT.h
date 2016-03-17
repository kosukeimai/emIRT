// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVX_POISIRT_H
#define GETVX_POISIRT_H

#include <RcppArmadillo.h>

void getVx(arma::mat &Vx,
				   const arma::mat &Ebeta,
				   const arma::mat &Vbeta,
                   const arma::mat &exi,
                   const arma::mat &i,
                   const double x_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) ;
#endif
