// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEX_POISIRT_H
#define GETEX_POISIRT_H

#include <RcppArmadillo.h>

void getEx(arma::mat &Ex,
				   const arma::mat &Vx,
				   const arma::mat &Ealpha,
				   const arma::mat &onesNJ,
                   const arma::mat &exi,
                   const arma::mat &xi,
                   const arma::mat &exixi,
                   const arma::mat &y,
                   const arma::mat &Epsi,
                   const arma::mat &Ebeta,
                   const arma::mat &Vbeta,
                   const arma::mat &i,
                   const double x_mu,
                   const double x_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) ;

#endif
