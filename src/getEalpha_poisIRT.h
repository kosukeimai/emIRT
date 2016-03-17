// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEA_POISIRT_H
#define GETEA_POISIRT_H

#include <RcppArmadillo.h>

void getEalpha(arma::mat &Ealpha,
				   arma::mat &Ex,
				   const arma::mat &onesNJ,
                   const arma::mat &exi,
                   const arma::mat &xi,
                   const arma::mat &exixi,
                   const arma::mat &y,
                   const arma::mat &i,
                   const arma::mat &Epsi,
                   const arma::mat &Ebeta,
                   const double alpha_mu,
                   const double alpha_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) ;

#endif
