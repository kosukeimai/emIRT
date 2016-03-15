// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEB_POISIRT_H
#define GETEB_POISIRT_H

#include <RcppArmadillo.h>

void getEbeta(arma::mat &Ebeta,
				   const arma::mat &Vbeta,
				   const arma::mat &Ealpha,
				   const arma::mat &onesNK,
                   const arma::mat &exi,
                   const arma::mat &xi,
                   const arma::mat &exixi,
                   const arma::mat &y,
                   const arma::mat &Epsi,
                   const arma::mat &Exfull,
                   const arma::mat &Vx,
                   const arma::mat &i,
                   const double beta_mu,
                   const double beta_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) ;

#endif
