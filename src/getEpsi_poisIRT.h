// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEPSI_POISIRT_H
#define GETEPSI_POISIRT_H

#include <RcppArmadillo.h>

void getEpsi(arma::mat &Epsi,
				   const arma::mat &Exfull,
				   const arma::mat &onesNK,
                   const arma::mat &exi,
                   const arma::mat &xi,
                   const arma::mat &exixi,
                   const arma::mat &y,
                   const arma::mat &i,
                   const arma::mat &Ealpha,
                   const arma::mat &Ebeta,
                   const double psi_mu,
                   const double psi_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) ;


#endif
