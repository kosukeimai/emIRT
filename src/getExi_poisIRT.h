// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEXI_POISIRT_H
#define GETEXI_POISIRT_H

#include <RcppArmadillo.h>

void getExi(arma::mat &exi,
					arma::mat &xi,
					arma::mat &exixi,
				   const arma::mat &Ealpha,
				   const arma::mat &Epsi,
                   const arma::mat &Ebeta,
                   const arma::mat &Ex,
                   const arma::mat &i,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) ;

#endif
