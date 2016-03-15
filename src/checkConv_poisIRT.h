// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef CHECK_POISIRT_HPP
#define CHECK_POISIRT_HPP

#include <RcppArmadillo.h>

int checkConv_poisIRT(const arma::mat &oldEalpha,
              const arma::mat &curEalpha,
			  const arma::mat &oldEpsi,
              const arma::mat &curEpsi,
			  const arma::mat &oldEbeta,
              const arma::mat &curEbeta,
			  const arma::mat &oldEx,
              const arma::mat &curEx,
              double thresh,
              int convtype
              ) ;

               
#endif
