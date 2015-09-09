// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef CHECK_HIERIRT_H
#define CHECK_HIERIRT_H

#include <RcppArmadillo.h>

int checkConv_hierIRT(const arma::mat &oldEa,
              const arma::mat &curEa,
			  const arma::mat &oldEb,
              const arma::mat &curEb,
			  const arma::mat &oldEgamma,
              const arma::mat &curEgamma,
			  const arma::mat &oldEx,
              const arma::mat &curEx,
			  int ND,
              double thresh,
              int convtype
               );

#endif
