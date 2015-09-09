// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef CHECK_DYNIRT_H
#define CHECK_DYNIRT_H

#include <RcppArmadillo.h>

int checkConv_dynIRT(const arma::mat &oldEx,
              const arma::mat &curEx,
			  const arma::mat &oldEb,
              const arma::mat &curEb,
			  const arma::mat &oldEa,
              const arma::mat &curEa,
              double thresh,
              int convtype
               );

#endif
