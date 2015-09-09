// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef CHECKCONV_ORDIRT_H
#define CHECKCONV_ORDIRT_H

#include <RcppArmadillo.h>

int checkConv_ordIRT(const arma::mat &oldEx,
              const arma::mat &curEx,
			  const arma::mat &oldEb,
              const arma::mat &curEb,
			  const arma::mat &oldEtau,
              const arma::mat &curEtau,
			  const arma::mat &oldEdd,
              const arma::mat &curEdd,
              double thresh,
              int convtype
               );

#endif
