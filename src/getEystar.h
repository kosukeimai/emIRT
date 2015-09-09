// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETYESTAR_H
#define GETYESTAR_H

#include <RcppArmadillo.h>

arma::mat getEystar(const arma::mat &alpha,
                    const arma::mat &beta,
                    const arma::mat &x,
                    const arma::mat &y,
                    const int D,
                    const int N,
                    const int J
                    ) ;

#endif
