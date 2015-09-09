// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEZSTAR_ORDIRT_H
#define GETEZSTAR_ORDIRT_H

#include <RcppArmadillo.h>

arma::mat getEzstar_ordIRT(const arma::mat &Edd,
                    const arma::mat &beta,
                    const arma::mat &x,
                    const arma::mat &Etau,
                    const arma::mat &y,
                    const int N,
                    const int J
                    );

#endif
