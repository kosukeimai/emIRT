// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVZSTAR_ORDIRT_HPP
#define GETVZSTAR_ORDIRT_HPP

#include <RcppArmadillo.h>

arma::mat getVzstar_ordIRT(const arma::mat &Edd,
                    const arma::mat &beta,
                    const arma::mat &x,
                    const arma::mat &Etau,
                    const arma::mat &y,
                    const int N,
                    const int J
                    );

#endif
