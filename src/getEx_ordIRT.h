// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEX_ORDIRT_H
#define GETEX_ORDIRT_H

#include <RcppArmadillo.h>

arma::mat getEx_ordIRT(const arma::mat &Ezstar,
                const arma::mat &Eb,
                const arma::mat &Etau,
                const arma::mat &Vx,
                const arma::mat &Edd,
                const arma::mat &xmu,
                const arma::mat &xsigma,
                const arma::cube &Vb2,
                const int N,
                const int J
                );
#endif
