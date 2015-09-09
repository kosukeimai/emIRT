// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETETT_ORDIRT_H
#define GETETT_ORDIRT_H

#include <RcppArmadillo.h>

arma::mat getEtt_ordIRT(const arma::mat &Etau,
                const arma::mat &Vtau,
                const int nJ
                );

#endif
