// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEBB_ORDIRT_H
#define GETEBB_ORDIRT_H

#include <RcppArmadillo.h>

arma::mat getEbb_ordIRT(const arma::mat &Eb,
                const arma::mat &Vb,
                const int nJ
                );

#endif
