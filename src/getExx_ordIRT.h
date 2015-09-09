// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEXX_ORDIRT_H
#define GETEXX_ORDIRT_H

#include <RcppArmadillo.h>

arma::mat getExx_ordIRT(const arma::mat &Ex,
                const arma::mat &Vx,
                const int nN
                );

#endif
