// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVX_ORDIRT_HPP
#define GETVX_ORDIRT_HPP

#include <RcppArmadillo.h>

arma::mat getVx_ordIRT(const arma::mat &Ebb,
                const arma::mat &Edd,
                const arma::mat &sigma2_x
                );
                
#endif
