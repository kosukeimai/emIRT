// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVBETA2_DYNIRT_HPP
#define GETVBETA2_DYNIRT_HPP

#include <RcppArmadillo.h>

void getVb2_dynIRT(arma::cube &Vb2,
				 const arma::cube &Ex2x2,
                 const arma::mat &sigma,
                 const int T
                 );

#endif
