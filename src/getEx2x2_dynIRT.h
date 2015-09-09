// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEX2X2_DYNIRT_H
#define GETEX2X2_DYNIRT_H

#include <RcppArmadillo.h>

void getEx2x2_dynIRT(arma::cube &Ex2x2,
				   const arma::mat &Ex,
                   const arma::mat &Vx,
                   const arma::mat &legis_by_session,
                   const arma::mat &Nlegis_session,
                   const int T
                   );

#endif
