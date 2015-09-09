// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEB2_DYNIRT_H
#define GETEB2_DYNIRT_H

#include <RcppArmadillo.h>

void getEb2_dynIRT(arma::mat &Eb2,
				 const arma::mat &Eystar,
                 const arma::mat &Ex,
                 const arma::cube &Vb2,
                 const arma::mat &bill_session,
                 const arma::mat &mubeta,
                 const arma::mat &sigmabeta,
                 const arma::mat &ones_col,
                 const int J
                 );

#endif
