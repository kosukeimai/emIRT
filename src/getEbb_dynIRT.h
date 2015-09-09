// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEBB_DYNIRT_H
#define GETEBB_DYNIRT_H

#include <RcppArmadillo.h>

arma::mat getEbb_dynIRT(const arma::mat &Eb,
				 const arma::cube &Vb2,
                 const arma::mat &bill_session,
                 const int nJ
                 );

#endif
