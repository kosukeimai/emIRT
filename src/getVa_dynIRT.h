// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVALPHA_DYNIRT_HPP
#define GETVALPHA_DYNIRT_HPP

#include <RcppArmadillo.h>

arma::mat getVa_dynIRT(const arma::cube &Vb2,
                 const arma::mat &bill_session,
                 const int J
                 );

#endif
