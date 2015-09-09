// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETLBS_DYNIRT_HPP
#define GETLBS_DYNIRT_HPP

#include <RcppArmadillo.h>

arma::mat getLBS_dynIRT(const arma::mat &startlegis,
                 const arma::mat &endlegis,
                 const int T,
                 const int N);


#endif
