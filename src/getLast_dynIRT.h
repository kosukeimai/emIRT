// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETLAST_DYNIRT_HPP
#define GETLAST_DYNIRT_HPP

#include <RcppArmadillo.h>

arma::mat getLast_dynIRT(const arma::mat &bill_session,
					const int T,
					const int J);


#endif
