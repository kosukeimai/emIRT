// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETNLEGIS_DYNIRT_HPP
#define GETNLEGIS_DYNIRT_HPP

#include <RcppArmadillo.h>

arma::mat getNlegis_dynIRT(const arma::mat &legis_by_session,
					const int T,
					const int N);


#endif
