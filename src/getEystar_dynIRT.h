// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETYESTAR_DYNIRT_H
#define GETYESTAR_DYNIRT_H

#include <RcppArmadillo.h>

void getEystar_dynIRT(arma::mat &Eystar,
					const arma::mat &alpha,
                    const arma::mat &beta,
                    const arma::mat &x,
                    const arma::mat &y,
                    const arma::mat &bill_session,
                    const arma::mat &startlegis,
                    const arma::mat &endlegis,
                    const int N,
                    const int J
                    );

#endif
