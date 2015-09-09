// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEB2_HIERIRT_H
#define GETEB2_HIERIRT_H

#include <RcppArmadillo.h>

void getEb2_hierIRT(arma::mat &Eb2,
					const arma::cube &Vb2,
					const arma::mat &betasigma,
					const arma::mat &betamu,
					const arma::mat &Eystar,
					const arma::mat &Ex2,
					const arma::mat &i,
                    const arma::mat &j,
                    const int NL,
                    const int NJ
                    );
#endif
