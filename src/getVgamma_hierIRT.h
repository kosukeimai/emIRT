// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVGAMMA_HIERIRT_HPP
#define GETVGAMMA_HIERIRT_HPP

#include <RcppArmadillo.h>

void getVgamma_hierIRT(arma::cube &Vgamma,
					arma::mat &gammasigma,
					arma::mat &Ebb,
					const arma::mat &g,
					const arma::mat &i,
                    const arma::mat &j,
                    const arma::mat &z,
                    const int NL,
                    const int NG
                    );
#endif
