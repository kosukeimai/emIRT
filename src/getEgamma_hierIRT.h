// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEGAMMA_HIERIRT_H
#define GETEGAMMA_HIERIRT_H

#include <RcppArmadillo.h>

void getEgamma_hierIRT(arma::mat &Egamma,
					const arma::cube &Vgamma,
					const arma::mat &gammasigma,
					const arma::mat &gammamu,
					const arma::mat &g,
					const arma::mat &i,
                    const arma::mat &j,
                    const arma::mat &z,
                    const arma::mat &Eb,
                    const arma::mat &Ebb,
                    const arma::mat &Eystar,
                    const arma::mat &Ea,
                    const arma::mat &Eta,
                    const int NL,
                    const int NG
                    );

#endif
