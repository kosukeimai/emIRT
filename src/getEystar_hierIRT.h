// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETYESTAR_HIERIRT_H
#define GETYESTAR_HIERIRT_H

#include <RcppArmadillo.h>

void getEystar_hierIRT(arma::mat &Eystar,
					const arma::mat &y,
					const arma::mat &z,
                    const arma::mat &g,
                    const arma::mat &i,
                    const arma::mat &j,
                    const arma::mat &Ea,
                    const arma::mat &Eb,
                    const arma::mat &Egamma,
                    const arma::mat &Eta,
                    const int ND,
                    const int NG,
                    const int NI,
                    const int NJ,
                    const int NL
                    );

#endif
