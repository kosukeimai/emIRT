// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETETA_HIERIRT_H
#define GETETA_HIERIRT_H

#include <RcppArmadillo.h>

void getEta_hierIRT(arma::mat &Eta,
					arma::mat &Eta2,
					const arma::mat &Veta,
					const arma::mat &Eystar,
					const arma::mat &Eb,
					const arma::mat &Eba,
					const arma::mat &Ebb,
					const arma::mat &Egamma,
                    const arma::mat &z,
                    const arma::mat &g,
					const arma::mat &i,
                    const arma::mat &j,
                    const int ND,
                    const int NI,
                    const int NL
                    );

#endif
