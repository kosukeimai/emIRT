// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEX2_HIERIRT_H
#define GETEX2_HIERIRT_H

#include <RcppArmadillo.h>

void getEx2_hierIRT(arma::mat &curEx2,
					const arma::mat &Egamma,
                    const arma::mat &Eta,
                    const arma::mat &g,
                    const arma::mat &z,
                    const int NI
                    );


#endif
