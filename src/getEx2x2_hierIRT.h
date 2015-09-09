// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEX2X2_HIERIRT_H
#define GETEX2X2_HIERIRT_H

#include <RcppArmadillo.h>

void getEx2x2_hierIRT(arma::cube &curEx2x2,
					const arma::mat &curEx2,
                    const arma::cube &Egg,
                    const arma::mat &Eg,
                    const arma::mat &Eta,
                    const arma::mat &Eta2,
                    const arma::mat &g,
                    const arma::mat &z,
                    const int NI
                    );


#endif
