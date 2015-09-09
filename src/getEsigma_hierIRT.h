// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETESIGMA_HIERIRT_H
#define GETESIGMA_HIERIRT_H

#include <RcppArmadillo.h>

void getEsigma_hierIRT(arma::mat &Esigma,
					const arma::mat &Eta2,
					const arma::mat &sigmav,
					const arma::mat &sigmas,
					const arma::mat &g,
                    const int NG,
                    const int NI
                    );

#endif
